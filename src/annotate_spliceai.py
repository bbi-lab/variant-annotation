"""Annotate variants with SpliceAI scores.

Supports two modes:
1. precomputed (default): query user-provided bgzipped/tabix-indexed SpliceAI VCFs
2. compute: run local ``spliceai`` CLI (requires SpliceAI + TensorFlow installation)

For rows with pipe-delimited genomic HGVS values (for example, protein reverse-translation
outputs), each output SpliceAI column is emitted as a pipe-delimited list aligned to
``mapped_hgvs_g`` positions.
"""

from __future__ import annotations

import argparse
import csv
from itertools import islice
import logging
import os
import re
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Any, Callable, Optional

logger = logging.getLogger(__name__)

SPLICEAI_COLS = [
    "spliceai.ds_ag",
    "spliceai.ds_al",
    "spliceai.ds_dg",
    "spliceai.ds_dl",
    "spliceai.dp_ag",
    "spliceai.dp_al",
    "spliceai.dp_dg",
    "spliceai.dp_dl",
    "spliceai.max_delta_score",
]

NC_TO_CHROM_GRCH38 = {
    "NC_000001.11": "1",
    "NC_000002.12": "2",
    "NC_000003.12": "3",
    "NC_000004.12": "4",
    "NC_000005.10": "5",
    "NC_000006.12": "6",
    "NC_000007.14": "7",
    "NC_000008.11": "8",
    "NC_000009.12": "9",
    "NC_000010.11": "10",
    "NC_000011.10": "11",
    "NC_000012.12": "12",
    "NC_000013.11": "13",
    "NC_000014.9": "14",
    "NC_000015.10": "15",
    "NC_000016.10": "16",
    "NC_000017.11": "17",
    "NC_000018.12": "18",
    "NC_000019.10": "19",
    "NC_000020.11": "20",
    "NC_000021.9": "21",
    "NC_000022.11": "22",
    "NC_000023.11": "X",
    "NC_000024.10": "Y",
    "NC_012920.1": "MT",
}

NC_TO_CHROM_GRCH37 = {
    "NC_000001.10": "1",
    "NC_000002.11": "2",
    "NC_000003.11": "3",
    "NC_000004.11": "4",
    "NC_000005.9": "5",
    "NC_000006.11": "6",
    "NC_000007.13": "7",
    "NC_000008.10": "8",
    "NC_000009.11": "9",
    "NC_000010.10": "10",
    "NC_000011.9": "11",
    "NC_000012.11": "12",
    "NC_000013.10": "13",
    "NC_000014.8": "14",
    "NC_000015.9": "15",
    "NC_000016.9": "16",
    "NC_000017.10": "17",
    "NC_000018.9": "18",
    "NC_000019.9": "19",
    "NC_000020.10": "20",
    "NC_000021.8": "21",
    "NC_000022.10": "22",
    "NC_000023.10": "X",
    "NC_000024.9": "Y",
    "NC_012920.1": "MT",
}


def get_nc_to_chrom_map(annotation: str) -> dict[str, str]:
    if annotation.lower() in ("grch37", "hg19"):
        return NC_TO_CHROM_GRCH37
    return NC_TO_CHROM_GRCH38


def _empty_score_dict() -> dict[str, str]:
    return {c: "" for c in SPLICEAI_COLS}


def _fmt_score(value: float) -> str:
    return f"{value:.4f}"


def split_pipe_preserve_positions(value: str) -> list[str]:
    raw = value or ""
    if "|" not in raw:
        return [raw.strip()]
    return [part.strip() for part in raw.split("|")]


def parse_spliceai_info(info: str, alt: Optional[str] = None) -> dict[str, str]:
    """Parse SpliceAI INFO and return output-column score values.

    SpliceAI INFO entry format:
    ALLELE|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL

    When multiple entries are present, this function takes the maximum per score field
    across entries (filtered by ALT allele when provided).
    """
    match = re.search(r"SpliceAI=([^;\s]+)", info)
    if not match:
        return _empty_score_dict()

    max_vals: dict[str, float] = {}
    found = False
    for entry in match.group(1).split(","):
        parts = entry.split("|")
        if len(parts) < 10:
            continue
        allele = parts[0].upper()
        if alt is not None and allele and allele != alt.upper():
            continue
        try:
            ds_ag = float(parts[2])
            ds_al = float(parts[3])
            ds_dg = float(parts[4])
            ds_dl = float(parts[5])
            dp_ag = float(parts[6])
            dp_al = float(parts[7])
            dp_dg = float(parts[8])
            dp_dl = float(parts[9])
        except ValueError:
            continue

        found = True
        current = {
            "spliceai.ds_ag": ds_ag,
            "spliceai.ds_al": ds_al,
            "spliceai.ds_dg": ds_dg,
            "spliceai.ds_dl": ds_dl,
            "spliceai.dp_ag": dp_ag,
            "spliceai.dp_al": dp_al,
            "spliceai.dp_dg": dp_dg,
            "spliceai.dp_dl": dp_dl,
        }
        for key, value in current.items():
            if key not in max_vals or value > max_vals[key]:
                max_vals[key] = value

    if not found:
        return _empty_score_dict()

    ds_max = max(
        max_vals["spliceai.ds_ag"],
        max_vals["spliceai.ds_al"],
        max_vals["spliceai.ds_dg"],
        max_vals["spliceai.ds_dl"],
    )
    out = {
        "spliceai.ds_ag": _fmt_score(max_vals["spliceai.ds_ag"]),
        "spliceai.ds_al": _fmt_score(max_vals["spliceai.ds_al"]),
        "spliceai.ds_dg": _fmt_score(max_vals["spliceai.ds_dg"]),
        "spliceai.ds_dl": _fmt_score(max_vals["spliceai.ds_dl"]),
        "spliceai.dp_ag": _fmt_score(max_vals["spliceai.dp_ag"]),
        "spliceai.dp_al": _fmt_score(max_vals["spliceai.dp_al"]),
        "spliceai.dp_dg": _fmt_score(max_vals["spliceai.dp_dg"]),
        "spliceai.dp_dl": _fmt_score(max_vals["spliceai.dp_dl"]),
        "spliceai.max_delta_score": _fmt_score(ds_max),
    }
    return out


def parse_hgvs_g_to_vcf(
    hgvs_g: str,
    nc_to_chrom: dict[str, str],
    fasta: Optional[Any] = None,
) -> Optional[tuple[str, int, str, str]]:
    m = re.match(r"^(NC_\d+\.\d+):g\.(.+)$", hgvs_g.strip())
    if not m:
        return None

    nc_acc, expr = m.group(1), m.group(2)
    chrom = nc_to_chrom.get(nc_acc)
    if not chrom:
        return None

    snv = re.match(r"^(\d+)([ACGTacgt]+)>([ACGTacgt]+)$", expr)
    if snv:
        return chrom, int(snv.group(1)), snv.group(2).upper(), snv.group(3).upper()

    del_single = re.match(r"^(\d+)del$", expr)
    if del_single and fasta is not None:
        pos = int(del_single.group(1))
        pad = str(fasta[chrom][pos - 2 : pos - 1]).upper()
        deleted = str(fasta[chrom][pos - 1 : pos]).upper()
        return chrom, pos - 1, pad + deleted, pad

    del_range = re.match(r"^(\d+)_(\d+)del$", expr)
    if del_range and fasta is not None:
        start = int(del_range.group(1))
        end = int(del_range.group(2))
        pad = str(fasta[chrom][start - 2 : start - 1]).upper()
        deleted = str(fasta[chrom][start - 1 : end]).upper()
        return chrom, start - 1, pad + deleted, pad

    ins = re.match(r"^(\d+)_(\d+)ins([ACGTacgt]+)$", expr)
    if ins and fasta is not None:
        pos = int(ins.group(1))
        inserted = ins.group(3).upper()
        ref_base = str(fasta[chrom][pos - 1 : pos]).upper()
        return chrom, pos, ref_base, ref_base + inserted

    delins = re.match(r"^(\d+)(?:_(\d+))?delins([ACGTacgt]+)$", expr)
    if delins and fasta is not None:
        start = int(delins.group(1))
        end = int(delins.group(2)) if delins.group(2) else start
        alt = delins.group(3).upper()
        ref = str(fasta[chrom][start - 1 : end]).upper()
        return chrom, start, ref, alt

    dup = re.match(r"^(\d+)(?:_(\d+))?dup$", expr)
    if dup and fasta is not None:
        start = int(dup.group(1))
        end = int(dup.group(2)) if dup.group(2) else start
        pad = str(fasta[chrom][start - 2 : start - 1]).upper()
        dup_seq = str(fasta[chrom][start - 1 : end]).upper()
        return chrom, start - 1, pad, pad + dup_seq

    return None


def tabix_query_for_hgvs_g(
    hgvs_g: str,
    nc_to_chrom: dict[str, str],
) -> Optional[tuple[str, int, Callable[[list[str]], bool]]]:
    """Create a tabix position query + line matcher without requiring FASTA."""
    m = re.match(r"^(NC_\d+\.\d+):g\.(.+)$", hgvs_g.strip())
    if not m:
        return None
    nc_acc, expr = m.group(1), m.group(2)
    chrom = nc_to_chrom.get(nc_acc)
    if not chrom:
        return None

    snv = re.match(r"^(\d+)([ACGTacgt]+)>([ACGTacgt]+)$", expr)
    if snv:
        pos, ref, alt = int(snv.group(1)), snv.group(2).upper(), snv.group(3).upper()

        def _match_snv(f: list[str], p: int = pos, r: str = ref, a: str = alt) -> bool:
            return int(f[1]) == p and f[3] == r and f[4] == a

        return chrom, pos, _match_snv

    del_single = re.match(r"^(\d+)del$", expr)
    if del_single:
        vcf_pos = int(del_single.group(1)) - 1

        def _match_del_single(f: list[str], p: int = vcf_pos) -> bool:
            return int(f[1]) == p and len(f[3]) == 2 and len(f[4]) == 1

        return chrom, vcf_pos, _match_del_single

    del_range = re.match(r"^(\d+)_(\d+)del$", expr)
    if del_range:
        start, end = int(del_range.group(1)), int(del_range.group(2))
        vcf_pos = start - 1
        del_len = end - start + 1

        def _match_del_range(f: list[str], p: int = vcf_pos, d: int = del_len) -> bool:
            return int(f[1]) == p and len(f[3]) - len(f[4]) == d and len(f[4]) == 1

        return (
            chrom,
            vcf_pos,
            _match_del_range,
        )

    ins = re.match(r"^(\d+)_(\d+)ins([ACGTacgt]+)$", expr)
    if ins:
        vcf_pos = int(ins.group(1))
        inserted = ins.group(3).upper()

        def _match_ins(f: list[str], p: int = vcf_pos, i: str = inserted) -> bool:
            return int(f[1]) == p and len(f[3]) == 1 and f[4] == f[3] + i

        return (
            chrom,
            vcf_pos,
            _match_ins,
        )

    delins = re.match(r"^(\d+)(?:_(\d+))?delins([ACGTacgt]+)$", expr)
    if delins:
        start = int(delins.group(1))
        end = int(delins.group(2)) if delins.group(2) else start
        del_len = end - start + 1
        alt = delins.group(3).upper()

        def _match_delins(f: list[str], p: int = start, d: int = del_len, a: str = alt) -> bool:
            return int(f[1]) == p and len(f[3]) == d and f[4] == a

        return (
            chrom,
            start,
            _match_delins,
        )

    dup = re.match(r"^(\d+)(?:_(\d+))?dup$", expr)
    if dup:
        start = int(dup.group(1))
        end = int(dup.group(2)) if dup.group(2) else start
        dup_len = end - start + 1
        vcf_pos = start - 1

        def _match_dup(f: list[str], p: int = vcf_pos, d: int = dup_len) -> bool:
            return int(f[1]) == p and len(f[4]) - len(f[3]) == d and len(f[3]) == 1

        return (
            chrom,
            vcf_pos,
            _match_dup,
        )

    return None


def _chrom_candidates(chrom: str) -> list[str]:
    candidates = [chrom]
    if chrom.startswith("chr"):
        candidates.append(chrom[3:])
    else:
        candidates.append(f"chr{chrom}")
    # Preserve order and uniqueness.
    return list(dict.fromkeys(candidates))


def _run_tabix_fetch_lines(vcf_path: Path, chrom: str, pos: int) -> list[str]:
    region = f"{chrom}:{pos}-{pos}"
    proc = subprocess.run(
        ["tabix", str(vcf_path), region],
        capture_output=True,
        text=True,
        check=False,
    )
    if proc.returncode != 0:
        return []
    return [line for line in proc.stdout.splitlines() if line and not line.startswith("#")]


def ensure_tabix_index(vcf_path: Path, force: bool = False) -> Path:
    tbi = vcf_path.with_suffix(vcf_path.suffix + ".tbi")
    if tbi.exists() and not force:
        return tbi

    logger.info("Indexing precomputed SpliceAI VCF: %s", vcf_path)
    proc = subprocess.run(
        ["tabix", "-f", "-p", "vcf", str(vcf_path)],
        capture_output=True,
        text=True,
        check=False,
    )
    if proc.returncode != 0:
        raise RuntimeError(
            f"tabix failed for {vcf_path}: {proc.stderr.strip() or proc.stdout.strip()}"
        )
    if not tbi.exists():
        raise RuntimeError(f"tabix did not produce index file for {vcf_path}")
    return tbi


def cache_precomputed_vcf(src_path: Path, cache_dir: Path, refresh: bool = False) -> Path:
    if not src_path.exists():
        raise FileNotFoundError(f"Precomputed VCF not found: {src_path}")
    cache_dir.mkdir(parents=True, exist_ok=True)

    dest = cache_dir / src_path.name
    if src_path.resolve() != dest.resolve():
        if refresh or not dest.exists():
            logger.info("Caching precomputed SpliceAI VCF %s -> %s", src_path, dest)
            shutil.copy2(src_path, dest)
    return dest


def lookup_precomputed_scores(
    hgvs_to_exact_vcf: dict[str, Optional[tuple[str, int, str, str]]],
    nc_to_chrom: dict[str, str],
    vcf_paths: list[Path],
) -> dict[str, dict[str, str]]:
    out: dict[str, dict[str, str]] = {}
    for hgvs_g, exact in hgvs_to_exact_vcf.items():
        matcher: Callable[[list[str]], bool]
        if exact is not None:
            chrom, pos, ref, alt = exact

            def _match_exact(f: list[str], p: int = pos, r: str = ref, a: str = alt) -> bool:
                return int(f[1]) == p and f[3] == r and f[4] == a

            matcher = _match_exact
            query_chrom = chrom
            query_pos = pos
            alt_for_parse = alt
        else:
            query = tabix_query_for_hgvs_g(hgvs_g, nc_to_chrom)
            if query is None:
                out[hgvs_g] = _empty_score_dict()
                continue
            query_chrom, query_pos, matcher = query
            alt_for_parse = None

        found_score = _empty_score_dict()
        found = False
        for vcf_path in vcf_paths:
            for chrom in _chrom_candidates(query_chrom):
                lines = _run_tabix_fetch_lines(vcf_path, chrom, query_pos)
                for line in lines:
                    fields = line.split("\t")
                    if len(fields) < 8:
                        continue
                    try:
                        if matcher(fields):
                            found_score = parse_spliceai_info(fields[7], alt=alt_for_parse)
                            found = True
                            break
                    except (ValueError, IndexError):
                        continue
                if found:
                    break
            if found:
                break

        out[hgvs_g] = found_score
    return out


def run_spliceai(
    vcf_variants: list[tuple[str, int, str, str]],
    genome_fasta: str,
    annotation: str,
) -> dict[tuple[str, int, str, str], dict[str, str]]:
    if not vcf_variants:
        return {}

    in_path = None
    out_path = None
    try:
        with tempfile.NamedTemporaryFile(suffix=".vcf", mode="w", delete=False) as in_vcf:
            in_path = in_vcf.name
            in_vcf.write("##fileformat=VCFv4.2\n")
            for chrom in sorted({chrom for chrom, *_ in vcf_variants}):
                in_vcf.write(f"##contig=<ID={chrom}>\n")
            in_vcf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
            for chrom, pos, ref, alt in vcf_variants:
                in_vcf.write(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\t.\t.\n")

        with tempfile.NamedTemporaryFile(suffix=".vcf", delete=False) as out_vcf:
            out_path = out_vcf.name

        proc = subprocess.run(
            [
                "spliceai",
                "-I",
                in_path,
                "-O",
                out_path,
                "-R",
                genome_fasta,
                "-A",
                annotation,
            ],
            capture_output=True,
            text=True,
            check=False,
        )
        if proc.returncode != 0:
            raise RuntimeError(f"spliceai failed: {proc.stderr.strip() or proc.stdout.strip()}")

        scores: dict[tuple[str, int, str, str], dict[str, str]] = {}
        with open(out_path, encoding="utf-8") as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 8:
                    continue
                key = (fields[0], int(fields[1]), fields[3], fields[4])
                scores[key] = parse_spliceai_info(fields[7], alt=fields[4])
        return scores
    finally:
        if in_path:
            Path(in_path).unlink(missing_ok=True)
        if out_path:
            Path(out_path).unlink(missing_ok=True)


def annotate_row_with_scores(
    row: dict[str, str],
    hgvs_g_col: str,
    hgvs_score_map: dict[str, dict[str, str]],
) -> dict[str, str]:
    hgvs_parts = split_pipe_preserve_positions(row.get(hgvs_g_col, ""))
    values: dict[str, list[str]] = {col: [] for col in SPLICEAI_COLS}
    for hgvs_g in hgvs_parts:
        score = hgvs_score_map.get(hgvs_g, _empty_score_dict()) if hgvs_g else _empty_score_dict()
        for col in SPLICEAI_COLS:
            values[col].append(score.get(col, ""))
    return {col: "|".join(parts) for col, parts in values.items()}


def _parse_args(argv: Optional[list[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=(
            "Annotate rows with SpliceAI delta scores using precomputed VCF lookups "
            "(default) or local SpliceAI computation."
        )
    )
    p.add_argument("input_file", help="Input CSV/TSV file")
    p.add_argument("output_file", help="Output CSV/TSV file")
    p.add_argument(
        "--mode",
        choices=["precomputed", "compute"],
        default="precomputed",
        help="Annotation mode (default: precomputed)",
    )
    p.add_argument(
        "--annotation",
        default="grch38",
        help="Genome build label for NC accession mapping and spliceai -A (default: grch38)",
    )
    p.add_argument(
        "--hgvs-g-col",
        default="mapped_hgvs_g",
        help="Column containing genomic HGVS values (pipe-delimited candidates allowed)",
    )
    p.add_argument("--delimiter", default="\t", help="Input/output delimiter (default TAB)")
    p.add_argument(
        "--cache-dir",
        default=os.environ.get("SPLICEAI_CACHE_DIR", "/tmp/spliceai_cache"),
        help="Cache directory for precomputed VCFs and indices",
    )
    p.add_argument(
        "--precomputed-snv-vcf",
        default=os.environ.get("SPLICEAI_PRECOMPUTED_SNV_VCF", ""),
        help="Path to precomputed SpliceAI SNV VCF (.vcf.gz)",
    )
    p.add_argument(
        "--precomputed-indel-vcf",
        default=os.environ.get("SPLICEAI_PRECOMPUTED_INDEL_VCF", ""),
        help="Path to precomputed SpliceAI indel VCF (.vcf.gz)",
    )
    p.add_argument(
        "--precomputed-vcf",
        action="append",
        default=[],
        help="Additional precomputed SpliceAI VCFs to search in order",
    )
    p.add_argument(
        "--genome",
        default="",
        help="Reference FASTA for compute mode and exact indel normalization in precomputed mode",
    )
    p.add_argument(
        "--prepare-cache-only",
        action="store_true",
        help="Only cache/index precomputed files; do not annotate",
    )
    p.add_argument(
        "--refresh-cache",
        action="store_true",
        help="Re-copy source precomputed VCFs into cache and rebuild .tbi indices",
    )
    p.add_argument(
        "--skip",
        type=int,
        default=0,
        help="Number of data rows to skip before annotation",
    )
    p.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Maximum number of data rows to annotate",
    )
    p.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        help="Logging verbosity",
    )
    return p.parse_args(argv)


def _precomputed_paths(args: argparse.Namespace) -> list[Path]:
    paths: list[str] = []
    if args.precomputed_snv_vcf:
        paths.append(args.precomputed_snv_vcf)
    if args.precomputed_indel_vcf:
        paths.append(args.precomputed_indel_vcf)
    paths.extend(args.precomputed_vcf or [])
    return [Path(p) for p in paths if p]


def _maybe_load_fasta(genome_path: str) -> Optional[object]:
    if not genome_path:
        return None
    try:
        import pyfaidx  # type: ignore
    except Exception:
        logger.warning(
            "pyfaidx is not installed; exact indel normalization from --genome will be unavailable"
        )
        return None
    return pyfaidx.Fasta(genome_path)


def _cache_and_index_precomputed_files(args: argparse.Namespace) -> list[Path]:
    precomputed_input_paths = _precomputed_paths(args)
    if not precomputed_input_paths:
        raise ValueError(
            "Precomputed mode requires --precomputed-snv-vcf/--precomputed-indel-vcf "
            "or one or more --precomputed-vcf arguments"
        )

    cache_dir = Path(args.cache_dir)
    cached_paths: list[Path] = []
    for src in precomputed_input_paths:
        cached = cache_precomputed_vcf(src, cache_dir, refresh=args.refresh_cache)
        ensure_tabix_index(cached, force=args.refresh_cache)
        cached_paths.append(cached)
    return cached_paths


def main(argv: Optional[list[str]] = None) -> None:
    args = _parse_args(argv)
    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )
    if args.skip < 0:
        raise ValueError(f"--skip must be >= 0, got: {args.skip}")
    if args.limit is not None and args.limit < 1:
        raise ValueError(f"--limit must be >= 1 when provided, got: {args.limit}")

    delim = "\t" if args.delimiter == "\\t" else args.delimiter
    nc_to_chrom = get_nc_to_chrom_map(args.annotation)

    if args.mode == "precomputed":
        precomputed_vcfs = _cache_and_index_precomputed_files(args)
        if args.prepare_cache_only:
            logger.info("SpliceAI precomputed cache is ready: %s", args.cache_dir)
            return
    else:
        if args.prepare_cache_only:
            raise ValueError("--prepare-cache-only can only be used with --mode precomputed")
        if not args.genome:
            raise ValueError("--genome is required in --mode compute")
        precomputed_vcfs = []

    fasta = _maybe_load_fasta(args.genome)

    input_path = Path(args.input_file)
    output_path = Path(args.output_file)

    with input_path.open("r", encoding="utf-8", newline="") as in_fh:
        reader = csv.DictReader(in_fh, delimiter=delim)
        if reader.fieldnames is None:
            raise ValueError(f"Input file appears empty: {input_path}")
        fieldnames = list(reader.fieldnames)
        rows = list(
            islice(
                reader,
                args.skip,
                None if args.limit is None else args.skip + args.limit,
            )
        )

    unique_hgvs: set[str] = set()
    for row in rows:
        for hgvs_g in split_pipe_preserve_positions(row.get(args.hgvs_g_col, "")):
            if hgvs_g:
                unique_hgvs.add(hgvs_g)

    hgvs_to_exact_vcf: dict[str, Optional[tuple[str, int, str, str]]] = {}
    for hgvs_g in unique_hgvs:
        hgvs_to_exact_vcf[hgvs_g] = parse_hgvs_g_to_vcf(hgvs_g, nc_to_chrom, fasta)

    if args.mode == "precomputed":
        logger.info(
            "Looking up SpliceAI annotations for %d unique variants from %d cached VCF(s)",
            len(hgvs_to_exact_vcf),
            len(precomputed_vcfs),
        )
        hgvs_score_map = lookup_precomputed_scores(hgvs_to_exact_vcf, nc_to_chrom, precomputed_vcfs)
    else:
        unique_vcf_keys = list({v for v in hgvs_to_exact_vcf.values() if v is not None})
        logger.info("Running local SpliceAI on %d unique variants", len(unique_vcf_keys))
        vcf_score_map = run_spliceai(unique_vcf_keys, args.genome, args.annotation)
        hgvs_score_map = {
            hgvs_g: vcf_score_map.get(vcf_key, _empty_score_dict()) if vcf_key is not None else _empty_score_dict()
            for hgvs_g, vcf_key in hgvs_to_exact_vcf.items()
        }

    out_fieldnames = fieldnames + [col for col in SPLICEAI_COLS if col not in fieldnames]
    annotated_rows = 0
    with output_path.open("w", encoding="utf-8", newline="") as out_fh:
        writer = csv.DictWriter(
            out_fh,
            fieldnames=out_fieldnames,
            delimiter=delim,
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        for row in rows:
            ann = annotate_row_with_scores(row, args.hgvs_g_col, hgvs_score_map)
            row.update(ann)
            writer.writerow(row)
            if ann["spliceai.max_delta_score"].replace("|", "").strip():
                annotated_rows += 1

    logger.info("Done. %d/%d rows received SpliceAI annotations", annotated_rows, len(rows))


if __name__ == "__main__":
    main()