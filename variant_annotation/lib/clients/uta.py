"""UTA-backed TranscriptSource implementation."""

from __future__ import annotations

from typing import Any
from urllib.parse import urlsplit
import psycopg2  # type: ignore[import-untyped]

from ..accessions import transcript_sort_key, looks_like_transcript_accession


def connect_uta(database_url: str) -> Any:
    """Open a psycopg2 connection to the UTA database at database_url."""

    parsed = urlsplit(database_url)
    path_parts = [p for p in parsed.path.split("/") if p]
    if not path_parts:
        raise RuntimeError(f"UTA_DB_URL is missing a database name: {database_url}")

    kwargs: dict[str, Any] = {"dbname": path_parts[0]}
    if parsed.hostname:
        kwargs["host"] = parsed.hostname
    if parsed.port:
        kwargs["port"] = parsed.port
    if parsed.username:
        kwargs["user"] = parsed.username
    if parsed.password:
        kwargs["password"] = parsed.password

    return psycopg2.connect(**kwargs)


class UtaClient:
    """Satisfies TranscriptSource via a psycopg2 connection to UTA."""

    def __init__(self, connection: Any, *, schema: str = "uta_20241220") -> None:
        self._conn = connection
        self._schema = schema

    def transcript_for_protein(self, protein_accession: str) -> str | None:
        """Return the preferred coding transcript for a RefSeq protein accession."""
        with self._conn.cursor() as cursor:
            cursor.execute(
                f"""
                SELECT tx_ac
                FROM {self._schema}.associated_accessions
                WHERE pro_ac = %s
                """,
                (protein_accession,),
            )
            transcripts = sorted(
                {row[0] for row in cursor.fetchall() if row and row[0] and looks_like_transcript_accession(row[0])},
                key=transcript_sort_key,
            )
        return transcripts[0] if transcripts else None

    def codon_at(self, transcript: str, aa_position: int) -> str | None:
        """Return the 3-letter uppercase codon at aa_position (1-based) in transcript."""
        cds_offset = (aa_position - 1) * 3
        with self._conn.cursor() as cursor:
            cursor.execute(
                f"""
                SELECT SUBSTRING(s.seq FROM t.cds_start_i + %s + 1 FOR 3)
                FROM {self._schema}.transcript t
                JOIN {self._schema}.seq_anno sa ON sa.ac = t.ac
                JOIN {self._schema}.seq s ON s.seq_id = sa.seq_id
                WHERE t.ac = %s
                LIMIT 1
                """,
                (cds_offset, transcript),
            )
            row = cursor.fetchone()
        if row and row[0] and len(row[0]) == 3:
            return row[0].upper()
        return None

    def coding_interval_reference(self, transcript: str, start: int, stop: int) -> str | None:
        """Return the reference bases of transcript over the inclusive 1-based c. interval [start, stop].

        Satisfies the VEP ReferenceSequence port used to recognise reference-identical input. Reads the
        stored transcript sequence at the CDS offset (c.1 is the first CDS base), the same offset scheme
        as :meth:`codon_at`. Returns None when the interval cannot be resolved — unknown transcript, or
        the stored sequence is shorter than the requested interval. Only plain CDS coordinates reach here:
        the caller parses only integer c. intervals, so UTR/intronic offset forms (c.-12, c.*5, c.123+4)
        are never passed.
        """
        if start < 1 or stop < start:
            return None
        length = stop - start + 1
        with self._conn.cursor() as cursor:
            cursor.execute(
                f"""
                SELECT SUBSTRING(s.seq FROM t.cds_start_i + %s + 1 FOR %s)
                FROM {self._schema}.transcript t
                JOIN {self._schema}.seq_anno sa ON sa.ac = t.ac
                JOIN {self._schema}.seq s ON s.seq_id = sa.seq_id
                WHERE t.ac = %s
                LIMIT 1
                """,
                (start - 1, length, transcript),
            )
            row = cursor.fetchone()
        if row and row[0] and len(row[0]) == length:
            return row[0].upper()
        return None
