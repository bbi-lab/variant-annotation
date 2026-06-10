"""HGVS-mapper-backed CoordinateTranslator implementation."""

from __future__ import annotations

from typing import Any
import hgvs.assemblymapper  # type: ignore[import-untyped]
import hgvs.parser  # type: ignore[import-untyped]
import hgvs.dataproviders.uta  # type: ignore[import-untyped]


class HgvsMapper:
    """Satisfies CoordinateTranslator via the biocommons hgvs AssemblyMapper."""

    def __init__(self, hdp: Any, *, assembly: str = "GRCh38") -> None:
        self._parser = hgvs.parser.Parser()
        self._mapper = hgvs.assemblymapper.AssemblyMapper(
            hdp,
            assembly_name=assembly,
            alt_aln_method="splign",
        )

    @classmethod
    def from_url(cls, uta_url: str, *, assembly: str = "GRCh38") -> "HgvsMapper":
        """Construct an HgvsMapper from a UTA database URL."""
        hdp = hgvs.dataproviders.uta.connect(uta_url)
        return cls(hdp, assembly=assembly)

    def c_to_p(self, c_hgvs: str) -> str:
        var = self._parser.parse(c_hgvs)
        return str(self._mapper.c_to_p(var))

    def g_to_c(self, g_hgvs: str, transcript: str) -> str:
        var = self._parser.parse(g_hgvs)
        c_var = self._mapper.g_to_t(var, transcript)
        return str(c_var)

    def c_to_g(self, c_hgvs: str) -> str:
        var = self._parser.parse(c_hgvs)
        return str(self._mapper.c_to_g(var))
