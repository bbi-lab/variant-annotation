"""Clear ClinGen Redis cache keys for this project namespace."""

from __future__ import annotations

import argparse
import logging
import sys
from typing import Optional

from src.lib.clingen import clear_clingen_cache


logger = logging.getLogger(__name__)


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Clear ClinGen Allele Registry cache keys stored in Redis.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument(
        "--prefix",
        default=None,
        help="Cache key prefix to clear (defaults to CLINGEN_CACHE_PREFIX or clingen:v1).",
    )
    p.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        help="Logging verbosity level.",
    )
    return p


def main(argv: Optional[list[str]] = None) -> None:
    args = _build_parser().parse_args(argv)

    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
        datefmt="%H:%M:%S",
    )

    try:
        deleted = clear_clingen_cache(prefix=args.prefix)
    except Exception as exc:
        logger.error("Failed to clear ClinGen cache: %s", exc)
        sys.exit(1)

    logger.info("Deleted %d ClinGen cache key(s)", deleted)


if __name__ == "__main__":
    main()
