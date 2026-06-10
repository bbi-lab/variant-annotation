"""Clear ClinGen Allele Registry cache keys from Redis.

The pipeline caches ClinGen API responses in Redis under a key namespace
prefix (default: ``clingen:v1``, configurable via the ``CLINGEN_CACHE_PREFIX``
environment variable). Use this script to invalidate the cache when you need
to force fresh API lookups — for example after a ClinGen data update, when
troubleshooting stale allele IDs, or when switching to a different ClinGen
environment.

All keys that begin with the chosen prefix are deleted atomically. The
number of deleted keys is reported.

Requires the Redis service to be running (see docker compose).
"""

from __future__ import annotations

import argparse
import logging
import sys
from typing import Optional

from variant_annotation.lib.clingen import clear_clingen_cache


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
