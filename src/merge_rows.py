"""Backward-compatible wrapper for ``src.replace_rows``.

Use ``python -m src.replace_rows`` for the renamed utility.
"""

from src.replace_rows import main, replace_rows as merge_rows


if __name__ == "__main__":
    main()
