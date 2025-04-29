"""CLI entry point for notebook."""
import sys

from notebook.app import main

sys.exit(main())  # type:ignore[no-untyped-call]
