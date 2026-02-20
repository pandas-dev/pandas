# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.

"""CLI entry point for jupyterlab server."""
import sys

from jupyterlab_server.app import main

sys.exit(main())  # type:ignore[no-untyped-call]
