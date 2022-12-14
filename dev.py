#!/usr/bin/env python
#
# Example stub for running `python -m dev.py`
#
# Copy this into your project root.

import os
import runpy
import sys

sys.path.remove(os.path.abspath(os.path.dirname(sys.argv[0])))
try:
    runpy.run_module("devpy", run_name="__main__")
except ImportError:
    print("Cannot import devpy; please install it using")
    print()
    print("  pip install git+https://github.com/scientific-python/devpy")
    print()
    sys.exit(1)
