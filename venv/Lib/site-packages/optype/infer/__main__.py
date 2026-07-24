"""``python -m optype.infer EXPR [PARAM ...]`` (also: the ``optype infer`` command)."""

import sys

from optype.infer._cli import run

if __name__ == "__main__":
    run(sys.argv[1:])
