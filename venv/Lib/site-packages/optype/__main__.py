"""Command-line entry point: ``optype infer EXPR [PARAM ...]``."""

import sys

from optype.infer._cli import run


def main() -> None:
    match sys.argv[1:]:
        case ["infer", *args]:
            run(args)
        case _:
            sys.exit("usage: optype infer EXPR [PARAM ...]")


if __name__ == "__main__":
    main()
