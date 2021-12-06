"""Mypy type checker command line tool."""

import sys
import os

from mypy.main import main


def console_entry() -> None:
    try:
        main(None, sys.stdout, sys.stderr)
        sys.stdout.flush()
        sys.stderr.flush()
    except BrokenPipeError:
        # Python flushes standard streams on exit; redirect remaining output
        # to devnull to avoid another BrokenPipeError at shutdown
        devnull = os.open(os.devnull, os.O_WRONLY)
        os.dup2(devnull, sys.stdout.fileno())
        sys.exit(2)


if __name__ == '__main__':
    console_entry()
