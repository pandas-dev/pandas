"""Mypy type checker command line tool."""

from __future__ import annotations

import os
import sys
import traceback

from mypy.fscache import FileSystemCache
from mypy.main import main, process_options
from mypy.util import FancyFormatter


def console_entry() -> None:
    try:
        main()
        sys.stdout.flush()
        sys.stderr.flush()
    except BrokenPipeError:
        # Python flushes standard streams on exit; redirect remaining output
        # to devnull to avoid another BrokenPipeError at shutdown
        devnull = os.open(os.devnull, os.O_WRONLY)
        os.dup2(devnull, sys.stdout.fileno())
        sys.exit(2)
    except KeyboardInterrupt:
        # Setting fscache prevents bogus errors when --package-root is set.
        _, options = process_options(args=sys.argv[1:], fscache=FileSystemCache())
        if options.show_traceback:
            sys.stdout.write(traceback.format_exc())
        formatter = FancyFormatter(sys.stdout, sys.stderr, False)
        msg = "Interrupted\n"
        sys.stdout.write(formatter.style(msg, color="red", bold=True))
        sys.stdout.flush()
        sys.stderr.flush()
        sys.exit(2)
    except Exception as e:
        # Try reporting any uncaught error canonically, otherwise just flush the traceback.
        try:
            import mypy.errors

            _, options = process_options(args=sys.argv[1:], fscache=FileSystemCache())
            mypy.errors.report_internal_error(e, None, 0, None, options)
        except Exception:
            pass
        sys.stdout.write(traceback.format_exc())
        sys.stdout.flush()
        sys.stderr.flush()
        sys.exit(2)


if __name__ == "__main__":
    console_entry()
