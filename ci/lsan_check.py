"""
pytest plugin running LeakSanitizer's check inside each pytest-xdist worker.

LSan normally checks at process exit, but xdist kills workers that take more
than 10 seconds to exit, which an ASan build does after a long run, so the
check never happens. Use with ``-p lsan_check`` and ``ci/`` on PYTHONPATH.
"""

import ctypes
import sys


def pytest_sessionfinish(session):
    workerinput = getattr(session.config, "workerinput", None)
    if workerinput is not None:
        # symbols of the ASan runtime linked into the python executable
        ctypes.CDLL(None).__lsan_do_leak_check()
        print(
            f"lsan_check: leak check ran in {workerinput['workerid']}", file=sys.stderr
        )
