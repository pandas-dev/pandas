# Licensed under the Apache License: http://www.apache.org/licenses/LICENSE-2.0
# For details: https://github.com/coveragepy/coveragepy/blob/main/NOTICE.txt

# pylint: disable=missing-module-docstring
# pragma: exclude file from coverage
# This will become the .pth file for subprocesses.

import os

if os.getenv("COVERAGE_PROCESS_START") or os.getenv("COVERAGE_PROCESS_CONFIG"):
    try:
        import coverage
    except:  # pylint: disable=bare-except
        pass
    else:
        coverage.process_startup(slug="pth")
