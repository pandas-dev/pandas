import re
import subprocess
import sys

import pytest

pat = re.compile(r"python3")


@pytest.fixture(autouse=True)
def clear_font_filehandles(request):
    # https://github.com/matplotlib/matplotlib/issues/22017#issuecomment-998241017
    yield
    import matplotlib

    matplotlib.font_manager._get_font.cache_clear()
    lsof = subprocess.run(["lsof", "-d", "0-25"], capture_output=True).stdout.decode(
        "utf-8"
    )
    print(f"{request.node.name}\n{lsof}", flush=True, file=sys.stderr)
    # for line in lsof.split("\n"):
    #     if re.search(pat, line):
    #         # sys.stderr for xdist
    #         # https://github.com/pytest-dev/pytest/issues/1693#issuecomment-233282644
    #         print(f"{request.node.name}: {line}", flush=True, file=sys.stderr)
