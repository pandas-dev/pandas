import io
import re
import subprocess
import sys

import pytest

pat = re.compile(r"python3")


@pytest.fixture(autouse=True)
def check_bufferedrandom_resourcewarning(request):
    yield
    lsof = subprocess.run(["lsof", "-d", "0-25"], capture_output=True).stdout.decode(
        "utf-8"
    )
    for line in lsof.split("\n"):
        if re.search(pat, line):
            # sys.stderr for xdist
            # https://github.com/pytest-dev/pytest/issues/1693#issuecomment-233282644
            print(f"{request.node.name}: {line}", flush=True, file=sys.stderr)


@pytest.fixture(autouse=True)
def no_font_loading(monkeypatch):
    import matplotlib

    class FT2FontCloser(matplotlib.ft2font.FT2Font):
        def __init__(self, path, *args, **kwargs):
            with open(path, "rb") as f:
                stream = io.BytesIO(f.read())
            super().__init__(stream, *args, **kwargs)

    monkeypatch.setattr("matplotlib.ft2font.FT2Font", FT2FontCloser)
