import io

import matplotlib
import pytest


class FT2FontCloser(matplotlib.ft2font.FT2Font):
    def __init__(self, path, *args, **kwargs):
        with open(path, "rb") as f:
            stream = io.BytesIO(f.read())
        super().__init__(stream, *args, **kwargs)


@pytest.fixture(autouse=True)
def no_font_loading(monkeypatch):
    """Remove requests.sessions.Session.request for all tests."""
    monkeypatch.setattr("matplotlib.ft2font.FT2Font", FT2FontCloser)
