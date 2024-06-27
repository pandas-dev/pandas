import pytest
import traitlets

from jupyter_lsp.session import LanguageServerSession


@pytest.mark.parametrize(
    "spec",
    [
        {},
        {"argv": [], "languages": []},
        {"languages": None},
        {"languages": 1},
        {"languages": [1, "two"]},
    ],
)
def test_bad_spec(spec):
    with pytest.raises(traitlets.TraitError):
        LanguageServerSession(spec=spec)
