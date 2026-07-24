import sys
from pathlib import Path

import pytest

ROOT_DIR = Path(__file__).parent.parent.parent


@pytest.mark.skipif(
    sys.version_info < (3, 11), reason="Requires Python 3.11+ for tomllib"
)
@pytest.mark.skipif(
    not (ROOT_DIR / "pixi.toml").exists(),
    reason="Test can only be run while the project is deployed from git",
)
def test_development_guidelines_matches_ci():
    """When pixi environments change, make sure to update the docs as well."""
    import tomllib

    with open(ROOT_DIR / "pixi.toml", "rb") as f:
        pixi_toml = tomllib.load(f)
    environments = set(pixi_toml["environments"])

    with open(ROOT_DIR / "docs/source/develop.rst", encoding="utf8") as f:
        docs = f.read()

    missing = environments - {"default", "build", "spark"} - set(docs.split())
    assert not missing, (
        "The following pixi environments are not listed in docs/source/develop.rst: "
        f"{missing}"
    )
