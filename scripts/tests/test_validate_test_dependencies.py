import pytest

from scripts.validate_test_dependencies import (
    ALLOWED_UNDECLARED,
    declared_packages,
    environments_running_tests,
    validate_test_dependencies,
)

PIXI = {
    "dependencies": {"python-dateutil": ">=2.9.0"},
    "feature": {
        "test-base": {"dependencies": {"pytest": ">=8.3.4"}},
        "optional-dependencies": {"dependencies": {"scipy": ">=1.16.1"}},
        "downstream": {"dependencies": {"dask": "*"}},
        "nightly": {"pypi-dependencies": {"numpy": "*"}},
        "documentation": {"dependencies": {"sphinx": "*"}},
    },
    "environments": {
        "py313": {"features": ["test-base", "optional-dependencies"]},
        "downstream": {"features": ["test-base", "downstream"]},
        "numpy-nightly": ["test-base", "nightly"],
        "docs": {"features": ["documentation"]},
    },
}


def test_declared_packages_unions_the_given_environments() -> None:
    result = declared_packages(PIXI, ["py313", "downstream"])
    assert result == {"python-dateutil", "pytest", "scipy", "dask"}


def test_declared_packages_excludes_other_environments() -> None:
    # sphinx is only installed in an environment that does not run pytest
    assert "sphinx" not in declared_packages(PIXI, ["py313", "downstream"])


def test_declared_packages_includes_pypi_dependencies() -> None:
    assert "numpy" in declared_packages(PIXI, ["numpy-nightly"])


def test_declared_packages_accepts_a_bare_feature_list() -> None:
    # pixi allows both `env = [...]` and `env = { features = [...] }`
    assert declared_packages(PIXI, ["numpy-nightly"]) == {
        "python-dateutil",
        "pytest",
        "numpy",
    }


def test_environments_running_tests_reads_the_workflow_matrix() -> None:
    environments = environments_running_tests()
    # the plain matrix entries, and one that only appears under `include`
    assert {"py311", "py312", "py313", "py314", "downstream"} <= environments


def test_repository_test_dependencies_are_declared(capsys) -> None:
    assert validate_test_dependencies() == 0
    assert capsys.readouterr().out == ""


@pytest.mark.parametrize("module", sorted(ALLOWED_UNDECLARED))
def test_allowed_undeclared_entries_have_a_reason(module) -> None:
    assert ALLOWED_UNDECLARED[module].strip()
