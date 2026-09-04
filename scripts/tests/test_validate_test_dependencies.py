import pytest

import scripts.validate_test_dependencies as vm
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
    "target": {
        "win": {"dependencies": {"tzdata": ">=2023.3"}},
    },
}


def test_declared_packages_unions_the_given_environments() -> None:
    everywhere, _ = declared_packages(PIXI, ["py313", "downstream"])
    assert everywhere == {"python-dateutil", "pytest", "scipy", "dask"}


def test_declared_packages_excludes_other_environments() -> None:
    # sphinx is only installed in an environment that does not run pytest
    everywhere, _ = declared_packages(PIXI, ["py313", "downstream"])
    assert "sphinx" not in everywhere


def test_declared_packages_includes_pypi_dependencies() -> None:
    everywhere, _ = declared_packages(PIXI, ["numpy-nightly"])
    assert "numpy" in everywhere


def test_declared_packages_accepts_a_bare_feature_list() -> None:
    # pixi allows both `env = [...]` and `env = { features = [...] }`
    everywhere, _ = declared_packages(PIXI, ["numpy-nightly"])
    assert everywhere == {"python-dateutil", "pytest", "numpy"}


def test_declared_packages_reports_platform_specific_targets() -> None:
    # [target.win.dependencies] installs tzdata only on Windows, so it is not
    # part of `everywhere` but is reported with the platform it runs on.
    everywhere, platform_only = declared_packages(PIXI, ["py313"])
    assert "tzdata" not in everywhere
    assert platform_only == {"tzdata": {"win"}}


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


def test_validate_reports_platform_only_gate_without_failing(
    monkeypatch, capsys, tmp_path
) -> None:
    # A gate on a package that is declared only for one platform (e.g. Windows)
    # runs there, so it is reported but must not fail the hook.
    pixi = tmp_path / "pixi.toml"
    pixi.write_text(
        "[dependencies]\n"
        'python-dateutil = ">=2.9.0"\n'
        "\n"
        "[environments]\n"
        "py313 = []\n"
        "\n"
        "[target.win.dependencies]\n"
        'tzdata = ">=2023.3"\n',
        encoding="utf-8",
    )
    monkeypatch.setattr(vm, "PIXI_PATH", pixi)
    monkeypatch.setattr(vm, "environments_running_tests", lambda: {"py313"})
    monkeypatch.setattr(vm, "gated_modules", lambda: {"tzdata": {"a.py"}})

    assert validate_test_dependencies() == 0
    out = capsys.readouterr().out
    assert "win" in out
    assert "never run" not in out


def test_validate_fails_when_gate_runs_nowhere(monkeypatch, capsys, tmp_path) -> None:
    pixi = tmp_path / "pixi.toml"
    pixi.write_text(
        '[dependencies]\npython-dateutil = ">=2.9.0"\n\n[environments]\npy313 = []\n',
        encoding="utf-8",
    )
    monkeypatch.setattr(vm, "PIXI_PATH", pixi)
    monkeypatch.setattr(vm, "environments_running_tests", lambda: {"py313"})
    monkeypatch.setattr(vm, "gated_modules", lambda: {"nonexistent": {"a.py"}})

    assert validate_test_dependencies() == 1
    out = capsys.readouterr().out
    assert "never run" in out
