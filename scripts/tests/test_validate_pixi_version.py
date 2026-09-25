from pathlib import Path

from scripts.validate_pixi_version import check_pixi_versions


def _write(path: Path, *lines: str) -> None:
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def test_consistent_pins_pass(tmp_path: Path) -> None:
    action = tmp_path / "action.yml"
    workflow = tmp_path / "workflow.yml"
    _write(
        action,
        "        pixi-version: v0.76.2",
        "        pixi-version: v0.76.2",
    )
    _write(workflow, "        pixi-version: v0.76.2")
    assert check_pixi_versions([action, workflow]) == 0


def test_mismatched_pins_fail(tmp_path: Path) -> None:
    action = tmp_path / "action.yml"
    workflow = tmp_path / "workflow.yml"
    _write(
        action,
        "        pixi-version: v0.76.2",
        "        pixi-version: v0.76.2",
    )
    _write(workflow, "        pixi-version: v0.76.3")
    assert check_pixi_versions([action, workflow]) == 1


def test_missing_pin_fails(tmp_path: Path) -> None:
    action = tmp_path / "action.yml"
    _write(action, "        run-install: false")
    assert check_pixi_versions([action]) == 1


def test_missing_file_fails(tmp_path: Path) -> None:
    assert check_pixi_versions([tmp_path / "nope.yml"]) == 1
