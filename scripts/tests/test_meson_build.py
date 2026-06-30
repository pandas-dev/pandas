from pathlib import Path


ROOT_MESON_BUILD = Path(__file__).parents[2] / "meson.build"


def test_cython_type_specs_enabled_for_c_and_cpp() -> None:
    meson_build = ROOT_MESON_BUILD.read_text(encoding="utf-8")

    expected = (
        "add_project_arguments("
        "'-DCYTHON_USE_TYPE_SPECS=1', language: ['c', 'cpp'])"
    )
    assert expected in meson_build
