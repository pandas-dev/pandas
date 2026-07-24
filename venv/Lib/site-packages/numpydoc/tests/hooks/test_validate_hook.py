"""Test the numpydoc validate pre-commit hook."""

import inspect
import re
from pathlib import Path

import pytest

from numpydoc.hooks.validate_docstrings import run_hook


@pytest.fixture
def example_module(request):
    fullpath = (
        Path(request.config.rootdir)
        / "numpydoc"
        / "tests"
        / "hooks"
        / "example_module.py"
    )
    return str(fullpath.relative_to(request.config.rootdir))


@pytest.mark.parametrize("config", [None, "fake_dir"])
def test_validate_hook(example_module, config, capsys):
    """Test that a file is correctly processed in the absence of config files."""

    expected = inspect.cleandoc(
        f"""
        {example_module!s}:4: ES01 No extended summary found

        {example_module!s}:4: PR01 Parameters {{'name'}} not documented

        {example_module!s}:4: SA01 See Also section not found

        {example_module!s}:4: EX01 No examples section found

        {example_module!s}:8: ES01 No extended summary found

        {example_module!s}:8: SA01 See Also section not found

        {example_module!s}:8: EX01 No examples section found

        {example_module!s}:11: GL08 The object does not have a docstring

        {example_module!s}:17: ES01 No extended summary found

        {example_module!s}:17: PR01 Parameters {{'**kwargs'}} not documented

        {example_module!s}:17: PR07 Parameter "*args" has no description

        {example_module!s}:17: SA01 See Also section not found

        {example_module!s}:17: EX01 No examples section found

        {example_module!s}:26: SS05 Summary must start with infinitive verb, not third person (e.g. use "Generate" instead of "Generates")

        {example_module!s}:26: ES01 No extended summary found

        {example_module!s}:26: SA01 See Also section not found

        {example_module!s}:26: EX01 No examples section found

        {example_module!s}:30: GL08 The object does not have a docstring

        {example_module!s}:31: SA01 See Also section not found

        {example_module!s}:31: EX01 No examples section found

        {example_module!s}:46: SA01 See Also section not found

        {example_module!s}:46: EX01 No examples section found

        {example_module!s}:58: ES01 No extended summary found

        {example_module!s}:58: PR01 Parameters {{'name'}} not documented

        {example_module!s}:58: SA01 See Also section not found

        {example_module!s}:58: EX01 No examples section found
        """
    )

    return_code = run_hook([example_module], config=config)
    assert return_code == 1
    assert capsys.readouterr().err.strip() == expected


def test_validate_hook_with_ignore(example_module, capsys):
    """
    Test that a file is correctly processed in the absence of config files
    with command line ignore options.
    """

    expected = inspect.cleandoc(
        f"""
        {example_module!s}:4: PR01 Parameters {{'name'}} not documented

        {example_module!s}:11: GL08 The object does not have a docstring

        {example_module!s}:17: PR01 Parameters {{'**kwargs'}} not documented

        {example_module!s}:17: PR07 Parameter "*args" has no description

        {example_module!s}:26: SS05 Summary must start with infinitive verb, not third person (e.g. use "Generate" instead of "Generates")

        {example_module!s}:30: GL08 The object does not have a docstring

        {example_module!s}:58: PR01 Parameters {{'name'}} not documented
        """
    )

    return_code = run_hook([example_module], ignore=["ES01", "SA01", "EX01"])

    assert return_code == 1
    assert capsys.readouterr().err.strip() == expected


def test_validate_hook_with_toml_config(example_module, tmp_path, capsys):
    """
    Test that a file is correctly processed with the config coming from
    a pyproject.toml file.
    """

    with open(tmp_path / "pyproject.toml", "w") as config_file:
        config_file.write(
            inspect.cleandoc(
                """
                [tool.numpydoc_validation]
                checks = [
                    "all",
                    "EX01",
                    "SA01",
                    "ES01",
                ]
                exclude = '\\.__init__$'
                override_SS05 = [
                    '^Creates',
                ]
                """
            )
        )

    expected = inspect.cleandoc(
        f"""
        {example_module!s}:4: PR01 Parameters {{'name'}} not documented

        {example_module!s}:17: PR01 Parameters {{'**kwargs'}} not documented

        {example_module!s}:17: PR07 Parameter "*args" has no description

        {example_module!s}:30: GL08 The object does not have a docstring
        """
    )

    return_code = run_hook([example_module], config=tmp_path)
    assert return_code == 1
    assert capsys.readouterr().err.strip() == expected


def test_validate_hook_with_setup_cfg(example_module, tmp_path, capsys):
    """
    Test that a file is correctly processed with the config coming from
    a setup.cfg file.
    """

    with open(tmp_path / "setup.cfg", "w") as config_file:
        config_file.write(
            inspect.cleandoc(
                """
                [tool:numpydoc_validation]
                checks = all,EX01,SA01,ES01
                exclude = \\.__init__$
                override_SS05 = ^Creates
                """
            )
        )

    expected = inspect.cleandoc(
        f"""
        {example_module!s}:4: PR01 Parameters {{'name'}} not documented

        {example_module!s}:17: PR01 Parameters {{'**kwargs'}} not documented

        {example_module!s}:17: PR07 Parameter "*args" has no description

        {example_module!s}:30: GL08 The object does not have a docstring
        """
    )

    return_code = run_hook([example_module], config=tmp_path)
    assert return_code == 1
    assert capsys.readouterr().err.strip() == expected


def test_validate_hook_exclude_option_pyproject(example_module, tmp_path, capsys):
    """
    Test that a file is correctly processed with the config coming from
    a pyproject.toml file and exclusions provided.
    """

    with open(tmp_path / "pyproject.toml", "w") as config_file:
        config_file.write(
            inspect.cleandoc(
                r"""
                [tool.numpydoc_validation]
                checks = [
                    "all",
                    "EX01",
                    "SA01",
                    "ES01",
                ]
                exclude = [
                    '\.do_something$',
                    '\.__init__$',
                ]
                override_SS05 = [
                    '^Creates',
                ]
                """
            )
        )

    expected = inspect.cleandoc(
        f"""
        {example_module!s}:4: PR01 Parameters {{'name'}} not documented

        {example_module!s}:30: GL08 The object does not have a docstring
        """
    )

    return_code = run_hook([example_module], config=tmp_path)
    assert return_code == 1
    assert capsys.readouterr().err.strip() == expected


def test_validate_hook_exclude_option_setup_cfg(example_module, tmp_path, capsys):
    """
    Test that a file is correctly processed with the config coming from
    a setup.cfg file and exclusions provided.
    """

    with open(tmp_path / "setup.cfg", "w") as config_file:
        config_file.write(
            inspect.cleandoc(
                """
                [tool:numpydoc_validation]
                checks = all,EX01,SA01,ES01
                exclude = \\.NewClass$,\\.__init__$
                override_SS05 = ^Creates
                """
            )
        )

    expected = inspect.cleandoc(
        f"""
        {example_module!s}:4: PR01 Parameters {{'name'}} not documented

        {example_module!s}:17: PR01 Parameters {{'**kwargs'}} not documented

        {example_module!s}:17: PR07 Parameter "*args" has no description
        """
    )

    return_code = run_hook([example_module], config=tmp_path)
    assert return_code == 1
    assert capsys.readouterr().err.strip() == expected


@pytest.mark.parametrize(
    "file_exists, expected_code",
    [(True, 0), (False, 1)],
)
def test_validate_hook_exclude_files_option_pyproject(
    example_module, file_exists, expected_code, tmp_path
):
    """
    Test that the hook correctly processes the toml config and either includes
    or excludes files based on the `exclude_files` option.
    """
    exclude = str(example_module) if file_exists else "does_not_exist.py"

    with open(tmp_path / "pyproject.toml", "w") as config_file:
        config_file.write(
            inspect.cleandoc(
                f"""
                [tool.numpydoc_validation]
                checks = [
                    "all",
                    "EX01",
                    "SA01",
                    "ES01",
                ]
                exclude = '\\.__init__$'
                override_SS05 = [
                    '^Creates',
                ]
                exclude_files = [
                    '{re.escape(exclude)}',
                ]"""
            )
        )

    return_code = run_hook([example_module], config=tmp_path)
    assert return_code == expected_code  # Should not-report/report findings.


@pytest.mark.parametrize(
    "file_exists, expected_code",
    [(True, 0), (False, 1)],
)
def test_validate_hook_exclude_files_option_setup_cfg(
    example_module, file_exists, expected_code, tmp_path
):
    """
    Test that the hook correctly processes the setup config and either includes
    or excludes files based on the `exclude_files` option.
    """
    exclude = str(example_module) if file_exists else "does_not_exist.py"

    with open(tmp_path / "setup.cfg", "w") as config_file:
        config_file.write(
            inspect.cleandoc(
                f"""
                [tool:numpydoc_validation]
                checks = all,EX01,SA01,ES01
                exclude = \\.NewClass$,\\.__init__$
                override_SS05 = ^Creates
                exclude_files = {re.escape(exclude)}
                """
            )
        )

    return_code = run_hook([example_module], config=tmp_path)
    assert return_code == expected_code  # Should not-report/report findings.
