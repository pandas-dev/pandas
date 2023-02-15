import pathlib
import sys

import pytest

if sys.version_info >= (3, 11):
    import tomllib
else:
    import tomli as tomllib

from scripts.validate_min_versions_in_sync import update_yaml_file_version


@pytest.mark.parametrize(
    "src_toml, src_yaml, expected_yaml",
    [
        (  # random
            pathlib.Path("scripts/tests/dummy_test_files/dummy_toml_file.toml"),
            pathlib.Path(
                "scripts/tests/dummy_test_files/dummy_yaml_unmodified_file_1.yaml"
            ),
            pathlib.Path(
                "scripts/tests/dummy_test_files/dummy_yaml_expected_file_1.yaml"
            ),
        ),
        (  # same version
            pathlib.Path("scripts/tests/dummy_test_files/dummy_toml_file.toml"),
            pathlib.Path(
                "scripts/tests/dummy_test_files/dummy_yaml_unmodified_file_2.yaml"
            ),
            pathlib.Path(
                "scripts/tests/dummy_test_files/dummy_yaml_expected_file_2.yaml"
            ),
        ),
        (  # duplicate package
            pathlib.Path("scripts/tests/dummy_test_files/dummy_toml_file.toml"),
            pathlib.Path(
                "scripts/tests/dummy_test_files/dummy_yaml_unmodified_file_3.yaml"
            ),
            pathlib.Path(
                "scripts/tests/dummy_test_files/dummy_yaml_expected_file_3.yaml"
            ),
        ),
        (  # empty version
            pathlib.Path("scripts/tests/dummy_test_files/dummy_toml_file.toml"),
            pathlib.Path(
                "scripts/tests/dummy_test_files/dummy_yaml_unmodified_file_4.yaml"
            ),
            pathlib.Path(
                "scripts/tests/dummy_test_files/dummy_yaml_expected_file_4.yaml"
            ),
        ),
        (  # range version
            pathlib.Path("scripts/tests/dummy_test_files/dummy_toml_file.toml"),
            pathlib.Path(
                "scripts/tests/dummy_test_files/dummy_yaml_unmodified_file_5.yaml"
            ),
            pathlib.Path(
                "scripts/tests/dummy_test_files/dummy_yaml_expected_file_5.yaml"
            ),
        ),
    ],
)
def test_update_yaml_file_version(src_toml, src_yaml, expected_yaml):
    with open(src_toml, "rb") as toml_f:
        toml_dic = tomllib.load(toml_f)
    with open(src_yaml) as yaml_f:
        result_yaml_file = update_yaml_file_version(yaml_f, toml_dic)
    with open(expected_yaml) as yaml_f:
        dummy_yaml_expected_file_1 = yaml_f.read()
    assert result_yaml_file == dummy_yaml_expected_file_1
