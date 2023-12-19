import pathlib
import sys

import pytest
import yaml

if sys.version_info >= (3, 11):
    import tomllib
else:
    import tomli as tomllib

from scripts.validate_min_versions_in_sync import (
    get_toml_map_from,
    get_yaml_map_from,
    pin_min_versions_to_yaml_file,
)


@pytest.mark.parametrize(
    "src_toml, src_yaml, expected_yaml",
    [
        (
            pathlib.Path("scripts/tests/data/deps_minimum.toml"),
            pathlib.Path("scripts/tests/data/deps_unmodified_random.yaml"),
            pathlib.Path("scripts/tests/data/deps_expected_random.yaml"),
        ),
        (
            pathlib.Path("scripts/tests/data/deps_minimum.toml"),
            pathlib.Path("scripts/tests/data/deps_unmodified_same_version.yaml"),
            pathlib.Path("scripts/tests/data/deps_expected_same_version.yaml"),
        ),
        (
            pathlib.Path("scripts/tests/data/deps_minimum.toml"),
            pathlib.Path("scripts/tests/data/deps_unmodified_duplicate_package.yaml"),
            pathlib.Path("scripts/tests/data/deps_expected_duplicate_package.yaml"),
        ),
        (
            pathlib.Path("scripts/tests/data/deps_minimum.toml"),
            pathlib.Path("scripts/tests/data/deps_unmodified_no_version.yaml"),
            pathlib.Path("scripts/tests/data/deps_expected_no_version.yaml"),
        ),
        (
            pathlib.Path("scripts/tests/data/deps_minimum.toml"),
            pathlib.Path("scripts/tests/data/deps_unmodified_range.yaml"),
            pathlib.Path("scripts/tests/data/deps_expected_range.yaml"),
        ),
    ],
)
def test_pin_min_versions_to_yaml_file(src_toml, src_yaml, expected_yaml) -> None:
    with open(src_toml, "rb") as toml_f:
        toml_map = tomllib.load(toml_f)
    with open(src_yaml, encoding="utf-8") as yaml_f:
        yaml_file_data = yaml_f.read()
    yaml_file = yaml.safe_load(yaml_file_data)
    yaml_dependencies = yaml_file["dependencies"]
    yaml_map = get_yaml_map_from(yaml_dependencies)
    toml_map = get_toml_map_from(toml_map)
    result_yaml_file = pin_min_versions_to_yaml_file(yaml_map, toml_map, yaml_file_data)
    with open(expected_yaml, encoding="utf-8") as yaml_f:
        dummy_yaml_expected_file_1 = yaml_f.read()
    assert result_yaml_file == dummy_yaml_expected_file_1
