import pathlib
from textwrap import dedent
import tomllib

import pytest
import yaml

from scripts.validate_min_versions_in_sync import (
    get_toml_map_from,
    get_versions_from_ci,
    get_yaml_map_from,
    pin_min_versions_to_yaml_file,
)

DATA_PATH = pathlib.Path(__file__).parents[2] / "scripts/tests/data/"


@pytest.mark.parametrize(
    "src_toml, src_yaml, expected_yaml",
    [
        (
            DATA_PATH / "deps_minimum.toml",
            DATA_PATH / "deps_unmodified_random.yaml",
            DATA_PATH / "deps_expected_random.yaml",
        ),
        (
            DATA_PATH / "deps_minimum.toml",
            DATA_PATH / "deps_unmodified_same_version.yaml",
            DATA_PATH / "deps_expected_same_version.yaml",
        ),
        (
            DATA_PATH / "deps_minimum.toml",
            DATA_PATH / "deps_unmodified_duplicate_package.yaml",
            DATA_PATH / "deps_expected_duplicate_package.yaml",
        ),
        (
            DATA_PATH / "deps_minimum.toml",
            DATA_PATH / "deps_unmodified_no_version.yaml",
            DATA_PATH / "deps_expected_no_version.yaml",
        ),
        (
            DATA_PATH / "deps_minimum.toml",
            DATA_PATH / "deps_unmodified_range.yaml",
            DATA_PATH / "deps_expected_range.yaml",
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


def test_get_versions_from_ci_parses_pixi_toml() -> None:
    content = dedent(
        """
        [dependencies]
        # Build dependencies
        meson = ">=1.2.3,<2"
        pip = ">=26"

        # Required dependencies
        python-dateutil = ">=2.8.2"

        [feature.numpy.dependencies]
        numpy = ">=1.26.0,<3"

        [feature.numpy21.dependencies]
        numpy = "2.1.*"

        [feature.test-base.dependencies]
        pytest = ">=8.3.4"
        pytest-xdist = ">=3.6.1"
        hypothesis = ">=6.116.0"

        [feature.test-network.dependencies]
        boto3 = "==1.40.46"

        [feature.test-clipboard.dependencies]
        pyqt = ">=5.15.9"
        qtpy = ">=2.4.2"

        [feature.pyarrow.dependencies]
        pyarrow = ">=16.0.0"

        [feature.pyarrow21.dependencies]
        pyarrow = "21.*"

        [feature.optional-dependencies.dependencies]
        beautifulsoup4 = ">=4.12.3"
        blosc = ">=1.21.3"
        pytables = ">=3.10.1"
        zstandard = ">=0.23.0"
        """
    ).splitlines()

    required, optional = get_versions_from_ci(content)

    assert required == {
        "numpy": "1.26.0",
        "python-dateutil": "2.8.2",
    }
    assert optional == {
        "beautifulsoup4": "4.12.3",
        "hypothesis": "6.116.0",
        "pyarrow": "16.0.0",
        "pytest": "8.3.4",
        "pytables": "3.10.1",
        "qtpy": "2.4.2",
        "zstandard": "0.23.0",
    }
