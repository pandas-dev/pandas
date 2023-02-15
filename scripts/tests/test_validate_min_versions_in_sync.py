import pathlib
import sys

if sys.version_info >= (3, 11):
    import tomllib
else:
    import tomli as tomllib

from scripts.validate_min_versions_in_sync import update_yaml_file_version


def test_update_yaml_file_version():
    # minimum dependencies file
    TOML_PATH = pathlib.Path("scripts/tests/dummy_toml_file.toml")
    # yaml file that needs to be corrected
    YAML_UNMODIFIED_PATH = pathlib.Path(
        "scripts/tests/dummy_yaml_unmodified_file_1.yaml"
    )
    # yaml file that is the expected result
    YAML_EXPECTED_PATH = pathlib.Path("scripts/tests/dummy_yaml_expected_file_1.yaml")
    toml_dic = {}
    with open(TOML_PATH, "rb") as toml_f:
        toml_dic = tomllib.load(toml_f)
    with open(YAML_UNMODIFIED_PATH) as yaml_f:
        result_yaml_file = update_yaml_file_version(yaml_f, toml_dic)
    with open(YAML_EXPECTED_PATH) as yaml_f:
        dummy_yaml_expected_file_1 = yaml_f.read()
    # if assert passes, then yaml file version update is correct
    assert result_yaml_file == dummy_yaml_expected_file_1
