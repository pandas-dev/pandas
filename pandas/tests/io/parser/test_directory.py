from csv import (
    DictWriter,
    reader as csv_reader,
)

import pytest


@pytest.fixture
def directory_data():
    return ["a", "b", "c"], [
        {"first": {"a": 1, "b": 2, "c": 3}},
        {"second": {"a": 4, "b": 5, "c": 6}},
        {"third": {"a": 7, "b": 8, "c": 9}},
    ]


@pytest.fixture
def directory_data_to_file(tmp_path, directory_data):
    field_names, data_list = directory_data
    for data in data_list:
        file_name = next(iter(data.keys()))
        path = tmp_path / f"{file_name}.csv"
        with path.open("w", newline="", encoding="utf-8") as file:
            writer = DictWriter(file, fieldnames=field_names)
            writer.writeheader()
            writer.writerow(data[file_name])
    return tmp_path


def test_directory_data(directory_data_to_file):
    assert len(list(directory_data_to_file.iterdir())) == 3
    for file in directory_data_to_file.iterdir():
        with file.open(encoding="utf-8") as f:
            reader = csv_reader(f)
            header = next(reader)
            assert header == ["a", "b", "c"]
