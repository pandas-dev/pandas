from csv import (
    DictWriter,
    reader as csv_reader,
)

import pytest

import pandas as pd
import pandas._testing as tm


@pytest.fixture(scope="module")
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


@pytest.fixture
def s3_bucket_public_with_directory_data(s3_bucket_public, directory_data_to_file):
    for file in directory_data_to_file.iterdir():
        s3_bucket_public.put_object(
            Key=file.name,
            Body=file.open("rb"),
            ContentType="text/csv",
        )
    return s3_bucket_public


def test_directory_data(directory_data_to_file):
    assert len(list(directory_data_to_file.iterdir())) == 3
    for file in directory_data_to_file.iterdir():
        with file.open(encoding="utf-8") as f:
            reader = csv_reader(f)
            header = next(reader)
            assert header == ["a", "b", "c"]


def test_read_directory_local(directory_data, directory_data_to_file):
    _, data_list = directory_data
    df_list = list(pd.read_csv(directory_data_to_file))
    assert len(df_list) == len(data_list)
    df_concat = pd.concat(df_list, ignore_index=True)
    df_concat = df_concat.sort_values(by=list(df_concat.columns)).reset_index(drop=True)
    expected = pd.DataFrame([value for data in data_list for value in data.values()])
    tm.assert_frame_equal(df_concat, expected)


def test_read_directory_s3(s3_bucket_public_with_directory_data, s3so, directory_data):
    _, data_list = directory_data
    df_list = list(
        pd.read_csv(
            f"s3://{s3_bucket_public_with_directory_data.name}/",
            storage_options=s3so,
        )
    )
    assert len(df_list) == len(data_list)
    df_concat = pd.concat(df_list, ignore_index=True)
    df_concat = df_concat.sort_values(by=list(df_concat.columns)).reset_index(drop=True)
    expected = pd.DataFrame([value for data in data_list for value in data.values()])
    tm.assert_frame_equal(df_concat, expected)
