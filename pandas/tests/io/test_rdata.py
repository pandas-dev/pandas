from io import BytesIO
import os
from urllib.error import HTTPError

import pytest

import pandas.util._test_decorators as td

from pandas import DataFrame
import pandas._testing as tm

from pandas.io.rdata import read_rdata

pyreadr = pytest.importorskip("pyreadr")


ghg_df = DataFrame(
    {
        "gas": {
            141: "Carbon dioxide",
            142: "Methane",
            143: "Nitrous oxide",
            144: "Fluorinated gases",
            145: "Total",
        },
        "year": {141: 2018, 142: 2018, 143: 2018, 144: 2018, 145: 2018},
        "emissions": {
            141: 5424.881502132882,
            142: 634.4571270782675,
            143: 434.52855537666636,
            144: 182.78243246177678,
            145: 6676.649617049592,
        },
    }
).rename_axis("rownames")

plants_df = DataFrame(
    {
        "plant_group": {
            16: "Pteridophytes",
            17: "Pteridophytes",
            18: "Pteridophytes",
            19: "Pteridophytes",
            20: "Pteridophytes",
        },
        "status": {
            16: "Data Deficient",
            17: "Extinct",
            18: "Not Threatened",
            19: "Possibly Threatened",
            20: "Threatened",
        },
        "count": {16: 398, 17: 65, 18: 1294, 19: 408, 20: 1275},
    }
).rename_axis("rownames")

sea_ice_df = DataFrame(
    {
        "year": {1012: 2016, 1013: 2017, 1014: 2018, 1015: 2019, 1016: 2020},
        "mo": {1012: 12, 1013: 12, 1014: 12, 1015: 12, 1016: 12},
        "data.type": {
            1012: "Goddard",
            1013: "Goddard",
            1014: "Goddard",
            1015: "Goddard",
            1016: "NRTSI-G",
        },
        "region": {1012: "S", 1013: "S", 1014: "S", 1015: "S", 1016: "S"},
        "extent": {1012: 8.28, 1013: 9.48, 1014: 9.19, 1015: 9.41, 1016: 10.44},
        "area": {1012: 5.51, 1013: 6.23, 1014: 5.59, 1015: 6.59, 1016: 6.5},
    }
).rename_axis("rownames")


@pytest.fixture(params=["rda", "rds"])
def rtype(request):
    return request.param


@pytest.fixture(params=[None, "gzip", "bz2", "xz"])
def comp(request):
    return request.param


def adj_int(df):
    """
    Convert int32 columns to int64.

    Since pyreadr engine reads ints int int32 and writes ints
    to floats  this method converts such columns for testing.
    """
    int_cols = df.select_dtypes("int32").columns
    df[int_cols] = df[int_cols].astype("int64")

    if "index" in df.columns:
        df["index"] = df["index"].astype("int64")

    if "year" in df.columns:
        df["year"] = df["year"].astype("int64")
    if "mo" in df.columns:
        df["mo"] = df["mo"].astype("int64")

    return df


# RDA READER

# PATH_OR_BUFFER


def test_read_rds_file(datapath):
    filename = datapath("io", "data", "rdata", "ghg_df.rds")
    r_df = read_rdata(filename)
    output = adj_int(r_df).tail()

    tm.assert_frame_equal(ghg_df, output)


def test_read_rda_file(datapath):
    filename = datapath("io", "data", "rdata", "env_data_dfs.rda")
    r_dfs = read_rdata(filename)

    r_dfs = {str(k): adj_int(v) for k, v in r_dfs.items()}

    assert list(r_dfs.keys()) == ["ghg_df", "plants_df", "sea_ice_df"]

    tm.assert_frame_equal(ghg_df, r_dfs["ghg_df"].tail())
    tm.assert_frame_equal(plants_df, r_dfs["plants_df"].tail())
    tm.assert_frame_equal(sea_ice_df, r_dfs["sea_ice_df"].tail())


def test_bytes_read_rds(datapath):
    filename = datapath("io", "data", "rdata", "sea_ice_df.rds")

    with open(filename, "rb") as f:
        r_df = read_rdata(f, file_format="rds")

    output = adj_int(r_df).tail()

    tm.assert_frame_equal(sea_ice_df, output)


def test_bytes_read_rda(datapath):
    filename = datapath("io", "data", "rdata", "env_data_dfs.rda")

    with open(filename, "rb") as f:
        r_dfs = read_rdata(f, file_format="rda")

    r_dfs = {str(k): adj_int(v) for k, v in r_dfs.items()}

    assert list(r_dfs.keys()) == ["ghg_df", "plants_df", "sea_ice_df"]

    tm.assert_frame_equal(ghg_df, r_dfs["ghg_df"].tail())
    tm.assert_frame_equal(plants_df, r_dfs["plants_df"].tail())
    tm.assert_frame_equal(sea_ice_df, r_dfs["sea_ice_df"].tail())


def test_bytesio_rds(datapath):
    filename = datapath("io", "data", "rdata", "sea_ice_df.rds")

    with open(filename, "rb") as f:
        with BytesIO(f.read()) as b_io:
            r_df = read_rdata(b_io, file_format="rds")

    output = adj_int(r_df).tail()

    tm.assert_frame_equal(sea_ice_df, output)


def test_bytesio_rda(datapath):
    filename = datapath("io", "data", "rdata", "env_data_dfs.rda")

    with open(filename, "rb") as f:
        with BytesIO(f.read()) as b_io:
            r_dfs = read_rdata(b_io, file_format="rda")

    r_dfs = {str(k): adj_int(v) for k, v in r_dfs.items()}

    assert list(r_dfs.keys()) == ["ghg_df", "plants_df", "sea_ice_df"]

    tm.assert_frame_equal(ghg_df, r_dfs["ghg_df"].tail())
    tm.assert_frame_equal(plants_df, r_dfs["plants_df"].tail())
    tm.assert_frame_equal(sea_ice_df, r_dfs["sea_ice_df"].tail())


# FILE FORMAT


def test_read_wrong_format(datapath):
    with pytest.raises(ValueError, match="not a valid value for file_format"):
        filename = datapath("io", "data", "rdata", "plants_df.rds")
        read_rdata(filename, file_format="r")


def test_read_wrong_file():
    with pytest.raises(FileNotFoundError, match="file cannot be found"):
        filename = os.path.join("data", "rdata", "plants_df.rda")
        read_rdata(filename)


def test_read_rds_non_df(datapath):
    from pyreadr import custom_errors

    with pytest.raises(
        custom_errors.LibrdataError,
        match="Invalid file, or file has unsupported features",
    ):
        filename = datapath("io", "data", "rdata", "ppm_ts.rds")
        read_rdata(filename)


def test_read_rda_non_dfs(datapath):
    from pyreadr import custom_errors

    with pytest.raises(
        custom_errors.LibrdataError,
        match="Invalid file, or file has unsupported features",
    ):
        filename = datapath("io", "data", "rdata", "env_data_non_dfs.rda")
        read_rdata(filename)


def test_read_not_rda_file(datapath):
    from pyreadr import custom_errors

    with pytest.raises(
        custom_errors.LibrdataError, match="The file contains an unrecognized object"
    ):
        filename = datapath("io", "data", "rdata", "ppm_df.csv")
        read_rdata(filename, file_format="rda")


def test_bytes_read_infer_rds(datapath):
    filename = datapath("io", "data", "rdata", "sea_ice_df.rds")

    with pytest.raises(ValueError, match="Unable to infer file format from file name"):
        with open(filename, "rb") as f:
            read_rdata(f)


def test_bytes_read_infer_rda(datapath):
    filename = datapath("io", "data", "rdata", "env_data_dfs.rda")

    with pytest.raises(ValueError, match="Unable to infer file format from file name"):
        with open(filename, "rb") as f:
            read_rdata(f)


# URL


@tm.network
def test_read_rda_url():
    url_df = DataFrame(
        {
            "carrier": {1: "9E", 2: "AA", 3: "AS", 4: "B6", 5: "DL"},
            "name": {
                1: "Endeavor Air Inc.",
                2: "American Airlines Inc.",
                3: "Alaska Airlines Inc.",
                4: "JetBlue Airways",
                5: "Delta Air Lines Inc.",
            },
        }
    ).rename_axis("rownames")

    url = (
        "https://github.com/hadley/nycflights13/blob/master/data/airlines.rda?raw=true"
    )
    r_dfs = read_rdata(url, file_format="rda")

    tm.assert_frame_equal(url_df, r_dfs["airlines"].head())


@tm.network
def test_read_unable_infer_format():
    with pytest.raises(ValueError, match="Unable to infer file format from file name"):
        url = (
            "https://github.com/hadley/nycflights13/"
            "blob/master/data/airlines.rda?raw=true"
        )
        read_rdata(url)


@tm.network
def test_read_wrong_url():
    with pytest.raises(HTTPError, match="HTTP Error 404: Not Found"):
        url = "https://example.com/data.rdata"
        read_rdata(url)


# S3


@tm.network
@td.skip_if_no("s3fs")
def test_read_rda_s3():
    s3 = "s3://assets.datacamp.com/production/course_1478/datasets/wine.RData"
    s3_df = DataFrame(
        {
            "Alcohol": {1: 13.2, 2: 13.16, 3: 14.37, 4: 13.24, 5: 14.2},
            "Malic.acid": {1: 1.78, 2: 2.36, 3: 1.95, 4: 2.59, 5: 1.76},
            "Ash": {1: 2.14, 2: 2.67, 3: 2.5, 4: 2.87, 5: 2.45},
            "Alcalinity.of.ash": {1: 11.2, 2: 18.6, 3: 16.8, 4: 21.0, 5: 15.2},
            "Magnesium": {1: 100, 2: 101, 3: 113, 4: 118, 5: 112},
            "Total.phenols": {1: 2.65, 2: 2.8, 3: 3.85, 4: 2.8, 5: 3.27},
            "Flavanoids": {1: 2.76, 2: 3.24, 3: 3.49, 4: 2.69, 5: 3.39},
            "Nonflavanoid.phenols": {1: 0.26, 2: 0.3, 3: 0.24, 4: 0.39, 5: 0.34},
            "Proanthocyanins": {1: 1.28, 2: 2.81, 3: 2.18, 4: 1.82, 5: 1.97},
            "Color.intensity": {1: 4.38, 2: 5.68, 3: 7.8, 4: 4.32, 5: 6.75},
            "Hue": {1: 3.4, 2: 3.17, 3: 3.45, 4: 2.93, 5: 2.85},
            "Proline": {1: 1050, 2: 1185, 3: 1480, 4: 735, 5: 1450},
        }
    ).rename_axis("rownames")
    r_dfs = read_rdata(s3)
    r_dfs["wine"] = adj_int(r_dfs["wine"])

    # pyreadr remove dots in colnames
    r_dfs["wine"].columns = r_dfs["wine"].columns.str.replace(" ", ".")

    tm.assert_frame_equal(s3_df, r_dfs["wine"].head())


# ENGINE


def test_read_rds_df_output(datapath):
    filename = datapath("io", "data", "rdata", "sea_ice_df.rds")
    r_df = read_rdata(filename)

    assert isinstance(r_df, DataFrame)


def test_read_rda_dict_output(datapath):
    filename = datapath("io", "data", "rdata", "env_data_dfs.rda")
    r_dfs = read_rdata(filename)

    assert isinstance(r_dfs, dict)
    assert list(r_dfs.keys()) == ["ghg_df", "plants_df", "sea_ice_df"]


# SELECT_FRAMES


def test_read_select_frames_rda_dfs(datapath):
    filename = datapath("io", "data", "rdata", "env_data_dfs.rda")
    r_dfs = read_rdata(filename, select_frames=["ghg_df", "sea_ice_df"])

    assert "plants_df" not in list(r_dfs.keys())
    assert "ghg_df" in list(r_dfs.keys())
    assert "sea_ice_df" in list(r_dfs.keys())


def test_read_wrong_select_frames(datapath):
    with pytest.raises(TypeError, match="not a valid type for select_frames"):
        filename = datapath("io", "data", "rdata", "env_data_dfs.rda")
        read_rdata(filename, select_frames="plants_df")


# ROWNAMES


def test_read_rownames_true_rds(datapath):
    filename = datapath("io", "data", "rdata", "sea_ice_df.rds")
    r_df = read_rdata(filename, rownames=True)

    if isinstance(r_df, DataFrame):
        assert r_df.index.name == "rownames"


def test_read_rownames_false_rds(datapath):
    filename = datapath("io", "data", "rdata", "sea_ice_df.rds")
    r_df = read_rdata(filename, rownames=False)

    if isinstance(r_df, DataFrame):
        assert r_df.index.name != "rownames"


def test_read_rownames_true_rda(datapath):
    filename = datapath("io", "data", "rdata", "env_data_dfs.rda")
    r_dfs = read_rdata(filename, rownames=True)

    assert r_dfs["ghg_df"].index.name == "rownames"
    assert r_dfs["plants_df"].index.name == "rownames"
    assert r_dfs["sea_ice_df"].index.name == "rownames"


def test_read_rownames_false_rda(datapath):
    filename = datapath("io", "data", "rdata", "env_data_dfs.rda")
    r_dfs = read_rdata(filename, rownames=False)

    assert r_dfs["ghg_df"].index.name != "rownames"
    assert r_dfs["plants_df"].index.name != "rownames"
    assert r_dfs["sea_ice_df"].index.name != "rownames"


# ENCODING


def test_non_utf8_data(datapath, rtype):
    filename = datapath("io", "data", "rdata", f"climate_non_utf8_df.{rtype}")
    with pytest.raises(UnicodeDecodeError, match=("'utf-8' codec can't decode byte")):
        read_rdata(filename)


# RDA WRITER

# PATH_OR_BUFFER


def test_write_read_file(rtype):
    with tm.ensure_clean("test.out") as path:
        ghg_df.to_rdata(path, file_format=rtype, index=False)
        r_dfs = read_rdata(path, file_format=rtype, rownames=False)

        expected = ghg_df.reset_index(drop=True)
        output = (
            adj_int(r_dfs["pandas_dataframe"]) if rtype == "rda" else adj_int(r_dfs)
        )

        tm.assert_frame_equal(output, expected)


def test_write_read_pathlib(rtype):
    from pathlib import Path

    with tm.ensure_clean_dir() as tmp_dir:
        tmp_file = Path(tmp_dir).joinpath("test.out")
        sea_ice_df.to_rdata(tmp_file, file_format=rtype, index=False)
        r_dfs = read_rdata(tmp_file, file_format=rtype, rownames=False)

        expected = sea_ice_df.reset_index(drop=True)
        output = (
            adj_int(r_dfs["pandas_dataframe"]) if rtype == "rda" else adj_int(r_dfs)
        )

        tm.assert_frame_equal(output, expected)


def test_write_read_filelike(rtype):
    with BytesIO() as b_io:
        sea_ice_df.to_rdata(b_io, file_format=rtype, index=False)
        r_dfs = read_rdata(
            b_io.getvalue(),
            file_format=rtype,
            rownames=False,
        )

        expected = sea_ice_df.reset_index(drop=True)
        output = (
            adj_int(r_dfs["pandas_dataframe"]) if rtype == "rda" else adj_int(r_dfs)
        )

        tm.assert_frame_equal(output, expected)


# FILE FORMAT


def test_write_wrong_format():
    with tm.ensure_clean("test.rda") as path:
        with pytest.raises(ValueError, match=("not a valid value for file_format")):
            ghg_df.to_rdata(path, file_format="csv")


def test_write_unable_to_infer():
    with tm.ensure_clean("test") as path:
        with pytest.raises(
            ValueError, match=("Unable to infer file format from file name")
        ):
            ghg_df.to_rdata(path)


# INDEX


def test_index_true(rtype):
    with tm.ensure_clean("test.out") as path:
        plants_df.rename_axis(None).to_rdata(path, file_format=rtype, index=True)
        r_dfs = read_rdata(path, file_format=rtype)

    r_df = r_dfs if rtype == "rds" else r_dfs["pandas_dataframe"]

    if isinstance(r_df, DataFrame):
        assert "index" in r_df.columns


def test_index_false(rtype):
    with tm.ensure_clean("test.out") as path:
        plants_df.rename_axis(None).to_rdata(path, file_format=rtype, index=False)
        r_dfs = read_rdata(path, file_format=rtype)

    r_df = r_dfs if rtype == "rds" else r_dfs["pandas_dataframe"]

    if isinstance(r_df, DataFrame):
        assert "index" not in r_df.columns


# COMPRESS


def test_compress_ok_comp(rtype, comp):
    with tm.ensure_clean("test.out") as path:
        ghg_df.to_rdata(path, file_format=rtype, compression=comp, index=False)
        r_dfs = read_rdata(path, file_format=rtype, rownames=False)

        expected = ghg_df.reset_index(drop=True)
        output = (
            adj_int(r_dfs["pandas_dataframe"]) if rtype == "rda" else adj_int(r_dfs)
        )

        tm.assert_frame_equal(output, expected)


def test_compress_zip(rtype):
    with tm.ensure_clean("test.out") as path:
        with pytest.raises(ValueError, match=("not a supported value for compression")):
            ghg_df.to_rdata(path, file_format=rtype, index=False, compression="zip")


# RDA_NAMES


def test_new_rda_name():
    with tm.ensure_clean("test.rda") as path:
        ghg_df.to_rdata(path, rda_name="py_df")
        r_dfs = read_rdata(path)

        assert "py_df" in list(r_dfs.keys())
