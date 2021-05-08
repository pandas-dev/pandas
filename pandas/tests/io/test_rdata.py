from io import BytesIO
import os
from urllib.error import HTTPError

import pytest

import pandas.util._test_decorators as td

from pandas import (
    DataFrame,
    Timestamp,
)
import pandas._testing as tm

from pandas.io.rdata._rdata import (
    LibrdataParserError,
    LibrdataWriterError,
)
from pandas.io.rdata.rdata_reader import read_rdata

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

ppm_df = DataFrame(
    {
        "date": {
            612: Timestamp("2020-12-16 23:42:25.920000256"),
            613: Timestamp("2021-01-16 11:17:31.199999744"),
            614: Timestamp("2021-02-15 21:00:00"),
            615: Timestamp("2021-03-18 06:42:28.800000256"),
            616: Timestamp("2021-04-17 17:17:31.199999744"),
        },
        "decimal_date": {
            612: 2020.9583,
            613: 2021.0417,
            614: 2021.125,
            615: 2021.2083,
            616: 2021.2917,
        },
        "monthly_average": {
            612: 414.25,
            613: 415.52,
            614: 416.75,
            615: 417.64,
            616: 419.05,
        },
        "deseasonalized": {
            612: 414.98,
            613: 415.26,
            614: 415.93,
            615: 416.18,
            616: 416.23,
        },
        "num_days": {612: 30, 613: 29, 614: 28, 615: 28, 616: 24},
        "std_dev_of_days": {612: 0.47, 613: 0.44, 614: 1.02, 615: 0.86, 616: 1.12},
        "unc_of_mon_mean": {612: 0.17, 613: 0.16, 614: 0.37, 615: 0.31, 616: 0.44},
    }
).rename_axis("rownames")


@pytest.fixture(params=["rda", "rds"])
def rtype(request):
    return request.param


@pytest.fixture(params=[None, "gzip", "bz2", "xz"])
def comp(request):
    return request.param


# RDATA READER


# PATH_OR_BUFFER


def test_read_rds_file(datapath):
    filename = datapath("io", "data", "rdata", "ghg_df.rds")
    r_dfs = read_rdata(filename)

    tm.assert_frame_equal(ghg_df, r_dfs["r_dataframe"].tail())


def test_read_rda_file(datapath):
    filename = datapath("io", "data", "rdata", "env_data_dfs.rda")
    r_dfs = read_rdata(filename)

    assert list(r_dfs.keys()) == ["ghg_df", "plants_df", "sea_ice_df"]

    tm.assert_frame_equal(ghg_df, r_dfs["ghg_df"].tail())
    tm.assert_frame_equal(plants_df, r_dfs["plants_df"].tail())
    tm.assert_frame_equal(sea_ice_df, r_dfs["sea_ice_df"].tail())


def test_read_rds_filelike(datapath):
    filename = datapath("io", "data", "rdata", "sea_ice_df.rds")

    with open(filename, "rb") as f:
        r_dfs = read_rdata(f, file_format="rds")

    tm.assert_frame_equal(sea_ice_df, r_dfs["r_dataframe"].tail())


def test_read_rda_filelike(datapath):
    filename = datapath("io", "data", "rdata", "env_data_dfs.rda")

    with open(filename, "rb") as f:
        r_dfs = read_rdata(f, file_format="rda")

    assert list(r_dfs.keys()) == ["ghg_df", "plants_df", "sea_ice_df"]

    tm.assert_frame_equal(ghg_df, r_dfs["ghg_df"].tail())
    tm.assert_frame_equal(plants_df, r_dfs["plants_df"].tail())
    tm.assert_frame_equal(sea_ice_df, r_dfs["sea_ice_df"].tail())


def test_bytesio_rds(datapath):
    filename = datapath("io", "data", "rdata", "sea_ice_df.rds")

    with open(filename, "rb") as f:
        with BytesIO(f.read()) as b_io:
            r_dfs = read_rdata(b_io, file_format="rds")

    tm.assert_frame_equal(sea_ice_df, r_dfs["r_dataframe"].tail())


def test_bytesio_rda(datapath):
    filename = datapath("io", "data", "rdata", "env_data_dfs.rda")

    with open(filename, "rb") as f:
        with BytesIO(f.read()) as b_io:
            r_dfs = read_rdata(b_io, file_format="rda")

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
    with pytest.raises(
        LibrdataParserError,
        match="Invalid file, or file has unsupported features",
    ):
        filename = datapath("io", "data", "rdata", "ppm_ts.rds")
        read_rdata(filename)


def test_read_rda_non_dfs(datapath):
    with pytest.raises(
        LibrdataParserError,
        match="Invalid file, or file has unsupported features",
    ):
        filename = datapath("io", "data", "rdata", "env_data_non_dfs.rda")
        read_rdata(filename)


def test_read_not_rda_file(datapath):
    with pytest.raises(
        LibrdataParserError, match="The file contains an unrecognized object"
    ):
        filename = datapath("io", "data", "rdata", "ppm_df.csv")
        read_rdata(filename, file_format="rda", compression=None)


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

    # librdata remove dots in colnames
    r_dfs["wine"].columns = r_dfs["wine"].columns.str.replace(" ", ".")

    tm.assert_frame_equal(s3_df, r_dfs["wine"].head())


# TYPE


def test_read_rds_df_output(datapath):
    filename = datapath("io", "data", "rdata", "sea_ice_df.rds")
    r_dfs = read_rdata(filename)

    assert isinstance(r_dfs, dict)
    assert list(r_dfs.keys()) == ["r_dataframe"]


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
    r_df = read_rdata(filename, rownames=True)["r_dataframe"]

    if isinstance(r_df, DataFrame):
        assert r_df.index.name == "rownames"


def test_read_rownames_false_rds(datapath):
    filename = datapath("io", "data", "rdata", "sea_ice_df.rds")
    r_df = read_rdata(filename, rownames=False)["r_dataframe"]

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
    with pytest.raises(SystemError, match=("returned a result with an error set")):
        read_rdata(filename)


# DATE / TIME


def test_utc_datetime_convert(datapath):
    filename = datapath("io", "data", "rdata", "ppm_df.rda")
    r_dfs = read_rdata(filename)

    assert str(r_dfs["ppm_df"]["date"].dtype) == "datetime64[ns]"

    tm.assert_frame_equal(ppm_df, r_dfs["ppm_df"].tail())


# RDATA WRITER

# PATH_OR_BUFFER


def test_write_read_file(rtype):
    with tm.ensure_clean("test.out") as path:
        ghg_df.to_rdata(path, file_format=rtype, index=False)
        r_dfs = read_rdata(path, file_format=rtype, rownames=False)

        expected = ghg_df.reset_index(drop=True)
        output = r_dfs["pandas_dataframe"] if rtype == "rda" else r_dfs["r_dataframe"]

        tm.assert_frame_equal(output, expected)


def test_write_read_pathlib(rtype):
    from pathlib import Path

    with tm.ensure_clean_dir() as tmp_dir:
        tmp_file = Path(tmp_dir).joinpath("test.out")
        sea_ice_df.to_rdata(tmp_file, file_format=rtype, index=False)
        r_dfs = read_rdata(tmp_file, file_format=rtype, rownames=False)

        expected = sea_ice_df.reset_index(drop=True)
        output = r_dfs["pandas_dataframe"] if rtype == "rda" else r_dfs["r_dataframe"]

        tm.assert_frame_equal(output, expected)


def test_write_read_filelike(rtype):
    with BytesIO() as b_io:
        sea_ice_df.to_rdata(b_io, file_format=rtype, compression=None, index=False)
        r_dfs = read_rdata(
            b_io.getvalue(),
            file_format=rtype,
            rownames=False,
            compression=None,
        )

        expected = sea_ice_df.reset_index(drop=True)
        output = r_dfs["pandas_dataframe"] if rtype == "rda" else r_dfs["r_dataframe"]

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


def test_write_index_true(rtype):
    with tm.ensure_clean("test.out") as path:
        plants_df.rename_axis(None).to_rdata(path, file_format=rtype, index=True)
        r_dfs = read_rdata(path, file_format=rtype)

    r_df = r_dfs if rtype == "rds" else r_dfs["pandas_dataframe"]

    if isinstance(r_df, DataFrame):
        assert "index" in r_df.columns


def test_write_index_false(rtype):
    with tm.ensure_clean("test.out") as path:
        plants_df.rename_axis(None).to_rdata(path, file_format=rtype, index=False)
        r_dfs = read_rdata(path, file_format=rtype)

    r_df = r_dfs if rtype == "rds" else r_dfs["pandas_dataframe"]

    if isinstance(r_df, DataFrame):
        assert "index" not in r_df.columns


# COMPRESS


def test_write_compress_all(rtype, comp):
    with tm.ensure_clean("test.out") as path:
        ghg_df.to_rdata(path, file_format=rtype, compression=comp, index=False)
        r_dfs = read_rdata(path, file_format=rtype, compression=comp, rownames=False)

        expected = ghg_df.reset_index(drop=True)
        output = r_dfs["pandas_dataframe"] if rtype == "rda" else r_dfs["r_dataframe"]

        tm.assert_frame_equal(output, expected)


def test_write_compress_zip(rtype):
    with tm.ensure_clean("test.out") as path:
        with pytest.raises(ValueError, match=("not a supported value for compression")):
            ghg_df.to_rdata(path, file_format=rtype, index=False, compression="zip")


# RDA_NAMES


def test_write_new_rda_name():
    with tm.ensure_clean("test.rda") as path:
        ghg_df.to_rdata(path, rda_name="py_df")
        r_dfs = read_rdata(path)

        assert "py_df" in list(r_dfs.keys())


# PROBLEM DATA


def test_write_nested_list(rtype, comp):
    plants_df["plants_dict"] = plants_df["plant_group"].apply(
        lambda x: plants_df["plant_group"].unique()
    )
    with tm.ensure_clean("test") as path:
        with pytest.raises(
            LibrdataWriterError,
            match=("DataFrame contains one more invalid types or data values"),
        ):
            plants_df.to_rdata(path, file_format=rtype, compression=comp)


# DATE / TIME


def test_write_read_utc_dateteime():
    with tm.ensure_clean("test.rda") as path:
        ppm_df.to_rdata(path, index=False)
        r_dfs = read_rdata(path, rownames=False)

        ppm_df["date"] = ppm_df["date"].dt.floor("S")

        tm.assert_frame_equal(ppm_df.reset_index(drop=True), r_dfs["pandas_dataframe"])
