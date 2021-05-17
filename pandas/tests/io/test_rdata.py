import gzip
from io import BytesIO
import os
import pickle
import shutil
from urllib.error import HTTPError

import numpy as np
import pytest

from pandas._libs.tslibs.np_datetime import OutOfBoundsDatetime
from pandas.compat import (
    IS64,
    PY38,
)
import pandas.util._test_decorators as td

from pandas import (
    Categorical,
    DataFrame,
    Timestamp,
    array,
    interval_range,
    period_range,
    to_datetime,
)
import pandas._testing as tm
from pandas.arrays import SparseArray

from pandas.io.rdata._rdata import (
    LibrdataReader,
    LibrdataReaderError,
    LibrdataWriter,
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
            754: Timestamp("2020-12-16 23:42:25.920000256"),
            755: Timestamp("2021-01-16 11:17:31.199999744"),
            756: Timestamp("2021-02-15 21:00:00"),
            757: Timestamp("2021-03-18 06:42:28.800000256"),
            758: Timestamp("2021-04-17 17:17:31.199999744"),
        },
        "decimal_date": {
            754: 2020.9583,
            755: 2021.0417,
            756: 2021.125,
            757: 2021.2083,
            758: 2021.2917,
        },
        "monthly_average": {
            754: 414.25,
            755: 415.52,
            756: 416.75,
            757: 417.64,
            758: 419.05,
        },
        "deseasonalized": {
            754: 414.98,
            755: 415.26,
            756: 415.93,
            757: 416.18,
            758: 416.23,
        },
        "num_days": {754: 30, 755: 29, 756: 28, 757: 28, 758: 24},
        "std_dev_of_days": {754: 0.47, 755: 0.44, 756: 1.02, 757: 0.86, 758: 1.12},
        "unc_of_mon_mean": {754: 0.17, 755: 0.16, 756: 0.37, 757: 0.31, 758: 0.44},
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
        LibrdataReaderError,
        match="Invalid file, or file has unsupported features",
    ):
        filename = datapath("io", "data", "rdata", "ppm_ts.rds")
        read_rdata(filename)


def test_read_rda_non_dfs(datapath):
    with pytest.raises(
        LibrdataReaderError,
        match="Invalid file, or file has unsupported features",
    ):
        filename = datapath("io", "data", "rdata", "env_data_non_dfs.rda")
        read_rdata(filename)


def test_read_not_rda_file(datapath):
    with pytest.raises(
        LibrdataReaderError, match="The file contains an unrecognized object"
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


@pytest.mark.slow
@tm.network
@td.skip_if_no("s3fs")
def test_read_rda_s3():
    # Public Data of CRAN Packages on GitHub
    rda_s3 = "s3://public-r-data/ghcran.Rdata"
    r_df = read_rdata(rda_s3, compression=None, rownames=False)

    # below needed to pass codespell on keyword
    r_df["ghcran"].columns.values[107] = "Repository"

    # test structure and not static data since data changes daily
    expected_cols = [
        "Package",
        "Type",
        "Title",
        "Version",
        "Date",
        "Author",
        "Maintainer",
        "Description",
        "License",
        "Depends",
        "Suggests",
        "NeedsCompilation",
        "Packaged",
        "Repository",
        "Date/Publication",
        "Contact",
        "Imports",
        "VignetteBuilder",
        "Encoding",
        "SystemRequirements",
        "RoxygenNote",
        "LazyLoad",
        "URL",
        "Authors@R",
        "Classification/ACM",
        "Classification/JEL",
        "LinkingTo",
        "BugReports",
        "LazyData",
        "Keywords",
        "Repository/R-Forge/Project",
        "Repository/R-Forge/Revision",
        "Repository/R-Forge/DateTimeStamp",
        "biocViews",
        "Collate",
        "Copyright",
        "ByteCompile",
        "ZipData",
        "BuildVignettes",
        "Additional_repositories",
        "Acknowledgements",
        "MailingList",
        "Enhances",
        "Classification/MSC",
        "OS_type",
        "BuildManual",
        "BuildResaveData",
        "References",
        "Note",
        "X-CRAN-Original-Maintainer",
        "RcppModules",
        "Data",
        "BioViews",
        "lazy-loading",
        "URLNote",
        "Reference",
        "KeepSource",
        "LazyDataCompression",
        "Language",
        "Requires",
        "Dependencies",
        "X-CRAN-Comment",
        "Citation",
        "Biarch",
        "Published",
        "RequiredLauncherGeneration",
        "SuggestsNote",
        "Priority",
        "Acknowledgments",
        "Revision",
        "License_is_FOSS",
        "License_restricts_use",
        "Archs",
        "LazyDataNote",
        "Affiliations",
        "LicenseDetails",
        "SCM",
        "Classification/ACM-2012",
        "X-CRAN-Original-Package",
        "Dialect",
        "Limitations",
        "Check",
        "Recommends",
        "LastChangedDate",
        "LastChangedRevision",
        "SVNRevision",
        "X-CRAN-Original-OS_type",
        "RcmdrModels",
        "Log-Exceptions",
        "Models",
        "DateNote",
        "SystemRequirementsNote",
        "Url",
        "Reverse depends",
        "Lazyload",
        "DependsNote",
        "VersionSplus",
        "MaintainerSplus",
        "VersionNote",
        "Disclaimer",
        "LicenseNote",
        "Namespace",
        "Address",
        "Keyword",
        "Contributors",
        "NOTE",
        "Acknowledgement",
        "Repository",
        "Lazydata",
        "RdMacros",
        "HowToCite",
        "Publication",
        "Reference Manual",
        "Special Acknowledgement",
        "SysDataCompression",
        "DisplayMode",
        "Nickname",
        "BuildKeepEmpty",
        "Twitter",
        "Remotes",
        "SystemRequirement",
        "Github",
    ]

    assert isinstance(r_df, dict)
    assert isinstance(r_df["ghcran"], DataFrame)
    assert r_df["ghcran"].columns.tolist() == expected_cols


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


def test_read_outbound_dates(datapath, rtype):
    filename = datapath("io", "data", "rdata", f"planetary_boundaries_df.{rtype}")
    with pytest.raises(
        OutOfBoundsDatetime, match=("cannot convert input with unit 's'")
    ):
        read_rdata(filename)


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


# COMPRESSION


def test_write_all_compression(rtype, comp):
    with tm.ensure_clean("test.out") as path:
        ghg_df.to_rdata(path, file_format=rtype, compression=comp, index=False)
        r_dfs = read_rdata(path, file_format=rtype, compression=comp, rownames=False)

        expected = ghg_df.reset_index(drop=True)
        output = r_dfs["pandas_dataframe"] if rtype == "rda" else r_dfs["r_dataframe"]

        tm.assert_frame_equal(output, expected)


def test_write_zip_compression(rtype):
    with tm.ensure_clean("test.out") as path:
        with pytest.raises(ValueError, match=("not a supported value for compression")):
            ghg_df.to_rdata(path, file_format=rtype, compression="zip")


@pytest.mark.skipif(
    not PY38,
    reason=("gzip.BadGzipFile exception added in 3.8"),
)
def test_write_read_mismatched_compression(rtype):
    with tm.ensure_clean("test.out") as path:
        with pytest.raises(gzip.BadGzipFile, match=("Not a gzipped file")):
            ghg_df.to_rdata(path, file_format=rtype, compression=None)
            read_rdata(path, file_format=rtype)


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


# DTYPES


@pytest.mark.skipif(
    not IS64,
    reason=("large dtypes not supported in 32-bit"),
)
def test_write_read_dtypes(rtype, comp):
    rda_name = "pandas_dataframe" if rtype == "rda" else "r_dataframe"

    dts = [
        Timestamp.min.ceil("S"),
        Timestamp(-(10 ** 18)),
        Timestamp(0),
        Timestamp(10 ** 18),
        Timestamp.now().floor("S"),
        Timestamp.max.floor("S"),
    ]

    arr = np.random.randn(6)
    arr[2:-2] = np.nan

    dtypes_df = DataFrame(
        {
            "categ": Categorical(
                ["ocean", "climate", "biosphere", "land", "freshwater", "atmosphere"]
            ),
            "interval": interval_range(start=10, periods=6, freq=10 * 2),
            "bool": [False, True, True, True, False, False],
            "int": [2 ** 31 - 1, 1, -(2 ** 31) + 1, -1, 0, 10 ** 9],
            "float": [0, np.pi, float("nan"), np.e, np.euler_gamma, 0],
            "string": array(
                ["acidification", "change", "loss", "use", "depletion", "aersols"],
                dtype="string",
            ),
            "sparse": SparseArray(arr),
            "period": period_range(
                start="2021-01-01 00:00:00", end="2021-06-01 00:00:00", freq="M"
            ),
            "datetime": to_datetime(dts),
            "datetime_tz": to_datetime(dts).tz_localize("utc"),
            "timedelta": [(dt - Timestamp(0)) for dt in dts],
        }
    )

    with tm.ensure_clean("test") as path:
        dtypes_df.to_rdata(path, file_format=rtype, index=False, compression=comp)
        r_df = read_rdata(path, file_format=rtype, rownames=False, compression=comp)[
            rda_name
        ]

    # convert non-primitive and non-datetimes to objects not supported in R
    excl_types = ["bool", "number", "object", "datetime", "datetimetz", "timedelta"]
    for col in dtypes_df.select_dtypes(exclude=excl_types).columns:
        dtypes_df[col] = dtypes_df[col].astype(str)

    # convert special types
    dtypes_df["sparse"] = np.array(dtypes_df["sparse"].values, dtype="float64")
    dtypes_df["datetime_tz"] = dtypes_df["datetime_tz"].dt.tz_localize(None)
    dtypes_df["timedelta"] = dtypes_df["timedelta"].dt.total_seconds()

    tm.assert_frame_equal(dtypes_df, r_df)


# CYTHON CLASSES


def test_reader_unpickled(datapath, rtype):
    if rtype == "rda":
        filename = datapath("io", "data", "rdata", "env_data_dfs.rda")
        rda_name = "sea_ice_df"
    elif rtype == "rds":
        filename = datapath("io", "data", "rdata", "plants_df.rds")
        rda_name = "r_dataframe"

    lbr1 = LibrdataReader()

    with tm.ensure_clean("test.pkl") as pklpath:
        with open(pklpath, "wb") as f_w:
            pickle.dump(lbr1, f_w)

        with open(pklpath, "rb") as f_r:
            lbr2 = pickle.load(f_r)

    with tm.ensure_clean("test") as r_temp:
        # need to decompress to temp file
        with gzip.open(filename, "rb") as f_r:
            with open(r_temp, "wb") as f_w:
                shutil.copyfileobj(f_r, f_w)

        df_output = read_rdata(
            r_temp, file_format=rtype, compression=None, rownames=False
        )[rda_name].to_dict()

        cy_output = lbr2.read_rdata(r_temp)

    lbr_output = {
        vcol: vdata
        for (kdata, vdata), (kcol, vcol) in zip(
            cy_output[rda_name]["data"].items(), cy_output[rda_name]["colnames"].items()
        )
    }

    assert lbr_output == df_output


def test_writer_unpickled(datapath, rtype):
    rda_name = "test_frame" if rtype == "rda" else "r_dataframe"

    lbw1 = LibrdataWriter()

    with tm.ensure_clean("test.pkl") as pklpath:
        with open(pklpath, "wb") as f_w:
            pickle.dump(lbw1, f_w)

        with open(pklpath, "rb") as f_r:
            lbw2 = pickle.load(f_r)

    rdict = {"dtypes": {k: str(v) for k, v in ghg_df.dtypes.to_dict().items()}}
    for k, v in rdict["dtypes"].items():
        if any(x in v for x in ("bool", "Boolean")):
            rdict["dtypes"][k] = "bool"

        elif any(x in v for x in ("int", "uint", "Int", "UInt")):
            rdict["dtypes"][k] = "int"

        elif any(x in v for x in ("float", "Float")):
            rdict["dtypes"][k] = "float"

        elif any(x in v for x in ("datetime", "Datetime")):
            rdict["dtypes"][k] = "datetime"

        elif any(x in v for x in ("object", "string", "String")):
            rdict["dtypes"][k] = "object"

    rdict["data"] = ghg_df.reset_index(drop=True).to_dict()

    expected = ghg_df.reset_index(drop=True)

    with tm.ensure_clean("test") as r_temp:
        lbw2.write_rdata(
            rfile=r_temp,
            rdict=rdict,
            rformat=rtype,
            tbl_name="test_frame",
        )

        output = read_rdata(
            r_temp,
            file_format=rtype,
            rownames=False,
            compression=None,
        )[rda_name]

    tm.assert_frame_equal(output, expected)
