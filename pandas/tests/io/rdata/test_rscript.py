from io import BytesIO
import os
import subprocess
from urllib.error import HTTPError

import pytest

from pandas.compat._optional import import_optional_dependency
import pandas.util._test_decorators as td

from pandas import DataFrame
import pandas._testing as tm

from pandas.io.rdata import (
    R_ARROW,
    R_RSQLITE,
    RSCRIPT_EXISTS,
    RScriptError,
    read_rdata,
)

pytestmark = pytest.mark.skipif(not RSCRIPT_EXISTS, reason="R is not installed.")

PYARROW = import_optional_dependency("pyarrow")

ghg_df = DataFrame(
    {
        "gas": {
            "141": "Carbon dioxide",
            "142": "Methane",
            "143": "Nitrous oxide",
            "144": "Fluorinated gases",
            "145": "Total",
        },
        "year": {"141": 2018, "142": 2018, "143": 2018, "144": 2018, "145": 2018},
        "emissions": {
            "141": 5424.88150213288,
            "142": 634.457127078267,
            "143": 434.528555376666,
            "144": 182.782432461777,
            "145": 6676.64961704959,
        },
    }
).rename_axis("rownames")

plants_df = DataFrame(
    {
        "plant_group": {
            "16": "Pteridophytes",
            "17": "Pteridophytes",
            "18": "Pteridophytes",
            "19": "Pteridophytes",
            "20": "Pteridophytes",
        },
        "status": {
            "16": "Data Deficient",
            "17": "Extinct",
            "18": "Not Threatened",
            "19": "Possibly Threatened",
            "20": "Threatened",
        },
        "count": {"16": 398, "17": 65, "18": 1294, "19": 408, "20": 1275},
    }
).rename_axis("rownames")

sea_ice_df = DataFrame(
    {
        "year": {"1012": 2016, "1013": 2017, "1014": 2018, "1015": 2019, "1016": 2020},
        "mo": {"1012": 12, "1013": 12, "1014": 12, "1015": 12, "1016": 12},
        "data.type": {
            "1012": "Goddard",
            "1013": "Goddard",
            "1014": "Goddard",
            "1015": "Goddard",
            "1016": "NRTSI-G",
        },
        "region": {"1012": "S", "1013": "S", "1014": "S", "1015": "S", "1016": "S"},
        "extent": {
            "1012": 8.28,
            "1013": 9.48,
            "1014": 9.19,
            "1015": 9.41,
            "1016": 10.44,
        },
        "area": {"1012": 5.51, "1013": 6.23, "1014": 5.59, "1015": 6.59, "1016": 6.5},
    }
).rename_axis("rownames")


def run_rscript(cmds) -> str:
    """
    Run R script at command line.

    This method will read write_rdata output and check
    console output.
    """

    r_batch = """
              args <- commandArgs(trailingOnly=TRUE)

              switch(args[2],
                  "rda" = load(args[1]),
                  "rds" = {
                     pandas_dataframe <- readRDS(args[1])
                  }
              )

              rm(args)
              mget(ls())
              """
    with open(cmds[1], "w") as f:
        f.write(r_batch)

    p = subprocess.Popen(
        cmds, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    output, error = p.communicate()
    if len(error) != 0:
        raise ValueError(error.decode("UTF-8"))

    return output.decode("UTF-8")


def adj_int(df):
    """
    Convert int32 columns to int64.

    Since parquet and feather modes parses ints int int32,
    this method converts for testing.
    """
    for col in df.select_dtypes("int32").columns:
        df[col] = df[col].astype("int64")

    return df


def handle_index_rownames(df):
    df = df.drop(["rownames"], axis=1).set_index("index").rename_axis(None)

    return df


@pytest.fixture(params=["rda", "rds"])
def rtype(request):
    return request.param


@pytest.fixture(
    params=[
        "csv",
        pytest.param(
            "parquet",
            marks=pytest.mark.skipif(
                not R_ARROW or not PYARROW,
                reason="R arrow or pyarrow not installed",
            ),
        ),
        pytest.param(
            "feather",
            marks=pytest.mark.skipif(
                not R_ARROW or not PYARROW,
                reason="R arrow or pyarrow not installed",
            ),
        ),
        pytest.param(
            "sqlite",
            marks=pytest.mark.skipif(not R_RSQLITE, reason="R RSQLite not installed"),
        ),
    ]
)
def mode(request):
    return request.param


@pytest.fixture(params=[True, False, None])
def ascii(request):
    return request.param


@pytest.fixture(params=[False, "gzip", "bzip2", "xz"])
def comp(request):
    return request.param


# RDA READER

# PATH_OR_BUFFER


def test_read_rds_file(datapath):
    filename = datapath("io", "data", "rdata", "ghg_df.rds")
    r_df = read_rdata(filename, engine="rscript")

    tm.assert_frame_equal(ghg_df, r_df.tail())


def test_read_rda_file(datapath):
    filename = datapath("io", "data", "rdata", "env_data_dfs.rda")
    r_dfs = read_rdata(filename, engine="rscript")

    assert list(r_dfs.keys()) == ["plants_df", "sea_ice_df", "ghg_df"]

    tm.assert_frame_equal(ghg_df, r_dfs["ghg_df"].tail())
    tm.assert_frame_equal(plants_df, r_dfs["plants_df"].tail())
    tm.assert_frame_equal(sea_ice_df, r_dfs["sea_ice_df"].tail())


def test_buffer_read_rds(datapath):
    filename = datapath("io", "data", "rdata", "sea_ice_df.rds")

    with open(filename, "rb") as f:
        r_df = read_rdata(f, file_format="rds", engine="rscript")

    r_df = adj_int(r_df)

    tm.assert_frame_equal(sea_ice_df, r_df.tail())


def test_bytes_read_rda(datapath):
    filename = datapath("io", "data", "rdata", "env_data_dfs.rda")

    with open(filename, "rb") as f:
        r_dfs = read_rdata(f.read(), file_format="rda", engine="rscript")

    r_dfs = {k: adj_int(v) for k, v in r_dfs.items()}

    assert list(r_dfs.keys()) == ["plants_df", "sea_ice_df", "ghg_df"]

    tm.assert_frame_equal(ghg_df, r_dfs["ghg_df"].tail())
    tm.assert_frame_equal(plants_df, r_dfs["plants_df"].tail())
    tm.assert_frame_equal(sea_ice_df, r_dfs["sea_ice_df"].tail())


def test_bytesio_rds(datapath):
    filename = datapath("io", "data", "rdata", "sea_ice_df.rds")

    with open(filename, "rb") as f:
        with BytesIO(f.read()) as b_io:
            r_df = read_rdata(b_io, file_format="rds", engine="rscript")

    r_df = adj_int(r_df)

    tm.assert_frame_equal(sea_ice_df, r_df.tail())


def test_bytesio_rda(datapath):
    filename = datapath("io", "data", "rdata", "env_data_dfs.rda")

    with open(filename, "rb") as f:
        with BytesIO(f.read()) as b_io:
            r_dfs = read_rdata(b_io, file_format="rda", engine="rscript")

    r_dfs = {k: adj_int(v) for k, v in r_dfs.items()}

    assert list(r_dfs.keys()) == ["plants_df", "sea_ice_df", "ghg_df"]

    tm.assert_frame_equal(ghg_df, r_dfs["ghg_df"].tail())
    tm.assert_frame_equal(plants_df, r_dfs["plants_df"].tail())
    tm.assert_frame_equal(sea_ice_df, r_dfs["sea_ice_df"].tail())


# FILE FORMAT


def test_read_wrong_format(datapath):
    with pytest.raises(ValueError, match="not a valid value for file_format"):
        filename = datapath("io", "data", "rdata", "plants_df.rds")
        read_rdata(filename, engine="rscript", file_format="r")


def test_read_wrong_file():
    with pytest.raises(FileNotFoundError, match="file cannot be found"):
        filename = os.path.join("data", "rdata", "plants_df.rda")
        read_rdata(filename, engine="rscript")


@pytest.mark.slow
def test_read_rds_non_dfs(datapath, mode):
    with pytest.raises(
        ValueError, match="No actual data frame or coercible data frames"
    ):
        filename = datapath("io", "data", "rdata", "ghg_t_tests.rds")
        read_rdata(filename, engine="rscript", mode=mode)


@pytest.mark.slow
def test_read_rda_non_dfs(datapath, mode):
    with pytest.raises(
        ValueError, match="No actual data frame or coercible data frames"
    ):
        filename = datapath("io", "data", "rdata", "env_data_non_dfs.rda")
        read_rdata(filename, engine="rscript", mode=mode)


def test_read_not_rda_file(datapath, mode):
    with pytest.raises(RScriptError, match="bad restore file magic number"):
        read_rdata(
            datapath("io", "data", "rdata", "ppm_df.csv"),
            file_format="rda",
            engine="rscript",
            mode=mode,
        )


def test_read_not_rds_file(datapath, mode):
    with pytest.raises(RScriptError, match="unknown input format"):
        read_rdata(
            datapath("io", "data", "rdata", "ppm_df.csv"),
            file_format="rds",
            engine="rscript",
            mode=mode,
        )


def test_bytes_read_infer_rds(datapath):
    filename = datapath("io", "data", "rdata", "sea_ice_df.rds")

    with pytest.raises(ValueError, match="Unable to infer file format from file name"):
        with open(filename, "rb") as f:
            read_rdata(f.read(), engine="rscript")


def test_bytes_read_infer_rda(datapath):
    filename = datapath("io", "data", "rdata", "env_data_dfs.rda")

    with pytest.raises(ValueError, match="Unable to infer file format from file name"):
        with open(filename, "rb") as f:
            read_rdata(f.read(), engine="rscript")


# URL


@tm.network
def test_read_rda_url():
    url_df = DataFrame(
        {
            "carrier": {"1": "9E", "2": "AA", "3": "AS", "4": "B6", "5": "DL"},
            "name": {
                "1": "Endeavor Air Inc.",
                "2": "American Airlines Inc.",
                "3": "Alaska Airlines Inc.",
                "4": "JetBlue Airways",
                "5": "Delta Air Lines Inc.",
            },
        }
    ).rename_axis("rownames")

    url = (
        "https://github.com/hadley/nycflights13/blob/master/data/airlines.rda?raw=true"
    )
    r_df = read_rdata(url, file_format="rda", engine="rscript")["airlines"]

    tm.assert_frame_equal(url_df, r_df.head())


@tm.network
def test_read_unable_infer_format():
    with pytest.raises(ValueError, match="Unable to infer file format from file name"):
        url = (
            "https://github.com/hadley/nycflights13/"
            "blob/master/data/airlines.rda?raw=true"
        )
        read_rdata(url, engine="rscript")


@tm.network
def test_read_wrong_url():
    with pytest.raises(HTTPError, match="HTTP Error 404: Not Found"):
        url = "https://example.com/data.rdata"
        read_rdata(url, engine="rscript")


# S3


@tm.network
@pytest.mark.slow
def test_read_rda_s3():
    s3 = "s3://assets.datacamp.com/production/course_1478/datasets/wine.RData"
    s3_df = DataFrame(
        {
            "Alcohol": {"1": 13.2, "2": 13.16, "3": 14.37, "4": 13.24, "5": 14.2},
            "Malic acid": {"1": 1.78, "2": 2.36, "3": 1.95, "4": 2.59, "5": 1.76},
            "Ash": {"1": 2.14, "2": 2.67, "3": 2.5, "4": 2.87, "5": 2.45},
            "Alcalinity of ash": {
                "1": 11.2,
                "2": 18.6,
                "3": 16.8,
                "4": 21.0,
                "5": 15.2,
            },
            "Magnesium": {"1": 100, "2": 101, "3": 113, "4": 118, "5": 112},
            "Total phenols": {"1": 2.65, "2": 2.8, "3": 3.85, "4": 2.8, "5": 3.27},
            "Flavanoids": {"1": 2.76, "2": 3.24, "3": 3.49, "4": 2.69, "5": 3.39},
            "Nonflavanoid phenols": {
                "1": 0.26,
                "2": 0.3,
                "3": 0.24,
                "4": 0.39,
                "5": 0.34,
            },
            "Proanthocyanins": {"1": 1.28, "2": 2.81, "3": 2.18, "4": 1.82, "5": 1.97},
            "Color intensity": {"1": 4.38, "2": 5.68, "3": 7.8, "4": 4.32, "5": 6.75},
            "Hue": {"1": 3.4, "2": 3.17, "3": 3.45, "4": 2.93, "5": 2.85},
            "Proline": {"1": 1050, "2": 1185, "3": 1480, "4": 735, "5": 1450},
        }
    ).rename_axis("rownames")
    r_dfs = read_rdata(s3, engine="rscript")

    tm.assert_frame_equal(s3_df, r_dfs["wine"].head())


# ENGINE


def test_read_rds_df_output(datapath):
    filename = datapath("io", "data", "rdata", "sea_ice_df.rds")
    r_dfs = read_rdata(filename, engine="rscript")

    assert isinstance(r_dfs, DataFrame)


def test_read_rda_dict_output(datapath):
    filename = datapath("io", "data", "rdata", "env_data_dfs.rda")
    r_dfs = read_rdata(filename, engine="rscript")

    assert isinstance(r_dfs, dict)
    assert list(r_dfs.keys()) == ["plants_df", "sea_ice_df", "ghg_df"]


def test_read_wrong_engine(datapath):
    with pytest.raises(ValueError, match="not a supported engine"):
        filename = datapath("io", "data", "rdata", "sea_ice_df.rds")
        read_rdata(filename, engine="rpy2")


# MODE


@pytest.mark.slow
def test_read_rds_mode_file(datapath, mode):
    filename = datapath("io", "data", "rdata", "ghg_df.rds")
    r_df = read_rdata(filename, engine="rscript", mode=mode)

    r_df = adj_int(r_df)

    tm.assert_frame_equal(ghg_df, r_df.tail())


@pytest.mark.slow
def test_read_rda_mode_file(datapath, mode):
    filename = datapath("io", "data", "rdata", "env_data_dfs.rda")
    r_dfs = read_rdata(filename, engine="rscript", mode=mode)

    if mode in ["parquet", "feather"]:
        (r_dfs["ghg_df"], r_dfs["plants_df"], r_dfs["sea_ice_df"]) = (
            adj_int(r_dfs["ghg_df"]),
            adj_int(r_dfs["plants_df"]),
            adj_int(r_dfs["sea_ice_df"]),
        )

    assert list(r_dfs.keys()) == ["plants_df", "sea_ice_df", "ghg_df"]

    tm.assert_frame_equal(ghg_df, r_dfs["ghg_df"].tail())
    tm.assert_frame_equal(plants_df, r_dfs["plants_df"].tail())
    tm.assert_frame_equal(sea_ice_df, r_dfs["sea_ice_df"].tail())


def test_read_wrong_mode(datapath):
    with pytest.raises(ValueError, match="not supported value for mode"):
        filename = datapath("io", "data", "rdata", "plants_df.rds")
        read_rdata(filename, engine="rscript", mode="pickle")


# USE_OBJECTS


def test_read_select_frames_rda_dfs(datapath):
    filename = datapath("io", "data", "rdata", "env_data_dfs.rda")
    r_dfs = read_rdata(
        filename, engine="rscript", select_frames=["ghg_df", "sea_ice_df"]
    )

    assert "plants_df" not in list(r_dfs.keys())
    assert "ghg_df" in list(r_dfs.keys())
    assert "sea_ice_df" in list(r_dfs.keys())


def test_read_select_frames_rda_objs(datapath):
    filename = datapath("io", "data", "rdata", "env_data_objs.rda")
    r_dfs = read_rdata(
        filename,
        engine="rscript",
        select_frames=["ppm_ts", "species_mtx", "plants_arry"],
    )

    assert "species_vec" not in list(r_dfs.keys())
    assert "ghg_df" not in list(r_dfs.keys())

    assert "ppm_ts" in list(r_dfs.keys())
    assert "species_mtx" in list(r_dfs.keys())
    assert "plants_arry" in list(r_dfs.keys())


def test_read_wrong_select_frames(datapath):
    with pytest.raises(TypeError, match="not a valid type for select_frames"):
        filename = datapath("io", "data", "rdata", "env_data_dfs.rda")
        read_rdata(filename, engine="rscript", select_frames="plants_df")


# ROWNAMES


def test_read_rownames_true_rds(datapath):
    filename = datapath("io", "data", "rdata", "sea_ice_df.rds")
    r_df = read_rdata(filename, engine="rscript", rownames=True)

    assert r_df.index.name == "rownames"


def test_read_rownames_false_rds(datapath):
    filename = datapath("io", "data", "rdata", "sea_ice_df.rds")
    r_df = read_rdata(filename, engine="rscript", rownames=False)

    assert r_df.index.name != "rownames"


def test_read_rownames_true_rda(datapath):
    filename = datapath("io", "data", "rdata", "env_data_dfs.rda")
    r_dfs = read_rdata(filename, engine="rscript", rownames=True)

    assert r_dfs["ghg_df"].index.name == "rownames"
    assert r_dfs["plants_df"].index.name == "rownames"
    assert r_dfs["sea_ice_df"].index.name == "rownames"


def test_read_rownames_false_rda(datapath):
    filename = datapath("io", "data", "rdata", "env_data_dfs.rda")
    r_dfs = read_rdata(filename, engine="rscript", rownames=False)

    assert r_dfs["ghg_df"].index.name != "rownames"
    assert r_dfs["plants_df"].index.name != "rownames"
    assert r_dfs["sea_ice_df"].index.name != "rownames"


# ENCODING


def test_non_utf8_data(datapath, rtype):
    filename = datapath("io", "data", "rdata", f"climate_non_utf8_df.{rtype}")

    expected = DataFrame(
        {
            "número": {
                "1": 1,
                "2": 2,
                "3": 3,
                "4": 4,
                "5": 5,
                "6": 6,
                "7": 7,
                "8": 8,
                "9": 9,
                "10": 10,
            },
            "punto central del climatismo": {
                "1": "Parada de la circulación de vuelco meridional del Atlántico",
                "2": "Desintegración de la capa de hielo de la Antártida occidental",
                "3": "Muerte de la selva amazónica",
                "4": "Cambio de monzón en África occidental",
                "5": "Permafrost e hidratos de metano",
                "6": "Muerte de los arrecifes de coral",
                "7": "Cambio de monzón de la India",
                "8": "Desintegración de la capa de hielo de Groenlandia",
                "9": "Desplazamiento del bosque boreal",
                "10": "Reducción del hielo marino del Ártico ",
            },
        },
        index=[str(i) for i in range(1, 11)],
    ).rename_axis("rownames")

    rdfs = read_rdata(filename, engine="rscript", encoding="iso-8859-1", mode="csv")

    output = rdfs["climate_df"] if rtype == "rda" else rdfs

    tm.assert_frame_equal(output, expected)


# RDA WRITER

# PATH_OR_BUFFER


@pytest.mark.slow
def test_write_read_file(datapath, rtype, mode):
    with tm.ensure_clean("test.out") as path:
        ghg_df.to_rdata(
            path, file_format=rtype, engine="rscript", mode=mode, index=False
        )
        r_dfs = read_rdata(
            path, file_format=rtype, engine="rscript", mode=mode, rownames=False
        )

        expected = ghg_df.reset_index(drop=True)
        output = r_dfs if rtype == "rds" else r_dfs["pandas_dataframe"]
        output["year"] = output["year"].astype("int64")

        tm.assert_frame_equal(output, expected)


@pytest.mark.slow
def test_write_read_bytes_io(datapath, rtype, mode):
    with BytesIO() as b_io:
        sea_ice_df.to_rdata(
            b_io, file_format=rtype, engine="rscript", mode=mode, index=False
        )
        r_dfs = read_rdata(
            b_io.getvalue(),
            file_format=rtype,
            engine="rscript",
            mode=mode,
            rownames=False,
        )

        expected = sea_ice_df.reset_index(drop=True)
        output = r_dfs if rtype == "rds" else r_dfs["pandas_dataframe"]
        output["year"] = output["year"].astype("int64")
        output["mo"] = output["mo"].astype("int64")

        tm.assert_frame_equal(output, expected)


# FILE_FORMAT


def test_write_rda_file(rtype):
    expected = """\
$pandas_dataframe
  rownames year mo data.type region extent area
1     1012 2016 12   Goddard      S   8.28 5.51
2     1013 2017 12   Goddard      S   9.48 6.23
3     1014 2018 12   Goddard      S   9.19 5.59
4     1015 2019 12   Goddard      S   9.41 6.59
5     1016 2020 12   NRTSI-G      S  10.44 6.50

"""
    with tm.ensure_clean_dir() as tmp_dir:
        out_file = os.path.join(tmp_dir, "rdata.out")
        r_code = os.path.join(tmp_dir, "r_test.R")

        sea_ice_df.to_rdata(out_file, file_format=rtype, engine="rscript")

        cmds = ["Rscript", r_code, out_file, rtype, "pandas_dataframe"]
        output = run_rscript(cmds)

        assert output == expected


def test_write_wrong_format():
    with tm.ensure_clean("test.rda") as path:
        with pytest.raises(ValueError, match=("not a valid value for file_format")):
            ghg_df.to_rdata(path, engine="rscript", file_format="csv")


def test_write_unable_to_infer():
    with tm.ensure_clean("test") as path:
        with pytest.raises(
            ValueError, match=("Unable to infer file format from file name")
        ):
            ghg_df.to_rdata(path, engine="rscript")


# ENGINE


@td.skip_if_no("pyreadr")
def test_write_engine_consistency(rtype):
    expected = """\
$pandas_dataframe
  rownames   plant_group              status count
1       16 Pteridophytes      Data Deficient   398
2       17 Pteridophytes             Extinct    65
3       18 Pteridophytes      Not Threatened  1294
4       19 Pteridophytes Possibly Threatened   408
5       20 Pteridophytes          Threatened  1275

"""
    with tm.ensure_clean_dir() as tmp_dir:
        out_file = os.path.join(tmp_dir, "rdata.out")
        r_code = os.path.join(tmp_dir, "r_test.R")

        plants_df.to_rdata(out_file, file_format=rtype, engine="pyreadr")
        cmds = ["Rscript", r_code, out_file, rtype, "pandas_dataframe"]
        pyr_output = run_rscript(cmds)

        plants_df.to_rdata(out_file, file_format=rtype, engine="rscript")
        cmds = ["Rscript", r_code, out_file, rtype, "pandas_dataframe"]
        rcomp_output = run_rscript(cmds)

        assert pyr_output == expected
        assert pyr_output == rcomp_output


def test_write_wrong_engine():
    with tm.ensure_clean("test.rda") as path:
        with pytest.raises(ValueError, match=("not a supported engine")):
            ghg_df.to_rdata(path, engine="rpy2")


# MODE


@pytest.mark.slow
def test_write_mode(rtype, mode):
    expected = """\
$pandas_dataframe
  rownames               gas year emissions
1      141    Carbon dioxide 2018 5424.8815
2      142           Methane 2018  634.4571
3      143     Nitrous oxide 2018  434.5286
4      144 Fluorinated gases 2018  182.7824
5      145             Total 2018 6676.6496

"""
    with tm.ensure_clean_dir() as tmp_dir:
        out_file = os.path.join(tmp_dir, "rdata.out")
        r_code = os.path.join(tmp_dir, "r_test.R")

        ghg_df.to_rdata(out_file, file_format=rtype, engine="rscript", mode=mode)
        cmds = ["Rscript", r_code, out_file, rtype, "pandas_dataframe"]
        output = run_rscript(cmds)

        assert output == expected


def test_write_wrong_mode():
    with tm.ensure_clean("test.rds") as path:
        with pytest.raises(ValueError, match=("not supported value for mode")):
            ghg_df.to_rdata(path, engine="rscript", mode="pickle")


# INDEX


@pytest.mark.slow
def test_write_index_false(rtype, mode):
    expected = """\
$pandas_dataframe
                gas year emissions
1    Carbon dioxide 2018 5424.8815
2           Methane 2018  634.4571
3     Nitrous oxide 2018  434.5286
4 Fluorinated gases 2018  182.7824
5             Total 2018 6676.6496

"""
    with tm.ensure_clean_dir() as tmp_dir:
        out_file = os.path.join(tmp_dir, "rdata.out")
        r_code = os.path.join(tmp_dir, "r_test.R")

        ghg_df.to_rdata(
            out_file, file_format=rtype, index=False, engine="rscript", mode=mode
        )

        cmds = ["Rscript", r_code, out_file, rtype, "pandas_dataframe"]
        output = run_rscript(cmds)

        assert output == expected


# ASCII


@pytest.mark.slow
def test_write_ascii_output(rtype, mode, ascii):
    expected = """\
$pandas_dataframe
  rownames               gas year emissions
1      141    Carbon dioxide 2018 5424.8815
2      142           Methane 2018  634.4571
3      143     Nitrous oxide 2018  434.5286
4      144 Fluorinated gases 2018  182.7824
5      145             Total 2018 6676.6496

"""
    with tm.ensure_clean_dir() as tmp_dir:
        out_file = os.path.join(tmp_dir, "rdata.out")
        r_code = os.path.join(tmp_dir, "r_test.R")

        ghg_df.to_rdata(
            out_file, file_format=rtype, engine="rscript", mode=mode, ascii=ascii
        )

        cmds = ["Rscript", r_code, out_file, rtype, "pandas_dataframe"]
        output = run_rscript(cmds)

        assert output == expected


def test_write_read_ascii(rtype):
    with tm.ensure_clean_dir() as tmp_dir:
        out_file = os.path.join(tmp_dir, "rdata.out")

        ghg_df.to_rdata(
            out_file,
            file_format=rtype,
            engine="rscript",
            index=False,
            ascii=True,
            compress=False,
        )

        with open(out_file) as f:
            r_dfs = read_rdata(f, file_format=rtype, engine="rscript", rownames=False)

        expected = ghg_df.reset_index(drop=True)
        output = r_dfs if rtype == "rds" else r_dfs["pandas_dataframe"]
        output["year"] = output["year"].astype("int64")

        tm.assert_frame_equal(output, expected)


# COMPRESS


@pytest.mark.slow
def test_write_compress_types(rtype, mode, comp):
    expected = """\
$pandas_dataframe
  rownames year mo data.type region extent area
1     1012 2016 12   Goddard      S   8.28 5.51
2     1013 2017 12   Goddard      S   9.48 6.23
3     1014 2018 12   Goddard      S   9.19 5.59
4     1015 2019 12   Goddard      S   9.41 6.59
5     1016 2020 12   NRTSI-G      S  10.44 6.50

"""
    with tm.ensure_clean_dir() as tmp_dir:
        out_file = os.path.join(tmp_dir, "rdata.out")
        r_code = os.path.join(tmp_dir, "r_test.R")

        sea_ice_df.to_rdata(
            out_file, file_format=rtype, engine="rscript", mode=mode, compress=comp
        )

        cmds = ["Rscript", r_code, out_file, rtype, "pandas_dataframe"]
        output = run_rscript(cmds)

        assert output == expected


def test_write_wrong_comp():
    with tm.ensure_clean("test.rds") as path:
        with pytest.raises(ValueError, match=("not a supported value for compress")):
            ghg_df.to_rdata(path, engine="rscript", compress="zip")


def test_write_none_comp():
    with tm.ensure_clean("test.rds") as path:
        with pytest.raises(RScriptError, match=("invalid 'compress' argument")):
            ghg_df.to_rdata(path, engine="rscript", compress=None)


# OTHER_FRAMES


@pytest.mark.slow
def test_write_other_frames(mode):
    expected = """\
$ghg_df
  rownames               gas year emissions
1      141    Carbon dioxide 2018 5424.8815
2      142           Methane 2018  634.4571
3      143     Nitrous oxide 2018  434.5286
4      144 Fluorinated gases 2018  182.7824
5      145             Total 2018 6676.6496

$plants_df
  rownames   plant_group              status count
1       16 Pteridophytes      Data Deficient   398
2       17 Pteridophytes             Extinct    65
3       18 Pteridophytes      Not Threatened  1294
4       19 Pteridophytes Possibly Threatened   408
5       20 Pteridophytes          Threatened  1275

$sea_ice_df
  rownames year mo data.type region extent area
1     1012 2016 12   Goddard      S   8.28 5.51
2     1013 2017 12   Goddard      S   9.48 6.23
3     1014 2018 12   Goddard      S   9.19 5.59
4     1015 2019 12   Goddard      S   9.41 6.59
5     1016 2020 12   NRTSI-G      S  10.44 6.50

"""
    with tm.ensure_clean_dir() as tmp_dir:
        out_file = os.path.join(tmp_dir, "rdata.rda")
        r_code = os.path.join(tmp_dir, "r_test.R")

        ghg_df.to_rdata(
            out_file,
            engine="rscript",
            other_frames=[plants_df, sea_ice_df],
            rda_names=["ghg_df", "plants_df", "sea_ice_df"],
            mode=mode,
        )

        cmds = ["Rscript", r_code, out_file, "rda", ""]
        output = run_rscript(cmds)

        assert output == expected


def test_write_other_frames_wrong_type():
    with tm.ensure_clean("test.rds") as path:
        with pytest.raises(
            TypeError, match=("objects in other_frames is not a DataFrame")
        ):
            ghg_df.to_rdata(
                path, engine="rscript", other_frames=plants_df, rda_names=["plants_df"]
            )


def test_write_read_other_frames(datapath):
    with tm.ensure_clean("test.rda") as path:
        ghg_df.to_rdata(
            path,
            engine="rscript",
            other_frames=[plants_df, sea_ice_df],
            rda_names=["ghg_df", "plants_df", "sea_ice_df"],
        )
        r_dfs = read_rdata(path, engine="rscript")

        assert list(r_dfs.keys()) == ["plants_df", "sea_ice_df", "ghg_df"]


# RDA NAMES


def test_write_mismatched_names_frames():
    with tm.ensure_clean("test.rds") as path:
        with pytest.raises(
            ValueError,
            match=("does not match number of current DataFrame and other_frames"),
        ):
            ghg_df.to_rdata(
                path,
                engine="rscript",
                other_frames=[plants_df, sea_ice_df],
                rda_names=["plants_df", "sea_ice_df"],
            )
