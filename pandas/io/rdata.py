from datetime import datetime
import io
import os
import platform
import subprocess
from tempfile import TemporaryDirectory
from typing import (
    Dict,
    List,
    Optional,
    Type,
    Union,
)

from pandas._typing import (
    Buffer,
    FilePathOrBuffer,
    StorageOptions,
)
from pandas.compat._optional import import_optional_dependency
from pandas.errors import (
    AbstractMethodError,
    ParserError,
)
from pandas.util._decorators import doc

from pandas.core.dtypes.common import is_list_like

from pandas.core.frame import DataFrame
from pandas.core.shared_docs import _shared_docs

from pandas.io.common import (
    file_exists,
    get_handle,
    is_fsspec_url,
    is_url,
    stringify_path,
)
from pandas.io.feather_format import read_feather
from pandas.io.parquet import read_parquet
from pandas.io.parsers import read_csv
from pandas.io.sql import read_sql


class RScriptError(Exception):
    """
    Exception raises when command line call to RScript throws a non-empty
    error message. Message will capture verbatim R output in console.
    """

    pass


def _executable_exists(name) -> bool:
    """
    Internal method to check if R exists on system.

    This method will return True if R is installed for Rscript command
    line call and if machine recognizes Rscript in Path env variable.
    """

    WHICH_CMD = "where" if platform.system() == "Windows" else "which"

    return (
        subprocess.call(
            [WHICH_CMD, name], stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )
        == 0
    )


RSCRIPT_EXISTS = _executable_exists("Rscript")


@doc(storage_options=_shared_docs["storage_options"])
def read_rdata(
    path_or_buffer: FilePathOrBuffer,
    file_format: str = "infer",
    engine: str = "pyreadr",
    mode: str = "csv",
    select_frames: Optional[List[str]] = None,
    rownames: bool = True,
    encoding: str = "utf-8",
    storage_options: StorageOptions = None,
) -> Union[DataFrame, Dict[str, DataFrame]]:
    r"""
    Read R data (.RData, .rda, .rds) into DataFrame or ``dict`` of DataFrames.

    .. versionadded:: 1.3.0

    Parameters
    ----------
    path_or_buffer : str, path object, or file-like object
        Any valid file path is acceptable. The string could be a URL.
        Valid URL schemes include http, ftp, s3, and file.

    file_format : {{'infer', 'rda', 'rdata', 'rds'}}, default 'infer'
        R serialization type as output from R's base::save or base::saveRDS
        commands. Default 'infer' will use extension in file name to
        to determine the format type.

    engine : {{'pyreadr'. 'rscript'}}, default 'pyreadr'
        Engine used to parse or read R data. Currently, two types are
        supported: ``pyreadr`` which requires the pyreadr package to be
        installed and ``rscript`` which requires R to be installed on machine.
        For ``rscript``, be sure the R bin installation folder is included in
        the system Path environment variable. The ``pyreadr`` is the faster
        parser to handle most needs but ``rscript`` engine provides fuller
        support of rda and rds formats since it calls native R commands.

    mode : {{'csv', 'parquet', 'feather', 'sqlite'}}, default 'csv'
        Python and R I/O transfer mode that only applies to ``rscript``
        engine (ignored for ``pyreadr``). Using ``csv`` (text approach), no
        additional packages are required. Using ``parquet`` or ``feather``
        (binary approach) requires pyarrow installed in Python and arrow
        package installed in R. Using ``sqlite`` (database approach) requires
        RSQLite package installed in R. Binary will usually be faster to process
        than text data. Database usually ensures data type integrity.

    select_frames : list, default None
        Selected names of DataFrames to return from R rda and RData types that
        can contain multiple objects.

    rownames : bool, default True
        Include original rownames in R data frames to map into a DataFrame index.

    encoding : str, optional, default 'utf-8'
        Encoding of R data. Currently, ``pyreadr`` engine only supports utf-8
        encoded data.

    {storage_options}

    Returns
    -------
    DataFrame or dict of DataFrames
        Depends on R data type where rds formats returns a single DataFrame and
        rda or RData formats return ``dict`` of DataFrames.

    See Also
    --------
    read_sas : Read SAS datasets into DataFrame.
    read_stata : Read Stata datasets into DataFrame.
    read_spss : Read SPSS datasets into DataFrame.

    Notes
    -----
    For ``pyreadr`` engine, any R data file that contains a non-data.frame object
    may raise parsing errors. For ``rscript`` engine, such objects will be
    ignored. Both methods will or attempt to return data.frame objects or any
    object that is coercible to R's data.frame such as matrix, tibble,
    and data.table. For arrays, method will attempt to convert to 2D
    structure and may not reproduce original R object representation.

    If object in rds types or all objects in rda or RData types are not  data
    frames, this method will raise an error and will not return None or an empty
    dictionary.

    For ``pyreadr`` engine, ``select_frames`` above is synonymous to ``use_objects``
    in package's `read_r` method. Also, ``timezone`` argument defaults to current
    system regional timezone in order to correspond to original date/times in R.

    Examples
    --------
    To read an .rds file which only contains a single object, below returns a
    DataFrame:

    >>> R_code = '''
    ... ghg_df <- data.frame(
    ...     gas = c('Carbon dioxide',
    ...             'Methane',
    ...             'Nitrous oxide',
    ...             'Fluorinated gases',
    ...             'Total'),
    ...     year = c(2018,
    ...              2018,
    ...              2018,
    ...              2018,
    ...              2018),
    ...     emissions = c(5424.88,
    ...                   634.46,
    ...                   434.53,
    ...                   182.78,
    ...                   6676.65)
    ... )
    ... saveRDS(ghg_df, file="ghg_df.rds")
    ... '''

    >>> ghg_df = pd.read_rdata("ghg_df.rds")  # doctest: +SKIP
    >>> ghg_df  # doctest: +SKIP
                            gas  year  emissions
    rownames
    1            Carbon dioxide  2018    5424.88
    2                   Methane  2018     634.46
    3             Nitrous oxide  2018     434.53
    4         Fluorinated gases  2018     182.78
    5                     Total  2018    6676.65

    To read an .rda or .RData file which can contain multiple objects, blue
    returns a ``dict`` of DataFrames:

    >>> R_code = '''
    ... plants_df <- pd.DataFrame(
    ...     plant_group = c('Pteridophytes',
    ...                     'Pteridophytes',
    ...                     'Pteridophytes',
    ...                     'Pteridophytes',
    ...                     'Pteridophytes'),
    ...     status = c('Data Deficient',
    ...                'Extinct',
    ...                'Not Threatened',
    ...                'Possibly Threatened',
    ...                'Threatened'),
    ...     count = c(398, 65, 1294, 408, 1275)
    ... )
    ... sea_ice_df <- pd.DataFrame(
    ...     year = c(2016, 2017, 2018, 2019, 2020),
    ...     mo = c(12, 12, 12, 12, 12],
    ...     data.type: c('Goddard',
    ...                  'Goddard',
    ...                  'Goddard',
    ...                  'Goddard',
    ...                  'NRTSI-G'),
    ...     region = c('S', 'S', 'S', 'S', 'S'),
    ...     extent = c(8.28, 9.48, 9.19, 9.41, 10.44),
    ...     area = c(5.51, 6.23, 5.59, 6.59, 6.5)
    ... )
    ... save(ghg_df, plants_df, sea_ice_df, file="env_data_dfs.rda")
    ... '''

    >>> env_dfs = pd.read_rdata("env_data_dfs.rda")  # doctest: +SKIP
    >>> env_dfs  # doctest: +SKIP
    {{'ghg_df':
                          gas  year  emissions
    rownames
    1          Carbon dioxide  2018    5424.88
    2                 Methane  2018     634.46
    3           Nitrous oxide  2018     434.53
    4       Fluorinated gases  2018     182.79
    5                   Total  2018    6676.65,
    'plants_df':
               plant_group               status  count
    rownames
    1        Pteridophytes       Data Deficient    398
    2        Pteridophytes              Extinct     65
    3        Pteridophytes       Not Threatened   1294
    4        Pteridophytes  Possibly Threatened    408
    5        Pteridophytes           Threatened   1275,
    'sea_ice_df':
           year  mo data.type region  extent  area
    rownames
    1      2016  12   Goddard      S    8.28  5.51
    2      2017  12   Goddard      S    9.48  6.23
    3      2018  12   Goddard      S    9.19  5.59
    4      2019  12   Goddard      S    9.41  6.59
    5      2020  12   NRTSI-G      S   10.44  6.50}}
    """

    return _parse(
        path_or_buffer=path_or_buffer,
        file_format=file_format,
        engine=engine,
        mode=mode,
        select_frames=select_frames,
        rownames=rownames,
        encoding=encoding,
        storage_options=storage_options,
    )


def _parse(
    path_or_buffer,
    file_format,
    engine,
    mode,
    select_frames,
    rownames,
    encoding,
    storage_options,
    **kwargs,
) -> Union[DataFrame, Dict[str, DataFrame]]:
    """
    Call internal parser classes.

    This method will conditionally call internal parsers:
    _PyReadrParser or _RscriptParser.

    Raises
    ------
    FileNotFoundError
        * If Rscript bin executable is not installed or found on machine.

    ImportError
        * If pyreadr for engine and pyarrow for mode is not installed.

    ValueError
        * If engine is neither pyreadr or rscript.
    """
    pyreadr = import_optional_dependency("pyreadr", errors="ignore")
    pyarrow = import_optional_dependency("pyarrow", errors="ignore")

    RDataReader: Union[Type[_PyReadrParser], Type[_RscriptParser]]

    if engine == "pyreadr":
        if pyreadr is None:
            raise ImportError("pyreadr not found, please install for this engine.")

        RDataReader = _PyReadrParser

    elif engine == "rscript":
        if RSCRIPT_EXISTS is None:
            raise FileNotFoundError(
                "R is either not installed on this system or its "
                "bin folder is not in Path environment variable."
            )

        if pyarrow is None and mode in ["parquet", "feather"]:
            raise ImportError("pyarrow not found, please install for this mode.")

        RDataReader = _RscriptParser
    else:
        raise ValueError(f"{engine} is not a supported engine.")

    rdr = RDataReader(
        path_or_buffer,
        file_format,
        engine,
        mode,
        select_frames,
        rownames,
        encoding,
        storage_options,
    )

    return rdr.parse_data()


def _get_data_from_filepath(
    filepath_or_buffer,
    encoding,
    compression,
    storage_options,
) -> Union[str, bytes, Buffer]:
    """
    Extract raw R data.

    The method accepts three input types:
        1. filepath (string-like)
        2. file-like object (e.g. open file object, BytesIO)
        3. R data file in ascii or binary content

    This method turns (1) into (2) to simplify the rest of the processing.
    It returns input types (2) and (3) unchanged.
    """
    filepath_or_buffer = stringify_path(filepath_or_buffer)

    if (
        not isinstance(filepath_or_buffer, str)
        or is_url(filepath_or_buffer)
        or is_fsspec_url(filepath_or_buffer)
        or file_exists(filepath_or_buffer)
    ):
        with get_handle(
            filepath_or_buffer,
            "rb",
            encoding=encoding,
            compression=compression,
            storage_options=storage_options,
            is_text=False,
        ) as handle_obj:
            filepath_or_buffer = (
                handle_obj.handle.read()
                if hasattr(handle_obj.handle, "read")
                else handle_obj.handle
            )
    else:
        raise FileNotFoundError(f"{filepath_or_buffer} file cannot be found.")

    return filepath_or_buffer


def _preprocess_data(data) -> Union[io.StringIO, io.BytesIO]:
    """
    Convert extracted raw data.

    This method will return underlying data of extracted R data formats.
    The data either has a `read` attribute (e.g. a file object or a
    StringIO/BytesIO) or is bytes that represents the R data.
    """

    if isinstance(data, str):
        data = io.StringIO(data)

    elif isinstance(data, bytes):
        data = io.BytesIO(data)

    return data


class _RDataReader:
    """
    Internal subclass to parse R data files into dict of DataFrames.

    Parameters
    ----------
    path_or_buffer : a valid str, path object or file-like object
        Any valid string path is acceptable. The string could be a URL. Valid
        URL schemes include http, ftp, s3, and file.

    file_format : {{'infer', 'rda', 'rdata', 'rds'}}, default 'infer'
        R serialization type.

    engine : {{'pyreadr', 'rscript'}}, default 'pyreadr'
        Engine used to parse or read R data.

    mode : {{'csv', 'parquet', 'feather', 'sqlite'}}, default 'csv'
        Python and R i/o transfer mode.

    select_frames : list, default None
        Selected names of DataFrames to return from R data.

    rownames : bool, default True
        Include original rownames in R data frames.

    encoding : str, optional, default 'utf-8'
        Encoding of R data.

    storage_options : dict, optional
        Extra options that make sense for a particular storage connection,
        e.g. host, port, username, password, etc.,

    See also
    --------
    pandas.io.rdata._PyReadrParser
    pandas.io.rdata._RscriptParser

    Notes
    -----
    To subclass this class effectively you must override the following methods:`
        * :func:`handle_rownames`
        * :func:`parse_data`

    See each method's respective documentation for details on their
    functionality.
    """

    def __init__(
        self,
        path_or_buffer,
        file_format,
        engine,
        mode,
        select_frames,
        rownames,
        encoding,
        storage_options,
    ) -> None:
        self.path_or_buffer = path_or_buffer
        self.file_format = file_format.lower()
        self.engine = engine
        self.mode = mode
        self.select_frames = select_frames
        self.rownames = rownames
        self.encoding = encoding
        self.storage_options = storage_options

    def verify_params(self) -> None:
        """
        Verify user entries of parameters.

        This method will check the values and types of select parameters
        and raise appropriate errors.
        """

        if self.file_format not in ["infer", "rda", "rdata", "rds"]:
            raise ValueError(
                f"'{self.file_format}' is not a valid value for file_format"
            )

        if (
            self.file_format == "infer"
            and isinstance(self.path_or_buffer, str)
            and not self.path_or_buffer.lower().endswith((".rda", ".rdata", ".rds"))
        ) or (self.file_format == "infer" and not isinstance(self.path_or_buffer, str)):
            raise ValueError(
                f"Unable to infer file format from file name: {self.path_or_buffer}. "
                "Please use known R data type (.rda, .rdata, .rds)."
            )

        if self.file_format == "infer":
            self.file_format = os.path.splitext(self.path_or_buffer.lower())[1][1:]

        if self.mode is not None and self.mode not in [
            "csv",
            "feather",
            "parquet",
            "sqlite",
        ]:
            raise ValueError(f"'{self.mode}' is not supported value for mode.")

        if self.select_frames is not None and not is_list_like(self.select_frames):
            raise TypeError(
                f"{type(self.select_frames).__name__} is "
                "not a valid type for select_frames"
            )

    def buffer_to_disk(self, tmp_dir: str) -> str:
        """
        Convert path or buffer to disk file.

        This method will convert path_or_buffer to temp file
        for pyreadr to parse and rscript to import.
        """

        r_temp = os.path.join(tmp_dir, "rdata.rda")

        handle_data = _get_data_from_filepath(
            filepath_or_buffer=self.path_or_buffer,
            encoding=self.encoding,
            compression=None,
            storage_options=self.storage_options,
        )

        with _preprocess_data(handle_data) as r_data:
            mode = "wb" if isinstance(r_data, io.BytesIO) else "w"
            with open(r_temp, mode) as f:
                f.write(r_data.read())

        return r_temp

    def handle_row_names(self) -> DataFrame:
        """
        Migrate R rownames to DataFrame index.

        This method will conditionally adjust index to reflect
        original R rownames.
        """

        raise AbstractMethodError(self)

    def parse_data(self) -> Union[DataFrame, Dict[str, DataFrame]]:
        """
        Parse R data files.

        This method will run engine methods to return a single DataFrame
        for rds type or dictionary of DataFrames for RData or rda types.
        """

        raise AbstractMethodError(self)


class _PyReadrParser(_RDataReader):
    """
    Internal class to parse R data types using third-party
    package, pyreadr.
    """

    def __init__(self, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs)
        self.verify_params()

    def handle_rownames(self, df) -> DataFrame:
        if not self.rownames:
            df = df.reset_index(drop=True)
            df.index.name = None

        if self.rownames and df.index.name != "rownames":
            df.index.name = "rownames"
            if df.index[0] == 0:
                df.index += 1

        return df

    def parse_data(self) -> Union[DataFrame, Dict[str, DataFrame]]:
        from pyreadr import read_r

        tz = datetime.now().astimezone().tzinfo
        with TemporaryDirectory() as tmp_dir:
            r_temp = self.buffer_to_disk(tmp_dir)
            rdata = read_r(r_temp, use_objects=self.select_frames, timezone=tz)

        rdata = {k: self.handle_rownames(df) for k, df in rdata.items()}
        rdata = rdata[None] if self.file_format == "rds" else dict(rdata)

        return rdata


class _RscriptParser(_RDataReader):
    """
    Internal class to parse R data types using temp script and data
    files and command line call to installed Rscript executable.
    """

    def __init__(self, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs)
        self.verify_params()

    def handle_rownames(self, df) -> DataFrame:
        if self.rownames:
            df = df.set_index("rownames")
        else:
            df = df.drop(["rownames"], axis=1)

        return df

    def run_rscript(self, tmp_dir, r_batch, cmds) -> str:
        """
        Run R script at command line.

        This method will call subprocess.Popen to run R script that
        saves temp data and meta files and returns R's console output.
        """

        with open(cmds[1], "w") as f:
            f.write(r_batch)

        p = subprocess.Popen(
            cmds,
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            cwd=tmp_dir,
        )
        output, error = p.communicate()
        if len(error) != 0:
            raise RScriptError(error.decode(self.encoding))

        return output.decode(self.encoding)

    def parse_data(self) -> Union[DataFrame, Dict[str, DataFrame]]:
        self.r_to_py_types = {
            "logical": "bool",
            "integer": "int64",
            "numeric": "float64",
            "factor": "category",
            "character": "str",
            "Date": "date",
            "POSIXct": "date",
        }

        switch_board = {
            "rda": {
                "csv": self.read_rdata_csv,
                "feather": self.read_rdata_feather,
                "parquet": self.read_rdata_parquet,
                "sqlite": self.read_rdata_sqlite,
            },
            "rdata": {
                "csv": self.read_rdata_csv,
                "feather": self.read_rdata_feather,
                "parquet": self.read_rdata_parquet,
                "sqlite": self.read_rdata_sqlite,
            },
            "rds": {
                "csv": self.read_rds_csv,
                "feather": self.read_rds_feather,
                "parquet": self.read_rds_parquet,
                "sqlite": self.read_rds_sqlite,
            },
        }

        rdata: Union[DataFrame, Dict[str, DataFrame], None]
        rdata = switch_board[self.file_format][self.mode]()

        rdata = (
            {k: v for k, v in rdata.items() if k in self.select_frames}
            if self.select_frames
            else rdata
        )
        rdata = {k: self.handle_rownames(df) for k, df in rdata.items()}

        rdata = rdata or None
        rdata = (
            rdata["r_df"]
            if (self.file_format == "rds" and rdata is not None)
            else rdata
        )

        if rdata is None:
            raise ValueError(
                "No actual data frame or coercible data frames found in R data file."
            )
        return rdata

    def read_rdata_csv(self) -> Dict[str, DataFrame]:
        """
        Read R rda data via IO csv.

        This method will call `load` and `write.csv` in R to export all
        data frames and metadata into temp csv files for pandas `read_csv`.     .
        """

        r_batch = """
                  args <- commandArgs(trailingOnly=TRUE)
                  load(args[1],  temp_env <- new.env())

                  env_list <- as.list.environment(temp_env)
                  rm(temp_env)

                  output_data_meta <- function(obj, nm) {
                    df <- tryCatch(data.frame(obj,
                                     check.names=FALSE,
                                     stringsAsFactors=FALSE
                                   ), error=function(e) NULL)

                    if (!is.null(df)) {
                        cat(nm, "\n", sep="")

                        df <- data.frame(rownames = row.names(df), df,
                                         check.names=FALSE,
                                         stringsAsFactors=FALSE)
                        writeLines(
                            c(paste0(colnames(df), collapse=","),
                              paste0(sapply(df,
                                     function(x) class(x)[1]), collapse=",")),
                            con=paste0("meta_", nm, ".txt")
                        )

                        write.csv(df, paste0("data_", nm, ".csv"),
                                  row.names=FALSE, na="")
                    }
                  }

                  output <- mapply(output_data_meta, env_list, names(env_list))
                  """

        with TemporaryDirectory() as tmp_dir:
            r_file = os.path.join(tmp_dir, "r_batch.R")
            rda_file = self.buffer_to_disk(tmp_dir)

            output = self.run_rscript(tmp_dir, r_batch, ["Rscript", r_file, rda_file])

            oline: str
            dfs: Dict[str, DataFrame] = {}
            for oline in filter(None, output.strip().split("\n")):
                with open(
                    os.path.join(tmp_dir, f"meta_{oline}.txt"),
                    encoding=self.encoding,
                ) as f:
                    flines = [fline.strip() for fline in f]

                r_hdrs: List[List[str]] = [h.split(",") for h in flines]
                py_types = {n: self.r_to_py_types[d] for n, d in zip(*r_hdrs)}

                dt_cols = [col for col, d in py_types.items() if d == "date"]
                py_types = {k: v for k, v in py_types.items() if v != "date"}

                try:
                    dfs[oline] = read_csv(
                        os.path.join(tmp_dir, f"data_{oline}.csv"),
                        dtype=py_types,  # type: ignore[arg-type]
                        parse_dates=dt_cols,
                        encoding=self.encoding,
                    )
                except (ParserError, ValueError):
                    dfs[oline] = read_csv(
                        os.path.join(tmp_dir, f"data_{oline}.csv"),
                        encoding=self.encoding,
                    )

        return dfs

    def read_rdata_feather(self) -> Dict[str, DataFrame]:
        """
        Read R rda data via IO feather.

        This method will call `readRDS` and `write_feather` in R to export all
        data frames into temp feather files for pandas `read_feather`.
        """

        r_batch = """
                  suppressPackageStartupMessages(library(arrow))

                  args <- commandArgs(trailingOnly=TRUE)

                  load(args[1],  temp_env <- new.env())
                  env_list <- as.list.environment(temp_env)
                  rm(temp_env)

                  output_data_meta <- function(obj, nm) {
                    df <- tryCatch(data.frame(obj,
                                    check.names=FALSE,
                                    stringsAsFactors=FALSE
                                   ), error=function(e) NULL)

                    if (!is.null(df)) {
                       cat(nm, "\n", sep="")
                       df <- data.frame(rownames = row.names(df), df,
                                        check.names=FALSE,
                                        stringsAsFactors=FALSE)
                       arrow::write_feather(df, paste0("data_", nm, ".feather"))
                    }
                  }

                  output <- mapply(output_data_meta, env_list, names(env_list))
                  """

        with TemporaryDirectory() as tmp_dir:
            r_file = os.path.join(tmp_dir, "r_batch.R")
            rda_file = self.buffer_to_disk(tmp_dir)

            output = self.run_rscript(tmp_dir, r_batch, ["Rscript", r_file, rda_file])

            oline: str
            dfs: Dict[str, DataFrame] = {
                oline: read_feather(os.path.join(tmp_dir, f"data_{oline}.feather"))
                for oline in filter(None, output.strip().split("\n"))
            }

        return dfs

    def read_rdata_parquet(self) -> Dict[str, DataFrame]:
        """
        Read R rda data via IO parquet.

        This method will call `load` and `write_parquet` in R to export all
        data frames into temp parquet files for pandas `read_parquet`.
        """

        r_batch = """
                  suppressPackageStartupMessages(library(arrow))

                  args <- commandArgs(trailingOnly=TRUE)

                  load(args[1],  temp_env <- new.env())
                  env_list <- as.list.environment(temp_env)
                  rm(temp_env)

                  output_data_meta <- function(obj, nm) {
                    df <- tryCatch(data.frame(obj,
                                    check.names=FALSE,
                                    stringsAsFactors=FALSE
                                   ), error=function(e) NULL)

                    if (!is.null(df)) {
                       cat(nm, "\n", sep="")
                       df <- data.frame(rownames = row.names(df), df,
                                        check.names=FALSE,
                                        stringsAsFactors=FALSE)
                       arrow::write_parquet(df, paste0("data_", nm, ".parquet"))
                    }
                  }

                  output <- mapply(output_data_meta, env_list, names(env_list))
                  """

        with TemporaryDirectory() as tmp_dir:
            r_file = os.path.join(tmp_dir, "r_batch.R")
            rda_file = self.buffer_to_disk(tmp_dir)

            output = self.run_rscript(tmp_dir, r_batch, ["Rscript", r_file, rda_file])

            oline: str
            dfs: Dict[str, DataFrame] = {
                oline: read_parquet(os.path.join(tmp_dir, f"data_{oline}.parquet"))
                for oline in filter(None, output.strip().split("\n"))
            }

        return dfs

    def read_rdata_sqlite(self) -> Dict[str, DataFrame]:
        """
        Read R rda data via IO sql.

        This method will call `load` and `dbWriteTable` in R to export all
        data frames into a temp SQLite database for pandas `read_sql`.
        """
        import sqlite3

        r_batch = """
                  suppressPackageStartupMessages(library(RSQLite))

                  args <- commandArgs(trailingOnly=TRUE)

                  load(args[1],  temp_env <- new.env())
                  env_list <- as.list.environment(temp_env)
                  rm(temp_env)

                  conn <- dbConnect(RSQLite::SQLite(), "r_data.db")
                  output_data_meta <- function(obj, nm) {
                    df <- tryCatch(data.frame(obj,
                                   check.names=FALSE,
                                   stringsAsFactors=FALSE
                                   ), error=function(e) NULL)

                    if (!is.null(df)) {
                       cat(nm, "\n", sep="")
                       df <- data.frame(rownames = row.names(df), df,
                                        check.names=FALSE,
                                        stringsAsFactors=FALSE)
                       dbWriteTable(conn, paste0("data_", nm), df, row.names=FALSE)
                    }
                  }

                  output <- mapply(output_data_meta, env_list, names(env_list))
                  dbDisconnect(conn)
                  """

        with TemporaryDirectory() as tmp_dir:
            r_db = os.path.join(tmp_dir, "r_data.db")
            r_file = os.path.join(tmp_dir, "r_batch.R")
            rda_file = self.buffer_to_disk(tmp_dir)

            output = self.run_rscript(tmp_dir, r_batch, ["Rscript", r_file, rda_file])

            oline: str
            conn = sqlite3.connect(r_db)
            dfs: Dict[str, DataFrame] = {
                oline: read_sql(f"SELECT * FROM data_{oline}", conn)
                for oline in filter(None, output.strip().split("\n"))
            }
            conn.close()

        return dfs

    def read_rds_csv(self) -> Dict[str, DataFrame]:
        """
        Read R rds data via IO csv.

        This method will call `readRDS` and `write.csv` in R to export single
        data frame and metadata into temp csv files for pandas `read_csv`.
        """

        r_batch = """
                  args <- commandArgs(trailingOnly=TRUE)

                  raw <- readRDS(args[1])
                  df <- tryCatch(data.frame(raw,
                                   check.names=FALSE,
                                   stringsAsFactors=FALSE
                                 ), error = function(e) NULL)

                  if(!is.null(df)) {
                      df <- data.frame(rownames = row.names(df), df,
                                       check.names=FALSE,
                                       stringsAsFactors=FALSE)
                      write.csv(df, file=args[2], row.names=FALSE)

                      cat(paste0(colnames(df), collapse=","),"|",
                          paste0(sapply(df, function(x)
                                       class(x)[1]), collapse=","),
                          sep="")
                  }
                  """

        dfs: Dict[str, DataFrame] = {}
        with TemporaryDirectory() as tmp_dir:
            r_data = os.path.join(tmp_dir, "r_data.csv")
            r_file = os.path.join(tmp_dir, "r_batch.R")

            rds_file = self.buffer_to_disk(tmp_dir)
            output = self.run_rscript(
                tmp_dir, r_batch, ["Rscript", r_file, rds_file, r_data]
            )

            if os.path.isfile(r_data):
                r_hdrs = [h.split(",") for h in output.split("|")]
                n: str
                py_types = {n: self.r_to_py_types[d] for n, d in zip(*r_hdrs)}

                dt_cols = [col for col, d in py_types.items() if d == "date"]
                py_types = {k: v for k, v in py_types.items() if v != "date"}

                try:
                    dfs["r_df"] = read_csv(
                        r_data,
                        dtype=py_types,  # type: ignore[arg-type]
                        parse_dates=dt_cols,
                        encoding=self.encoding,
                    )
                except (ParserError, ValueError):
                    dfs["r_df"] = read_csv(r_data)

        return dfs

    def read_rds_feather(self) -> Dict[str, DataFrame]:
        """
        Read R rds data via IO feather.

        This method will call `readRDS` and `write_feather` in R to export single
        data frame into a temp feather file for pandas `read_feather`.
        """

        r_batch = """
                  suppressPackageStartupMessages(library(arrow))
                  args <- commandArgs(trailingOnly=TRUE)

                  raw <- readRDS(args[1])
                  df <- tryCatch(data.frame(raw,
                                   check.names=FALSE,
                                   stringsAsFactors=FALSE
                                 ), error = function(e) NULL)

                  if(!is.null(df)) {
                      df <- data.frame(rownames = row.names(df), df,
                                       check.names=FALSE,
                                       stringsAsFactors=FALSE)
                      arrow::write_feather(df, args[2])
                  }
                  """

        with TemporaryDirectory() as tmp_dir:
            r_data = os.path.join(tmp_dir, "r_data.feather")
            r_file = os.path.join(tmp_dir, "r_batch.R")

            rds_file = self.buffer_to_disk(tmp_dir)
            self.run_rscript(tmp_dir, r_batch, ["Rscript", r_file, rds_file, r_data])

            dfs: Dict[str, DataFrame] = (
                {"r_df": read_feather(r_data)} if os.path.isfile(r_data) else {}
            )

        return dfs

    def read_rds_parquet(self) -> Dict[str, DataFrame]:
        """
        Read R rds data via IO parquet.

        This method will call `readRDS` and `write_parquet` in R to export
        single data frame into a temp parquet file for pandas `read_parquet`.
        """

        r_batch = """
                  suppressPackageStartupMessages(library(arrow))
                  args <- commandArgs(trailingOnly=TRUE)

                  raw <- readRDS(args[1])
                  df <- tryCatch(data.frame(raw,
                                   check.names=FALSE,
                                   stringsAsFactors=FALSE
                                 ), error = function(e) NULL)

                  if(!is.null(df)) {
                      df <- data.frame(rownames = row.names(df), df,
                                       check.names=FALSE,
                                       stringsAsFactors=FALSE)
                      arrow::write_parquet(df, args[2])
                  }
                  """

        with TemporaryDirectory() as tmp_dir:
            r_data = os.path.join(tmp_dir, "r_data.parquet")
            r_file = os.path.join(tmp_dir, "r_batch.R")

            rds_file = self.buffer_to_disk(tmp_dir)
            self.run_rscript(tmp_dir, r_batch, ["Rscript", r_file, rds_file, r_data])

            dfs: Dict[str, DataFrame] = (
                {"r_df": read_parquet(r_data, engine="pyarrow")}
                if os.path.isfile(r_data)
                else {}
            )

        return dfs

    def read_rds_sqlite(self) -> Dict[str, DataFrame]:
        """
        Read R rds data via IO sql.

        This method will call `readRDS` and `dbWriteTable` in R to export
        single data frame into a temp SQLite database for pandas `read_sql`.
        """
        import sqlite3

        r_batch = """
                  suppressPackageStartupMessages(library(RSQLite))
                  args <- commandArgs(trailingOnly=TRUE)

                  raw <- readRDS(args[1])
                  df <- tryCatch(data.frame(raw,
                                   check.names=FALSE,
                                   stringsAsFactors=FALSE
                                 ), error = function(e) NULL)

                  if(!is.null(df)) {
                      conn <- dbConnect(RSQLite::SQLite(), args[2])
                      df <- data.frame(rownames = row.names(df), df,
                                       check.names=FALSE,
                                       stringsAsFactors=FALSE)
                      dbWriteTable(conn, "rdata", df, row.names=FALSE)
                      dbDisconnect(conn)
                  }
                  """

        dfs: Dict[str, DataFrame] = {}
        with TemporaryDirectory() as tmp_dir:
            r_data = os.path.join(tmp_dir, "r_data.db")
            r_file = os.path.join(tmp_dir, "r_batch.R")

            rds_file = self.buffer_to_disk(tmp_dir)
            self.run_rscript(tmp_dir, r_batch, ["Rscript", r_file, rds_file, r_data])

            if os.path.isfile(r_data):
                conn = sqlite3.connect(r_data)
                dfs["r_df"] = read_sql("SELECT * FROM rdata", conn)
                conn.close()

        return dfs


class RDataWriter:
    """
    Subclass to write pandas DataFrames into R data files.

    Parameters
    ----------
    path_or_buffer : a valid str, path object or file-like object
        Any valid string path is acceptable.

    file_format : {{'infer', 'rda', 'rdata', 'rds'}}, default 'infer'
        R serialization type.

    engine : {{'rscript','pyreadr'}}, default 'utf-8'
        Engine used to write R data.

    mode : {{'csv', 'parquet', 'feather'}}, default 'csv'
        Python and R i/o transfer mode.

    other_frames : list, optional
        Other DataFrames to be included in rda (not rds) files
        that can contain multiple objects.

    rda_names : list, default ["pandas_dataframe"]
        Names for all exported objects in rda file.

    index : bool, default True
        Include index or MultiIndex in output as separate columns.

    ascii : bool, default False
        Write data in ASCII representation.

    compress : bool or {{'gzip', 'bzip2', 'xz'}}, default 'gzip'
        Compression types for R data. For pyreadr engine, only gzip
        is supported. Use False for uncompressed files.

    encoding : str, optional, default 'utf-8'
        Encoding of R data.

    storage_options : dict, optional
        Extra options that make sense for a particular storage connection,
        e.g. host, port, username, password, etc.

    See also
    --------
    pandas.io.rdata.PyReadrWriter
    pandas.io.rdata.RscriptWriter

    Notes
    -----
    To subclass this class effectively you must override the following methods:`
        * :func:`write_data`

    See each method's respective documentation for details on their
    functionality.
    """

    def __init__(
        self,
        frame: DataFrame,
        path_or_buffer: FilePathOrBuffer,
        file_format: str = "infer",
        engine: str = "rscript",
        mode: str = "csv",
        other_frames: Optional[List[DataFrame]] = None,
        rda_names: List[str] = ["pandas_dataframe"],
        index: bool = True,
        ascii: bool = False,
        compress: Union[bool, str] = "gzip",
        encoding: str = "utf-8",
        storage_options: StorageOptions = None,
    ) -> None:
        self.frame = frame
        self.path_or_buffer = path_or_buffer
        self.file_format = file_format.lower()
        self.engine = engine
        self.mode = mode
        self.other_frames = other_frames
        self.rda_names = rda_names
        self.index = index
        self.ascii = ascii
        self.compress = compress
        self.encoding = encoding
        self.storage_options = storage_options

    def verify_params(self) -> None:
        """
        Verify user entries of parameters.

        This method will check the values and types of select parameters
        and raise appropriate errors.
        """

        if self.file_format not in ["infer", "rda", "rdata", "rds"]:
            raise ValueError(
                f"{self.file_format} is not a valid value for file_format."
            )

        if (
            self.file_format == "infer"
            and isinstance(self.path_or_buffer, str)
            and not self.path_or_buffer.lower().endswith((".rda", ".rdata", ".rds"))
        ):
            raise ValueError(
                f"Unable to infer file format from file name: {self.path_or_buffer}"
                "Please use known R data type (.rda, .rdata, .rds)."
            )

        if self.file_format == "infer" and isinstance(self.path_or_buffer, str):
            self.file_format = os.path.splitext(self.path_or_buffer.lower())[1][1:]

        if self.mode is not None and self.mode not in [
            "csv",
            "feather",
            "parquet",
            "sqlite",
        ]:
            raise ValueError(f"{self.mode} is not supported value for mode.")

        if self.other_frames is not None and not is_list_like(self.other_frames):
            raise TypeError(
                f"{type(self.other_frames).__name__} is not "
                " a valid type for other_frames."
            )
        elif self.other_frames is not None:
            for df in self.other_frames:
                if not isinstance(df, DataFrame):
                    raise TypeError(
                        "One or more of the objects in "
                        "other_frames is not a DataFrame."
                    )

        if self.rda_names is not None and not is_list_like(self.rda_names):
            raise TypeError(
                f"{type(self.rda_names).__name__} is not a valid type for rda_names."
            )

        if self.compress is not None and self.compress not in [
            True,
            False,
            "gzip",
            "bzip2",
            "xz",
        ]:
            raise ValueError(f"{self.compress} is not a supported value for compress.")

    def disk_to_buffer(self, r_file: str) -> None:
        """
        Save temp file to path or buffer.

        This method will convert written R data to path_or_buffer.
        """

        with open(r_file, "rb") as rdata:
            with get_handle(
                self.path_or_buffer,
                "wb",
                compression=None,
                storage_options=self.storage_options,
                is_text=False,
            ) as handles:
                handles.handle.write(rdata.read())  # type: ignore[arg-type]

        return None

    def write_data(self) -> None:
        """
        Write DataFrames to R data files.

        This method will run engine methods to export DataFrames
        to R data files.
        """

        raise AbstractMethodError(self)


class PyReadrWriter(RDataWriter):
    """
    Main class called in `pandas.core.frame` to write DataFrame to R
    data types using third-party package, pyreadr.
    """

    def __init__(self, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs)
        self.verify_params()

    def write_data(self) -> None:
        from pyreadr import (
            write_rdata,
            write_rds,
        )

        self.frame = (
            self.frame.reset_index()
            if self.index
            else self.frame.reset_index(drop=True)
        )

        with TemporaryDirectory() as tmp_dir:
            r_temp = os.path.join(tmp_dir, "rdata.rda")

            if self.file_format in ["rda", "rdata"]:
                write_rdata(
                    path=r_temp,
                    df=self.frame,
                    df_name=self.rda_names[0],
                    compress=self.compress,
                )
            elif self.file_format == "rds":
                write_rds(path=r_temp, df=self.frame, compress=self.compress)

            self.disk_to_buffer(r_temp)

        return None


class RscriptWriter(RDataWriter):
    """
    Main class called in `pandas.core.frame` to write DataFrame(s) to R
    data types using command line to Rscript.
    """

    def __init__(self, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs)
        self.verify_params()
        self.handle_objects()

    def handle_objects(self) -> None:

        self.all_frames = (
            [self.frame] + self.other_frames if self.other_frames else [self.frame]
        )

        if len(self.rda_names) != len(self.all_frames):
            raise ValueError(
                f"Length of {self.rda_names} does not match number "
                "of current DataFrame and other_frames"
            )

        return None

    def run_rscript(self, tmp_dir, r_batch, cmds) -> None:
        """
        Run R script at command line.

        This method will call subprocess.Popen to run R script
        and return only non-empty error R output in console.
        """

        with open(cmds[1], "w") as f:
            f.write(r_batch)

        a = subprocess.Popen(
            cmds,
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            cwd=tmp_dir,
        )
        output, error = a.communicate()
        if len(error) != 0:
            raise RScriptError(error.decode(self.encoding))

        return None

    def write_data(self) -> None:
        self.py_to_r_types = {
            "int32": "integer",
            "int64": "integer",
            "float64": "numeric",
            "category": "factor",
            "object": "character",
            "bool": "logical",
            "datetime64[ns]": "POSIXct",
        }

        switch_board = {
            "rda": {
                "csv": self.write_rdata_csv,
                "feather": self.write_rdata_feather,
                "parquet": self.write_rdata_parquet,
                "sqlite": self.write_rdata_sqlite,
            },
            "rdata": {
                "csv": self.write_rdata_csv,
                "feather": self.write_rdata_feather,
                "parquet": self.write_rdata_parquet,
                "sqlite": self.write_rdata_sqlite,
            },
            "rds": {
                "csv": self.write_rds_csv,
                "feather": self.write_rds_feather,
                "parquet": self.write_rds_parquet,
                "sqlite": self.write_rds_sqlite,
            },
        }

        switch_board[self.file_format][self.mode]()

        return None

    def write_rdata_csv(self) -> None:
        """
        Write R rda data via IO csv.

        This method will export one or more DataFrames into temp data
        and metadata csv files and call `read.csv` and `save` in R.
        """

        r_batch = """
                  args <- commandArgs(trailingOnly=TRUE)

                  py_names <- strsplit(args[1], ",")[[1]]

                  for(obj in py_names) {
                      meta <- paste0("meta_", obj, ".txt")
                      r_types <- strsplit(readLines(meta, n=-1,
                                          warn=FALSE), ",")[[1]]

                      data <- paste0("data_", obj, ".csv")
                      df <- tryCatch(
                                read.csv(data, colClasses=r_types),
                                error = function(e) read.csv(data)
                            )
                      assign(obj, df)
                      rm(df)
                  }

                  r_ascii <- as.logical(args[3])
                  r_compress <- ifelse(args[4] %in% c("True", "False"),
                                       as.logical(args[4]),
                                       args[4])

                  dfs <- names(Filter(is.data.frame, mget(ls())))
                  save(list=dfs, file=args[2],
                       ascii=r_ascii, compress=r_compress)
                  """

        with TemporaryDirectory() as tmp_dir:
            for nm, df in zip(self.rda_names, self.all_frames):

                data_file = os.path.join(tmp_dir, f"data_{nm}.csv")
                meta_file = os.path.join(tmp_dir, f"meta_{nm}.txt")
                r_code = os.path.join(tmp_dir, "rbatch.R")
                r_temp = os.path.join(tmp_dir, "rdata.rda")

                df = df.reset_index() if self.index else df
                df.to_csv(data_file, index=False)

                with open(meta_file, "w") as f:
                    f.write(
                        ",".join(
                            self.py_to_r_types[p]
                            for p in df.dtypes.astype(str).tolist()
                        )
                    )

            cmds = [
                "Rscript",
                r_code,
                ",".join(self.rda_names),
                r_temp,
                str(self.ascii),
                str(self.compress),
            ]
            self.run_rscript(tmp_dir, r_batch, cmds)

            self.disk_to_buffer(r_temp)

        return None

    def write_rdata_feather(self) -> None:
        """
        Write R rda data via IO feather.

        This method will export one or more DataFrames into temp
        feather files and call `read_feather` and `save` in R.
        """

        r_batch = """
                  suppressPackageStartupMessages(library(arrow))
                  args <- commandArgs(trailingOnly=TRUE)

                  py_names <- strsplit(args[1], ",")[[1]]

                  for(obj in py_names) {
                      data <- paste0("data_", obj, ".feather")
                      df <- arrow::read_feather(data)
                      assign(obj, df)
                      rm(df)
                  }

                  r_ascii <- as.logical(args[3])
                  r_compress <- ifelse(args[4] %in% c("True", "False"),
                                       as.logical(args[4]),
                                       args[4])

                  dfs <- names(Filter(is.data.frame, mget(ls())))
                  save(list=dfs, file=args[2],
                       ascii=r_ascii, compress=r_compress)
                  """

        with TemporaryDirectory() as tmp_dir:
            for nm, df in zip(self.rda_names, self.all_frames):

                data_file = os.path.join(tmp_dir, f"data_{nm}.feather")
                r_code = os.path.join(tmp_dir, "rbatch.R")
                r_temp = os.path.join(tmp_dir, "rdata.rda")

                df = df.reset_index() if self.index else df.reset_index(drop=True)
                df.to_feather(data_file)

            cmds = [
                "Rscript",
                r_code,
                ",".join(self.rda_names),
                r_temp,
                str(self.ascii),
                str(self.compress),
            ]
            self.run_rscript(tmp_dir, r_batch, cmds)

            self.disk_to_buffer(r_temp)

    def write_rdata_parquet(self) -> None:
        """
        Write R rda data via IO parquet.

        This method will export one or more DataFrames into temp
        parquet files and call `read_parquet` and `save` in R.
        """

        r_batch = """
                  suppressPackageStartupMessages(library(arrow))
                  args <- commandArgs(trailingOnly=TRUE)

                  py_names <- strsplit(args[1], ",")[[1]]

                  for(obj in py_names) {
                      data <- paste0("data_", obj, ".parquet")
                      df <- arrow::read_parquet(data)
                      assign(obj, df)
                      rm(df)
                  }

                  r_ascii <- as.logical(args[3])
                  r_compress <- ifelse(args[4] %in% c("True", "False"),
                                       as.logical(args[4]),
                                       args[4])

                  dfs <- names(Filter(is.data.frame, mget(ls())))
                  save(list=dfs, file=args[2],
                       ascii=r_ascii, compress=r_compress)
                  """

        with TemporaryDirectory() as tmp_dir:
            for nm, df in zip(self.rda_names, self.all_frames):

                data_file = os.path.join(tmp_dir, f"data_{nm}.parquet")
                r_code = os.path.join(tmp_dir, "rbatch.R")
                r_temp = os.path.join(tmp_dir, "rdata.rda")

                df = df.reset_index() if self.index else df
                df.to_parquet(data_file, index=False)

            cmds = [
                "Rscript",
                r_code,
                ",".join(self.rda_names),
                r_temp,
                str(self.ascii),
                str(self.compress),
            ]
            self.run_rscript(tmp_dir, r_batch, cmds)

            self.disk_to_buffer(r_temp)

    def write_rdata_sqlite(self) -> None:
        """
        Write R rda data via IO sql.

        This method will export one or more DataFrames into a temp
        SQLite database and call `dbReadTable` and `save` in R.
        """
        import sqlite3

        r_batch = """
                  suppressPackageStartupMessages(library(RSQLite))
                  args <- commandArgs(trailingOnly=TRUE)

                  conn <- dbConnect(RSQLite::SQLite(), args[1])
                  py_names <- strsplit(args[2], ",")[[1]]

                  for(obj in py_names) {
                      data <- paste0("data_", obj)
                      df <- dbReadTable(conn, data)
                      assign(obj, df)
                      rm(df)
                  }
                  dbDisconnect(conn)

                  r_ascii <- as.logical(args[4])
                  r_compress <- ifelse(args[5] %in% c("True", "False"),
                                       as.logical(args[5]),
                                       args[5])

                  dfs <- names(Filter(is.data.frame, mget(ls())))
                  save(list=dfs, file=args[3],
                       ascii=r_ascii, compress=r_compress)
                  """

        with TemporaryDirectory() as tmp_dir:
            r_db = os.path.join(tmp_dir, "rdata.db")
            conn = sqlite3.connect(r_db)

            for nm, df in zip(self.rda_names, self.all_frames):
                r_code = os.path.join(tmp_dir, "rbatch.R")
                r_temp = os.path.join(tmp_dir, "rdata.rda")

                df = df.reset_index() if self.index else df
                df.to_sql(f"data_{nm}", conn, index=False)

            conn.close()
            cmds = [
                "Rscript",
                r_code,
                r_db,
                ",".join(self.rda_names),
                r_temp,
                str(self.ascii),
                str(self.compress),
            ]
            self.run_rscript(tmp_dir, r_batch, cmds)

            self.disk_to_buffer(r_temp)

    def write_rds_csv(self) -> None:
        """
        Write R rds data via IO csv.

        This method will export a single DataFrame into temp csv
        data and call `read.csv` and `saveRDS` in R.
        """

        r_batch = """
                  args <- commandArgs(trailingOnly=TRUE)
                  py_data <- args[1]
                  r_types <- strsplit(args[2], ",")[[1]]

                  df <- tryCatch(
                           read.csv(py_data, colClasses=r_types),
                           error = function(e) read.csv(py_data)
                  )

                  r_ascii <- as.logical(args[4])
                  r_compress <- ifelse(args[5] %in% c("True", "False"),
                                       as.logical(args[5]),
                                       args[5])

                  saveRDS(df, file=args[3],
                          ascii=r_ascii, compress=r_compress)
                  """

        with TemporaryDirectory() as tmp_dir:
            r_code = os.path.join(tmp_dir, "rbatch.R")
            py_data = os.path.join(tmp_dir, "pydata.csv")
            r_temp = os.path.join(tmp_dir, "rdata.rds")

            py_df = self.frame.reset_index() if self.index else self.frame
            r_types = ",".join(py_df.dtypes.astype(str).replace(self.py_to_r_types))

            py_df.to_csv(py_data, index=False)

            cmds = [
                "Rscript",
                r_code,
                py_data,
                r_types,
                r_temp,
                str(self.ascii),
                str(self.compress),
            ]
            self.run_rscript(tmp_dir, r_batch, cmds)

            self.disk_to_buffer(r_temp)

        return None

    def write_rds_feather(self) -> None:
        """
        Write R rds data via IO feather.

        This method will export a single DataFrame into a temp
        feather file to call `read_feather` and `saveRDS` in R.
        """

        r_batch = """
                  suppressPackageStartupMessages(library(arrow))
                  args <- commandArgs(trailingOnly=TRUE)

                  df <- arrow::read_feather(args[1])

                  r_ascii <- as.logical(args[3])
                  r_compress <- ifelse(args[4] %in% c("True", "False"),
                                       as.logical(args[4]),
                                       args[4])

                  saveRDS(df, file=args[2],
                          ascii=r_ascii, compress=r_compress)
                  """

        with TemporaryDirectory() as tmp_dir:
            r_code = os.path.join(tmp_dir, "rbatch.R")
            py_data = os.path.join(tmp_dir, "pydata.feather")
            r_temp = os.path.join(tmp_dir, "rdata.rds")

            py_df = (
                self.frame.reset_index()
                if self.index
                else self.frame.reset_index(drop=True)
            )

            py_df.to_feather(py_data)

            cmds = [
                "Rscript",
                r_code,
                py_data,
                r_temp,
                str(self.ascii),
                str(self.compress),
            ]
            self.run_rscript(tmp_dir, r_batch, cmds)

            self.disk_to_buffer(r_temp)

    def write_rds_parquet(self) -> None:
        """
        Write R rds data via IO parquet.

        This method will export a single DataFrame into a temp
        parquet file for `read_parquet` and `saveRDS` in R.
        """

        r_batch = """
                  suppressPackageStartupMessages(library(arrow))
                  args <- commandArgs(trailingOnly=TRUE)

                  df <- arrow::read_parquet(args[1])

                  r_ascii <- as.logical(args[3])
                  r_compress <- ifelse(args[4] %in% c("True", "False"),
                                       as.logical(args[4]),
                                       args[4])

                  saveRDS(df, file=args[2],
                          ascii=r_ascii, compress=r_compress)
                  """

        with TemporaryDirectory() as tmp_dir:
            r_code = os.path.join(tmp_dir, "rbatch.R")
            py_data = os.path.join(tmp_dir, "pydata.parquet")
            r_temp = os.path.join(tmp_dir, "rdata.rds")

            py_df = self.frame.reset_index() if self.index else self.frame

            py_df.to_parquet(py_data, index=False)

            cmds = [
                "Rscript",
                r_code,
                py_data,
                r_temp,
                str(self.ascii),
                str(self.compress),
            ]
            self.run_rscript(tmp_dir, r_batch, cmds)

            self.disk_to_buffer(r_temp)

    def write_rds_sqlite(self) -> None:
        """
        Write R rds data via IO sql.

        This method will export a single DataFrame into a temp
        parquet file for `dbReadTable` and `saveRDS` in R.
        """
        import sqlite3

        r_batch = """
                  suppressPackageStartupMessages(library(RSQLite))
                  args <- commandArgs(trailingOnly=TRUE)

                  conn <- dbConnect(RSQLite::SQLite(), args[1])
                  df <- dbReadTable(conn, "pydata")

                  r_ascii <- as.logical(args[3])
                  r_compress <- ifelse(args[4] %in% c("True", "False"),
                                       as.logical(args[4]),
                                       args[4])

                  saveRDS(df, file=args[2],
                          ascii=r_ascii, compress=r_compress)
                  dbDisconnect(conn)
                  """

        with TemporaryDirectory() as tmp_dir:
            r_code = os.path.join(tmp_dir, "rbatch.R")
            py_data = os.path.join(tmp_dir, "pydata.db")
            r_temp = os.path.join(tmp_dir, "rdata.rds")

            py_df = self.frame.reset_index() if self.index else self.frame

            conn = sqlite3.connect(py_data)
            py_df.to_sql("pydata", conn, index=False)
            conn.close()

            cmds = [
                "Rscript",
                r_code,
                py_data,
                r_temp,
                str(self.ascii),
                str(self.compress),
            ]
            self.run_rscript(tmp_dir, r_batch, cmds)

            self.disk_to_buffer(r_temp)
