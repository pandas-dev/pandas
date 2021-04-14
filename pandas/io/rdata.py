from datetime import datetime
import io
import os
from tempfile import TemporaryDirectory
from typing import (
    Dict,
    List,
    Optional,
    Union,
)

from pandas._typing import (
    Buffer,
    CompressionOptions,
    FilePathOrBuffer,
    StorageOptions,
)
from pandas.compat._optional import import_optional_dependency
from pandas.errors import AbstractMethodError
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


@doc(storage_options=_shared_docs["storage_options"])
def read_rdata(
    path_or_buffer: FilePathOrBuffer,
    file_format: str = "infer",
    select_frames: Optional[List[str]] = None,
    rownames: bool = True,
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

    select_frames : list, default None
        Selected names of DataFrames to return from R rda and RData types that
        can contain multiple objects.

    rownames : bool, default True
        Include original rownames in R data frames to map into a DataFrame index.

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
    Any R data file that contains a non-data.frame object may raise parsing errors.
    Method will return data.frame, matrix, and data.frame like object such as
    tibbles and data.tables.

    For ``pyreadr`` engine, ``select_frames`` above is synonymous to ``use_objects``
    in package's `read_r` method. Also, ``timezone`` argument defaults to current
    system regional timezone in order to correspond to original date/times in R.

    Examples
    --------
    For an .rds file which only contains a single R object, method returns a
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

    For an .RData or .rda file which can contain multiple R objects, method
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

    import_optional_dependency("pyreadr")

    rdr = _PyReadrParser(
        path_or_buffer,
        file_format,
        select_frames,
        rownames,
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

    select_frames : list, default None
        Selected names of DataFrames to return from R data.

    rownames : bool, default True
        Include original rownames in R data frames.

    storage_options : dict, optional
        Extra options that make sense for a particular storage connection,
        e.g. host, port, username, password, etc.,

    See also
    --------
    pandas.io.rdata._PyReadrParser

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
        select_frames,
        rownames,
        storage_options,
    ) -> None:
        self.path_or_buffer = path_or_buffer
        self.file_format = file_format.lower()
        self.select_frames = select_frames
        self.rownames = rownames
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

        if self.select_frames is not None and not is_list_like(self.select_frames):
            raise TypeError(
                f"{type(self.select_frames).__name__} is "
                "not a valid type for select_frames"
            )

    def buffer_to_disk(self, tmp_dir: str) -> str:
        """
        Convert path or buffer to disk file.

        This method will convert path_or_buffer to temp file
        for pyreadr to parse from disk.
        """

        r_temp = os.path.join(tmp_dir, "rdata.rda")

        handle_data = _get_data_from_filepath(
            filepath_or_buffer=self.path_or_buffer,
            encoding="utf-8",
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


class RDataWriter:
    """
    Subclass to write pandas DataFrames into R data files.

    Parameters
    ----------
    path_or_buffer : a valid str, path object or file-like object
        Any valid string path is acceptable.

    file_format : {{'infer', 'rda', 'rdata', 'rds'}}, default 'infer'
        R serialization type.

    rda_name : str, default "pandas_dataframe"
        Name for exported DataFrame in rda file.

    index : bool, default True
        Include index or MultiIndex in output as separate columns.

    compression : {'gzip', 'bz2', 'xz', None}, default 'gzip'
        Compression type for on-the-fly decompression of on-disk data.

    storage_options : dict, optional
        Extra options that make sense for a particular storage connection,
        e.g. host, port, username, password, etc.

    See also
    --------
    pandas.io.rdata.PyReadrWriter

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
        rda_name: str = "pandas_dataframe",
        index: bool = True,
        compression: CompressionOptions = "gzip",
        storage_options: StorageOptions = None,
    ) -> None:
        self.frame = frame
        self.path_or_buffer = path_or_buffer
        self.file_format = file_format.lower()
        self.rda_name = rda_name
        self.index = index
        self.compression = compression
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

        if self.compression is not None and self.compression not in [
            "gzip",
            "bz2",
            "xz",
        ]:
            raise ValueError(
                f"{self.compression} is not a supported value for compression."
            )

    def disk_to_buffer(self, r_file: str) -> None:
        """
        Save temp file to path or buffer.

        This method will convert written R data to path_or_buffer.
        """

        with open(r_file, "rb") as rdata:
            with get_handle(
                self.path_or_buffer,
                "wb",
                compression=self.compression,
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
                    df_name=self.rda_name,
                    compress=None,
                )
            elif self.file_format == "rds":
                write_rds(
                    path=r_temp,
                    df=self.frame,
                    compress=None,
                )

            self.disk_to_buffer(r_temp)

        return None
