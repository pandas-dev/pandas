"""
Read R data files (RData, rda, rds).

This IO module interfaces with the librdata C library by Evan Miller:
  https://github.com/WizardMac/librdata
"""
from __future__ import annotations

import io
import os
from tempfile import TemporaryDirectory

from pandas._typing import (
    Buffer,
    CompressionOptions,
    FilePathOrBuffer,
    StorageOptions,
)
from pandas.util._decorators import doc

from pandas.core.dtypes.common import is_list_like

from pandas.core.api import to_datetime
from pandas.core.arrays import Categorical
from pandas.core.frame import (
    DataFrame,
    Index,
    Series,
)
from pandas.core.shared_docs import _shared_docs

from pandas.io.common import (
    file_exists,
    get_handle,
    is_fsspec_url,
    is_url,
    stringify_path,
)
from pandas.io.rdata._rdata import LibrdataReader


@doc(storage_options=_shared_docs["storage_options"])
def read_rdata(
    path_or_buffer: FilePathOrBuffer,
    file_format: str = "infer",
    select_frames: list[str] | None = None,
    rownames: bool = True,
    compression: CompressionOptions = "gzip",
    storage_options: StorageOptions = None,
) -> dict[str, DataFrame]:
    r"""
    Read R data (.RData, .rda, .rds) into DataFrame or ``dict`` of DataFrames.

    .. versionadded:: 1.3.0

    Parameters
    ----------
    path_or_buffer : str, path object, or file-like object
        Any valid file path is acceptable. The string could be a URL.
        Valid URL schemes include http, ftp, s3, and file.

    file_format : {{'infer', 'rdata', 'rda', 'rds'}}, default 'infer'
        R serialization type as output from R's base::save or base::saveRDS
        commands. Default 'infer' will use extension in file name to
        to determine the format type.

    select_frames : list, default returns all DataFrames
        Selected names of DataFrames to return from R RData and rdata types that
        can contain multiple objects.

    rownames : bool, default True
        Include original rownames in R data frames to map into a DataFrame index.

    compression : {{'infer', 'gzip', 'bz2', 'zip', 'xz', None}}, default 'gzip'
        For on-the-fly decompression of on-disk data. If 'infer', then use
        gzip, bz2, zip or xz if path_or_buffer is a string ending in
        '.gz', '.bz2', '.zip', or 'xz', respectively, and no decompression
        otherwise. If using 'zip', the ZIP file must contain only one data
        file to be read in. Set to None for no decompression. This method will
        default to 'gzip' since 'gzip2` is the default compression in R for
        RData and rds types.

    {storage_options}

    Returns
    -------
    Dict of DataFrames
        Depends on R data type where rds formats returns a ``dict`` of a single
        DataFrame and RData or rda formats can return ``dict`` of one or more
        DataFrames.

    See Also
    --------
    read_sas : Read SAS datasets into DataFrame.
    read_stata : Read Stata datasets into DataFrame.
    read_spss : Read SPSS datasets into DataFrame.

    Notes
    -----
    Any R data file that contains a non-data.frame object may raise parsing errors.
    Method will return data.frame and data.frame like objects such as tibbles and
    data.tables. For more information of R serialization data types, see docs on
    `rds`<https://www.rdocumentation.org/packages/base/versions/3.6.2/topics/readRDS>__
    and `rda`<https://www.rdocumentation.org/packages/base/versions/3.6.2/topics/save>__
    formats.

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
    {{'r_dataframe':
                          gas  year  emissions
    rownames
    1          Carbon dioxide  2018    5424.88
    2                 Methane  2018     634.46
    3           Nitrous oxide  2018     434.53
    4       Fluorinated gases  2018     182.79
    5                   Total  2018    6676.65}}

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

    rdr = _RDataReader(
        path_or_buffer,
        file_format,
        select_frames,
        rownames,
        compression,
        storage_options,
    )

    r_dfs = rdr.parse_data()

    return r_dfs


def get_data_from_filepath(
    filepath_or_buffer,
    encoding,
    compression,
    storage_options,
) -> str | bytes | Buffer:
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


def preprocess_data(data) -> io.StringIO | io.BytesIO:
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

    file_format : {{'infer', 'rdata', 'rda', 'rds'}}, default 'infer'
        R serialization type.

    select_frames : list, default None
        Selected names of DataFrames to return from R data.

    rownames : bool, default True
        Include original rownames in R data frames.

    compression : {'infer', 'gzip', 'bz2', 'zip', 'xz', None}, default 'infer'
        Compression type for on-the-fly decompression of on-disk data.
        If 'infer', then use extension for gzip, bz2, zip or xz.

    storage_options : dict, optional
        Extra options that make sense for a particular storage connection,
        e.g. host, port, username, password, etc.
    """

    def __init__(
        self,
        path_or_buffer,
        file_format,
        select_frames,
        rownames,
        compression,
        storage_options,
    ) -> None:
        self.path_or_buffer = path_or_buffer
        self.file_format = file_format.lower()
        self.select_frames = select_frames
        self.rownames = rownames
        self.compression = compression
        self.storage_options = storage_options
        self.verify_params()

    def verify_params(self) -> None:
        """
        Verify user entries of parameters.

        This method will check the values and types of select parameters
        and raise appropriate errors.
        """

        path_ext: str | None = (
            os.path.splitext(self.path_or_buffer.lower())[1][1:]
            if isinstance(self.path_or_buffer, str)
            else None
        )

        if self.file_format not in ["infer", "rdata", "rda", "rds"]:
            raise ValueError(
                f"'{self.file_format}' is not a valid value for file_format"
            )

        if (
            self.file_format == "infer"
            and isinstance(self.path_or_buffer, str)
            and path_ext not in ["rdata", "rda", "rds"]
        ) or (self.file_format == "infer" and not isinstance(self.path_or_buffer, str)):
            raise ValueError(
                f"Unable to infer file format from file name: {self.path_or_buffer}. "
                "Please use known R data type (rdata, rda, rds)."
            )

        if self.file_format == "infer" and isinstance(path_ext, str):
            self.file_format = path_ext

        if self.select_frames is not None and not is_list_like(self.select_frames):
            raise TypeError(
                f"{type(self.select_frames).__name__} is "
                "not a valid type for select_frames"
            )

    def buffer_to_disk(self, tmp_dir: str) -> str:
        """
        Convert path or buffer to disk file.

        This method will convert path_or_buffer to temp file
        to parse RData from disk.
        """

        r_temp = os.path.join(tmp_dir, "rdata.rda")

        handle_data = get_data_from_filepath(
            filepath_or_buffer=self.path_or_buffer,
            encoding="utf-8",
            compression=self.compression,
            storage_options=self.storage_options,
        )

        with preprocess_data(handle_data) as r_data:
            if isinstance(r_data, io.BytesIO):
                with open(r_temp, "wb") as f:
                    f.write(r_data.read())

        return r_temp

    def build_frame(self, data_dict: dict) -> DataFrame:
        """
        Builds DataFrame from raw, nested parsed RData dict.

        Converts special class variables (bools, factors, dates, datetimes),
        then binds all columns together with DataFrame constructor.
        """

        final_dict = {
            k: Series(v)
            for k, v in data_dict["data"].items()
            if k not in ["dtypes", "colnames", "rownames"]
        }

        rdf = DataFrame(data=final_dict)

        for col, dtype in data_dict["dtypes"].items():
            if dtype == "bool":
                rdf[col] = rdf[col].astype(bool)

            if dtype == "factor":
                rdf[col] = Categorical(rdf[col])

            if dtype == "date":
                rdf[col] = to_datetime(rdf[col], unit="d")

            if dtype == "datetime":
                rdf[col] = to_datetime(rdf[col], unit="s")

        colnames = (
            None
            if data_dict["colnames"] is None
            else list(data_dict["colnames"].values())
        )
        if colnames is not None:
            rdf.columns = Index(colnames)

        rownames = (
            None
            if data_dict["rownames"] is None
            else list(data_dict["rownames"].values())
        )
        if self.rownames:
            if rownames is not None:
                rdf.index = Index(rownames)
            else:
                rdf.index += 1
            rdf.index.name = "rownames"

        return rdf

    def parse_data(self) -> dict[str, DataFrame]:
        """
        Parse R data files into DataFrames

        This method will retrieve dictionary of R data and build
        DataFrame for each item in data file
        """

        lbr = LibrdataReader()

        with TemporaryDirectory() as tmp_dir:
            r_temp = self.buffer_to_disk(tmp_dir)
            rdict = lbr.read_rdata(r_temp)

        r_dfs = {k: self.build_frame(v) for k, v in rdict.items()}

        if self.select_frames:
            r_dfs = {k: v for k, v in r_dfs.items() if k in self.select_frames}

        return r_dfs
