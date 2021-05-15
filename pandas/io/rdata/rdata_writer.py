"""
write R data files (RData, rda, rds).

This IO module interfaces with the librdata C library by Evan Miller:
  https://github.com/WizardMac/librdata
"""
from __future__ import annotations

import os
from tempfile import TemporaryDirectory

from pandas._typing import (
    CompressionOptions,
    FilePathOrBuffer,
    StorageOptions,
)

from pandas.core.frame import DataFrame

from pandas.io.common import get_handle
from pandas.io.rdata._rdata import LibrdataWriter


class RDataWriter:
    """
    Subclass to write pandas DataFrames into R data files.

    Parameters
    ----------
    path_or_buffer : a valid str, path object or file-like object
        Any valid string path is acceptable.

    file_format : {{'infer', 'rdata', 'rda', 'rds'}}, default 'infer'
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
                f"{self.file_format} is not a valid value for file_format."
            )

        if (
            self.file_format == "infer"
            and isinstance(self.path_or_buffer, str)
            and path_ext not in ["rdata", "rda", "rds"]
        ):
            raise ValueError(
                f"Unable to infer file format from file name: {self.path_or_buffer}"
                "Please use known R data type (rdata, rda, rds)."
            )

        if self.file_format == "infer" and isinstance(path_ext, str):
            self.file_format = path_ext

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

        Converts non-primitive and non-datetimes to object to align to R
        atomic types, then exports dictionaries of each column with meta data.
        """

        self.frame = (
            self.frame.reset_index()
            if self.index
            else self.frame.reset_index(drop=True)
        )

        excl_types = ["bool", "number", "object", "datetime", "datetimetz", "timedelta"]
        for col in self.frame.select_dtypes(exclude=excl_types).columns:
            self.frame[col] = self.frame[col].astype(str)

        for col in self.frame.select_dtypes(include=["datetimetz"]).columns:
            self.frame[col] = self.frame[col].dt.tz_localize(None)

        for col in self.frame.select_dtypes(include=["timedelta"]).columns:
            self.frame[col] = self.frame[col].dt.total_seconds()

        rdict = {"dtypes": {k: str(v) for k, v in self.frame.dtypes.to_dict().items()}}

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

        for col in self.frame.select_dtypes(include=["datetime"]).columns:
            self.frame[col] = self.frame[col].values.view("int64") / (10 ** 9)

        rdict["data"] = self.frame.to_dict()

        lbw = LibrdataWriter()

        with TemporaryDirectory() as tmp_dir:
            r_temp = os.path.join(tmp_dir, "rdata.rda")
            lbw.write_rdata(
                rfile=r_temp,
                rdict=rdict,
                rformat=self.file_format,
                tbl_name=self.rda_name,
            )

            self.disk_to_buffer(r_temp)

        return None
