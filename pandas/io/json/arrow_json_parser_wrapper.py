from __future__ import annotations

from typing import TYPE_CHECKING

from pandas._typing import ReadBuffer
from pandas.compat._optional import import_optional_dependency

from pandas.core.dtypes.inference import is_integer

if TYPE_CHECKING:
    from pandas import DataFrame


class ArrowJsonParserWrapper:
    """
    Wrapper for the pyarrow engine for read_json()
    """

    def __init__(self, src: ReadBuffer[bytes]) -> None:
        super().__init__()
        self.src = src

    def _parse_kwd(self) -> None:
        """
        Validates keywords before passing to pyarrow
        """
        ...

    def _get_pyarrow_options(self) -> None:
        ...

    def read(self) -> DataFrame:
        """
        Reads the contents of a JSON file into a DataFrame and
        processes it according to the kwargs passed in the
        constructor.

        Returns
        -------
        DataFrame
            The DataFrame created from the JSON file.
        """
        pyarrow_json = import_optional_dependency("pyarrow.json")
        table = pyarrow_json.read_json(self.src)

        frame = table.to_pandas()
        return frame

    def _finalize_output(self, frame: DataFrame) -> DataFrame:
        """
        Processes data read in based on kwargs.

        Parameters
        ----------
        frame: DataFrame
                The DataFrame to process.

        Returns
        -------
        DataFrame
                The processed DataFrame.
        """
        num_cols = len(frame.columns)
        multi_index_named = True
        if self.header is None:
            if self.names is None:
                if self.prefix is not None:
                    self.names = [f"{self.prefix}{i}" for i in range(num_cols)]
                elif self.header is None:
                    self.names = range(num_cols)
            if len(self.names) != num_cols:
                # usecols is passed through to pyarrow, we only handle index col here
                # The only way self.names is not the same length as number of cols is
                # if we have int index_col. We should just pad the names(they will get
                # removed anyways) to expected length then.
                self.names = list(range(num_cols - len(self.names))) + self.names
                multi_index_named = False
            frame.columns = self.names
        # we only need the frame not the names
        frame.columns, frame = self._do_date_conversions(frame.columns, frame)
        if self.index_col is not None:
            for i, item in enumerate(self.index_col):
                if is_integer(item):
                    self.index_col[i] = frame.columns[item]
                else:
                    # String case
                    if item not in frame.columns:
                        raise ValueError(f"Index {item} invalid")
            frame.set_index(self.index_col, drop=True, inplace=True)
            # Clear names if headerless and no name given
            if self.header is None and not multi_index_named:
                frame.index.names = [None] * len(frame.index.names)

        if self.kwds.get("dtype") is not None:
            try:
                frame = frame.astype(self.kwds.get("dtype"))
            except TypeError as e:
                # GH#44901 reraise to keep api consistent
                raise ValueError(e)
        return frame
