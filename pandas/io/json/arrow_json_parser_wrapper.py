from __future__ import annotations

from typing import TYPE_CHECKING

from pandas._typing import ReadBuffer
from pandas.compat._optional import import_optional_dependency

if TYPE_CHECKING:
    from pandas import DataFrame


class ArrowJsonParserWrapper:
    """
    Wrapper for the pyarrow engine for read_json()
    """

    def __init__(self, src: ReadBuffer[bytes]) -> None:
        self.src = src

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
