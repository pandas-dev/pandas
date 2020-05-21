from pandas.core.dtypes.base import ExtensionDtype
from typing import Type
from pandas._libs import missing as libmissing
from pandas.core.arrays.boolean import BooleanDtype

import numpy as np


class DateDtype(ExtensionDtype):
    """
    An ExtensionDtype to hold a single date.

    The attributes name & type are set when subclasses are created.
    """

    base = None
    type: Type = np.ndarray
    na_value = np.nan

    def __init__(self, day: int, month: int, year: int):

        assert month in range(13)

        self.date_array = np.array(
            [np.int8(date_piece) for date_piece in [day, month, year]]
        )
        print(type(self.date_array[0]))

    @property
    def name(self) -> str:
        """
        The alias for DateDtype is ``'string'``.
        """
        return "date"

    @property
    def na_value(self) -> "Scalar":
        """
        BooleanDtype uses :attr:`pandas.NA` as the missing NA value.

        .. warning::

           `na_value` may change in a future release.
        """
        return libmissing.NA

    def __repr__(self):
        return "DateDtype"

    # TODO Construct from string with pattern


if __name__ == "__main__":
    date = DateDtype(2019, 10, 21)
    print(date)
    print(date.date_array)
    print(date.type)
    print(date.na_value)
    # test = BooleanDtype.construct_from_string("boolean")
