"""Accessors for arrow-backed data."""

from __future__ import annotations

from typing import TYPE_CHECKING

from pandas.compat import pa_version_under7p0

if not pa_version_under7p0:
    import pyarrow as pa
    import pyarrow.compute as pc

    from pandas.core.dtypes.dtypes import ArrowDtype

if TYPE_CHECKING:
    from pandas import (
        DataFrame,
        Series,
    )


class StructAccessor:
    """
    Accessor object for structured data properties of the Series values.

    Parameters
    ----------
    data : Series
        Series containing Arrow struct data.
    """

    _validation_msg = (
        "Can only use the '.struct' accessor with 'struct[pyarrow]' dtype, not {dtype}."
    )

    def __init__(self, data=None) -> None:
        self._parent = data
        self._validate(data)

    def _validate(self, data):
        dtype = data.dtype
        if not isinstance(dtype, ArrowDtype):
            # Raise AttributeError so that inspect can handle non-struct Series.
            raise AttributeError(self._validation_msg.format(dtype=dtype))

        if not pa.types.is_struct(dtype.pyarrow_dtype):
            # Raise AttributeError so that inspect can handle non-struct Series.
            raise AttributeError(self._validation_msg.format(dtype=dtype))

    @property
    def dtypes(self) -> Series:
        """
        Return the dtype object of each child field of the struct.

        Returns
        -------
        pandas.Series
            The data type of each child field.

        Examples
        --------
        >>> import pyarrow as pa
        >>> s = pd.Series(
        ...     [
        ...         {"version": 1, "project": "pandas"},
        ...         {"version": 2, "project": "pandas"},
        ...         {"version": 1, "project": "numpy"},
        ...     ],
        ...     dtype=pd.ArrowDtype(pa.struct(
        ...         [("version", pa.int64()), ("project", pa.string())]
        ...     ))
        ... )
        >>> s.struct.dtypes
        version     int64[pyarrow]
        project    string[pyarrow]
        dtype: object
        """
        from pandas import (
            Index,
            Series,
        )

        pa_type = self._parent.dtype.pyarrow_dtype
        types = [ArrowDtype(struct.type) for struct in pa_type]
        names = [struct.name for struct in pa_type]
        return Series(types, index=Index(names))

    def field(self, name_or_index: str | int) -> Series:
        """
        Extract a child field of a struct as a Series.

        Parameters
        ----------
        name_or_index : str | int
            Name or index of the child field to extract.

        Returns
        -------
        pandas.Series
            The data corresponding to the selected child field.

        See Also
        --------
        Series.struct.explode : Return all child fields as a DataFrame.

        Examples
        --------
        >>> import pyarrow as pa
        >>> s = pd.Series(
        ...     [
        ...         {"version": 1, "project": "pandas"},
        ...         {"version": 2, "project": "pandas"},
        ...         {"version": 1, "project": "numpy"},
        ...     ],
        ...     dtype=pd.ArrowDtype(pa.struct(
        ...         [("version", pa.int64()), ("project", pa.string())]
        ...     ))
        ... )

        Extract by field name.

        >>> s.struct.field("project")
        0    pandas
        1    pandas
        2     numpy
        Name: project, dtype: string[pyarrow]

        Extract by field index.

        >>> s.struct.field(0)
        0    1
        1    2
        2    1
        Name: version, dtype: int64[pyarrow]
        """
        from pandas import Series

        pa_arr = self._parent.array._pa_array
        if isinstance(name_or_index, int):
            index = name_or_index
        elif isinstance(name_or_index, str):
            index = pa_arr.type.get_field_index(name_or_index)
        else:
            raise ValueError(
                "name_or_index must be an int or str, "
                f"got {type(name_or_index).__name__}"
            )

        pa_field = pa_arr.type[index]
        field_arr = pc.struct_field(pa_arr, [index])
        return Series(
            field_arr,
            dtype=ArrowDtype(field_arr.type),
            index=self._parent.index,
            name=pa_field.name,
        )

    def explode(self) -> DataFrame:
        """
        Extract all child fields of a struct as a DataFrame.

        Returns
        -------
        pandas.DataFrame
            The data corresponding to all child fields.

        See Also
        --------
        Series.struct.field : Return a single child field as a Series.

        Examples
        --------
        >>> import pyarrow as pa
        >>> s = pd.Series(
        ...     [
        ...         {"version": 1, "project": "pandas"},
        ...         {"version": 2, "project": "pandas"},
        ...         {"version": 1, "project": "numpy"},
        ...     ],
        ...     dtype=pd.ArrowDtype(pa.struct(
        ...         [("version", pa.int64()), ("project", pa.string())]
        ...     ))
        ... )

        >>> s.struct.explode()
           version project
        0        1  pandas
        1        2  pandas
        2        1   numpy
        """
        from pandas import concat

        pa_type = self._parent.dtype.pyarrow_dtype
        return concat(
            [self.field(i) for i in range(pa_type.num_fields)], axis="columns"
        )
