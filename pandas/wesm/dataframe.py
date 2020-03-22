# MIT License
#
# Copyright (c) 2020 Wes McKinney

from abc import ABC, abstractmethod
from collections import abc
from typing import Any, Hashable, Iterable, Optional, Sequence

# ----------------------------------------------------------------------
# A simple data type class hierarchy for illustration


class DataType(ABC):
    """
    A metadata object representing the logical value type of a cell in a data
    frame column. This metadata does not guarantee an specific underlying data
    representation
    """

    def __eq__(self, other: "DataType"):  # type: ignore
        return self.equals(other)

    def __str__(self):
        return self.to_string()

    def __repr__(self):
        return str(self)

    @abstractmethod
    def to_string(self) -> str:
        """
        Return human-readable representation of the data type
        """

    @abstractmethod
    def equals(self, other: "DataType") -> bool:
        """
        Return true if other DataType contains the same metadata as this
        DataType
        """
        pass


class PrimitiveType(DataType):
    def equals(self, other: DataType) -> bool:
        return type(self) == type(other)


class NullType(PrimitiveType):
    """
    A data type whose values are always null
    """

    def to_string(self):
        return "null"


class Boolean(PrimitiveType):
    def to_string(self):
        return "bool"


class NumberType(PrimitiveType):
    pass


class IntegerType(NumberType):
    pass


class SignedIntegerType(IntegerType):
    pass


class Int8(SignedIntegerType):
    def to_string(self):
        return "int8"


class Int16(SignedIntegerType):
    def to_string(self):
        return "int16"


class Int32(SignedIntegerType):
    def to_string(self):
        return "int32"


class Int64(SignedIntegerType):
    def to_string(self):
        return "int64"


class Binary(PrimitiveType):
    """
    A variable-size binary (bytes) value
    """

    def to_string(self):
        return "binary"


class String(PrimitiveType):
    """
    A UTF8-encoded string value
    """

    def to_string(self):
        return "string"


class Object(PrimitiveType):
    """
    Any PyObject value
    """

    def to_string(self):
        return "object"


class Categorical(DataType):
    """
    A categorical value is an ordinal (integer) value that references a
    sequence of category values of an arbitrary data type
    """

    def __init__(
        self, index_type: IntegerType, category_type: DataType, ordered: bool = False
    ):
        self.index_type = index_type
        self.category_type = category_type
        self.ordered = ordered

    def equals(self, other: DataType) -> bool:
        return (
            isinstance(other, Categorical)
            and self.index_type == other.index_type
            and self.category_type == other.category_type
            and self.ordered == other.ordered
        )

    def to_string(self):
        return "categorical(indices={}, categories={}, ordered={})".format(
            str(self.index_type), str(self.category_type), self.ordered
        )


# ----------------------------------------------------------------------
# Classes representing a column in a DataFrame


class Column(ABC):
    @property
    @abstractmethod
    def name(self) -> Optional[Hashable]:
        pass

    @property
    @abstractmethod
    def type(self) -> DataType:
        """
        Return the logical type of each column cell value
        """
        pass

    def to_numpy(self):
        """
        Access column's data as a NumPy array. Recommended to return a view if
        able but not required
        """
        raise NotImplementedError("Conversion to NumPy not available")

    def to_arrow(self, **kwargs):
        """
        Access column's data in the Apache Arrow format as pyarrow.Array or
        ChunkedArray. Recommended to return a view if able but not required
        """
        raise NotImplementedError("Conversion to Arrow not available")


# ----------------------------------------------------------------------
# DataFrame: the main public API


class DataFrame(ABC, abc.Mapping):
    """
    An abstract data frame base class.

    A "data frame" represents an ordered collection of named columns. A
    column's "name" is permitted to be any hashable Python value, but strings
    are common. Names are not required to be unique. Columns may be accessed by
    name (when the name is unique) or by position.
    """

    @property
    def __dataframe__(self):
        """
        Idempotence of data frame protocol
        """
        return self

    def __iter__(self):
        # TBD: Decide what iterating should return
        return iter(self.column_names)

    def __len__(self):
        return self.num_rows

    @property
    @abstractmethod
    def num_columns(self) -> int:
        """
        Return the number of columns in the DataFrame
        """
        pass

    @property
    @abstractmethod
    def num_rows(self) -> Optional[int]:
        """
        Return the number of rows in the DataFrame (if known)
        """
        pass

    @abstractmethod
    def iter_column_names(self) -> Iterable[Any]:
        """
        Return the column names as an iterable
        """
        pass

    # TODO: Should this be a method or property?
    @property
    @abstractmethod
    def column_names(self) -> Sequence[Any]:
        """
        Return the column names as a materialized sequence
        """
        pass

    # TODO: Should this be a method or property?
    @property
    def row_names(self) -> Sequence[Any]:
        """
        Return the row names (if any) as a materialized sequence. It is not
        necessary to implement this method
        """
        raise NotImplementedError("row_names")

    def __getitem__(self, key: Hashable) -> Column:
        return self.column_by_name(key)

    @abstractmethod
    def column_by_name(self, key: Hashable) -> Column:
        """
        Return the column whose name is the indicated key
        """
        pass

    @abstractmethod
    def column_by_index(self, i: int) -> Column:
        """
        Return the column at the indicated position
        """
        pass


class MutableDataFrame(DataFrame, abc.MutableMapping):
    # TODO: Mutable data frames are fraught at this interface level and
    # need more discussion
    pass
