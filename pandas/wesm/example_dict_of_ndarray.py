# MIT License
#
# Copyright (c) 2020 Wes McKinney

from typing import Dict, Hashable, Sequence
import pandas.wesm.dataframe as dataframe

import numpy as np


_numeric_types = {
    "int8": dataframe.Int8(),
    "int16": dataframe.Int16(),
    "int32": dataframe.Int32(),
    "int64": dataframe.Int64(),
}


def _integer_factory(dtype):
    return _numeric_types[dtype.name]


def _constant_factory(type_instance):
    def factory(*unused):
        return type_instance

    return factory


_type_factories = {
    "b": _constant_factory(dataframe.Boolean()),
    "i": _integer_factory,
    "O": _constant_factory(dataframe.Object()),
    "S": _constant_factory(dataframe.Binary()),
    "U": _constant_factory(dataframe.String()),
}


class NumPyColumn(dataframe.Column):
    def __init__(self, name, data):
        self._name = name
        self._data = data

    @property
    def name(self) -> Hashable:
        return self._name

    @property
    def type(self) -> dataframe.DataType:
        factory = _type_factories.get(self._data.dtype.kind)
        if factory is None:
            raise NotImplementedError(
                "Data frame type for NumPy Type {} "
                "not known".format(str(self._data.dtype))
            )
        return factory(self._data.dtype)

    def to_numpy(self):
        return self._data


class DictDataFrame(dataframe.DataFrame):
    """
    Construct data frame from dict of NumPy arrays

    Parameters
    ----------
    data : dict
    names : sequence, default None
        If not passed, the names will be determined by the data's keys
    num_rows : int, default None
        If not passed, determined from the data
    """

    def __init__(
        self,
        columns: Dict[Hashable, np.ndarray],
        names: Sequence[Hashable] = None,
        num_rows: int = None,
    ):
        if names is None:
            names = list(columns.keys())

        assert len(columns) == len(names)

        self._columns = columns.copy()
        self._names = list(names)
        # self._name_to_index = {i: k for i, k in enumerate(self._names)}

        if len(columns) > 0:
            assert num_rows is None
            self._num_rows = len(next(iter(columns.values())))
        else:
            self._num_rows = num_rows

    @property
    def num_columns(self):
        return len(self._columns)

    @property
    def num_rows(self):
        return self._num_rows

    def iter_column_names(self):
        return iter(self._names)

    @property
    def column_names(self):
        return self._names

    def column_by_name(self, key: Hashable) -> NumPyColumn:
        return NumPyColumn(key, self._columns[key])

    def column_by_index(self, i: int) -> NumPyColumn:
        return NumPyColumn(self._names[i], self._columns[self._names[i]])


def get_example():
    data = {
        "a": np.array([1, 2, 3, 4, 5], dtype="int64"),
        "b": np.array(["a", "b", "c", "d", "e"]),
        "c": np.array([True, False, True, False, True]),
    }
    names = ["a", "b", "c"]
    return data, names, DictDataFrame(data, names=names)


def test_basic_behavior():
    raw_data, names, df = get_example()

    assert len(df) == 5
    assert df.num_columns == 3
    assert df.num_rows == 5

    for i, name in enumerate(df.column_names):
        assert name == names[i]

    for i, name in enumerate(df.iter_column_names()):
        assert name == names[i]

    expected_types = {
        "a": dataframe.Int64(),
        "b": dataframe.String(),
        "c": dataframe.Boolean(),
    }

    for i, name in enumerate(names):
        col = df[name]
        assert col.name == name
        assert col.type == expected_types[name]
        assert col.to_numpy() is raw_data[name]
        assert df.column_by_index(i).name == col.name
