# pylint: disable=E1101,E1103
# pylint: disable=W0703,W0622,W0613,W0201
import numpy as np

from pandas.core.dtypes.common import is_list_like
from pandas import compat

from pandas.core.frame import DataFrame
from pandas.core.index import MultiIndex

from pandas.core.frame import _shared_docs
from pandas.util._decorators import Appender


@Appender(_shared_docs['melt'] %
          dict(caller='pd.melt(df, ',
               versionadded="",
               other='DataFrame.melt'))
def melt(frame, id_vars=None, value_vars=None, var_name=None,
         value_name='value', col_level=None):
    # TODO: what about the existing index?
    if id_vars is not None:
        if not is_list_like(id_vars):
            id_vars = [id_vars]
        elif (isinstance(frame.columns, MultiIndex) and
              not isinstance(id_vars, list)):
            raise ValueError('id_vars must be a list of tuples when columns'
                             ' are a MultiIndex')
        else:
            id_vars = list(id_vars)
    else:
        id_vars = []

    if value_vars is not None:
        if not is_list_like(value_vars):
            value_vars = [value_vars]
        elif (isinstance(frame.columns, MultiIndex) and
              not isinstance(value_vars, list)):
            raise ValueError('value_vars must be a list of tuples when'
                             ' columns are a MultiIndex')
        else:
            value_vars = list(value_vars)
        frame = frame.loc[:, id_vars + value_vars]
    else:
        frame = frame.copy()

    if col_level is not None:  # allow list or other?
        # frame is a copy
        frame.columns = frame.columns.get_level_values(col_level)

    if var_name is None:
        if isinstance(frame.columns, MultiIndex):
            if len(frame.columns.names) == len(set(frame.columns.names)):
                var_name = frame.columns.names
            else:
                var_name = ['variable_{i}'.format(i=i)
                            for i in range(len(frame.columns.names))]
        else:
            var_name = [frame.columns.name if frame.columns.name is not None
                        else 'variable']
    if isinstance(var_name, compat.string_types):
        var_name = [var_name]

    N, K = frame.shape
    K -= len(id_vars)

    mdata = {}
    for col in id_vars:
        mdata[col] = np.tile(frame.pop(col).values, K)

    mcolumns = id_vars + var_name + [value_name]

    mdata[value_name] = frame.values.ravel('F')
    for i, col in enumerate(var_name):
        # asanyarray will keep the columns as an Index
        mdata[col] = np.asanyarray(frame.columns
                                   ._get_level_values(i)).repeat(N)

    return DataFrame(mdata, columns=mcolumns)
