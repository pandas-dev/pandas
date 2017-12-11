# pylint: disable=E1101,E1103
# pylint: disable=W0703,W0622,W0613,W0201
import numpy as np

from pandas.core.dtypes.common import is_list_like
from pandas import compat
from pandas.core.categorical import Categorical

from pandas.core.frame import DataFrame
from pandas.core.index import MultiIndex

from pandas.core.frame import _shared_docs
from pandas.util._decorators import Appender

import re
from pandas.core.dtypes.missing import notna
from pandas.core.tools.numeric import to_numeric


# ensure that the the correct tuple/list is used when selecting
# from a MultiIndex/Index frame
def check_vars(frame, var, var_string, num_col_levels):
    for v in var:
        if num_col_levels > 1:
            if not isinstance(v, tuple):
                raise ValueError('{} must be a list of tuples'
                                 ' when columns are a MultiIndex'
                                 .format(var_string))
            elif len(v) != num_col_levels:
                raise ValueError('all tuples in {} must be length {}'
                                 .format(var_string,
                                         frame.columns.nlevels))
        else:
            if is_list_like(v) and len(v) > 1:
                raise ValueError('DataFrame has only a single level of '
                                 'columns. {} is not a column'.format(v))


def melt_one(frame, id_vars=None, value_vars=None, var_name=None,
             value_name='value', col_level=None, extra_group=0,
             var_end=None):
    """
    melts exactly one group. Parameters are a single list not
    a list of lists like the main melt function
    """

    # TODO: what about the existing index?
    # Assume all column levels used when none given
    if len(col_level) == 0:
        num_col_levels = frame.columns.nlevels
    else:
        num_col_levels = len(col_level)

    check_vars(frame, id_vars, 'id_vars', num_col_levels)
    check_vars(frame, value_vars, 'value_vars', num_col_levels)

    if var_name != [] and len(var_name) != num_col_levels:
        raise ValueError('Length of var_name must match effective number of '
                         'column levels.')

    # allow both integer location and label for column levels
    if col_level != []:
        droplevels = list(range(frame.columns.nlevels))
        for level in col_level:
            if isinstance(level, int):
                droplevels.remove(level)
            else:
                droplevels.remove(frame.columns.names.index(level))
        if droplevels != []:
            frame = frame.copy()
            frame.columns = frame.columns.droplevel(droplevels)

    for iv in id_vars:
        if iv not in frame.columns:
            raise KeyError('{} not in columns'.format(iv))

    if value_vars != []:
        for vv in value_vars:
            if vv not in frame.columns:
                raise KeyError('{} not in columns'.format(vv))

    # use column level names if available, if not auto-name them
    if var_name == []:
        names = list(frame.columns.names)
        if len(names) == 1:
            if names[0] is None:
                var_name.append('variable')
            else:
                var_name.append(names[0])
        elif names.count(None) == 1:
            names[names.index(None)] = 'variable'
            var_name = names
        else:
            # small API break - use column level names when avaialable
            missing_name_count = 0
            for name in names:
                if name is None:
                    var_name.append('variable_{}'.format(missing_name_count))
                    missing_name_count += 1
                else:
                    var_name.append(name)

    # when using default var_name, append int to make unique
    if var_end is not None:
        var_name = [vn + '_' + str(var_end) for vn in var_name]

    N = len(frame)

    # find integer location of all the melted columns
    non_id_ilocs = []
    if value_vars != []:
        for v in value_vars:
            for i, v1 in enumerate(frame.columns):
                if v == v1:
                    non_id_ilocs.append(i)
    else:
        if id_vars == []:
            non_id_ilocs = list(range(frame.shape[1]))
        else:
            for i, v in enumerate(frame.columns):
                if v not in id_vars:
                    non_id_ilocs.append(i)

    K = len(non_id_ilocs)

    mdata = {}
    mcolumns = []

    # id_vars do not get melted, but need to repeat for each
    # column in melted group. extra_group is used for cases
    # when first group is not the longest
    for col in id_vars:
        pandas_obj = frame[col]
        if isinstance(pandas_obj, DataFrame):
            for i in range(pandas_obj.shape[1]):
                col_name = col + '_id_' + str(i)
                mdata[col_name] = np.tile(pandas_obj.iloc[:, i].values,
                                          K + extra_group)
                mcolumns.append(col_name)
        else:
            mdata[col] = np.tile(pandas_obj, K + extra_group)
            mcolumns.append(col)

    # melt all the columns into one long array
    values = np.concatenate([frame.iloc[:, i] for i in non_id_ilocs])
    if extra_group > 0:
        values = np.concatenate((values, np.full([N * extra_group], np.nan)))
    mdata[value_name[0]] = values

    # the column names of the melted groups need to repeat
    for i, col in enumerate(var_name):
        values = frame.columns[non_id_ilocs]._get_level_values(i)

        if isinstance(values, MultiIndex):
            # asanyarray will keep the columns as an Index
            values = np.asanyarray(values).repeat(N)
        else:
            # faster to use lists than np.repeat?
            data_list = []
            for v in values.tolist():
                data_list.extend([v] * N)
            values = data_list

        # Append missing values for any groups shorter than largest
        if extra_group > 0:
            values = np.concatenate((values,
                                     np.full([N * extra_group], np.nan)))
        mdata[col] = values
    mcolumns += var_name + value_name

    return mdata, mcolumns


def convert_to_list(val):
        if val is None:
            return []
        elif isinstance(val, np.ndarray):
            return val.tolist()
        elif not is_list_like(val):
            return [val]
        else:
            return list(val)


@Appender(_shared_docs['melt'] %
          dict(caller='pd.melt(df, ',
               versionadded="",
               other='DataFrame.melt'))
def melt(frame, id_vars=None, value_vars=None, var_name=None,
         value_name='value', col_level=None):

    # much easier to handle parameters when they are all lists
    id_vars = convert_to_list(id_vars)
    value_vars = convert_to_list(value_vars)
    var_name = convert_to_list(var_name)
    value_name = convert_to_list(value_name)
    col_level = convert_to_list(col_level)

    # Whan a list of list is passed, assume multiple melt groups
    if value_vars != [] and isinstance(value_vars[0], list):
        if var_name != []:
            if len(value_vars) != len(var_name):
                raise ValueError('Number of inner lists of value_vars must '
                                 'equal length of var_name '
                                 '{} != {}'.format(len(value_vars),
                                                   len(var_name)))
        else:
            # for consistency, when the default var_name is used
            var_name = [[]] * len(value_vars)

        if len(value_name) > 1:
            if len(value_vars) != len(value_name):
                raise ValueError('Number of inner lists of value_vars must '
                                 'equal length of value_name '
                                 '{} != {}'.format(len(value_vars),
                                                   len(value_name)))
        # allow for value_name to be a single item and attach int to it
        else:
            value_name = [value_name[0] + '_' + str(i)
                          for i in range(len(value_vars))]

        # get the total number of columns in each melt group
        # This is not just the length of the list because this function
        # handles columns with the same names, which is commong when
        # using multiindex frames
        value_vars_length = []
        for vv in value_vars:
            count = 0
            for col in frame.columns.values:
                if col in vv:
                    count += 1
            value_vars_length.append(count)

        # Need the max number of columns for all the melt groups to
        # correctly append NaNs to end of unbalanced melt groups
        max_group_len = max(value_vars_length)

        # store each melted group as a dictionary in a list
        mdata_list = []

        # store columns from each melted group in a list of lists
        mcolumns_list = []

        # individually melt each group
        vars_zipped = zip(value_vars, var_name, value_name, value_vars_length)
        for i, (val_v, var_n, val_n, vvl) in enumerate(vars_zipped):
            var_n = convert_to_list(var_n)
            val_n = convert_to_list(val_n)

            # only melt the id_vars for the first group
            id_vars_ = [] if i > 0 else id_vars

            # append int at end of var_name to make unique
            var_end = i if var_n == [] else None

            md, mc = melt_one(frame, id_vars=id_vars_, value_vars=val_v,
                              var_name=var_n, value_name=val_n,
                              col_level=col_level,
                              extra_group=max_group_len - vvl,
                              var_end=var_end)

            mdata_list.append(md)
            mcolumns_list.append(mc)

        # make one large dictionary with all data for constructor
        mdata = {}
        for d in mdata_list:
            mdata.update(d)

        mcolumns = [e for lst in mcolumns_list for e in lst]
        return DataFrame(mdata, columns=mcolumns)

    else:
        mdata, mcolumns = melt_one(frame, id_vars=id_vars,
                                   value_vars=value_vars, var_name=var_name,
                                   value_name=value_name, col_level=col_level)
        return DataFrame(mdata, columns=mcolumns)


def lreshape(data, groups, dropna=True, label=None):
    """
    Reshape long-format data to wide. Generalized inverse of DataFrame.pivot

    Parameters
    ----------
    data : DataFrame
    groups : dict
        {new_name : list_of_columns}
    dropna : boolean, default True

    Examples
    --------
    >>> import pandas as pd
    >>> data = pd.DataFrame({'hr1': [514, 573], 'hr2': [545, 526],
    ...                      'team': ['Red Sox', 'Yankees'],
    ...                      'year1': [2007, 2007], 'year2': [2008, 2008]})
    >>> data
       hr1  hr2     team  year1  year2
    0  514  545  Red Sox   2007   2008
    1  573  526  Yankees   2007   2008

    >>> pd.lreshape(data, {'year': ['year1', 'year2'], 'hr': ['hr1', 'hr2']})
          team  year   hr
    0  Red Sox  2007  514
    1  Yankees  2007  573
    2  Red Sox  2008  545
    3  Yankees  2008  526

    Returns
    -------
    reshaped : DataFrame
    """
    if isinstance(groups, dict):
        keys = list(groups.keys())
        values = list(groups.values())
    else:
        keys, values = zip(*groups)

    all_cols = list(set.union(*[set(x) for x in values]))
    id_cols = list(data.columns.difference(all_cols))

    K = len(values[0])

    for seq in values:
        if len(seq) != K:
            raise ValueError('All column lists must be same length')

    mdata = {}
    pivot_cols = []

    for target, names in zip(keys, values):
        to_concat = [data[col].values for col in names]

        import pandas.core.dtypes.concat as _concat
        mdata[target] = _concat._concat_compat(to_concat)
        pivot_cols.append(target)

    for col in id_cols:
        mdata[col] = np.tile(data[col].values, K)

    if dropna:
        mask = np.ones(len(mdata[pivot_cols[0]]), dtype=bool)
        for c in pivot_cols:
            mask &= notna(mdata[c])
        if not mask.all():
            mdata = {k: v[mask] for k, v in compat.iteritems(mdata)}

    from pandas import DataFrame
    return DataFrame(mdata, columns=id_cols + pivot_cols)


def wide_to_long(df, stubnames, i, j, sep="", suffix=r'\d+'):
    r"""
    Wide panel to long format. Less flexible but more user-friendly than melt.

    With stubnames ['A', 'B'], this function expects to find one or more
    group of columns with format Asuffix1, Asuffix2,..., Bsuffix1, Bsuffix2,...
    You specify what you want to call this suffix in the resulting long format
    with `j` (for example `j='year'`)

    Each row of these wide variables are assumed to be uniquely identified by
    `i` (can be a single column name or a list of column names)

    All remaining variables in the data frame are left intact.

    Parameters
    ----------
    df : DataFrame
        The wide-format DataFrame
    stubnames : str or list-like
        The stub name(s). The wide format variables are assumed to
        start with the stub names.
    i : str or list-like
        Column(s) to use as id variable(s)
    j : str
        The name of the subobservation variable. What you wish to name your
        suffix in the long format.
    sep : str, default ""
        A character indicating the separation of the variable names
        in the wide format, to be stripped from the names in the long format.
        For example, if your column names are A-suffix1, A-suffix2, you
        can strip the hypen by specifying `sep='-'`

        .. versionadded:: 0.20.0

    suffix : str, default '\\d+'
        A regular expression capturing the wanted suffixes. '\\d+' captures
        numeric suffixes. Suffixes with no numbers could be specified with the
        negated character class '\\D+'. You can also further disambiguate
        suffixes, for example, if your wide variables are of the form
        Aone, Btwo,.., and you have an unrelated column Arating, you can
        ignore the last one by specifying `suffix='(!?one|two)'`

        .. versionadded:: 0.20.0

        .. versionchanged:: 0.22.0
            When all suffixes are numeric, they are cast to int64/float64.

    Returns
    -------
    DataFrame
        A DataFrame that contains each stub name as a variable, with new index
        (i, j)

    Examples
    --------
    >>> import pandas as pd
    >>> import numpy as np
    >>> np.random.seed(123)
    >>> df = pd.DataFrame({"A1970" : {0 : "a", 1 : "b", 2 : "c"},
    ...                    "A1980" : {0 : "d", 1 : "e", 2 : "f"},
    ...                    "B1970" : {0 : 2.5, 1 : 1.2, 2 : .7},
    ...                    "B1980" : {0 : 3.2, 1 : 1.3, 2 : .1},
    ...                    "X"     : dict(zip(range(3), np.random.randn(3)))
    ...                   })
    >>> df["id"] = df.index
    >>> df
      A1970 A1980  B1970  B1980         X  id
    0     a     d    2.5    3.2 -1.085631   0
    1     b     e    1.2    1.3  0.997345   1
    2     c     f    0.7    0.1  0.282978   2
    >>> pd.wide_to_long(df, ["A", "B"], i="id", j="year")
    ... # doctest: +NORMALIZE_WHITESPACE
                    X  A    B
    id year
    0  1970 -1.085631  a  2.5
    1  1970  0.997345  b  1.2
    2  1970  0.282978  c  0.7
    0  1980 -1.085631  d  3.2
    1  1980  0.997345  e  1.3
    2  1980  0.282978  f  0.1

    With multuple id columns

    >>> df = pd.DataFrame({
    ...     'famid': [1, 1, 1, 2, 2, 2, 3, 3, 3],
    ...     'birth': [1, 2, 3, 1, 2, 3, 1, 2, 3],
    ...     'ht1': [2.8, 2.9, 2.2, 2, 1.8, 1.9, 2.2, 2.3, 2.1],
    ...     'ht2': [3.4, 3.8, 2.9, 3.2, 2.8, 2.4, 3.3, 3.4, 2.9]
    ... })
    >>> df
       birth  famid  ht1  ht2
    0      1      1  2.8  3.4
    1      2      1  2.9  3.8
    2      3      1  2.2  2.9
    3      1      2  2.0  3.2
    4      2      2  1.8  2.8
    5      3      2  1.9  2.4
    6      1      3  2.2  3.3
    7      2      3  2.3  3.4
    8      3      3  2.1  2.9
    >>> l = pd.wide_to_long(df, stubnames='ht', i=['famid', 'birth'], j='age')
    >>> l
    ... # doctest: +NORMALIZE_WHITESPACE
                      ht
    famid birth age
    1     1     1    2.8
                2    3.4
          2     1    2.9
                2    3.8
          3     1    2.2
                2    2.9
    2     1     1    2.0
                2    3.2
          2     1    1.8
                2    2.8
          3     1    1.9
                2    2.4
    3     1     1    2.2
                2    3.3
          2     1    2.3
                2    3.4
          3     1    2.1
                2    2.9

    Going from long back to wide just takes some creative use of `unstack`

    >>> w = l.unstack()
    >>> w.columns = w.columns.map('{0[0]}{0[1]}'.format)
    >>> w.reset_index()
       famid  birth  ht1  ht2
    0      1      1  2.8  3.4
    1      1      2  2.9  3.8
    2      1      3  2.2  2.9
    3      2      1  2.0  3.2
    4      2      2  1.8  2.8
    5      2      3  1.9  2.4
    6      3      1  2.2  3.3
    7      3      2  2.3  3.4
    8      3      3  2.1  2.9

    Less wieldy column names are also handled

    >>> np.random.seed(0)
    >>> df = pd.DataFrame({'A(quarterly)-2010': np.random.rand(3),
    ...                    'A(quarterly)-2011': np.random.rand(3),
    ...                    'B(quarterly)-2010': np.random.rand(3),
    ...                    'B(quarterly)-2011': np.random.rand(3),
    ...                    'X' : np.random.randint(3, size=3)})
    >>> df['id'] = df.index
    >>> df # doctest: +NORMALIZE_WHITESPACE, +ELLIPSIS
       A(quarterly)-2010  A(quarterly)-2011  B(quarterly)-2010  ...
    0           0.548814           0.544883           0.437587  ...
    1           0.715189           0.423655           0.891773  ...
    2           0.602763           0.645894           0.963663  ...
       X  id
    0  0   0
    1  1   1
    2  1   2

    >>> pd.wide_to_long(df, ['A(quarterly)', 'B(quarterly)'], i='id',
    ...                 j='year', sep='-')
    ... # doctest: +NORMALIZE_WHITESPACE
             X  A(quarterly)  B(quarterly)
    id year
    0  2010  0      0.548814     0.437587
    1  2010  1      0.715189     0.891773
    2  2010  1      0.602763     0.963663
    0  2011  0      0.544883     0.383442
    1  2011  1      0.423655     0.791725
    2  2011  1      0.645894     0.528895

    If we have many columns, we could also use a regex to find our
    stubnames and pass that list on to wide_to_long

    >>> stubnames = sorted(
    ...     set([match[0] for match in df.columns.str.findall(
    ...         r'[A-B]\(.*\)').values if match != [] ])
    ... )
    >>> list(stubnames)
    ['A(quarterly)', 'B(quarterly)']

    All of the above examples have integers as suffixes. It is possible to
    have non-integers as suffixes.

    >>> df = pd.DataFrame({
    ...     'famid': [1, 1, 1, 2, 2, 2, 3, 3, 3],
    ...     'birth': [1, 2, 3, 1, 2, 3, 1, 2, 3],
    ...     'ht_one': [2.8, 2.9, 2.2, 2, 1.8, 1.9, 2.2, 2.3, 2.1],
    ...     'ht_two': [3.4, 3.8, 2.9, 3.2, 2.8, 2.4, 3.3, 3.4, 2.9]
    ... })
    >>> df
       birth  famid  ht_one  ht_two
    0      1      1     2.8     3.4
    1      2      1     2.9     3.8
    2      3      1     2.2     2.9
    3      1      2     2.0     3.2
    4      2      2     1.8     2.8
    5      3      2     1.9     2.4
    6      1      3     2.2     3.3
    7      2      3     2.3     3.4
    8      3      3     2.1     2.9

    >>> l = pd.wide_to_long(df, stubnames='ht', i=['famid', 'birth'], j='age',
                            sep='_', suffix='\w')
    >>> l
    ... # doctest: +NORMALIZE_WHITESPACE
                      ht
    famid birth age
    1     1     one  2.8
                two  3.4
          2     one  2.9
                two  3.8
          3     one  2.2
                two  2.9
    2     1     one  2.0
                two  3.2
          2     one  1.8
                two  2.8
          3     one  1.9
                two  2.4
    3     1     one  2.2
                two  3.3
          2     one  2.3
                two  3.4
          3     one  2.1
                two  2.9

    Notes
    -----
    All extra variables are left untouched. This simply uses
    `pandas.melt` under the hood, but is hard-coded to "do the right thing"
    in a typical case.
    """
    def get_var_names(df, stub, sep, suffix):
        regex = r'^{stub}{sep}{suffix}$'.format(
            stub=re.escape(stub), sep=re.escape(sep), suffix=suffix)
        pattern = re.compile(regex)
        return [col for col in df.columns if pattern.match(col)]

    def melt_stub(df, stub, i, j, value_vars, sep):
        newdf = melt(df, id_vars=i, value_vars=value_vars,
                     value_name=stub.rstrip(sep), var_name=j)
        newdf[j] = Categorical(newdf[j])
        newdf[j] = newdf[j].str.replace(re.escape(stub + sep), "")

        # GH17627 Cast numerics suffixes to int/float
        newdf[j] = to_numeric(newdf[j], errors='ignore')

        return newdf.set_index(i + [j])

    if any(col in stubnames for col in df.columns):
        raise ValueError("stubname can't be identical to a column name")

    if not is_list_like(stubnames):
        stubnames = [stubnames]
    else:
        stubnames = list(stubnames)

    if not is_list_like(i):
        i = [i]
    else:
        i = list(i)

    if df[i].duplicated().any():
        raise ValueError("the id variables need to uniquely identify each row")

    value_vars = [get_var_names(df, stub, sep, suffix) for stub in stubnames]

    value_vars_flattened = [e for sublist in value_vars for e in sublist]
    id_vars = list(set(df.columns.tolist()).difference(value_vars_flattened))

    melted = []
    for s, v in zip(stubnames, value_vars):
        melted.append(melt_stub(df, s, i, j, v, sep))
    melted = melted[0].join(melted[1:], how='outer')

    if len(i) == 1:
        new = df[id_vars].set_index(i).join(melted)
        return new

    new = df[id_vars].merge(melted.reset_index(), on=i).set_index(i + [j])

    return new
