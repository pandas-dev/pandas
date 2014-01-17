"""
Contains data structures designed for manipulating panel (3-dimensional) data
"""
# pylint: disable=E1103,W0231,W0212,W0621
from __future__ import division
from pandas.compat import (map, zip, range, lrange, lmap, u, OrderedDict,
                           OrderedDefaultdict)
from pandas import compat
import sys
import numpy as np
from pandas.core.common import (PandasError, _try_sort, _default_index,
                                _infer_dtype_from_scalar, notnull)
from pandas.core.categorical import Categorical
from pandas.core.index import (Index, MultiIndex, _ensure_index,
                               _get_combined_index)
from pandas.core.indexing import _maybe_droplevels, _is_list_like
from pandas.core.internals import (BlockManager,
                                   create_block_manager_from_arrays,
                                   create_block_manager_from_blocks)
from pandas.core.series import Series
from pandas.core.frame import DataFrame
from pandas.core.generic import NDFrame, _shared_docs
from pandas.tools.util import cartesian_product
from pandas import compat
from pandas.util.decorators import deprecate, Appender, Substitution
import pandas.core.common as com
import pandas.core.ops as ops
import pandas.core.nanops as nanops
import pandas.computation.expressions as expressions


_shared_doc_kwargs = dict(
    axes='items, major_axis, minor_axis',
    klass="Panel",
    axes_single_arg="{0,1,2,'items','major_axis','minor_axis'}")
_shared_doc_kwargs['args_transpose'] = ("three positional arguments: each one"
                                        "of\n        %s" %
                                        _shared_doc_kwargs['axes_single_arg'])


def _ensure_like_indices(time, panels):
    """
    Makes sure that time and panels are conformable
    """
    n_time = len(time)
    n_panel = len(panels)
    u_panels = np.unique(panels)  # this sorts!
    u_time = np.unique(time)
    if len(u_time) == n_time:
        time = np.tile(u_time, len(u_panels))
    if len(u_panels) == n_panel:
        panels = np.repeat(u_panels, len(u_time))
    return time, panels


def panel_index(time, panels, names=['time', 'panel']):
    """
    Returns a multi-index suitable for a panel-like DataFrame

    Parameters
    ----------
    time : array-like
        Time index, does not have to repeat
    panels : array-like
        Panel index, does not have to repeat
    names : list, optional
        List containing the names of the indices

    Returns
    -------
    multi_index : MultiIndex
        Time index is the first level, the panels are the second level.

    Examples
    --------
    >>> years = range(1960,1963)
    >>> panels = ['A', 'B', 'C']
    >>> panel_idx = panel_index(years, panels)
    >>> panel_idx
    MultiIndex([(1960, 'A'), (1961, 'A'), (1962, 'A'), (1960, 'B'),
                (1961, 'B'), (1962, 'B'), (1960, 'C'), (1961, 'C'),
                (1962, 'C')], dtype=object)

    or

    >>> import numpy as np
    >>> years = np.repeat(range(1960,1963), 3)
    >>> panels = np.tile(['A', 'B', 'C'], 3)
    >>> panel_idx = panel_index(years, panels)
    >>> panel_idx
    MultiIndex([(1960, 'A'), (1960, 'B'), (1960, 'C'), (1961, 'A'),
                (1961, 'B'), (1961, 'C'), (1962, 'A'), (1962, 'B'),
                (1962, 'C')], dtype=object)
    """
    time, panels = _ensure_like_indices(time, panels)
    time_factor = Categorical.from_array(time)
    panel_factor = Categorical.from_array(panels)

    labels = [time_factor.labels, panel_factor.labels]
    levels = [time_factor.levels, panel_factor.levels]
    return MultiIndex(levels, labels, sortorder=None, names=names,
                      verify_integrity=False)


class Panel(NDFrame):

    """
    Represents wide format panel data, stored as 3-dimensional array

    Parameters
    ----------
    data : ndarray (items x major x minor), or dict of DataFrames
    items : Index or array-like
        axis=0
    major_axis : Index or array-like
        axis=1
    minor_axis : Index or array-like
        axis=2
    dtype : dtype, default None
        Data type to force, otherwise infer
    copy : boolean, default False
        Copy data from inputs. Only affects DataFrame / 2d ndarray input
    """

    @property
    def _constructor(self):
        return type(self)

    _constructor_sliced = DataFrame

    def __init__(self, data=None, items=None, major_axis=None, minor_axis=None,
                 copy=False, dtype=None):
        self._init_data(data=data, items=items, major_axis=major_axis,
                        minor_axis=minor_axis, copy=copy, dtype=dtype)

    def _init_data(self, data, copy, dtype, **kwargs):
        """
        Generate ND initialization; axes are passed
        as required objects to __init__
        """
        if data is None:
            data = {}
        if dtype is not None:
            dtype = self._validate_dtype(dtype)

        passed_axes = [kwargs.get(a) for a in self._AXIS_ORDERS]
        axes = None
        if isinstance(data, BlockManager):
            if any(x is not None for x in passed_axes):
                axes = [x if x is not None else y
                        for x, y in zip(passed_axes, data.axes)]
            mgr = data
        elif isinstance(data, dict):
            mgr = self._init_dict(data, passed_axes, dtype=dtype)
            copy = False
            dtype = None
        elif isinstance(data, (np.ndarray, list)):
            mgr = self._init_matrix(data, passed_axes, dtype=dtype, copy=copy)
            copy = False
            dtype = None
        else:  # pragma: no cover
            raise PandasError('Panel constructor not properly called!')

        NDFrame.__init__(self, mgr, axes=axes, copy=copy, dtype=dtype)

    def _init_dict(self, data, axes, dtype=None):
        haxis = axes.pop(self._info_axis_number)

        # prefilter if haxis passed
        if haxis is not None:
            haxis = _ensure_index(haxis)
            data = OrderedDict((k, v) for k, v
                               in compat.iteritems(data) if k in haxis)
        else:
            ks = list(data.keys())
            if not isinstance(data, OrderedDict):
                ks = _try_sort(ks)
            haxis = Index(ks)

        for k, v in compat.iteritems(data):
            if isinstance(v, dict):
                data[k] = self._constructor_sliced(v)

        # extract axis for remaining axes & create the slicemap
        raxes = [self._extract_axis(self, data, axis=i)
                 if a is None else a for i, a in enumerate(axes)]
        raxes_sm = self._extract_axes_for_slice(self, raxes)

        # shallow copy
        arrays = []
        haxis_shape = [len(a) for a in raxes]
        for h in haxis:
            v = values = data.get(h)
            if v is None:
                values = np.empty(haxis_shape, dtype=dtype)
                values.fill(np.nan)
            elif isinstance(v, self._constructor_sliced):
                d = raxes_sm.copy()
                d['copy'] = False
                v = v.reindex(**d)
                if dtype is not None:
                    v = v.astype(dtype)
                values = v.values
            arrays.append(values)

        return self._init_arrays(arrays, haxis, [haxis] + raxes)

    def _init_arrays(self, arrays, arr_names, axes):
        return create_block_manager_from_arrays(arrays, arr_names, axes)

    @classmethod
    def from_dict(cls, data, intersect=False, orient='items', dtype=None):
        """
        Construct Panel from dict of DataFrame objects

        Parameters
        ----------
        data : dict
            {field : DataFrame}
        intersect : boolean
            Intersect indexes of input DataFrames
        orient : {'items', 'minor'}, default 'items'
            The "orientation" of the data. If the keys of the passed dict
            should be the items of the result panel, pass 'items'
            (default). Otherwise if the columns of the values of the passed
            DataFrame objects should be the items (which in the case of
            mixed-dtype data you should do), instead pass 'minor'


        Returns
        -------
        Panel
        """
        orient = orient.lower()
        if orient == 'minor':
            new_data = OrderedDefaultdict(dict)
            for col, df in compat.iteritems(data):
                for item, s in compat.iteritems(df):
                    new_data[item][col] = s
            data = new_data
        elif orient != 'items':  # pragma: no cover
            raise ValueError('Orientation must be one of {items, minor}.')

        d = cls._homogenize_dict(cls, data, intersect=intersect, dtype=dtype)
        ks = list(d['data'].keys())
        if not isinstance(d['data'], OrderedDict):
            ks = list(sorted(ks))
        d[cls._info_axis_name] = Index(ks)
        return cls(**d)

    def __getitem__(self, key):
        if isinstance(self._info_axis, MultiIndex):
            return self._getitem_multilevel(key)
        return super(Panel, self).__getitem__(key)

    def _getitem_multilevel(self, key):
        info = self._info_axis
        loc = info.get_loc(key)
        if isinstance(loc, (slice, np.ndarray)):
            new_index = info[loc]
            result_index = _maybe_droplevels(new_index, key)
            slices = [loc] + [slice(None) for x in range(
                self._AXIS_LEN - 1)]
            new_values = self.values[slices]

            d = self._construct_axes_dict(self._AXIS_ORDERS[1:])
            d[self._info_axis_name] = result_index
            result = self._constructor(new_values, **d)
            return result
        else:
            return self._get_item_cache(key)

    def _init_matrix(self, data, axes, dtype=None, copy=False):
        values = self._prep_ndarray(self, data, copy=copy)

        if dtype is not None:
            try:
                values = values.astype(dtype)
            except Exception:
                raise ValueError('failed to cast to %s' % dtype)

        shape = values.shape
        fixed_axes = []
        for i, ax in enumerate(axes):
            if ax is None:
                ax = _default_index(shape[i])
            else:
                ax = _ensure_index(ax)
            fixed_axes.append(ax)

        return create_block_manager_from_blocks([values], fixed_axes)

    #----------------------------------------------------------------------
    # Comparison methods

    def _compare_constructor(self, other, func):
        if not self._indexed_same(other):
            raise Exception('Can only compare identically-labeled '
                            'same type objects')

        new_data = {}
        for col in self._info_axis:
            new_data[col] = func(self[col], other[col])

        d = self._construct_axes_dict(copy=False)
        return self._constructor(data=new_data, **d)

    #----------------------------------------------------------------------
    # Magic methods

    def __unicode__(self):
        """
        Return a string representation for a particular Panel

        Invoked by unicode(df) in py2 only.
        Yields a Unicode String in both py2/py3.
        """

        class_name = str(self.__class__)

        shape = self.shape
        dims = u('Dimensions: %s') % ' x '.join(
            ["%d (%s)" % (s, a) for a, s in zip(self._AXIS_ORDERS, shape)])

        def axis_pretty(a):
            v = getattr(self, a)
            if len(v) > 0:
                return u('%s axis: %s to %s') % (a.capitalize(),
                                                 com.pprint_thing(v[0]),
                                                 com.pprint_thing(v[-1]))
            else:
                return u('%s axis: None') % a.capitalize()

        output = '\n'.join(
            [class_name, dims] + [axis_pretty(a) for a in self._AXIS_ORDERS])
        return output

    def _get_plane_axes_index(self, axis):
        """
        Get my plane axes indexes: these are already
        (as compared with higher level planes),
        as we are returning a DataFrame axes indexes
        """
        axis_name = self._get_axis_name(axis)

        if axis_name == 'major_axis':
            index = 'minor_axis'
            columns = 'items'
        if axis_name == 'minor_axis':
            index = 'major_axis'
            columns = 'items'
        elif axis_name == 'items':
            index = 'major_axis'
            columns = 'minor_axis'

        return index, columns

    def _get_plane_axes(self, axis):
        """
        Get my plane axes indexes: these are already
        (as compared with higher level planes),
        as we are returning a DataFrame axes
        """
        return [ self._get_axis(axi) for axi in self._get_plane_axes_index(axis) ]

    fromDict = from_dict

    def to_sparse(self, fill_value=None, kind='block'):
        """
        Convert to SparsePanel

        Parameters
        ----------
        fill_value : float, default NaN
        kind : {'block', 'integer'}

        Returns
        -------
        y : SparseDataFrame
        """
        from pandas.core.sparse import SparsePanel
        frames = dict(compat.iteritems(self))
        return SparsePanel(frames, items=self.items,
                           major_axis=self.major_axis,
                           minor_axis=self.minor_axis,
                           default_kind=kind,
                           default_fill_value=fill_value)

    def to_excel(self, path, na_rep='', engine=None, **kwargs):
        """
        Write each DataFrame in Panel to a separate excel sheet

        Parameters
        ----------
        path : string or ExcelWriter object
            File path or existing ExcelWriter
        na_rep : string, default ''
            Missing data representation
        engine : string, default None
            write engine to use - you can also set this via the options
            ``io.excel.xlsx.writer``, ``io.excel.xls.writer``, and
            ``io.excel.xlsm.writer``.

        Other Parameters
        ----------------
        float_format : string, default None
            Format string for floating point numbers
        cols : sequence, optional
            Columns to write
        header : boolean or list of string, default True
            Write out column names. If a list of string is given it is
            assumed to be aliases for the column names
        index : boolean, default True
            Write row names (index)
        index_label : string or sequence, default None
            Column label for index column(s) if desired. If None is given, and
            `header` and `index` are True, then the index names are used. A
            sequence should be given if the DataFrame uses MultiIndex.
        startow : upper left cell row to dump data frame
        startcol : upper left cell column to dump data frame

        Notes
        -----
        Keyword arguments (and na_rep) are passed to the ``to_excel`` method
        for each DataFrame written.
        """
        from pandas.io.excel import ExcelWriter

        if isinstance(path, compat.string_types):
            writer = ExcelWriter(path, engine=engine)
        else:
            writer = path
        kwargs['na_rep'] = na_rep

        for item, df in compat.iteritems(self):
            name = str(item)
            df.to_excel(writer, name, **kwargs)
        writer.save()

    def as_matrix(self):
        self._consolidate_inplace()
        return self._data.as_matrix()

    #----------------------------------------------------------------------
    # Getting and setting elements

    def get_value(self, *args):
        """
        Quickly retrieve single value at (item, major, minor) location

        Parameters
        ----------
        item : item label (panel item)
        major : major axis label (panel item row)
        minor : minor axis label (panel item column)

        Returns
        -------
        value : scalar value
        """
        nargs = len(args)
        nreq = self._AXIS_LEN

        # require an arg for each axis
        if nargs != nreq:
            raise TypeError('There must be an argument for each axis, you gave'
                            ' {0} args, but {1} are required'.format(nargs,
                                                                     nreq))

        # hm, two layers to the onion
        frame = self._get_item_cache(args[0])
        return frame.get_value(*args[1:])

    def set_value(self, *args):
        """
        Quickly set single value at (item, major, minor) location

        Parameters
        ----------
        item : item label (panel item)
        major : major axis label (panel item row)
        minor : minor axis label (panel item column)
        value : scalar

        Returns
        -------
        panel : Panel
            If label combo is contained, will be reference to calling Panel,
            otherwise a new object
        """
        # require an arg for each axis and the value
        nargs = len(args)
        nreq = self._AXIS_LEN + 1

        if nargs != nreq:
            raise TypeError('There must be an argument for each axis plus the '
                            'value provided, you gave {0} args, but {1} are '
                            'required'.format(nargs, nreq))

        try:
            frame = self._get_item_cache(args[0])
            frame.set_value(*args[1:])
            return self
        except KeyError:
            axes = self._expand_axes(args)
            d = self._construct_axes_dict_from(self, axes, copy=False)
            result = self.reindex(**d)
            args = list(args)
            likely_dtype, args[-1] = _infer_dtype_from_scalar(args[-1])
            made_bigger = not np.array_equal(
                axes[0], self._info_axis)
            # how to make this logic simpler?
            if made_bigger:
                com._possibly_cast_item(result, args[0], likely_dtype)

            return result.set_value(*args)

    def _box_item_values(self, key, values):
        if self.ndim == values.ndim:
            result = self._constructor(values)

            # a dup selection will yield a full ndim
            if result._get_axis(0).is_unique:
                result = result[key]

            return result

        d = self._construct_axes_dict_for_slice(self._AXIS_ORDERS[1:])
        return self._constructor_sliced(values, **d)

    def _slice(self, slobj, axis=0, raise_on_error=False, typ=None):
        new_data = self._data.get_slice(slobj,
                                        axis=axis,
                                        raise_on_error=raise_on_error)
        return self._constructor(new_data)

    def __setitem__(self, key, value):
        shape = tuple(self.shape)
        if isinstance(value, self._constructor_sliced):
            value = value.reindex(
                **self._construct_axes_dict_for_slice(self._AXIS_ORDERS[1:]))
            mat = value.values
        elif isinstance(value, np.ndarray):
            if value.shape != shape[1:]:
                raise ValueError(
                    'shape of value must be {0}, shape of given object was '
                    '{1}'.format(shape[1:], tuple(map(int, value.shape))))
            mat = np.asarray(value)
        elif np.isscalar(value):
            dtype, value = _infer_dtype_from_scalar(value)
            mat = np.empty(shape[1:], dtype=dtype)
            mat.fill(value)
        else:
            raise TypeError('Cannot set item of type: %s' % str(type(value)))

        mat = mat.reshape(tuple([1]) + shape[1:])
        NDFrame._set_item(self, key, mat)

    def _unpickle_panel_compat(self, state):  # pragma: no cover
        "Unpickle the panel"
        _unpickle = com._unpickle_array
        vals, items, major, minor = state

        items = _unpickle(items)
        major = _unpickle(major)
        minor = _unpickle(minor)
        values = _unpickle(vals)
        wp = Panel(values, items, major, minor)
        self._data = wp._data

    def conform(self, frame, axis='items'):
        """
        Conform input DataFrame to align with chosen axis pair.

        Parameters
        ----------
        frame : DataFrame
        axis : {'items', 'major', 'minor'}

            Axis the input corresponds to. E.g., if axis='major', then
            the frame's columns would be items, and the index would be
            values of the minor axis

        Returns
        -------
        DataFrame
        """
        axes = self._get_plane_axes(axis)
        return frame.reindex(**self._extract_axes_for_slice(self, axes))

    def head(self, n=5):
        raise NotImplementedError

    def tail(self, n=5):
        raise NotImplementedError

    def _needs_reindex_multi(self, axes, method, level):
        """ don't allow a multi reindex on Panel or above ndim """
        return False

    def dropna(self, axis=0, how='any', inplace=False, **kwargs):
        """
        Drop 2D from panel, holding passed axis constant

        Parameters
        ----------
        axis : int, default 0
            Axis to hold constant. E.g. axis=1 will drop major_axis entries
            having a certain amount of NA data
        how : {'all', 'any'}, default 'any'
            'any': one or more values are NA in the DataFrame along the
            axis. For 'all' they all must be.
        inplace : bool, default False
            If True, do operation inplace and return None.

        Returns
        -------
        dropped : Panel
        """
        axis = self._get_axis_number(axis)

        values = self.values
        mask = com.notnull(values)

        for ax in reversed(sorted(set(range(self._AXIS_LEN)) - set([axis]))):
            mask = mask.sum(ax)

        per_slice = np.prod(values.shape[:axis] + values.shape[axis + 1:])

        if how == 'all':
            cond = mask > 0
        else:
            cond = mask == per_slice

        new_ax = self._get_axis(axis)[cond]
        result = self.reindex_axis(new_ax, axis=axis)
        if inplace:
            self._update_inplace(result)
        else:
            return result

    def _combine(self, other, func, axis=0):
        if isinstance(other, Panel):
            return self._combine_panel(other, func)
        elif isinstance(other, DataFrame):
            return self._combine_frame(other, func, axis=axis)
        elif np.isscalar(other):
            return self._combine_const(other, func)

    def _combine_const(self, other, func):
        new_values = func(self.values, other)
        d = self._construct_axes_dict()
        return self._constructor(new_values, **d)

    def _combine_frame(self, other, func, axis=0):
        index, columns = self._get_plane_axes(axis)
        axis = self._get_axis_number(axis)

        other = other.reindex(index=index, columns=columns)

        if axis == 0:
            new_values = func(self.values, other.values)
        elif axis == 1:
            new_values = func(self.values.swapaxes(0, 1), other.values.T)
            new_values = new_values.swapaxes(0, 1)
        elif axis == 2:
            new_values = func(self.values.swapaxes(0, 2), other.values)
            new_values = new_values.swapaxes(0, 2)

        return self._constructor(new_values, self.items, self.major_axis,
                                 self.minor_axis)

    def _combine_panel(self, other, func):
        items = self.items + other.items
        major = self.major_axis + other.major_axis
        minor = self.minor_axis + other.minor_axis

        # could check that everything's the same size, but forget it
        this = self.reindex(items=items, major=major, minor=minor)
        other = other.reindex(items=items, major=major, minor=minor)

        result_values = func(this.values, other.values)

        return self._constructor(result_values, items, major, minor)

    def major_xs(self, key, copy=True):
        """
        Return slice of panel along major axis

        Parameters
        ----------
        key : object
            Major axis label
        copy : boolean, default True
            Copy data

        Returns
        -------
        y : DataFrame
            index -> minor axis, columns -> items
        """
        return self.xs(key, axis=self._AXIS_LEN - 2, copy=copy)

    def minor_xs(self, key, copy=True):
        """
        Return slice of panel along minor axis

        Parameters
        ----------
        key : object
            Minor axis label
        copy : boolean, default True
            Copy data

        Returns
        -------
        y : DataFrame
            index -> major axis, columns -> items
        """
        return self.xs(key, axis=self._AXIS_LEN - 1, copy=copy)

    def xs(self, key, axis=1, copy=True):
        """
        Return slice of panel along selected axis

        Parameters
        ----------
        key : object
            Label
        axis : {'items', 'major', 'minor}, default 1/'major'
        copy : boolean, default True
            Copy data

        Returns
        -------
        y : ndim(self)-1
        """
        axis = self._get_axis_number(axis)
        if axis == 0:
            data = self[key]
            if copy:
                data = data.copy()
            return data

        self._consolidate_inplace()
        axis_number = self._get_axis_number(axis)
        new_data = self._data.xs(key, axis=axis_number, copy=copy)
        return self._construct_return_type(new_data)

    _xs = xs

    def _ixs(self, i, axis=0):
        """
        i : int, slice, or sequence of integers
        axis : int
        """

        key = self._get_axis(axis)[i]

        # xs cannot handle a non-scalar key, so just reindex here
        if _is_list_like(key):
            indexer = {self._get_axis_name(axis): key}
            return self.reindex(**indexer)

        # a reduction
        if axis == 0:
            values = self._data.iget(i)
            return self._box_item_values(key, values)

        # xs by position
        self._consolidate_inplace()
        new_data = self._data.xs(i, axis=axis, copy=True, takeable=True)
        return self._construct_return_type(new_data)

    def groupby(self, function, axis='major'):
        """
        Group data on given axis, returning GroupBy object

        Parameters
        ----------
        function : callable
            Mapping function for chosen access
        axis : {'major', 'minor', 'items'}, default 'major'

        Returns
        -------
        grouped : PanelGroupBy
        """
        from pandas.core.groupby import PanelGroupBy
        axis = self._get_axis_number(axis)
        return PanelGroupBy(self, function, axis=axis)

    def to_frame(self, filter_observations=True):
        """
        Transform wide format into long (stacked) format as DataFrame whose
        columns are the Panel's items and whose index is a MultiIndex formed
        of the Panel's major and minor axes.

        Parameters
        ----------
        filter_observations : boolean, default True
            Drop (major, minor) pairs without a complete set of observations
            across all the items

        Returns
        -------
        y : DataFrame
        """
        _, N, K = self.shape

        if filter_observations:
            # shaped like the return DataFrame
            mask = com.notnull(self.values).all(axis=0)
            # size = mask.sum()
            selector = mask.ravel()
        else:
            # size = N * K
            selector = slice(None, None)

        data = {}
        for item in self.items:
            data[item] = self[item].values.ravel()[selector]

        def construct_multi_parts(idx, n_repeat, n_shuffle=1):
            axis_idx = idx.to_hierarchical(n_repeat, n_shuffle)
            labels = [x[selector] for x in axis_idx.labels]
            levels = axis_idx.levels
            names = axis_idx.names
            return labels, levels, names

        def construct_index_parts(idx, major=True):
            levels = [idx]
            if major:
                labels = [np.arange(N).repeat(K)[selector]]
                names = idx.name or 'major'
            else:
                labels = np.arange(K).reshape(1, K)[np.zeros(N, dtype=int)]
                labels = [labels.ravel()[selector]]
                names = idx.name or 'minor'
            names = [names]
            return labels, levels, names

        if isinstance(self.major_axis, MultiIndex):
            major_labels, major_levels, major_names = construct_multi_parts(
                self.major_axis, n_repeat=K)
        else:
            major_labels, major_levels, major_names = construct_index_parts(
                self.major_axis)

        if isinstance(self.minor_axis, MultiIndex):
            minor_labels, minor_levels, minor_names = construct_multi_parts(
                self.minor_axis, n_repeat=N, n_shuffle=K)
        else:
            minor_labels, minor_levels, minor_names = construct_index_parts(
                self.minor_axis, major=False)

        levels = major_levels + minor_levels
        labels = major_labels + minor_labels
        names = major_names + minor_names

        index = MultiIndex(levels=levels, labels=labels,
                           names=names, verify_integrity=False)

        return DataFrame(data, index=index, columns=self.items)

    to_long = deprecate('to_long', to_frame)
    toLong = deprecate('toLong', to_frame)

    def apply(self, func, axis='major', **kwargs):
        """
        Applies function along input axis of the Panel

        Parameters
        ----------
        func : function
            Function to apply to each combination of 'other' axes
            e.g. if axis = 'items', then the combination of major_axis/minor_axis
            will be passed a Series
        axis : {'major', 'minor', 'items'}
        Additional keyword arguments will be passed as keywords to the function

        Examples
        --------
        >>> p.apply(numpy.sqrt) # returns a Panel
        >>> p.apply(lambda x: x.sum(), axis=0) # equiv to p.sum(0)
        >>> p.apply(lambda x: x.sum(), axis=1) # equiv to p.sum(1)
        >>> p.apply(lambda x: x.sum(), axis=2) # equiv to p.sum(2)

        Returns
        -------
        result : Pandas Object
        """

        if kwargs and not isinstance(func, np.ufunc):
            f = lambda x: func(x, **kwargs)
        else:
            f = func

        # 2d-slabs
        if isinstance(axis, (tuple,list)) and len(axis) == 2:
            return self._apply_2d(f, axis=axis)

        axis = self._get_axis_number(axis)

        # try ufunc like
        if isinstance(f, np.ufunc):
            try:
                result = np.apply_along_axis(func, axis, self.values)
                return self._wrap_result(result, axis=axis)
            except (AttributeError):
                pass

        # 1d
        return self._apply_1d(f, axis=axis)

    def _apply_1d(self, func, axis):

        axis_name = self._get_axis_name(axis)
        ax = self._get_axis(axis)
        ndim = self.ndim
        values = self.values

        # iter thru the axes
        slice_axis = self._get_axis(axis)
        slice_indexer = [0]*(ndim-1)
        indexer = np.zeros(ndim, 'O')
        indlist = list(range(ndim))
        indlist.remove(axis)
        indexer[axis] = slice(None, None)
        indexer.put(indlist, slice_indexer)
        planes = [ self._get_axis(axi) for axi in indlist ]
        shape = np.array(self.shape).take(indlist)

        # all the iteration points
        points = cartesian_product(planes)

        results = []
        for i in range(np.prod(shape)):

            # construct the object
            pts = tuple([ p[i] for p in points ])
            indexer.put(indlist, slice_indexer)

            obj = Series(values[tuple(indexer)],index=slice_axis,name=pts)
            result = func(obj)

            results.append(result)

            # increment the indexer
            slice_indexer[-1] += 1
            n = -1
            while (slice_indexer[n] >= shape[n]) and (n > (1-ndim)):
                slice_indexer[n-1] += 1
                slice_indexer[n] = 0
                n -= 1

        # empty object
        if not len(results):
            return self._constructor(**self._construct_axes_dict())

        # same ndim as current
        if isinstance(results[0],Series):
            arr = np.vstack([ r.values for r in results ])
            arr = arr.T.reshape(tuple([len(slice_axis)] + list(shape)))
            tranp = np.array([axis]+indlist).argsort()
            arr = arr.transpose(tuple(list(tranp)))
            return self._constructor(arr,**self._construct_axes_dict())

        # ndim-1 shape
        results = np.array(results).reshape(shape)
        if results.ndim == 2 and axis_name != self._info_axis_name:
            results = results.T
            planes = planes[::-1]
        return self._construct_return_type(results,planes)

    def _apply_2d(self, func, axis):
        """ handle 2-d slices, equiv to iterating over the other axis """

        ndim = self.ndim
        axis = [ self._get_axis_number(a) for a in axis ]

        # construct slabs, in 2-d this is a DataFrame result
        indexer_axis = list(range(ndim))
        for a in axis:
            indexer_axis.remove(a)
        indexer_axis = indexer_axis[0]

        slicer = [ slice(None,None) ] * ndim
        ax = self._get_axis(indexer_axis)

        results = []
        for i, e in enumerate(ax):

            slicer[indexer_axis] = i
            sliced = self.iloc[tuple(slicer)]

            obj = func(sliced)
            results.append((e,obj))

        return self._construct_return_type(dict(results))

    def _reduce(self, op, axis=0, skipna=True, numeric_only=None,
                filter_type=None, **kwds):
        axis_name = self._get_axis_name(axis)
        axis_number = self._get_axis_number(axis_name)
        f = lambda x: op(x, axis=axis_number, skipna=skipna, **kwds)

        result = f(self.values)

        axes = self._get_plane_axes(axis_name)
        if result.ndim == 2 and axis_name != self._info_axis_name:
            result = result.T

        return self._construct_return_type(result, axes)

    def _construct_return_type(self, result, axes=None, **kwargs):
        """ return the type for the ndim of the result """
        ndim = getattr(result,'ndim',None)

        # need to assume they are the same
        if ndim is None:
            if isinstance(result,dict):
                ndim = getattr(list(compat.itervalues(result))[0],'ndim',None)

                # a saclar result
                if ndim is None:
                    ndim = 0

                # have a dict, so top-level is +1 dim
                else:
                    ndim += 1

        # scalar
        if ndim == 0:
            return Series(result)

        # same as self
        elif self.ndim == ndim:
            """ return the construction dictionary for these axes """
            if axes is None:
                return self._constructor(result)
            return self._constructor(result, **self._construct_axes_dict())

        # sliced
        elif self.ndim == ndim + 1:
            if axes is None:
                return self._constructor_sliced(result)
            return self._constructor_sliced(
                result, **self._extract_axes_for_slice(self, axes))

        raise PandasError('invalid _construct_return_type [self->%s] '
                          '[result->%s]' % (self, result))

    def _wrap_result(self, result, axis):
        axis = self._get_axis_name(axis)
        axes = self._get_plane_axes(axis)
        if result.ndim == 2 and axis != self._info_axis_name:
            result = result.T

        return self._construct_return_type(result, axes)

    @Appender(_shared_docs['reindex'] % _shared_doc_kwargs)
    def reindex(self, items=None, major_axis=None, minor_axis=None, **kwargs):
        major_axis = (major_axis if major_axis is not None
                      else kwargs.pop('major', None))
        minor_axis = (minor_axis if minor_axis is not None
                      else kwargs.pop('minor', None))
        return super(Panel, self).reindex(items=items, major_axis=major_axis,
                                          minor_axis=minor_axis, **kwargs)

    @Appender(_shared_docs['rename'] % _shared_doc_kwargs)
    def rename(self, items=None, major_axis=None, minor_axis=None, **kwargs):
        major_axis = (major_axis if major_axis is not None
                      else kwargs.pop('major', None))
        minor_axis = (minor_axis if minor_axis is not None
                      else kwargs.pop('minor', None))
        return super(Panel, self).rename(items=items, major_axis=major_axis,
                                         minor_axis=minor_axis, **kwargs)

    @Appender(_shared_docs['reindex_axis'] % _shared_doc_kwargs)
    def reindex_axis(self, labels, axis=0, method=None, level=None, copy=True,
                     limit=None, fill_value=np.nan):
        return super(Panel, self).reindex_axis(labels=labels, axis=axis,
                                               method=method, level=level,
                                               copy=copy, limit=limit,
                                               fill_value=fill_value)

    @Appender(_shared_docs['transpose'] % _shared_doc_kwargs)
    def transpose(self, *args, **kwargs):
        return super(Panel, self).transpose(*args, **kwargs)

    def count(self, axis='major'):
        """
        Return number of observations over requested axis.

        Parameters
        ----------
        axis : {'items', 'major', 'minor'} or {0, 1, 2}

        Returns
        -------
        count : DataFrame
        """
        i = self._get_axis_number(axis)

        values = self.values
        mask = np.isfinite(values)
        result = mask.sum(axis=i)

        return self._wrap_result(result, axis)

    def shift(self, lags, freq=None, axis='major'):
        """
        Shift major or minor axis by specified number of leads/lags. Drops
        periods right now compared with DataFrame.shift

        Parameters
        ----------
        lags : int
        axis : {'major', 'minor'}

        Returns
        -------
        shifted : Panel
        """
        values = self.values
        items = self.items
        major_axis = self.major_axis
        minor_axis = self.minor_axis

        if freq:
            return self.tshift(lags, freq, axis=axis)

        if lags > 0:
            vslicer = slice(None, -lags)
            islicer = slice(lags, None)
        elif lags == 0:
            vslicer = islicer = slice(None)
        else:
            vslicer = slice(-lags, None)
            islicer = slice(None, lags)

        axis = self._get_axis_name(axis)
        if axis == 'major_axis':
            values = values[:, vslicer, :]
            major_axis = major_axis[islicer]
        elif axis == 'minor_axis':
            values = values[:, :, vslicer]
            minor_axis = minor_axis[islicer]
        else:
            raise ValueError('Invalid axis')

        return self._constructor(values, items=items, major_axis=major_axis,
                                 minor_axis=minor_axis)

    def tshift(self, periods=1, freq=None, axis='major', **kwds):
        return super(Panel, self).tshift(periods, freq, axis, **kwds)

    def join(self, other, how='left', lsuffix='', rsuffix=''):
        """
        Join items with other Panel either on major and minor axes column

        Parameters
        ----------
        other : Panel or list of Panels
            Index should be similar to one of the columns in this one
        how : {'left', 'right', 'outer', 'inner'}
            How to handle indexes of the two objects. Default: 'left'
            for joining on index, None otherwise
            * left: use calling frame's index
            * right: use input frame's index
            * outer: form union of indexes
            * inner: use intersection of indexes
        lsuffix : string
            Suffix to use from left frame's overlapping columns
        rsuffix : string
            Suffix to use from right frame's overlapping columns

        Returns
        -------
        joined : Panel
        """
        from pandas.tools.merge import concat

        if isinstance(other, Panel):
            join_major, join_minor = self._get_join_index(other, how)
            this = self.reindex(major=join_major, minor=join_minor)
            other = other.reindex(major=join_major, minor=join_minor)
            merged_data = this._data.merge(other._data, lsuffix, rsuffix)
            return self._constructor(merged_data)
        else:
            if lsuffix or rsuffix:
                raise ValueError('Suffixes not supported when passing '
                                 'multiple panels')

            if how == 'left':
                how = 'outer'
                join_axes = [self.major_axis, self.minor_axis]
            elif how == 'right':
                raise ValueError('Right join not supported with multiple '
                                 'panels')
            else:
                join_axes = None

            return concat([self] + list(other), axis=0, join=how,
                          join_axes=join_axes, verify_integrity=True)

    def update(self, other, join='left', overwrite=True, filter_func=None,
               raise_conflict=False):
        """
        Modify Panel in place using non-NA values from passed
        Panel, or object coercible to Panel. Aligns on items

        Parameters
        ----------
        other : Panel, or object coercible to Panel
        join : How to join individual DataFrames
            {'left', 'right', 'outer', 'inner'}, default 'left'
        overwrite : boolean, default True
            If True then overwrite values for common keys in the calling panel
        filter_func : callable(1d-array) -> 1d-array<boolean>, default None
            Can choose to replace values other than NA. Return True for values
            that should be updated
        raise_conflict : bool
            If True, will raise an error if a DataFrame and other both
            contain data in the same place.
        """

        if not isinstance(other, self._constructor):
            other = self._constructor(other)

        axis_name = self._info_axis_name
        axis_values = self._info_axis
        other = other.reindex(**{axis_name: axis_values})

        for frame in axis_values:
            self[frame].update(other[frame], join, overwrite, filter_func,
                               raise_conflict)

    def _get_join_index(self, other, how):
        if how == 'left':
            join_major, join_minor = self.major_axis, self.minor_axis
        elif how == 'right':
            join_major, join_minor = other.major_axis, other.minor_axis
        elif how == 'inner':
            join_major = self.major_axis.intersection(other.major_axis)
            join_minor = self.minor_axis.intersection(other.minor_axis)
        elif how == 'outer':
            join_major = self.major_axis.union(other.major_axis)
            join_minor = self.minor_axis.union(other.minor_axis)
        return join_major, join_minor

    # miscellaneous data creation
    @staticmethod
    def _extract_axes(self, data, axes, **kwargs):
        """ return a list of the axis indicies """
        return [self._extract_axis(self, data, axis=i, **kwargs) for i, a
                in enumerate(axes)]

    @staticmethod
    def _extract_axes_for_slice(self, axes):
        """ return the slice dictionary for these axes """
        return dict([(self._AXIS_SLICEMAP[i], a)
                     for i, a in zip(self._AXIS_ORDERS[self._AXIS_LEN -
                                                       len(axes):], axes)])

    @staticmethod
    def _prep_ndarray(self, values, copy=True):
        if not isinstance(values, np.ndarray):
            values = np.asarray(values)
            # NumPy strings are a pain, convert to object
            if issubclass(values.dtype.type, compat.string_types):
                values = np.array(values, dtype=object, copy=True)
        else:
            if copy:
                values = values.copy()
        if values.ndim != self._AXIS_LEN:
            raise ValueError("The number of dimensions required is {0}, "
                             "but the number of dimensions of the "
                             "ndarray given was {1}".format(self._AXIS_LEN,
                                                            values.ndim))
        return values

    @staticmethod
    def _homogenize_dict(self, frames, intersect=True, dtype=None):
        """
        Conform set of _constructor_sliced-like objects to either
        an intersection of indices / columns or a union.

        Parameters
        ----------
        frames : dict
        intersect : boolean, default True

        Returns
        -------
        dict of aligned results & indicies
        """

        result = dict()
        # caller differs dict/ODict, presered type
        if isinstance(frames, OrderedDict):
            result = OrderedDict()

        adj_frames = OrderedDict()
        for k, v in compat.iteritems(frames):
            if isinstance(v, dict):
                adj_frames[k] = self._constructor_sliced(v)
            else:
                adj_frames[k] = v

        axes = self._AXIS_ORDERS[1:]
        axes_dict = dict([(a, ax) for a, ax in zip(axes, self._extract_axes(
            self, adj_frames, axes, intersect=intersect))])

        reindex_dict = dict(
            [(self._AXIS_SLICEMAP[a], axes_dict[a]) for a in axes])
        reindex_dict['copy'] = False
        for key, frame in compat.iteritems(adj_frames):
            if frame is not None:
                result[key] = frame.reindex(**reindex_dict)
            else:
                result[key] = None

        axes_dict['data'] = result
        return axes_dict

    @staticmethod
    def _extract_axis(self, data, axis=0, intersect=False):

        index = None
        if len(data) == 0:
            index = Index([])
        elif len(data) > 0:
            raw_lengths = []
            indexes = []

        have_raw_arrays = False
        have_frames = False

        for v in data.values():
            if isinstance(v, self._constructor_sliced):
                have_frames = True
                indexes.append(v._get_axis(axis))
            elif v is not None:
                have_raw_arrays = True
                raw_lengths.append(v.shape[axis])

        if have_frames:
            index = _get_combined_index(indexes, intersect=intersect)

        if have_raw_arrays:
            lengths = list(set(raw_lengths))
            if len(lengths) > 1:
                raise ValueError('ndarrays must match shape on axis %d' % axis)

            if have_frames:
                if lengths[0] != len(index):
                    raise AssertionError('Length of data and index must match')
            else:
                index = Index(np.arange(lengths[0]))

        if index is None:
            index = Index([])

        return _ensure_index(index)

    @classmethod
    def _add_aggregate_operations(cls, use_numexpr=True):
        """ add the operations to the cls; evaluate the doc strings again """

        # doc strings substitors
        _agg_doc = """
Wrapper method for %%s

Parameters
----------
other : %s or %s""" % (cls._constructor_sliced.__name__, cls.__name__) + """
axis : {""" + ', '.join(cls._AXIS_ORDERS) + "}" + """
Axis to broadcast over

Returns
-------
""" + cls.__name__ + "\n"

        def _panel_arith_method(op, name, str_rep=None, default_axis=None,
                                fill_zeros=None, **eval_kwargs):
            def na_op(x, y):
                try:
                    result = expressions.evaluate(op, str_rep, x, y,
                                                  raise_on_error=True,
                                                  **eval_kwargs)
                except TypeError:
                    result = op(x, y)

                # handles discrepancy between numpy and numexpr on division/mod
                # by 0 though, given that these are generally (always?)
                # non-scalars, I'm not sure whether it's worth it at the moment
                result = com._fill_zeros(result, y, fill_zeros)
                return result

            @Substitution(name)
            @Appender(_agg_doc)
            def f(self, other, axis=0):
                return self._combine(other, na_op, axis=axis)
            f.__name__ = name
            return f

        # add `div`, `mul`, `pow`, etc..
        ops.add_flex_arithmetic_methods(
            cls, _panel_arith_method, use_numexpr=use_numexpr,
            flex_comp_method=ops._comp_method_PANEL)

Panel._setup_axes(axes=['items', 'major_axis', 'minor_axis'],
                  info_axis=0,
                  stat_axis=1,
                  aliases={'major': 'major_axis',
                           'minor': 'minor_axis'},
                  slicers={'major_axis': 'index',
                           'minor_axis': 'columns'})

ops.add_special_arithmetic_methods(Panel, **ops.panel_special_funcs)
Panel._add_aggregate_operations()
Panel._add_numeric_operations()

WidePanel = Panel
LongPanel = DataFrame
