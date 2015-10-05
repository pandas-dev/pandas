# pylint: disable=W0231,E1101
import warnings
import operator
import weakref
import gc

import numpy as np
import pandas.lib as lib

import pandas as pd
from pandas.core.base import PandasObject
from pandas.core.index import (Index, MultiIndex, _ensure_index,
                               InvalidIndexError)
import pandas.core.indexing as indexing
from pandas.tseries.index import DatetimeIndex
from pandas.tseries.period import PeriodIndex
from pandas.core.internals import BlockManager
import pandas.core.common as com
import pandas.core.datetools as datetools
from pandas import compat
from pandas.compat import map, zip, lrange, string_types, isidentifier
from pandas.core.common import (isnull, notnull, is_list_like,
                                _values_from_object, _maybe_promote,
                                _maybe_box_datetimelike, ABCSeries,
                                SettingWithCopyError, SettingWithCopyWarning,
                                AbstractMethodError)
import pandas.core.nanops as nanops
from pandas.util.decorators import Appender, Substitution, deprecate_kwarg
from pandas.core import config


# goal is to be able to define the docs close to function, while still being
# able to share
_shared_docs = dict()
_shared_doc_kwargs = dict(axes='keywords for axes',
                          klass='NDFrame',
                          axes_single_arg='int or labels for object',
                          args_transpose='axes to permute (int or label for'
                                         ' object)')


def is_dictlike(x):
    return isinstance(x, (dict, com.ABCSeries))


def _single_replace(self, to_replace, method, inplace, limit):
    if self.ndim != 1:
        raise TypeError('cannot replace {0} with method {1} on a {2}'
                        .format(to_replace, method, type(self).__name__))

    orig_dtype = self.dtype
    result = self if inplace else self.copy()
    fill_f = com._get_fill_func(method)

    mask = com.mask_missing(result.values, to_replace)
    values = fill_f(result.values, limit=limit, mask=mask)

    if values.dtype == orig_dtype and inplace:
        return

    result = pd.Series(values, index=self.index,
                       dtype=self.dtype).__finalize__(self)

    if inplace:
        self._update_inplace(result._data)
        return

    return result


class NDFrame(PandasObject):

    """
    N-dimensional analogue of DataFrame. Store multi-dimensional in a
    size-mutable, labeled data structure

    Parameters
    ----------
    data : BlockManager
    axes : list
    copy : boolean, default False
    """
    _internal_names = ['_data', '_cacher', '_item_cache', '_cache',
                       'is_copy', '_subtyp', '_index',
                       '_default_kind', '_default_fill_value', '_metadata',
                       '__array_struct__', '__array_interface__']
    _internal_names_set = set(_internal_names)
    _accessors = frozenset([])
    _metadata = []
    is_copy = None

    def __init__(self, data, axes=None, copy=False, dtype=None,
                 fastpath=False):

        if not fastpath:
            if dtype is not None:
                data = data.astype(dtype)
            elif copy:
                data = data.copy()

            if axes is not None:
                for i, ax in enumerate(axes):
                    data = data.reindex_axis(ax, axis=i)

        object.__setattr__(self, 'is_copy', None)
        object.__setattr__(self, '_data', data)
        object.__setattr__(self, '_item_cache', {})

    def _validate_dtype(self, dtype):
        """ validate the passed dtype """

        if dtype is not None:
            dtype = com._coerce_to_dtype(dtype)

            # a compound dtype
            if dtype.kind == 'V':
                raise NotImplementedError("compound dtypes are not implemented"
                                          "in the {0} constructor"
                                          .format(self.__class__.__name__))
        return dtype

    def _init_mgr(self, mgr, axes=None, dtype=None, copy=False):
        """ passed a manager and a axes dict """
        for a, axe in axes.items():
            if axe is not None:
                mgr = mgr.reindex_axis(
                    axe, axis=self._get_block_manager_axis(a), copy=False)

        # make a copy if explicitly requested
        if copy:
            mgr = mgr.copy()
        if dtype is not None:
            # avoid further copies if we can
            if len(mgr.blocks) > 1 or mgr.blocks[0].values.dtype != dtype:
                mgr = mgr.astype(dtype=dtype)
        return mgr

    #----------------------------------------------------------------------
    # Construction

    @property
    def _constructor(self):
        """Used when a manipulation result has the same dimesions as the
        original.
        """
        raise AbstractMethodError(self)

    def __unicode__(self):
        # unicode representation based upon iterating over self
        # (since, by definition, `PandasContainers` are iterable)
        prepr = '[%s]' % ','.join(map(com.pprint_thing, self))
        return '%s(%s)' % (self.__class__.__name__, prepr)

    def _dir_additions(self):
        """ add the string-like attributes from the info_axis """
        return set([c for c in self._info_axis
                if isinstance(c, string_types) and isidentifier(c)])

    @property
    def _constructor_sliced(self):
        """Used when a manipulation result has one lower dimension(s) as the
        original, such as DataFrame single columns slicing.
        """
        raise AbstractMethodError(self)

    @property
    def _constructor_expanddim(self):
        """Used when a manipulation result has one higher dimension as the
        original, such as Series.to_frame() and DataFrame.to_panel()
        """
        raise NotImplementedError

    #----------------------------------------------------------------------
    # Axis

    @classmethod
    def _setup_axes(
        cls, axes, info_axis=None, stat_axis=None, aliases=None, slicers=None,
            axes_are_reversed=False, build_axes=True, ns=None):
        """ provide axes setup for the major PandasObjects

            axes : the names of the axes in order (lowest to highest)
            info_axis_num : the axis of the selector dimension (int)
            stat_axis_num : the number of axis for the default stats (int)
            aliases : other names for a single axis (dict)
            slicers : how axes slice to others (dict)
            axes_are_reversed : boolean whether to treat passed axes as
                reversed (DataFrame)
            build_axes : setup the axis properties (default True)
            """

        cls._AXIS_ORDERS = axes
        cls._AXIS_NUMBERS = dict((a, i) for i, a in enumerate(axes))
        cls._AXIS_LEN = len(axes)
        cls._AXIS_ALIASES = aliases or dict()
        cls._AXIS_IALIASES = dict((v, k)
                                  for k, v in cls._AXIS_ALIASES.items())
        cls._AXIS_NAMES = dict(enumerate(axes))
        cls._AXIS_SLICEMAP = slicers or None
        cls._AXIS_REVERSED = axes_are_reversed

        # typ
        setattr(cls, '_typ', cls.__name__.lower())

        # indexing support
        cls._ix = None

        if info_axis is not None:
            cls._info_axis_number = info_axis
            cls._info_axis_name = axes[info_axis]

        if stat_axis is not None:
            cls._stat_axis_number = stat_axis
            cls._stat_axis_name = axes[stat_axis]

        # setup the actual axis
        if build_axes:

            def set_axis(a, i):
                setattr(cls, a, lib.AxisProperty(i))
                cls._internal_names_set.add(a)

            if axes_are_reversed:
                m = cls._AXIS_LEN - 1
                for i, a in cls._AXIS_NAMES.items():
                    set_axis(a, m - i)
            else:
                for i, a in cls._AXIS_NAMES.items():
                    set_axis(a, i)

        # addtl parms
        if isinstance(ns, dict):
            for k, v in ns.items():
                setattr(cls, k, v)

    def _construct_axes_dict(self, axes=None, **kwargs):
        """ return an axes dictionary for myself """
        d = dict([(a, self._get_axis(a)) for a in (axes or self._AXIS_ORDERS)])
        d.update(kwargs)
        return d

    @staticmethod
    def _construct_axes_dict_from(self, axes, **kwargs):
        """ return an axes dictionary for the passed axes """
        d = dict([(a, ax) for a, ax in zip(self._AXIS_ORDERS, axes)])
        d.update(kwargs)
        return d

    def _construct_axes_dict_for_slice(self, axes=None, **kwargs):
        """ return an axes dictionary for myself """
        d = dict([(self._AXIS_SLICEMAP[a], self._get_axis(a))
                 for a in (axes or self._AXIS_ORDERS)])
        d.update(kwargs)
        return d

    def _construct_axes_from_arguments(self, args, kwargs, require_all=False):
        """ construct and returns axes if supplied in args/kwargs
            if require_all, raise if all axis arguments are not supplied
            return a tuple of (axes, kwargs) """

        # construct the args
        args = list(args)
        for a in self._AXIS_ORDERS:

            # if we have an alias for this axis
            alias = self._AXIS_IALIASES.get(a)
            if alias is not None:
                if a in kwargs:
                    if alias in kwargs:
                        raise TypeError(
                            "arguments are mutually exclusive for [%s,%s]" %
                            (a, alias)
                        )
                    continue
                if alias in kwargs:
                    kwargs[a] = kwargs.pop(alias)
                    continue

            # look for a argument by position
            if a not in kwargs:
                try:
                    kwargs[a] = args.pop(0)
                except (IndexError):
                    if require_all:
                        raise TypeError(
                            "not enough/duplicate arguments specified!")

        axes = dict([(a, kwargs.pop(a, None)) for a in self._AXIS_ORDERS])
        return axes, kwargs

    @classmethod
    def _from_axes(cls, data, axes, **kwargs):
        # for construction from BlockManager
        if isinstance(data, BlockManager):
            return cls(data, **kwargs)
        else:
            if cls._AXIS_REVERSED:
                axes = axes[::-1]
            d = cls._construct_axes_dict_from(cls, axes, copy=False)
            d.update(kwargs)
            return cls(data, **d)

    def _get_axis_number(self, axis):
        axis = self._AXIS_ALIASES.get(axis, axis)
        if com.is_integer(axis):
            if axis in self._AXIS_NAMES:
                return axis
        else:
            try:
                return self._AXIS_NUMBERS[axis]
            except:
                pass
        raise ValueError('No axis named {0} for object type {1}'
                         .format(axis, type(self)))

    def _get_axis_name(self, axis):
        axis = self._AXIS_ALIASES.get(axis, axis)
        if isinstance(axis, string_types):
            if axis in self._AXIS_NUMBERS:
                return axis
        else:
            try:
                return self._AXIS_NAMES[axis]
            except:
                pass
        raise ValueError('No axis named {0} for object type {1}'
                         .format(axis, type(self)))

    def _get_axis(self, axis):
        name = self._get_axis_name(axis)
        return getattr(self, name)

    def _get_block_manager_axis(self, axis):
        """ map the axis to the block_manager axis """
        axis = self._get_axis_number(axis)
        if self._AXIS_REVERSED:
            m = self._AXIS_LEN - 1
            return m - axis
        return axis

    def _get_axis_resolvers(self, axis):
        # index or columns
        axis_index = getattr(self, axis)
        d = dict()
        prefix = axis[0]

        for i, name in enumerate(axis_index.names):
            if name is not None:
                key = level = name
            else:
                # prefix with 'i' or 'c' depending on the input axis
                # e.g., you must do ilevel_0 for the 0th level of an unnamed
                # multiiindex
                key = '{prefix}level_{i}'.format(prefix=prefix, i=i)
                level = i

            level_values = axis_index.get_level_values(level)
            s = level_values.to_series()
            s.index = axis_index
            d[key] = s

        # put the index/columns itself in the dict
        if isinstance(axis_index, MultiIndex):
            dindex = axis_index
        else:
            dindex = axis_index.to_series()

        d[axis] = dindex
        return d

    def _get_index_resolvers(self):
        d = {}
        for axis_name in self._AXIS_ORDERS:
            d.update(self._get_axis_resolvers(axis_name))
        return d

    @property
    def _info_axis(self):
        return getattr(self, self._info_axis_name)

    @property
    def _stat_axis(self):
        return getattr(self, self._stat_axis_name)

    @property
    def shape(self):
        "Return a tuple of axis dimensions"
        return tuple(len(self._get_axis(a)) for a in self._AXIS_ORDERS)

    @property
    def axes(self):
        "Return index label(s) of the internal NDFrame"
        # we do it this way because if we have reversed axes, then
        # the block manager shows then reversed
        return [self._get_axis(a) for a in self._AXIS_ORDERS]

    @property
    def ndim(self):
        "Number of axes / array dimensions"
        return self._data.ndim

    @property
    def size(self):
        "number of elements in the NDFrame"
        return np.prod(self.shape)

    def _expand_axes(self, key):
        new_axes = []
        for k, ax in zip(key, self.axes):
            if k not in ax:
                if type(k) != ax.dtype.type:
                    ax = ax.astype('O')
                new_axes.append(ax.insert(len(ax), k))
            else:
                new_axes.append(ax)

        return new_axes

    def set_axis(self, axis, labels):
        """ public verson of axis assignment """
        setattr(self,self._get_axis_name(axis),labels)

    def _set_axis(self, axis, labels):
        self._data.set_axis(axis, labels)
        self._clear_item_cache()

    _shared_docs['transpose'] = """
        Permute the dimensions of the %(klass)s

        Parameters
        ----------
        args : %(args_transpose)s
        copy : boolean, default False
            Make a copy of the underlying data. Mixed-dtype data will
            always result in a copy

        Examples
        --------
        >>> p.transpose(2, 0, 1)
        >>> p.transpose(2, 0, 1, copy=True)

        Returns
        -------
        y : same as input
        """

    @Appender(_shared_docs['transpose'] % _shared_doc_kwargs)
    def transpose(self, *args, **kwargs):

        # construct the args
        axes, kwargs = self._construct_axes_from_arguments(
            args, kwargs, require_all=True)
        axes_names = tuple([self._get_axis_name(axes[a])
                            for a in self._AXIS_ORDERS])
        axes_numbers = tuple([self._get_axis_number(axes[a])
                             for a in self._AXIS_ORDERS])

        # we must have unique axes
        if len(axes) != len(set(axes)):
            raise ValueError('Must specify %s unique axes' % self._AXIS_LEN)

        new_axes = self._construct_axes_dict_from(
            self, [self._get_axis(x) for x in axes_names])
        new_values = self.values.transpose(axes_numbers)
        if kwargs.pop('copy', None) or (len(args) and args[-1]):
            new_values = new_values.copy()

        if kwargs:
            raise TypeError('transpose() got an unexpected keyword '
                    'argument "{0}"'.format(list(kwargs.keys())[0]))

        return self._constructor(new_values, **new_axes).__finalize__(self)

    def swapaxes(self, axis1, axis2, copy=True):
        """
        Interchange axes and swap values axes appropriately

        Returns
        -------
        y : same as input
        """
        i = self._get_axis_number(axis1)
        j = self._get_axis_number(axis2)

        if i == j:
            if copy:
                return self.copy()
            return self

        mapping = {i: j, j: i}

        new_axes = (self._get_axis(mapping.get(k, k))
                    for k in range(self._AXIS_LEN))
        new_values = self.values.swapaxes(i, j)
        if copy:
            new_values = new_values.copy()

        return self._constructor(new_values, *new_axes).__finalize__(self)

    def pop(self, item):
        """
        Return item and drop from frame. Raise KeyError if not found.
        """
        result = self[item]
        del self[item]
        try:
            result._reset_cacher()
        except AttributeError:
            pass

        return result

    def squeeze(self):
        """ squeeze length 1 dimensions """
        try:
            return self.ix[tuple([slice(None) if len(a) > 1 else a[0]
                                  for a in self.axes])]
        except:
            return self

    def swaplevel(self, i, j, axis=0):
        """
        Swap levels i and j in a MultiIndex on a particular axis

        Parameters
        ----------
        i, j : int, string (can be mixed)
            Level of index to be swapped. Can pass level name as string.

        Returns
        -------
        swapped : type of caller (new object)
        """
        axis = self._get_axis_number(axis)
        result = self.copy()
        labels = result._data.axes[axis]
        result._data.set_axis(axis, labels.swaplevel(i, j))
        return result

    #----------------------------------------------------------------------
    # Rename

    # TODO: define separate funcs for DataFrame, Series and Panel so you can
    # get completion on keyword arguments.
    _shared_docs['rename'] = """
        Alter axes input function or functions. Function / dict values must be
        unique (1-to-1). Labels not contained in a dict / Series will be left
        as-is.

        Parameters
        ----------
        %(axes)s : dict-like or function, optional
            Transformation to apply to that axis values

        copy : boolean, default True
            Also copy underlying data
        inplace : boolean, default False
            Whether to return a new %(klass)s. If True then value of copy is
            ignored.

        Returns
        -------
        renamed : %(klass)s (new object)
        """

    @Appender(_shared_docs['rename'] % dict(axes='axes keywords for this'
                                            ' object', klass='NDFrame'))
    def rename(self, *args, **kwargs):

        axes, kwargs = self._construct_axes_from_arguments(args, kwargs)
        copy = kwargs.pop('copy', True)
        inplace = kwargs.pop('inplace', False)

        if kwargs:
            raise TypeError('rename() got an unexpected keyword '
                    'argument "{0}"'.format(list(kwargs.keys())[0]))

        if (com._count_not_none(*axes.values()) == 0):
            raise TypeError('must pass an index to rename')

        # renamer function if passed a dict
        def _get_rename_function(mapper):
            if isinstance(mapper, (dict, ABCSeries)):
                def f(x):
                    if x in mapper:
                        return mapper[x]
                    else:
                        return x
            else:
                f = mapper

            return f

        self._consolidate_inplace()
        result = self if inplace else self.copy(deep=copy)

        # start in the axis order to eliminate too many copies
        for axis in lrange(self._AXIS_LEN):
            v = axes.get(self._AXIS_NAMES[axis])
            if v is None:
                continue
            f = _get_rename_function(v)

            baxis = self._get_block_manager_axis(axis)
            result._data = result._data.rename_axis(f, axis=baxis, copy=copy)
            result._clear_item_cache()

        if inplace:
            self._update_inplace(result._data)
        else:
            return result.__finalize__(self)

    rename.__doc__ = _shared_docs['rename']

    def rename_axis(self, mapper, axis=0, copy=True, inplace=False):
        """
        Alter index and / or columns using input function or functions.
        Function / dict values must be unique (1-to-1). Labels not contained in
        a dict / Series will be left as-is.

        Parameters
        ----------
        mapper : dict-like or function, optional
        axis : int or string, default 0
        copy : boolean, default True
            Also copy underlying data
        inplace : boolean, default False

        Returns
        -------
        renamed : type of caller
        """
        axis = self._get_axis_name(axis)
        d = {'copy': copy, 'inplace': inplace}
        d[axis] = mapper
        return self.rename(**d)

    #----------------------------------------------------------------------
    # Comparisons

    def _indexed_same(self, other):
        return all([self._get_axis(a).equals(other._get_axis(a))
                    for a in self._AXIS_ORDERS])

    def __neg__(self):
        values = _values_from_object(self)
        if values.dtype == np.bool_:
            arr = operator.inv(values)
        else:
            arr = operator.neg(values)
        return self.__array_wrap__(arr)

    def __invert__(self):
        try:
            arr = operator.inv(_values_from_object(self))
            return self.__array_wrap__(arr)
        except:

            # inv fails with 0 len
            if not np.prod(self.shape):
                return self

            raise

    def equals(self, other):
        """
        Determines if two NDFrame objects contain the same elements. NaNs in the
        same location are considered equal.
        """
        if not isinstance(other, self._constructor):
            return False
        return self._data.equals(other._data)

    #----------------------------------------------------------------------
    # Iteration

    def __hash__(self):
        raise TypeError('{0!r} objects are mutable, thus they cannot be'
                        ' hashed'.format(self.__class__.__name__))

    def __iter__(self):
        """
        Iterate over infor axis
        """
        return iter(self._info_axis)

    # can we get a better explanation of this?
    def keys(self):
        """Get the 'info axis' (see Indexing for more)

        This is index for Series, columns for DataFrame and major_axis for
        Panel."""
        return self._info_axis

    def iteritems(self):
        """Iterate over (label, values) on info axis

        This is index for Series, columns for DataFrame, major_axis for Panel,
        and so on.
        """
        for h in self._info_axis:
            yield h, self[h]

    # originally used to get around 2to3's changes to iteritems.
    # Now unnecessary. Sidenote: don't want to deprecate this for a while,
    # otherwise libraries that use 2to3 will have issues.
    def iterkv(self, *args, **kwargs):
        "iteritems alias used to get around 2to3. Deprecated"
        warnings.warn("iterkv is deprecated and will be removed in a future "
                      "release, use ``iteritems`` instead.",
                      FutureWarning, stacklevel=2)
        return self.iteritems(*args, **kwargs)

    def __len__(self):
        """Returns length of info axis """
        return len(self._info_axis)

    def __contains__(self, key):
        """True if the key is in the info axis """
        return key in self._info_axis

    @property
    def empty(self):
        "True if NDFrame is entirely empty [no items]"
        return not all(len(self._get_axis(a)) > 0 for a in self._AXIS_ORDERS)

    def __nonzero__(self):
        raise ValueError("The truth value of a {0} is ambiguous. "
                         "Use a.empty, a.bool(), a.item(), a.any() or a.all()."
                         .format(self.__class__.__name__))

    __bool__ = __nonzero__

    def bool(self):
        """ Return the bool of a single element PandasObject
            This must be a boolean scalar value, either True or False

            Raise a ValueError if the PandasObject does not have exactly
            1 element, or that element is not boolean """
        v = self.squeeze()
        if isinstance(v, (bool, np.bool_)):
            return bool(v)
        elif np.isscalar(v):
            raise ValueError("bool cannot act on a non-boolean single element "
                             "{0}".format(self.__class__.__name__))

        self.__nonzero__()

    def __abs__(self):
        return self.abs()

    #----------------------------------------------------------------------
    # Array Interface

    def __array__(self, dtype=None):
        return _values_from_object(self)

    def __array_wrap__(self, result, context=None):
        d = self._construct_axes_dict(self._AXIS_ORDERS, copy=False)
        return self._constructor(result, **d).__finalize__(self)

    # ideally we would define this to avoid the getattr checks, but
    # is slower
    #@property
    #def __array_interface__(self):
    #    """ provide numpy array interface method """
    #    values = self.values
    #    return dict(typestr=values.dtype.str,shape=values.shape,data=values)

    def to_dense(self):
        "Return dense representation of NDFrame (as opposed to sparse)"
        # compat
        return self

    #----------------------------------------------------------------------
    # Picklability

    def __getstate__(self):
        meta = dict((k, getattr(self, k, None)) for k in self._metadata)
        return dict(_data=self._data, _typ=self._typ,
                    _metadata=self._metadata, **meta)

    def __setstate__(self, state):

        if isinstance(state, BlockManager):
            self._data = state
        elif isinstance(state, dict):
            typ = state.get('_typ')
            if typ is not None:

                # set in the order of internal names
                # to avoid definitional recursion
                # e.g. say fill_value needing _data to be
                # defined
                meta = set(self._internal_names + self._metadata)
                for k in list(meta):
                    if k in state:
                        v = state[k]
                        object.__setattr__(self, k, v)

                for k, v in state.items():
                    if k not in meta:
                        object.__setattr__(self, k, v)

            else:
                self._unpickle_series_compat(state)
        elif isinstance(state[0], dict):
            if len(state) == 5:
                self._unpickle_sparse_frame_compat(state)
            else:
                self._unpickle_frame_compat(state)
        elif len(state) == 4:
            self._unpickle_panel_compat(state)
        elif len(state) == 2:
            self._unpickle_series_compat(state)
        else:  # pragma: no cover
            # old pickling format, for compatibility
            self._unpickle_matrix_compat(state)

        self._item_cache = {}

    #----------------------------------------------------------------------
    # IO

    #----------------------------------------------------------------------
    # I/O Methods

    def to_json(self, path_or_buf=None, orient=None, date_format='epoch',
                double_precision=10, force_ascii=True, date_unit='ms',
                default_handler=None):
        """
        Convert the object to a JSON string.

        Note NaN's and None will be converted to null and datetime objects
        will be converted to UNIX timestamps.

        Parameters
        ----------
        path_or_buf : the path or buffer to write the result string
            if this is None, return a StringIO of the converted string
        orient : string

            * Series

              - default is 'index'
              - allowed values are: {'split','records','index'}

            * DataFrame

              - default is 'columns'
              - allowed values are:
                {'split','records','index','columns','values'}

            * The format of the JSON string

              - split : dict like
                {index -> [index], columns -> [columns], data -> [values]}
              - records : list like
                [{column -> value}, ... , {column -> value}]
              - index : dict like {index -> {column -> value}}
              - columns : dict like {column -> {index -> value}}
              - values : just the values array

        date_format : {'epoch', 'iso'}
            Type of date conversion. `epoch` = epoch milliseconds,
            `iso`` = ISO8601, default is epoch.
        double_precision : The number of decimal places to use when encoding
            floating point values, default 10.
        force_ascii : force encoded string to be ASCII, default True.
        date_unit : string, default 'ms' (milliseconds)
            The time unit to encode to, governs timestamp and ISO8601
            precision.  One of 's', 'ms', 'us', 'ns' for second, millisecond,
            microsecond, and nanosecond respectively.
        default_handler : callable, default None
            Handler to call if object cannot otherwise be converted to a
            suitable format for JSON. Should receive a single argument which is
            the object to convert and return a serialisable object.

        Returns
        -------
        same type as input object with filtered info axis

        """

        from pandas.io import json
        return json.to_json(
            path_or_buf=path_or_buf,
            obj=self, orient=orient,
            date_format=date_format,
            double_precision=double_precision,
            force_ascii=force_ascii,
            date_unit=date_unit,
            default_handler=default_handler)

    def to_hdf(self, path_or_buf, key, **kwargs):
        """ activate the HDFStore

        Parameters
        ----------
        path_or_buf : the path (string) or HDFStore object
        key : string
            indentifier for the group in the store
        mode : optional, {'a', 'w', 'r', 'r+'}, default 'a'

          ``'r'``
              Read-only; no data can be modified.
          ``'w'``
              Write; a new file is created (an existing file with the same
              name would be deleted).
          ``'a'``
              Append; an existing file is opened for reading and writing,
              and if the file does not exist it is created.
          ``'r+'``
              It is similar to ``'a'``, but the file must already exist.
        format   : 'fixed(f)|table(t)', default is 'fixed'
            fixed(f) : Fixed format
                       Fast writing/reading. Not-appendable, nor searchable
            table(t) : Table format
                       Write as a PyTables Table structure which may perform
                       worse but allow more flexible operations like searching
                       / selecting subsets of the data
        append   : boolean, default False
            For Table formats, append the input data to the existing
        complevel : int, 1-9, default 0
            If a complib is specified compression will be applied
            where possible
        complib : {'zlib', 'bzip2', 'lzo', 'blosc', None}, default None
            If complevel is > 0 apply compression to objects written
            in the store wherever possible
        fletcher32 : bool, default False
            If applying compression use the fletcher32 checksum
        dropna : boolean, default False.
            If true, ALL nan rows will not be written to store.

        """

        from pandas.io import pytables
        return pytables.to_hdf(path_or_buf, key, self, **kwargs)

    def to_msgpack(self, path_or_buf=None, **kwargs):
        """
        msgpack (serialize) object to input file path

        THIS IS AN EXPERIMENTAL LIBRARY and the storage format
        may not be stable until a future release.

        Parameters
        ----------
        path : string File path, buffer-like, or None
            if None, return generated string
        append : boolean whether to append to an existing msgpack
            (default is False)
        compress : type of compressor (zlib or blosc), default to None (no
            compression)
        """

        from pandas.io import packers
        return packers.to_msgpack(path_or_buf, self, **kwargs)

    def to_sql(self, name, con, flavor='sqlite', schema=None, if_exists='fail',
               index=True, index_label=None, chunksize=None, dtype=None):
        """
        Write records stored in a DataFrame to a SQL database.

        Parameters
        ----------
        name : string
            Name of SQL table
        con : SQLAlchemy engine or DBAPI2 connection (legacy mode)
            Using SQLAlchemy makes it possible to use any DB supported by that
            library.
            If a DBAPI2 object, only sqlite3 is supported.
        flavor : {'sqlite', 'mysql'}, default 'sqlite'
            The flavor of SQL to use. Ignored when using SQLAlchemy engine.
            'mysql' is deprecated and will be removed in future versions, but it
            will be further supported through SQLAlchemy engines.
        schema : string, default None
            Specify the schema (if database flavor supports this). If None, use
            default schema.
        if_exists : {'fail', 'replace', 'append'}, default 'fail'
            - fail: If table exists, do nothing.
            - replace: If table exists, drop it, recreate it, and insert data.
            - append: If table exists, insert data. Create if does not exist.
        index : boolean, default True
            Write DataFrame index as a column.
        index_label : string or sequence, default None
            Column label for index column(s). If None is given (default) and
            `index` is True, then the index names are used.
            A sequence should be given if the DataFrame uses MultiIndex.
        chunksize : int, default None
            If not None, then rows will be written in batches of this size at a
            time.  If None, all rows will be written at once.
        dtype : dict of column name to SQL type, default None
            Optional specifying the datatype for columns. The SQL type should
            be a SQLAlchemy type, or a string for sqlite3 fallback connection.

        """
        from pandas.io import sql
        sql.to_sql(
            self, name, con, flavor=flavor, schema=schema, if_exists=if_exists,
            index=index, index_label=index_label, chunksize=chunksize,
            dtype=dtype)

    def to_pickle(self, path):
        """
        Pickle (serialize) object to input file path

        Parameters
        ----------
        path : string
            File path
        """
        from pandas.io.pickle import to_pickle
        return to_pickle(self, path)

    def to_clipboard(self, excel=None, sep=None, **kwargs):
        """
        Attempt to write text representation of object to the system clipboard
        This can be pasted into Excel, for example.

        Parameters
        ----------
        excel : boolean, defaults to True
                if True, use the provided separator, writing in a csv
                format for allowing easy pasting into excel.
                if False, write a string representation of the object
                to the clipboard
        sep : optional, defaults to tab
        other keywords are passed to to_csv

        Notes
        -----
        Requirements for your platform
          - Linux: xclip, or xsel (with gtk or PyQt4 modules)
          - Windows: none
          - OS X: none
        """
        from pandas.io import clipboard
        clipboard.to_clipboard(self, excel=excel, sep=sep, **kwargs)

    #----------------------------------------------------------------------
    # Fancy Indexing

    @classmethod
    def _create_indexer(cls, name, indexer):
        """ create an indexer like _name in the class """

        if getattr(cls, name, None) is None:
            iname = '_%s' % name
            setattr(cls, iname, None)

            def _indexer(self):
                i = getattr(self, iname)
                if i is None:
                    i = indexer(self, name)
                    setattr(self, iname, i)
                return i

            setattr(cls, name, property(_indexer, doc=indexer.__doc__))

            # add to our internal names set
            cls._internal_names_set.add(iname)

    def get(self, key, default=None):
        """
        Get item from object for given key (DataFrame column, Panel slice,
        etc.). Returns default value if not found

        Parameters
        ----------
        key : object

        Returns
        -------
        value : type of items contained in object
        """
        try:
            return self[key]
        except (KeyError, ValueError, IndexError):
            return default

    def __getitem__(self, item):
        return self._get_item_cache(item)

    def _get_item_cache(self, item):
        """ return the cached item, item represents a label indexer """
        cache = self._item_cache
        res = cache.get(item)
        if res is None:
            values = self._data.get(item)
            res = self._box_item_values(item, values)
            cache[item] = res
            res._set_as_cached(item, self)

            # for a chain
            res.is_copy = self.is_copy
        return res

    def _set_as_cached(self, item, cacher):
        """ set the _cacher attribute on the calling object with
            a weakref to cacher """
        self._cacher = (item, weakref.ref(cacher))

    def _reset_cacher(self):
        """ reset the cacher """
        if hasattr(self,'_cacher'):
            del self._cacher

    def _iget_item_cache(self, item):
        """ return the cached item, item represents a positional indexer """
        ax = self._info_axis
        if ax.is_unique:
            lower = self._get_item_cache(ax[item])
        else:
            lower = self.take(item, axis=self._info_axis_number, convert=True)
        return lower

    def _box_item_values(self, key, values):
        raise AbstractMethodError(self)

    def _maybe_cache_changed(self, item, value):
        """
        the object has called back to us saying
        maybe it has changed

        numpy < 1.8 has an issue with object arrays and aliasing
        GH6026
        """
        self._data.set(item, value, check=pd._np_version_under1p8)

    @property
    def _is_cached(self):
        """ boolean : return if I am cached """
        return getattr(self, '_cacher', None) is not None

    def _get_cacher(self):
        """ return my cacher or None """
        cacher = getattr(self, '_cacher', None)
        if cacher is not None:
            cacher = cacher[1]()
        return cacher

    @property
    def _is_view(self):
        """ boolean : return if I am a view of another array """
        return self._data.is_view

    def _maybe_update_cacher(self, clear=False, verify_is_copy=True):
        """

        see if we need to update our parent cacher
        if clear, then clear our cache

        Parameters
        ----------
        clear : boolean, default False
            clear the item cache
        verify_is_copy : boolean, default True
            provide is_copy checks

        """

        cacher = getattr(self, '_cacher', None)
        if cacher is not None:
            ref = cacher[1]()

            # we are trying to reference a dead referant, hence
            # a copy
            if ref is None:
                del self._cacher
            else:
                try:
                    ref._maybe_cache_changed(cacher[0], self)
                except:
                    pass

        if verify_is_copy:
            self._check_setitem_copy(stacklevel=5, t='referant')

        if clear:
            self._clear_item_cache()

    def _clear_item_cache(self, i=None):
        if i is not None:
            self._item_cache.pop(i, None)
        else:
            self._item_cache.clear()

    def _slice(self, slobj, axis=0, kind=None):
        """
        Construct a slice of this container.

        kind parameter is maintained for compatibility with Series slicing.

        """
        axis = self._get_block_manager_axis(axis)
        result = self._constructor(self._data.get_slice(slobj, axis=axis))
        result = result.__finalize__(self)

        # this could be a view
        # but only in a single-dtyped view slicable case
        is_copy = axis!=0 or result._is_view
        result._set_is_copy(self, copy=is_copy)
        return result

    def _set_item(self, key, value):
        self._data.set(key, value)
        self._clear_item_cache()

    def _set_is_copy(self, ref=None, copy=True):
        if not copy:
            self.is_copy = None
        else:
            if ref is not None:
                self.is_copy = weakref.ref(ref)
            else:
                self.is_copy = None

    def _check_is_chained_assignment_possible(self):
        """
        check if we are a view, have a cacher, and are of mixed type
        if so, then force a setitem_copy check

        should be called just near setting a value

        will return a boolean if it we are a view and are cached, but a single-dtype
        meaning that the cacher should be updated following setting
        """
        if self._is_view and self._is_cached:
            ref = self._get_cacher()
            if ref is not None and ref._is_mixed_type:
                self._check_setitem_copy(stacklevel=4, t='referant', force=True)
            return True
        elif self.is_copy:
            self._check_setitem_copy(stacklevel=4, t='referant')
        return False

    def _check_setitem_copy(self, stacklevel=4, t='setting', force=False):
        """

        Parameters
        ----------
        stacklevel : integer, default 4
           the level to show of the stack when the error is output
        t : string, the type of setting error
        force : boolean, default False
           if True, then force showing an error

        validate if we are doing a settitem on a chained copy.

        If you call this function, be sure to set the stacklevel such that the
        user will see the error *at the level of setting*

        It is technically possible to figure out that we are setting on
        a copy even WITH a multi-dtyped pandas object. In other words, some blocks
        may be views while other are not. Currently _is_view will ALWAYS return False
        for multi-blocks to avoid having to handle this case.

        df = DataFrame(np.arange(0,9), columns=['count'])
        df['group'] = 'b'

        # this technically need not raise SettingWithCopy if both are view (which is not
        # generally guaranteed but is usually True
        # however, this is in general not a good practice and we recommend using .loc
        df.iloc[0:5]['group'] = 'a'

        """

        if force or self.is_copy:

            value = config.get_option('mode.chained_assignment')
            if value is None:
                return

            # see if the copy is not actually refererd; if so, then disolve
            # the copy weakref
            try:
                gc.collect(2)
                if not gc.get_referents(self.is_copy()):
                    self.is_copy = None
                    return
            except:
                pass

            # we might be a false positive
            try:
                if self.is_copy().shape == self.shape:
                    self.is_copy = None
                    return
            except:
                pass

            # a custom message
            if isinstance(self.is_copy, string_types):
                t = self.is_copy

            elif t == 'referant':
                t = ("\n"
                     "A value is trying to be set on a copy of a slice from a "
                     "DataFrame\n\n"
                     "See the caveats in the documentation: "
                     "http://pandas.pydata.org/pandas-docs/stable/indexing.html#indexing-view-versus-copy")

            else:
                t = ("\n"
                     "A value is trying to be set on a copy of a slice from a "
                     "DataFrame.\n"
                     "Try using .loc[row_indexer,col_indexer] = value instead\n\n"
                     "See the caveats in the documentation: "
                     "http://pandas.pydata.org/pandas-docs/stable/indexing.html#indexing-view-versus-copy")

            if value == 'raise':
                raise SettingWithCopyError(t)
            elif value == 'warn':
                warnings.warn(t, SettingWithCopyWarning, stacklevel=stacklevel)

    def __delitem__(self, key):
        """
        Delete item
        """
        deleted = False

        maybe_shortcut = False
        if hasattr(self, 'columns') and isinstance(self.columns, MultiIndex):
            try:
                maybe_shortcut = key not in self.columns._engine
            except TypeError:
                pass

        if maybe_shortcut:
            # Allow shorthand to delete all columns whose first len(key)
            # elements match key:
            if not isinstance(key, tuple):
                key = (key,)
            for col in self.columns:
                if isinstance(col, tuple) and col[:len(key)] == key:
                    del self[col]
                    deleted = True
        if not deleted:
            # If the above loop ran and didn't delete anything because
            # there was no match, this call should raise the appropriate
            # exception:
            self._data.delete(key)

        # delete from the caches
        try:
            del self._item_cache[key]
        except KeyError:
            pass

    def take(self, indices, axis=0, convert=True, is_copy=True):
        """
        Analogous to ndarray.take

        Parameters
        ----------
        indices : list / array of ints
        axis : int, default 0
        convert : translate neg to pos indices (default)
        is_copy : mark the returned frame as a copy

        Returns
        -------
        taken : type of caller
        """

        self._consolidate_inplace()
        new_data = self._data.take(indices,
                                   axis=self._get_block_manager_axis(axis),
                                   convert=True, verify=True)
        result = self._constructor(new_data).__finalize__(self)

        # maybe set copy if we didn't actually change the index
        if is_copy:
            if not result._get_axis(axis).equals(self._get_axis(axis)):
                result._set_is_copy(self)

        return result

    def xs(self, key, axis=0, level=None, copy=None, drop_level=True):
        """
        Returns a cross-section (row(s) or column(s)) from the Series/DataFrame.
        Defaults to cross-section on the rows (axis=0).

        Parameters
        ----------
        key : object
            Some label contained in the index, or partially in a MultiIndex
        axis : int, default 0
            Axis to retrieve cross-section on
        level : object, defaults to first n levels (n=1 or len(key))
            In case of a key partially contained in a MultiIndex, indicate
            which levels are used. Levels can be referred by label or position.
        copy : boolean [deprecated]
            Whether to make a copy of the data
        drop_level : boolean, default True
            If False, returns object with same levels as self.

        Examples
        --------
        >>> df
           A  B  C
        a  4  5  2
        b  4  0  9
        c  9  7  3
        >>> df.xs('a')
        A    4
        B    5
        C    2
        Name: a
        >>> df.xs('C', axis=1)
        a    2
        b    9
        c    3
        Name: C

        >>> df
                            A  B  C  D
        first second third
        bar   one    1      4  1  8  9
              two    1      7  5  5  0
        baz   one    1      6  6  8  0
              three  2      5  3  5  3
        >>> df.xs(('baz', 'three'))
               A  B  C  D
        third
        2      5  3  5  3
        >>> df.xs('one', level=1)
                     A  B  C  D
        first third
        bar   1      4  1  8  9
        baz   1      6  6  8  0
        >>> df.xs(('baz', 2), level=[0, 'third'])
                A  B  C  D
        second
        three   5  3  5  3

        Returns
        -------
        xs : Series or DataFrame

        Notes
        -----
        xs is only for getting, not setting values.

        MultiIndex Slicers is a generic way to get/set values on any level or levels
        it is a superset of xs functionality, see :ref:`MultiIndex Slicers <advanced.mi_slicers>`

        """
        if copy is not None:
            warnings.warn("copy keyword is deprecated, "
                          "default is to return a copy or a view if possible")

        axis = self._get_axis_number(axis)
        labels = self._get_axis(axis)
        if level is not None:
            loc, new_ax = labels.get_loc_level(key, level=level,
                                               drop_level=drop_level)

            # convert to a label indexer if needed
            if isinstance(loc, slice):
                lev_num = labels._get_level_number(level)
                if labels.levels[lev_num].inferred_type == 'integer':
                    loc = labels[loc]

            # create the tuple of the indexer
            indexer = [slice(None)] * self.ndim
            indexer[axis] = loc
            indexer = tuple(indexer)

            result = self.ix[indexer]
            setattr(result, result._get_axis_name(axis), new_ax)
            return result

        if axis == 1:
            return self[key]

        self._consolidate_inplace()

        index = self.index
        if isinstance(index, MultiIndex):
            loc, new_index = self.index.get_loc_level(key,
                                                      drop_level=drop_level)
        else:
            loc = self.index.get_loc(key)

            if isinstance(loc, np.ndarray):
                if loc.dtype == np.bool_:
                    inds, = loc.nonzero()
                    return self.take(inds, axis=axis, convert=False)
                else:
                    return self.take(loc, axis=axis, convert=True)

            if not np.isscalar(loc):
                new_index = self.index[loc]

        if np.isscalar(loc):
            from pandas import Series
            new_values = self._data.fast_xs(loc)

            # may need to box a datelike-scalar
            #
            # if we encounter an array-like and we only have 1 dim
            # that means that their are list/ndarrays inside the Series!
            # so just return them (GH 6394)
            if not is_list_like(new_values) or self.ndim == 1:
                return _maybe_box_datetimelike(new_values)

            result = Series(new_values,
                            index=self.columns,
                            name=self.index[loc],
                            copy=copy,
                            dtype=new_values.dtype)

        else:
            result = self.iloc[loc]
            result.index = new_index

        # this could be a view
        # but only in a single-dtyped view slicable case
        result._set_is_copy(self, copy=not result._is_view)
        return result

    _xs = xs

    # TODO: Check if this was clearer in 0.12
    def select(self, crit, axis=0):
        """
        Return data corresponding to axis labels matching criteria

        Parameters
        ----------
        crit : function
            To be called on each index (label). Should return True or False
        axis : int

        Returns
        -------
        selection : type of caller
        """
        axis = self._get_axis_number(axis)
        axis_name = self._get_axis_name(axis)
        axis_values = self._get_axis(axis)

        if len(axis_values) > 0:
            new_axis = axis_values[
                np.asarray([bool(crit(label)) for label in axis_values])]
        else:
            new_axis = axis_values

        return self.reindex(**{axis_name: new_axis})

    def reindex_like(self, other, method=None, copy=True, limit=None,
                     tolerance=None):
        """ return an object with matching indicies to myself

        Parameters
        ----------
        other : Object
        method : string or None
        copy : boolean, default True
        limit : int, default None
            Maximum number of consecutive labels to fill for inexact matches.
        tolerance : optional
            Maximum distance between labels of the other object and this
            object for inexact matches.

            .. versionadded:: 0.17.0

        Notes
        -----
        Like calling s.reindex(index=other.index, columns=other.columns,
                               method=...)

        Returns
        -------
        reindexed : same as input
        """
        d = other._construct_axes_dict(axes=self._AXIS_ORDERS,
                method=method, copy=copy, limit=limit,
                tolerance=tolerance)

        return self.reindex(**d)

    def drop(self, labels, axis=0, level=None, inplace=False, errors='raise'):
        """
        Return new object with labels in requested axis removed

        Parameters
        ----------
        labels : single label or list-like
        axis : int or axis name
        level : int or level name, default None
            For MultiIndex
        inplace : bool, default False
            If True, do operation inplace and return None.
        errors : {'ignore', 'raise'}, default 'raise'
            If 'ignore', suppress error and existing labels are dropped.

            .. versionadded:: 0.16.1

        Returns
        -------
        dropped : type of caller
        """
        axis = self._get_axis_number(axis)
        axis_name = self._get_axis_name(axis)
        axis, axis_ = self._get_axis(axis), axis

        if axis.is_unique:
            if level is not None:
                if not isinstance(axis, MultiIndex):
                    raise AssertionError('axis must be a MultiIndex')
                new_axis = axis.drop(labels, level=level, errors=errors)
            else:
                new_axis = axis.drop(labels, errors=errors)
            dropped = self.reindex(**{axis_name: new_axis})
            try:
                dropped.axes[axis_].set_names(axis.names, inplace=True)
            except AttributeError:
                pass
            result = dropped

        else:
            labels = com._index_labels_to_array(labels)
            if level is not None:
                if not isinstance(axis, MultiIndex):
                    raise AssertionError('axis must be a MultiIndex')
                indexer = ~lib.ismember(axis.get_level_values(level).values,
                                        set(labels))
            else:
                indexer = ~axis.isin(labels)

            slicer = [slice(None)] * self.ndim
            slicer[self._get_axis_number(axis_name)] = indexer

            result = self.ix[tuple(slicer)]

        if inplace:
            self._update_inplace(result)
        else:
            return result

    def _update_inplace(self, result, verify_is_copy=True):
        """
        replace self internals with result.

        Parameters
        ----------
        verify_is_copy : boolean, default True
            provide is_copy checks

        """
        # NOTE: This does *not* call __finalize__ and that's an explicit
        # decision that we may revisit in the future.

        self._reset_cache()
        self._clear_item_cache()
        self._data = getattr(result,'_data',result)
        self._maybe_update_cacher(verify_is_copy=verify_is_copy)

    def add_prefix(self, prefix):
        """
        Concatenate prefix string with panel items names.

        Parameters
        ----------
        prefix : string

        Returns
        -------
        with_prefix : type of caller
        """
        new_data = self._data.add_prefix(prefix)
        return self._constructor(new_data).__finalize__(self)

    def add_suffix(self, suffix):
        """
        Concatenate suffix string with panel items names

        Parameters
        ----------
        suffix : string

        Returns
        -------
        with_suffix : type of caller
        """
        new_data = self._data.add_suffix(suffix)
        return self._constructor(new_data).__finalize__(self)

    _shared_docs['sort_values'] = """
        Sort by the values along either axis

        .. versionadded:: 0.17.0

        Parameters
        ----------
        by : string name or list of names which refer to the axis items
        axis : %(axes)s to direct sorting
        ascending : bool or list of bool
             Sort ascending vs. descending. Specify list for multiple sort orders.
             If this is a list of bools, must match the length of the by
        inplace : bool
             if True, perform operation in-place
        kind : {`quicksort`, `mergesort`, `heapsort`}
             Choice of sorting algorithm. See also ndarray.np.sort for more information.
             `mergesort` is the only stable algorithm. For DataFrames, this option is
             only applied when sorting on a single column or label.
        na_position : {'first', 'last'}
             `first` puts NaNs at the beginning, `last` puts NaNs at the end

        Returns
        -------
        sorted_obj : %(klass)s
        """
    def sort_values(self, by, axis=0, ascending=True, inplace=False,
                    kind='quicksort', na_position='last'):
        raise AbstractMethodError(self)

    _shared_docs['sort_index'] = """
        Sort object by labels (along an axis)

        Parameters
        ----------
        axis : %(axes)s to direct sorting
        level : int or level name or list of ints or list of level names
            if not None, sort on values in specified index level(s)
        ascending : boolean, default True
            Sort ascending vs. descending
        inplace : bool
            if True, perform operation in-place
        kind : {`quicksort`, `mergesort`, `heapsort`}
             Choice of sorting algorithm. See also ndarray.np.sort for more information.
             `mergesort` is the only stable algorithm. For DataFrames, this option is
             only applied when sorting on a single column or label.
        na_position : {'first', 'last'}
             `first` puts NaNs at the beginning, `last` puts NaNs at the end
        sort_remaining : bool
            if true and sorting by level and index is multilevel, sort by other levels
            too (in order) after sorting by specified level

        Returns
        -------
        sorted_obj : %(klass)s
        """

    @Appender(_shared_docs['sort_index'] % dict(axes="axes", klass="NDFrame"))
    def sort_index(self, axis=0, level=None, ascending=True, inplace=False,
                   kind='quicksort', na_position='last', sort_remaining=True):
        axis = self._get_axis_number(axis)
        axis_name = self._get_axis_name(axis)
        labels = self._get_axis(axis)

        if level is not None:
            raise NotImplementedError("level is not implemented")
        if inplace:
            raise NotImplementedError("inplace is not implemented")

        sort_index = labels.argsort()
        if not ascending:
            sort_index = sort_index[::-1]

        new_axis = labels.take(sort_index)
        return self.reindex(**{axis_name: new_axis})

    _shared_docs['reindex'] = """
        Conform %(klass)s to new index with optional filling logic, placing
        NA/NaN in locations having no value in the previous index. A new object
        is produced unless the new index is equivalent to the current one and
        copy=False

        Parameters
        ----------
        %(axes)s : array-like, optional (can be specified in order, or as
            keywords)
            New labels / index to conform to. Preferably an Index object to
            avoid duplicating data
        method : {None, 'backfill'/'bfill', 'pad'/'ffill', 'nearest'}, optional
            Method to use for filling holes in reindexed DataFrame:
              * default: don't fill gaps
              * pad / ffill: propagate last valid observation forward to next valid
              * backfill / bfill: use next valid observation to fill gap
              * nearest: use nearest valid observations to fill gap
        copy : boolean, default True
            Return a new object, even if the passed indexes are the same
        level : int or name
            Broadcast across a level, matching Index values on the
            passed MultiIndex level
        fill_value : scalar, default np.NaN
            Value to use for missing values. Defaults to NaN, but can be any
            "compatible" value
        limit : int, default None
            Maximum number of consecutive elements to forward or backward fill
        tolerance : optional
            Maximum distance between original and new labels for inexact
            matches. The values of the index at the matching locations most
            satisfy the equation ``abs(index[indexer] - target) <= tolerance``.

            .. versionadded:: 0.17.0

        Examples
        --------
        >>> df.reindex(index=[date1, date2, date3], columns=['A', 'B', 'C'])

        Returns
        -------
        reindexed : %(klass)s
        """
    # TODO: Decide if we care about having different examples for different
    #       kinds

    @Appender(_shared_docs['reindex'] % dict(axes="axes", klass="NDFrame"))
    def reindex(self, *args, **kwargs):

        # construct the args
        axes, kwargs = self._construct_axes_from_arguments(args, kwargs)
        method = com._clean_reindex_fill_method(kwargs.pop('method', None))
        level = kwargs.pop('level', None)
        copy = kwargs.pop('copy', True)
        limit = kwargs.pop('limit', None)
        tolerance = kwargs.pop('tolerance', None)
        fill_value = kwargs.pop('fill_value', np.nan)

        if kwargs:
            raise TypeError('reindex() got an unexpected keyword '
                    'argument "{0}"'.format(list(kwargs.keys())[0]))

        self._consolidate_inplace()

        # if all axes that are requested to reindex are equal, then only copy
        # if indicated must have index names equal here as well as values
        if all([self._get_axis(axis).identical(ax)
                for axis, ax in axes.items() if ax is not None]):
            if copy:
                return self.copy()
            return self

        # check if we are a multi reindex
        if self._needs_reindex_multi(axes, method, level):
            try:
                return self._reindex_multi(axes, copy, fill_value)
            except:
                pass

        # perform the reindex on the axes
        return self._reindex_axes(axes, level, limit, tolerance,
                                  method, fill_value, copy).__finalize__(self)

    def _reindex_axes(self, axes, level, limit, tolerance, method,
                      fill_value, copy):
        """ perform the reinxed for all the axes """
        obj = self
        for a in self._AXIS_ORDERS:
            labels = axes[a]
            if labels is None:
                continue

            ax = self._get_axis(a)
            new_index, indexer = ax.reindex(
                labels, level=level, limit=limit, tolerance=tolerance,
                method=method)

            axis = self._get_axis_number(a)
            obj = obj._reindex_with_indexers(
                {axis: [new_index, indexer]},
                fill_value=fill_value, copy=copy, allow_dups=False)

        return obj

    def _needs_reindex_multi(self, axes, method, level):
        """ check if we do need a multi reindex """
        return ((com._count_not_none(*axes.values()) == self._AXIS_LEN) and
                method is None and level is None and not self._is_mixed_type)

    def _reindex_multi(self, axes, copy, fill_value):
        return NotImplemented

    _shared_docs['reindex_axis'] = (
        """Conform input object to new index with optional filling logic,
        placing NA/NaN in locations having no value in the previous index. A
        new object is produced unless the new index is equivalent to the
        current one and copy=False

        Parameters
        ----------
        labels : array-like
            New labels / index to conform to. Preferably an Index object to
            avoid duplicating data
        axis : %(axes_single_arg)s
        method : {None, 'backfill'/'bfill', 'pad'/'ffill', 'nearest'}, optional
            Method to use for filling holes in reindexed DataFrame:
              * default: don't fill gaps
              * pad / ffill: propagate last valid observation forward to next valid
              * backfill / bfill: use next valid observation to fill gap
              * nearest: use nearest valid observations to fill gap
        copy : boolean, default True
            Return a new object, even if the passed indexes are the same
        level : int or name
            Broadcast across a level, matching Index values on the
            passed MultiIndex level
        limit : int, default None
            Maximum number of consecutive elements to forward or backward fill
        tolerance : optional
            Maximum distance between original and new labels for inexact
            matches. The values of the index at the matching locations most
            satisfy the equation ``abs(index[indexer] - target) <= tolerance``.

            .. versionadded:: 0.17.0

        Examples
        --------
        >>> df.reindex_axis(['A', 'B', 'C'], axis=1)

        See also
        --------
        reindex, reindex_like

        Returns
        -------
        reindexed : %(klass)s
        """)

    @Appender(_shared_docs['reindex_axis'] % _shared_doc_kwargs)
    def reindex_axis(self, labels, axis=0, method=None, level=None, copy=True,
                     limit=None, fill_value=np.nan):
        self._consolidate_inplace()

        axis_name = self._get_axis_name(axis)
        axis_values = self._get_axis(axis_name)
        method = com._clean_reindex_fill_method(method)
        new_index, indexer = axis_values.reindex(labels, method, level,
                                                 limit=limit)
        return self._reindex_with_indexers(
            {axis: [new_index, indexer]}, fill_value=fill_value, copy=copy)

    def _reindex_with_indexers(self, reindexers,
                               fill_value=np.nan, copy=False,
                               allow_dups=False):
        """ allow_dups indicates an internal call here """

        # reindex doing multiple operations on different axes if indiciated
        new_data = self._data
        for axis in sorted(reindexers.keys()):
            index, indexer = reindexers[axis]
            baxis = self._get_block_manager_axis(axis)

            if index is None:
                continue

            index = _ensure_index(index)
            if indexer is not None:
                indexer = com._ensure_int64(indexer)

            # TODO: speed up on homogeneous DataFrame objects
            new_data = new_data.reindex_indexer(index, indexer, axis=baxis,
                                                fill_value=fill_value,
                                                allow_dups=allow_dups,
                                                copy=copy)

        if copy and new_data is self._data:
            new_data = new_data.copy()

        return self._constructor(new_data).__finalize__(self)

    def _reindex_axis(self, new_index, fill_method, axis, copy):
        new_data = self._data.reindex_axis(new_index, axis=axis,
                                           method=fill_method, copy=copy)

        if new_data is self._data and not copy:
            return self
        else:
            return self._constructor(new_data).__finalize__(self)

    def filter(self, items=None, like=None, regex=None, axis=None):
        """
        Restrict the info axis to set of items or wildcard

        Parameters
        ----------
        items : list-like
            List of info axis to restrict to (must not all be present)
        like : string
            Keep info axis where "arg in col == True"
        regex : string (regular expression)
            Keep info axis with re.search(regex, col) == True
        axis : int or None
            The axis to filter on. By default this is the info axis. The "info
            axis" is the axis that is used when indexing with ``[]``. For
            example, ``df = DataFrame({'a': [1, 2, 3, 4]]}); df['a']``. So,
            the ``DataFrame`` columns are the info axis.

        Notes
        -----
        Arguments are mutually exclusive, but this is not checked for

        """
        import re

        if axis is None:
            axis = self._info_axis_name
        axis_name = self._get_axis_name(axis)
        axis_values = self._get_axis(axis_name)

        if items is not None:
            return self.reindex(**{axis_name: [r for r in items
                                               if r in axis_values]})
        elif like:
            matchf = lambda x: (like in x if isinstance(x, string_types)
                                else like in str(x))
            return self.select(matchf, axis=axis_name)
        elif regex:
            matcher = re.compile(regex)
            return self.select(lambda x: matcher.search(str(x)) is not None,
                               axis=axis_name)
        else:
            raise TypeError('Must pass either `items`, `like`, or `regex`')

    def head(self, n=5):
        """
        Returns first n rows
        """
        l = len(self)
        if l == 0 or n==0:
            return self
        return self.iloc[:n]

    def tail(self, n=5):
        """
        Returns last n rows
        """
        l = len(self)
        if l == 0 or n == 0:
            return self
        return self.iloc[-n:]


    def sample(self, n=None, frac=None, replace=False, weights=None, random_state=None, axis=None):
        """
        Returns a random sample of items from an axis of object.

        .. versionadded:: 0.16.1

        Parameters
        ----------
        n : int, optional
            Number of items from axis to return. Cannot be used with `frac`.
            Default = 1 if `frac` = None.
        frac : float, optional
            Fraction of axis items to return. Cannot be used with `n`.
        replace : boolean, optional
            Sample with or without replacement. Default = False.
        weights : str or ndarray-like, optional
            Default 'None' results in equal probability weighting.
            If passed a Series, will align with target object on index. Index
            values in weights not found in sampled object will be ignored and
            index values in sampled object not in weights will be assigned
            weights of zero.
            If called on a DataFrame, will accept the name of a column
            when axis = 0.
            Unless weights are a Series, weights must be same length as axis
            being sampled.
            If weights do not sum to 1, they will be normalized to sum to 1.
            Missing values in the weights column will be treated as zero.
            inf and -inf values not allowed.
        random_state : int or numpy.random.RandomState, optional
            Seed for the random number generator (if int), or numpy RandomState
            object.
        axis : int or string, optional
            Axis to sample. Accepts axis number or name. Default is stat axis
            for given data type (0 for Series and DataFrames, 1 for Panels).

        Returns
        -------
        A new object of same type as caller.
        """

        if axis is None:
            axis = self._stat_axis_number

        axis = self._get_axis_number(axis)
        axis_length = self.shape[axis]

        # Process random_state argument
        rs = com._random_state(random_state)

        # Check weights for compliance
        if weights is not None:

            # If a series, align with frame
            if isinstance(weights, pd.Series):
                weights = weights.reindex(self.axes[axis])

            # Strings acceptable if a dataframe and axis = 0
            if isinstance(weights, string_types):
                if isinstance(self, pd.DataFrame):
                    if axis == 0:
                        try:
                            weights = self[weights]
                        except KeyError:
                            raise KeyError("String passed to weights not a valid column")
                    else:
                        raise ValueError("Strings can only be passed to weights when sampling from rows on a DataFrame")
                else:
                    raise ValueError("Strings cannot be passed as weights when sampling from a Series or Panel.")

            weights = pd.Series(weights, dtype='float64')

            if len(weights) != axis_length:
                raise ValueError("Weights and axis to be sampled must be of same length")

            if (weights == np.inf).any() or (weights == -np.inf).any():
                raise ValueError("weight vector may not include `inf` values")

            if (weights < 0).any():
                raise ValueError("weight vector many not include negative values")

            # If has nan, set to zero.
            weights = weights.fillna(0)

            # Renormalize if don't sum to 1
            if weights.sum() != 1:
                if weights.sum() != 0:
                    weights = weights / weights.sum()
                else:
                    raise ValueError("Invalid weights: weights sum to zero")

            weights = weights.values

        # If no frac or n, default to n=1.
        if n is None and frac is None:
            n = 1
        elif n is not None and frac is None and n % 1 != 0:
            raise ValueError("Only integers accepted as `n` values")
        elif n is None and frac is not None:
            n = int(round(frac * axis_length))
        elif n is not None and frac is not None:
            raise ValueError('Please enter a value for `frac` OR `n`, not both')

        # Check for negative sizes
        if n < 0:
            raise ValueError("A negative number of rows requested. Please provide positive value.")

        locs = rs.choice(axis_length, size=n, replace=replace, p=weights)
        return self.take(locs, axis=axis, is_copy=False)


    _shared_docs['pipe'] = ("""
        Apply func(self, \*args, \*\*kwargs)

        .. versionadded:: 0.16.2

        Parameters
        ----------
        func : function
            function to apply to the %(klass)s.
            ``args``, and ``kwargs`` are passed into ``func``.
            Alternatively a ``(callable, data_keyword)`` tuple where
            ``data_keyword`` is a string indicating the keyword of
            ``callable`` that expects the %(klass)s.
        args : positional arguments passed into ``func``.
        kwargs : a dictionary of keyword arguments passed into ``func``.

        Returns
        -------
        object : the return type of ``func``.

        Notes
        -----

        Use ``.pipe`` when chaining together functions that expect
        on Series or DataFrames. Instead of writing

        >>> f(g(h(df), arg1=a), arg2=b, arg3=c)

        You can write

        >>> (df.pipe(h)
        ...    .pipe(g, arg1=a)
        ...    .pipe(f, arg2=b, arg3=c)
        ... )

        If you have a function that takes the data as (say) the second
        argument, pass a tuple indicating which keyword expects the
        data. For example, suppose ``f`` takes its data as ``arg2``:

        >>> (df.pipe(h)
        ...    .pipe(g, arg1=a)
        ...    .pipe((f, 'arg2'), arg1=a, arg3=c)
        ...  )

        See Also
        --------
        pandas.DataFrame.apply
        pandas.DataFrame.applymap
        pandas.Series.map
    """
    )
    @Appender(_shared_docs['pipe'] % _shared_doc_kwargs)
    def pipe(self, func, *args, **kwargs):
        if isinstance(func, tuple):
            func, target = func
            if target in kwargs:
                msg = '%s is both the pipe target and a keyword argument' % target
                raise ValueError(msg)
            kwargs[target] = self
            return func(*args, **kwargs)
        else:
            return func(self, *args, **kwargs)

    #----------------------------------------------------------------------
    # Attribute access

    def __finalize__(self, other, method=None, **kwargs):
        """
        propagate metadata from other to self

        Parameters
        ----------
        other : the object from which to get the attributes that we are going
            to propagate
        method : optional, a passed method name ; possibly to take different
            types of propagation actions based on this

        """
        if isinstance(other, NDFrame):
            for name in self._metadata:
                object.__setattr__(self, name, getattr(other, name, None))
        return self

    def __getattr__(self, name):
        """After regular attribute access, try looking up the name
        This allows simpler access to columns for interactive use.
        """
        # Note: obj.x will always call obj.__getattribute__('x') prior to
        # calling obj.__getattr__('x').

        if (name in self._internal_names_set
                or name in self._metadata
                or name in self._accessors):
            return object.__getattribute__(self, name)
        else:
            if name in self._info_axis:
                return self[name]
            raise AttributeError("'%s' object has no attribute '%s'" %
                                 (type(self).__name__, name))

    def __setattr__(self, name, value):
        """After regular attribute access, try setting the name
        This allows simpler access to columns for interactive use."""
        # first try regular attribute access via __getattribute__, so that
        # e.g. ``obj.x`` and ``obj.x = 4`` will always reference/modify
        # the same attribute.

        try:
            object.__getattribute__(self, name)
            return object.__setattr__(self, name, value)
        except AttributeError:
            pass

        # if this fails, go on to more involved attribute setting
        # (note that this matches __getattr__, above).
        if name in self._internal_names_set:
            object.__setattr__(self, name, value)
        elif name in self._metadata:
            object.__setattr__(self, name, value)
        else:
            try:
                existing = getattr(self, name)
                if isinstance(existing, Index):
                    object.__setattr__(self, name, value)
                elif name in self._info_axis:
                    self[name] = value
                else:
                    object.__setattr__(self, name, value)
            except (AttributeError, TypeError):
                object.__setattr__(self, name, value)

    #----------------------------------------------------------------------
    # Getting and setting elements

    #----------------------------------------------------------------------
    # Consolidation of internals

    def _protect_consolidate(self, f):
        """ consolidate _data. if the blocks have changed, then clear the cache """
        blocks_before = len(self._data.blocks)
        result = f()
        if len(self._data.blocks) != blocks_before:
            self._clear_item_cache()
        return result

    def _consolidate_inplace(self):
        """ we are inplace consolidating; return None """
        def f():
            self._data = self._data.consolidate()
        self._protect_consolidate(f)

    def consolidate(self, inplace=False):
        """
        Compute NDFrame with "consolidated" internals (data of each dtype
        grouped together in a single ndarray). Mainly an internal API function,
        but available here to the savvy user

        Parameters
        ----------
        inplace : boolean, default False
            If False return new object, otherwise modify existing object

        Returns
        -------
        consolidated : type of caller
        """
        if inplace:
            self._consolidate_inplace()
        else:
            f = lambda: self._data.consolidate()
            cons_data = self._protect_consolidate(f)
            return self._constructor(cons_data).__finalize__(self)

    @property
    def _is_mixed_type(self):
        f = lambda: self._data.is_mixed_type
        return self._protect_consolidate(f)

    @property
    def _is_numeric_mixed_type(self):
        f = lambda: self._data.is_numeric_mixed_type
        return self._protect_consolidate(f)

    @property
    def _is_datelike_mixed_type(self):
        f = lambda: self._data.is_datelike_mixed_type
        return self._protect_consolidate(f)

    def _check_inplace_setting(self, value):
        """ check whether we allow in-place setting with this type of value """

        if self._is_mixed_type:
            if not self._is_numeric_mixed_type:

                # allow an actual np.nan thru
                try:
                    if np.isnan(value):
                        return True
                except:
                    pass

                raise TypeError(
                    'Cannot do inplace boolean setting on mixed-types with a non np.nan value')

        return True

    def _get_numeric_data(self):
        return self._constructor(
            self._data.get_numeric_data()).__finalize__(self)

    def _get_bool_data(self):
        return self._constructor(self._data.get_bool_data()).__finalize__(self)

    #----------------------------------------------------------------------
    # Internal Interface Methods

    def as_matrix(self, columns=None):
        """
        Convert the frame to its Numpy-array representation.

        Parameters
        ----------
        columns: list, optional, default:None
            If None, return all columns, otherwise, returns specified columns.

        Returns
        -------
        values : ndarray
            If the caller is heterogeneous and contains booleans or objects,
            the result will be of dtype=object. See Notes.


        Notes
        -----
        Return is NOT a Numpy-matrix, rather, a Numpy-array.

        The dtype will be a lower-common-denominator dtype (implicit
        upcasting); that is to say if the dtypes (even of numeric types)
        are mixed, the one that accommodates all will be chosen. Use this
        with care if you are not dealing with the blocks.

        e.g. If the dtypes are float16 and float32, dtype will be upcast to
        float32.  If dtypes are int32 and uint8, dtype will be upcase to
        int32.

        This method is provided for backwards compatibility. Generally,
        it is recommended to use '.values'.

        See Also
        --------
        pandas.DataFrame.values
        """
        self._consolidate_inplace()
        if self._AXIS_REVERSED:
            return self._data.as_matrix(columns).T
        return self._data.as_matrix(columns)

    @property
    def values(self):
        """Numpy representation of NDFrame

        Notes
        -----
        The dtype will be a lower-common-denominator dtype (implicit
        upcasting); that is to say if the dtypes (even of numeric types)
        are mixed, the one that accommodates all will be chosen. Use this
        with care if you are not dealing with the blocks.

        e.g. If the dtypes are float16 and float32, dtype will be upcast to
        float32.  If dtypes are int32 and uint8, dtype will be upcase to
        int32.
        """
        return self.as_matrix()

    @property
    def _values(self):
        """ internal implementation """
        return self.values

    @property
    def _get_values(self):
        # compat
        return self.as_matrix()

    def get_values(self):
        """ same as values (but handles sparseness conversions) """
        return self.as_matrix()

    def get_dtype_counts(self):
        """ Return the counts of dtypes in this object """
        from pandas import Series
        return Series(self._data.get_dtype_counts())

    def get_ftype_counts(self):
        """ Return the counts of ftypes in this object """
        from pandas import Series
        return Series(self._data.get_ftype_counts())

    @property
    def dtypes(self):
        """ Return the dtypes in this object """
        from pandas import Series
        return Series(self._data.get_dtypes(), index=self._info_axis,
                      dtype=np.object_)

    @property
    def ftypes(self):
        """
        Return the ftypes (indication of sparse/dense and dtype)
        in this object.
        """
        from pandas import Series
        return Series(self._data.get_ftypes(), index=self._info_axis,
                      dtype=np.object_)

    def as_blocks(self, copy=True):
        """
        Convert the frame to a dict of dtype -> Constructor Types that each has
        a homogeneous dtype.

        NOTE: the dtypes of the blocks WILL BE PRESERVED HERE (unlike in
              as_matrix)

        Parameters
        ----------
        copy : boolean, default True

               .. versionadded: 0.16.1

        Returns
        -------
        values : a dict of dtype -> Constructor Types
        """
        self._consolidate_inplace()

        bd = {}
        for b in self._data.blocks:
            bd.setdefault(str(b.dtype), []).append(b)

        result = {}
        for dtype, blocks in bd.items():
            # Must combine even after consolidation, because there may be
            # sparse items which are never consolidated into one block.
            combined = self._data.combine(blocks, copy=copy)
            result[dtype] = self._constructor(combined).__finalize__(self)

        return result

    @property
    def blocks(self):
        "Internal property, property synonym for as_blocks()"
        return self.as_blocks()

    def astype(self, dtype, copy=True, raise_on_error=True, **kwargs):
        """
        Cast object to input numpy.dtype
        Return a copy when copy = True (be really careful with this!)

        Parameters
        ----------
        dtype : numpy.dtype or Python type
        raise_on_error : raise on invalid input
        kwargs : keyword arguments to pass on to the constructor

        Returns
        -------
        casted : type of caller
        """

        mgr = self._data.astype(
            dtype=dtype, copy=copy, raise_on_error=raise_on_error, **kwargs)
        return self._constructor(mgr).__finalize__(self)

    def copy(self, deep=True):
        """
        Make a copy of this object

        Parameters
        ----------
        deep : boolean or string, default True
            Make a deep copy, i.e. also copy data

        Returns
        -------
        copy : type of caller
        """
        data = self._data.copy(deep=deep)
        return self._constructor(data).__finalize__(self)

    def _convert(self, datetime=False, numeric=False, timedelta=False,
                 coerce=False, copy=True):
        """
        Attempt to infer better dtype for object columns

        Parameters
        ----------
        datetime : boolean, default False
            If True, convert to date where possible.
        numeric : boolean, default False
            If True, attempt to convert to numbers (including strings), with
            unconvertible values becoming NaN.
        timedelta : boolean, default False
            If True, convert to timedelta where possible.
        coerce : boolean, default False
            If True, force conversion with unconvertible values converted to
            nulls (NaN or NaT)
        copy : boolean, default True
            If True, return a copy even if no copy is necessary (e.g. no
            conversion was done). Note: This is meant for internal use, and
            should not be confused with inplace.

        Returns
        -------
        converted : same as input object
        """
        return self._constructor(
            self._data.convert(datetime=datetime,
                                numeric=numeric,
                                timedelta=timedelta,
                                coerce=coerce,
                                copy=copy)).__finalize__(self)

    # TODO: Remove in 0.18 or 2017, which ever is sooner
    def convert_objects(self, convert_dates=True, convert_numeric=False,
                        convert_timedeltas=True, copy=True):
        """
        Attempt to infer better dtype for object columns

        Parameters
        ----------
        convert_dates : boolean, default True
            If True, convert to date where possible. If 'coerce', force
            conversion, with unconvertible values becoming NaT.
        convert_numeric : boolean, default False
            If True, attempt to coerce to numbers (including strings), with
            unconvertible values becoming NaN.
        convert_timedeltas : boolean, default True
            If True, convert to timedelta where possible. If 'coerce', force
            conversion, with unconvertible values becoming NaT.
        copy : boolean, default True
            If True, return a copy even if no copy is necessary (e.g. no
            conversion was done). Note: This is meant for internal use, and
            should not be confused with inplace.

        Returns
        -------
        converted : same as input object
        """
        from warnings import warn
        warn("convert_objects is deprecated.  Use the data-type specific "
             "converters pd.to_datetime, pd.to_timedelta and pd.to_numeric.",
             FutureWarning, stacklevel=2)

        return self._constructor(
            self._data.convert(convert_dates=convert_dates,
                               convert_numeric=convert_numeric,
                               convert_timedeltas=convert_timedeltas,
                               copy=copy)).__finalize__(self)

    #----------------------------------------------------------------------
    # Filling NA's

    _shared_docs['fillna'] = (
        """
        Fill NA/NaN values using the specified method

        Parameters
        ----------
        value : scalar, dict, Series, or DataFrame
            Value to use to fill holes (e.g. 0), alternately a dict/Series/DataFrame of
            values specifying which value to use for each index (for a Series) or
            column (for a DataFrame). (values not in the dict/Series/DataFrame will not be
            filled). This value cannot be a list.
        method : {'backfill', 'bfill', 'pad', 'ffill', None}, default None
            Method to use for filling holes in reindexed Series
            pad / ffill: propagate last valid observation forward to next valid
            backfill / bfill: use NEXT valid observation to fill gap
        axis : %(axes_single_arg)s
        inplace : boolean, default False
            If True, fill in place. Note: this will modify any
            other views on this object, (e.g. a no-copy slice for a column in a
            DataFrame).
        limit : int, default None
            If method is specified, this is the maximum number of consecutive
            NaN values to forward/backward fill. In other words, if there is
            a gap with more than this number of consecutive NaNs, it will only
            be partially filled. If method is not specified, this is the
            maximum number of entries along the entire axis where NaNs will be
            filled.
        downcast : dict, default is None
            a dict of item->dtype of what to downcast if possible,
            or the string 'infer' which will try to downcast to an appropriate
            equal type (e.g. float64 to int64 if possible)

        See also
        --------
        reindex, asfreq

        Returns
        -------
        filled : %(klass)s
        """
    )

    @Appender(_shared_docs['fillna'] % _shared_doc_kwargs)
    def fillna(self, value=None, method=None, axis=None, inplace=False,
               limit=None, downcast=None):
        if isinstance(value, (list, tuple)):
            raise TypeError('"value" parameter must be a scalar or dict, but '
                            'you passed a "{0}"'.format(type(value).__name__))
        self._consolidate_inplace()

        # set the default here, so functions examining the signaure
        # can detect if something was set (e.g. in groupby) (GH9221)
        if axis is None:
            axis = 0
        axis = self._get_axis_number(axis)
        method = com._clean_fill_method(method)

        from pandas import DataFrame
        if value is None:
            if method is None:
                raise ValueError('must specify a fill method or value')
            if self._is_mixed_type and axis == 1:
                if inplace:
                    raise NotImplementedError()
                result = self.T.fillna(method=method, limit=limit).T

                # need to downcast here because of all of the transposes
                result._data = result._data.downcast()

                return result

            # > 3d
            if self.ndim > 3:
                raise NotImplementedError(
                    'Cannot fillna with a method for > 3dims'
                )

            # 3d
            elif self.ndim == 3:

                # fill in 2d chunks
                result = dict([(col, s.fillna(method=method, value=value))
                               for col, s in compat.iteritems(self)])
                return self._constructor.from_dict(result).__finalize__(self)

            # 2d or less
            method = com._clean_fill_method(method)
            new_data = self._data.interpolate(method=method,
                                              axis=axis,
                                              limit=limit,
                                              inplace=inplace,
                                              coerce=True,
                                              downcast=downcast)
        else:
            if method is not None:
                raise ValueError('cannot specify both a fill method and value')

            if len(self._get_axis(axis)) == 0:
                return self

            if self.ndim == 1:
                if isinstance(value, (dict, com.ABCSeries)):
                    from pandas import Series
                    value = Series(value)
                elif not com.is_list_like(value):
                    pass
                else:
                    raise ValueError("invalid fill value with a %s" % type(value))

                new_data = self._data.fillna(value=value,
                                             limit=limit,
                                             inplace=inplace,
                                             downcast=downcast)

            elif isinstance(value, (dict, com.ABCSeries)):
                if axis == 1:
                    raise NotImplementedError('Currently only can fill '
                                              'with dict/Series column '
                                              'by column')

                result = self if inplace else self.copy()
                for k, v in compat.iteritems(value):
                    if k not in result:
                        continue
                    obj = result[k]
                    obj.fillna(v, limit=limit, inplace=True)
                return result
            elif not com.is_list_like(value):
                new_data = self._data.fillna(value=value,
                                             limit=limit,
                                             inplace=inplace,
                                             downcast=downcast)
            elif isinstance(value, DataFrame) and self.ndim == 2:
                new_data = self.where(self.notnull(), value)
            else:
                raise ValueError("invalid fill value with a %s" % type(value))

        if inplace:
            self._update_inplace(new_data)
        else:
            return self._constructor(new_data).__finalize__(self)

    def ffill(self, axis=None, inplace=False, limit=None, downcast=None):
        "Synonym for NDFrame.fillna(method='ffill')"
        return self.fillna(method='ffill', axis=axis, inplace=inplace,
                           limit=limit, downcast=downcast)

    def bfill(self, axis=None, inplace=False, limit=None, downcast=None):
        "Synonym for NDFrame.fillna(method='bfill')"
        return self.fillna(method='bfill', axis=axis, inplace=inplace,
                           limit=limit, downcast=downcast)

    def replace(self, to_replace=None, value=None, inplace=False, limit=None,
                regex=False, method='pad', axis=None):
        """
        Replace values given in 'to_replace' with 'value'.

        Parameters
        ----------
        to_replace : str, regex, list, dict, Series, numeric, or None

            * str or regex:

                - str: string exactly matching `to_replace` will be replaced
                  with `value`
                - regex: regexs matching `to_replace` will be replaced with
                  `value`

            * list of str, regex, or numeric:

                - First, if `to_replace` and `value` are both lists, they
                  **must** be the same length.
                - Second, if ``regex=True`` then all of the strings in **both**
                  lists will be interpreted as regexs otherwise they will match
                  directly. This doesn't matter much for `value` since there
                  are only a few possible substitution regexes you can use.
                - str and regex rules apply as above.

            * dict:

                - Nested dictionaries, e.g., {'a': {'b': nan}}, are read as
                  follows: look in column 'a' for the value 'b' and replace it
                  with nan. You can nest regular expressions as well. Note that
                  column names (the top-level dictionary keys in a nested
                  dictionary) **cannot** be regular expressions.
                - Keys map to column names and values map to substitution
                  values. You can treat this as a special case of passing two
                  lists except that you are specifying the column to search in.

            * None:

                - This means that the ``regex`` argument must be a string,
                  compiled regular expression, or list, dict, ndarray or Series
                  of such elements. If `value` is also ``None`` then this
                  **must** be a nested dictionary or ``Series``.

            See the examples section for examples of each of these.
        value : scalar, dict, list, str, regex, default None
            Value to use to fill holes (e.g. 0), alternately a dict of values
            specifying which value to use for each column (columns not in the
            dict will not be filled). Regular expressions, strings and lists or
            dicts of such objects are also allowed.
        inplace : boolean, default False
            If True, in place. Note: this will modify any
            other views on this object (e.g. a column form a DataFrame).
            Returns the caller if this is True.
        limit : int, default None
            Maximum size gap to forward or backward fill
        regex : bool or same types as `to_replace`, default False
            Whether to interpret `to_replace` and/or `value` as regular
            expressions. If this is ``True`` then `to_replace` *must* be a
            string. Otherwise, `to_replace` must be ``None`` because this
            parameter will be interpreted as a regular expression or a list,
            dict, or array of regular expressions.
        method : string, optional, {'pad', 'ffill', 'bfill'}
            The method to use when for replacement, when ``to_replace`` is a
            ``list``.

        See also
        --------
        NDFrame.reindex
        NDFrame.asfreq
        NDFrame.fillna

        Returns
        -------
        filled : NDFrame

        Raises
        ------
        AssertionError
            * If `regex` is not a ``bool`` and `to_replace` is not ``None``.
        TypeError
            * If `to_replace` is a ``dict`` and `value` is not a ``list``,
              ``dict``, ``ndarray``, or ``Series``
            * If `to_replace` is ``None`` and `regex` is not compilable into a
              regular expression or is a list, dict, ndarray, or Series.
        ValueError
            * If `to_replace` and `value` are ``list`` s or ``ndarray`` s, but
              they are not the same length.

        Notes
        -----
        * Regex substitution is performed under the hood with ``re.sub``. The
          rules for substitution for ``re.sub`` are the same.
        * Regular expressions will only substitute on strings, meaning you
          cannot provide, for example, a regular expression matching floating
          point numbers and expect the columns in your frame that have a
          numeric dtype to be matched. However, if those floating point numbers
          *are* strings, then you can do this.
        * This method has *a lot* of options. You are encouraged to experiment
          and play with this method to gain intuition about how it works.

        """
        if not com.is_bool(regex) and to_replace is not None:
            raise AssertionError("'to_replace' must be 'None' if 'regex' is "
                                 "not a bool")
        if axis is not None:
            from warnings import warn
            warn('the "axis" argument is deprecated and will be removed in'
                 'v0.13; this argument has no effect')

        self._consolidate_inplace()

        if value is None:
            # passing a single value that is scalar like
            # when value is None (GH5319), for compat
            if not is_dictlike(to_replace) and not is_dictlike(regex):
                to_replace = [to_replace]

            if isinstance(to_replace, (tuple, list)):
                return _single_replace(self, to_replace, method, inplace,
                                       limit)

            if not is_dictlike(to_replace):
                if not is_dictlike(regex):
                    raise TypeError('If "to_replace" and "value" are both None'
                                    ' and "to_replace" is not a list, then '
                                    'regex must be a mapping')
                to_replace = regex
                regex = True

            items = list(compat.iteritems(to_replace))
            keys, values = zip(*items)

            are_mappings = [is_dictlike(v) for v in values]

            if any(are_mappings):
                if not all(are_mappings):
                    raise TypeError("If a nested mapping is passed, all values"
                                    " of the top level mapping must be "
                                    "mappings")
                # passed a nested dict/Series
                to_rep_dict = {}
                value_dict = {}

                for k, v in items:
                    keys, values = zip(*v.items())
                    if set(keys) & set(values):
                        raise ValueError("Replacement not allowed with "
                                         "overlapping keys and values")
                    to_rep_dict[k] = list(keys)
                    value_dict[k] = list(values)

                to_replace, value = to_rep_dict, value_dict
            else:
                to_replace, value = keys, values

            return self.replace(to_replace, value, inplace=inplace,
                                limit=limit, regex=regex)
        else:

            # need a non-zero len on all axes
            for a in self._AXIS_ORDERS:
                if not len(self._get_axis(a)):
                    return self

            new_data = self._data
            if is_dictlike(to_replace):
                if is_dictlike(value):  # {'A' : NA} -> {'A' : 0}
                    res = self if inplace else self.copy()
                    for c, src in compat.iteritems(to_replace):
                        if c in value and c in self:
                            res[c] = res[c].replace(to_replace=src,
                                                    value=value[c],
                                                    inplace=False,
                                                    regex=regex)
                    return None if inplace else res

                # {'A': NA} -> 0
                elif not com.is_list_like(value):
                    for k, src in compat.iteritems(to_replace):
                        if k in self:
                            new_data = new_data.replace(to_replace=src,
                                                        value=value,
                                                        filter=[k],
                                                        inplace=inplace,
                                                        regex=regex)
                else:
                    raise TypeError('value argument must be scalar, dict, or '
                                    'Series')

            elif com.is_list_like(to_replace):  # [NA, ''] -> [0, 'missing']
                if com.is_list_like(value):
                    if len(to_replace) != len(value):
                        raise ValueError('Replacement lists must match '
                                         'in length. Expecting %d got %d ' %
                                         (len(to_replace), len(value)))

                    new_data = self._data.replace_list(src_list=to_replace,
                                                       dest_list=value,
                                                       inplace=inplace,
                                                       regex=regex)

                else:  # [NA, ''] -> 0
                    new_data = self._data.replace(to_replace=to_replace,
                                                  value=value,
                                                  inplace=inplace,
                                                  regex=regex)
            elif to_replace is None:
                if not (com.is_re_compilable(regex) or
                        com.is_list_like(regex) or
                        is_dictlike(regex)):
                    raise TypeError("'regex' must be a string or a compiled "
                                    "regular expression or a list or dict of "
                                    "strings or regular expressions, you "
                                    "passed a"
                                    " {0!r}".format(type(regex).__name__))
                return self.replace(regex, value, inplace=inplace, limit=limit,
                                    regex=True)
            else:

                # dest iterable dict-like
                if is_dictlike(value):  # NA -> {'A' : 0, 'B' : -1}
                    new_data = self._data

                    for k, v in compat.iteritems(value):
                        if k in self:
                            new_data = new_data.replace(to_replace=to_replace,
                                                        value=v,
                                                        filter=[k],
                                                        inplace=inplace,
                                                        regex=regex)

                elif not com.is_list_like(value):  # NA -> 0
                    new_data = self._data.replace(to_replace=to_replace, value=value,
                                                  inplace=inplace, regex=regex)
                else:
                    msg = ('Invalid "to_replace" type: '
                           '{0!r}').format(type(to_replace).__name__)
                    raise TypeError(msg)  # pragma: no cover

        new_data = new_data.convert(copy=not inplace, numeric=False)

        if inplace:
            self._update_inplace(new_data)
        else:
            return self._constructor(new_data).__finalize__(self)

    def interpolate(self, method='linear', axis=0, limit=None, inplace=False,
                    limit_direction='forward', downcast=None, **kwargs):
        """
        Interpolate values according to different methods.

        Please note that only ``method='linear'`` is supported for DataFrames/Series
        with a MultiIndex.

        Parameters
        ----------
        method : {'linear', 'time', 'index', 'values', 'nearest', 'zero',
                  'slinear', 'quadratic', 'cubic', 'barycentric', 'krogh',
                  'polynomial', 'spline' 'piecewise_polynomial', 'pchip'}

            * 'linear': ignore the index and treat the values as equally
              spaced. This is the only method supported on MultiIndexes.
              default
            * 'time': interpolation works on daily and higher resolution
              data to interpolate given length of interval
            * 'index', 'values': use the actual numerical values of the index
            * 'nearest', 'zero', 'slinear', 'quadratic', 'cubic',
              'barycentric', 'polynomial' is passed to
              ``scipy.interpolate.interp1d``. Both 'polynomial' and 'spline'
              require that you also specify an `order` (int),
              e.g. df.interpolate(method='polynomial', order=4).
              These use the actual numerical values of the index.
            * 'krogh', 'piecewise_polynomial', 'spline', and 'pchip' are all
              wrappers around the scipy interpolation methods of similar
              names. These use the actual numerical values of the index. See
              the scipy documentation for more on their behavior
              `here <http://docs.scipy.org/doc/scipy/reference/interpolate.html#univariate-interpolation>`__
              `and here <http://docs.scipy.org/doc/scipy/reference/tutorial/interpolate.html>`__

        axis : {0, 1}, default 0
            * 0: fill column-by-column
            * 1: fill row-by-row
        limit : int, default None.
            Maximum number of consecutive NaNs to fill.
        limit_direction : {'forward', 'backward', 'both'}, defaults to 'forward'
            If limit is specified, consecutive NaNs will be filled in this
            direction.

            .. versionadded:: 0.17.0

        inplace : bool, default False
            Update the NDFrame in place if possible.
        downcast : optional, 'infer' or None, defaults to None
            Downcast dtypes if possible.
        kwargs : keyword arguments to pass on to the interpolating function.

        Returns
        -------
        Series or DataFrame of same shape interpolated at the NaNs

        See Also
        --------
        reindex, replace, fillna

        Examples
        --------

        Filling in NaNs

        >>> s = pd.Series([0, 1, np.nan, 3])
        >>> s.interpolate()
        0    0
        1    1
        2    2
        3    3
        dtype: float64

        """
        if self.ndim > 2:
            raise NotImplementedError("Interpolate has not been implemented "
                                      "on Panel and Panel 4D objects.")

        if axis == 0:
            ax = self._info_axis_name
            _maybe_transposed_self = self
        elif axis == 1:
            _maybe_transposed_self = self.T
            ax = 1
        else:
            _maybe_transposed_self = self
        ax = _maybe_transposed_self._get_axis_number(ax)

        if _maybe_transposed_self.ndim == 2:
            alt_ax = 1 - ax
        else:
            alt_ax = ax

        if isinstance(_maybe_transposed_self.index, MultiIndex) and method != 'linear':
            raise ValueError("Only `method=linear` interpolation is supported "
                             "on MultiIndexes.")

        if _maybe_transposed_self._data.get_dtype_counts().get('object') == len(_maybe_transposed_self.T):
            raise TypeError("Cannot interpolate with all NaNs.")

        # create/use the index
        if method == 'linear':
            index = np.arange(len(_maybe_transposed_self._get_axis(alt_ax)))  # prior default
        else:
            index = _maybe_transposed_self._get_axis(alt_ax)

        if pd.isnull(index).any():
            raise NotImplementedError("Interpolation with NaNs in the index "
                                      "has not been implemented. Try filling "
                                      "those NaNs before interpolating.")
        new_data = _maybe_transposed_self._data.interpolate(
            method=method,
            axis=ax,
            index=index,
            values=_maybe_transposed_self,
            limit=limit,
            limit_direction=limit_direction,
            inplace=inplace,
            downcast=downcast,
            **kwargs
        )
        if inplace:
            if axis == 1:
                new_data = self._constructor(new_data).T._data
            self._update_inplace(new_data)
        else:
            res = self._constructor(new_data).__finalize__(self)
            if axis == 1:
                res = res.T
            return res

    #----------------------------------------------------------------------
    # Action Methods

    def isnull(self):
        """
        Return a boolean same-sized object indicating if the values are null

        See also
        --------
        notnull : boolean inverse of isnull
        """
        return isnull(self).__finalize__(self)

    def notnull(self):
        """Return a boolean same-sized object indicating if the values are
        not null

        See also
        --------
        isnull : boolean inverse of notnull
        """
        return notnull(self).__finalize__(self)

    def clip(self, lower=None, upper=None, out=None, axis=None):
        """
        Trim values at input threshold(s)

        Parameters
        ----------
        lower : float or array_like, default None
        upper : float or array_like, default None
        axis : int or string axis name, optional
            Align object with lower and upper along the given axis.

        Returns
        -------
        clipped : Series

        Examples
        --------
        >>> df
          0         1
        0  0.335232 -1.256177
        1 -1.367855  0.746646
        2  0.027753 -1.176076
        3  0.230930 -0.679613
        4  1.261967  0.570967
        >>> df.clip(-1.0, 0.5)
                  0         1
        0  0.335232 -1.000000
        1 -1.000000  0.500000
        2  0.027753 -1.000000
        3  0.230930 -0.679613
        4  0.500000  0.500000
        >>> t
        0   -0.3
        1   -0.2
        2   -0.1
        3    0.0
        4    0.1
        dtype: float64
        >>> df.clip(t, t + 1, axis=0)
                  0         1
        0  0.335232 -0.300000
        1 -0.200000  0.746646
        2  0.027753 -0.100000
        3  0.230930  0.000000
        4  1.100000  0.570967
        """
        if out is not None:  # pragma: no cover
            raise Exception('out argument is not supported yet')

        # GH 2747 (arguments were reversed)
        if lower is not None and upper is not None:
            if lib.isscalar(lower) and lib.isscalar(upper):
                lower, upper = min(lower, upper), max(lower, upper)

        result = self
        if lower is not None:
            result = result.clip_lower(lower, axis)
        if upper is not None:
            result = result.clip_upper(upper, axis)

        return result

    def clip_upper(self, threshold, axis=None):
        """
        Return copy of input with values above given value(s) truncated

        Parameters
        ----------
        threshold : float or array_like
        axis : int or string axis name, optional
            Align object with threshold along the given axis.

        See also
        --------
        clip

        Returns
        -------
        clipped : same type as input
        """
        if np.any(isnull(threshold)):
            raise ValueError("Cannot use an NA value as a clip threshold")

        subset = self.le(threshold, axis=axis) | isnull(self)
        return self.where(subset, threshold, axis=axis)

    def clip_lower(self, threshold, axis=None):
        """
        Return copy of the input with values below given value(s) truncated

        Parameters
        ----------
        threshold : float or array_like
        axis : int or string axis name, optional
            Align object with threshold along the given axis.

        See also
        --------
        clip

        Returns
        -------
        clipped : same type as input
        """
        if np.any(isnull(threshold)):
            raise ValueError("Cannot use an NA value as a clip threshold")

        subset = self.ge(threshold, axis=axis) | isnull(self)
        return self.where(subset, threshold, axis=axis)

    def groupby(self, by=None, axis=0, level=None, as_index=True, sort=True,
                group_keys=True, squeeze=False):
        """
        Group series using mapper (dict or key function, apply given function
        to group, return result as series) or by a series of columns

        Parameters
        ----------
        by : mapping function / list of functions, dict, Series, or tuple /
            list of column names.
            Called on each element of the object index to determine the groups.
            If a dict or Series is passed, the Series or dict VALUES will be
            used to determine the groups
        axis : int, default 0
        level : int, level name, or sequence of such, default None
            If the axis is a MultiIndex (hierarchical), group by a particular
            level or levels
        as_index : boolean, default True
            For aggregated output, return object with group labels as the
            index. Only relevant for DataFrame input. as_index=False is
            effectively "SQL-style" grouped output
        sort : boolean, default True
            Sort group keys. Get better performance by turning this off.
            Note this does not influence the order of observations within each group.
            groupby preserves the order of rows within each group.
        group_keys : boolean, default True
            When calling apply, add group keys to index to identify pieces
        squeeze : boolean, default False
            reduce the dimensionality of the return type if possible,
            otherwise return a consistent type

        Examples
        --------
        DataFrame results

        >>> data.groupby(func, axis=0).mean()
        >>> data.groupby(['col1', 'col2'])['col3'].mean()

        DataFrame with hierarchical index

        >>> data.groupby(['col1', 'col2']).mean()

        Returns
        -------
        GroupBy object

        """
        from pandas.core.groupby import groupby

        if level is None and by is None:
            raise TypeError("You have to supply one of 'by' and 'level'")
        axis = self._get_axis_number(axis)
        return groupby(self, by=by, axis=axis, level=level, as_index=as_index,
                       sort=sort, group_keys=group_keys, squeeze=squeeze)

    def asfreq(self, freq, method=None, how=None, normalize=False):
        """
        Convert all TimeSeries inside to specified frequency using DateOffset
        objects. Optionally provide fill method to pad/backfill missing values.

        Parameters
        ----------
        freq : DateOffset object, or string
        method : {'backfill', 'bfill', 'pad', 'ffill', None}
            Method to use for filling holes in reindexed Series
            pad / ffill: propagate last valid observation forward to next valid
            backfill / bfill: use NEXT valid observation to fill method
        how : {'start', 'end'}, default end
            For PeriodIndex only, see PeriodIndex.asfreq
        normalize : bool, default False
            Whether to reset output index to midnight

        Returns
        -------
        converted : type of caller
        """
        from pandas.tseries.resample import asfreq
        return asfreq(self, freq, method=method, how=how,
                      normalize=normalize)

    def at_time(self, time, asof=False):
        """
        Select values at particular time of day (e.g. 9:30AM)

        Parameters
        ----------
        time : datetime.time or string

        Returns
        -------
        values_at_time : type of caller
        """
        try:
            indexer = self.index.indexer_at_time(time, asof=asof)
            return self.take(indexer, convert=False)
        except AttributeError:
            raise TypeError('Index must be DatetimeIndex')

    def between_time(self, start_time, end_time, include_start=True,
                     include_end=True):
        """
        Select values between particular times of the day (e.g., 9:00-9:30 AM)

        Parameters
        ----------
        start_time : datetime.time or string
        end_time : datetime.time or string
        include_start : boolean, default True
        include_end : boolean, default True

        Returns
        -------
        values_between_time : type of caller
        """
        try:
            indexer = self.index.indexer_between_time(
                start_time, end_time, include_start=include_start,
                include_end=include_end)
            return self.take(indexer, convert=False)
        except AttributeError:
            raise TypeError('Index must be DatetimeIndex')

    def resample(self, rule, how=None, axis=0, fill_method=None,
                 closed=None, label=None, convention='start',
                 kind=None, loffset=None, limit=None, base=0):
        """
        Convenience method for frequency conversion and resampling of regular
        time-series data.

        Parameters
        ----------
        rule : string
            the offset string or object representing target conversion
        how : string
            method for down- or re-sampling, default to 'mean' for
            downsampling
        axis : int, optional, default 0
        fill_method : string, default None
            fill_method for upsampling
        closed : {'right', 'left'}
            Which side of bin interval is closed
        label : {'right', 'left'}
            Which bin edge label to label bucket with
        convention : {'start', 'end', 's', 'e'}
        kind : "period"/"timestamp"
        loffset : timedelta
            Adjust the resampled time labels
        limit : int, default None
            Maximum size gap to when reindexing with fill_method
        base : int, default 0
            For frequencies that evenly subdivide 1 day, the "origin" of the
            aggregated intervals. For example, for '5min' frequency, base could
            range from 0 through 4. Defaults to 0


        Examples
        --------

        Start by creating a series with 9 one minute timestamps.

        >>> index = pd.date_range('1/1/2000', periods=9, freq='T')
        >>> series = pd.Series(range(9), index=index)
        >>> series
        2000-01-01 00:00:00    0
        2000-01-01 00:01:00    1
        2000-01-01 00:02:00    2
        2000-01-01 00:03:00    3
        2000-01-01 00:04:00    4
        2000-01-01 00:05:00    5
        2000-01-01 00:06:00    6
        2000-01-01 00:07:00    7
        2000-01-01 00:08:00    8
        Freq: T, dtype: int64

        Downsample the series into 3 minute bins and sum the values
        of the timestamps falling into a bin.

        >>> series.resample('3T', how='sum')
        2000-01-01 00:00:00     3
        2000-01-01 00:03:00    12
        2000-01-01 00:06:00    21
        Freq: 3T, dtype: int64

        Downsample the series into 3 minute bins as above, but label each
        bin using the right edge instead of the left. Please note that the
        value in the bucket used as the label is not included in the bucket,
        which it labels. For example, in the original series the
        bucket ``2000-01-01 00:03:00`` contains the value 3, but the summed
        value in the resampled bucket with the label``2000-01-01 00:03:00``
        does not include 3 (if it did, the summed value would be 6, not 3).
        To include this value close the right side of the bin interval as
        illustrated in the example below this one.

        >>> series.resample('3T', how='sum', label='right')
        2000-01-01 00:03:00     3
        2000-01-01 00:06:00    12
        2000-01-01 00:09:00    21
        Freq: 3T, dtype: int64

        Downsample the series into 3 minute bins as above, but close the right
        side of the bin interval.

        >>> series.resample('3T', how='sum', label='right', closed='right')
        2000-01-01 00:00:00     0
        2000-01-01 00:03:00     6
        2000-01-01 00:06:00    15
        2000-01-01 00:09:00    15
        Freq: 3T, dtype: int64

        Upsample the series into 30 second bins.

        >>> series.resample('30S')[0:5] #select first 5 rows
        2000-01-01 00:00:00     0
        2000-01-01 00:00:30   NaN
        2000-01-01 00:01:00     1
        2000-01-01 00:01:30   NaN
        2000-01-01 00:02:00     2
        Freq: 30S, dtype: float64

        Upsample the series into 30 second bins and fill the ``NaN``
        values using the ``pad`` method.

        >>> series.resample('30S', fill_method='pad')[0:5]
        2000-01-01 00:00:00    0
        2000-01-01 00:00:30    0
        2000-01-01 00:01:00    1
        2000-01-01 00:01:30    1
        2000-01-01 00:02:00    2
        Freq: 30S, dtype: int64

        Upsample the series into 30 second bins and fill the
        ``NaN`` values using the ``bfill`` method.

        >>> series.resample('30S', fill_method='bfill')[0:5]
        2000-01-01 00:00:00    0
        2000-01-01 00:00:30    1
        2000-01-01 00:01:00    1
        2000-01-01 00:01:30    2
        2000-01-01 00:02:00    2
        Freq: 30S, dtype: int64

        Pass a custom function to ``how``.

        >>> def custom_resampler(array_like):
        ...     return np.sum(array_like)+5

        >>> series.resample('3T', how=custom_resampler)
        2000-01-01 00:00:00     8
        2000-01-01 00:03:00    17
        2000-01-01 00:06:00    26
        Freq: 3T, dtype: int64

        """

        from pandas.tseries.resample import TimeGrouper
        axis = self._get_axis_number(axis)
        sampler = TimeGrouper(rule, label=label, closed=closed, how=how,
                              axis=axis, kind=kind, loffset=loffset,
                              fill_method=fill_method, convention=convention,
                              limit=limit, base=base)
        return sampler.resample(self).__finalize__(self)

    def first(self, offset):
        """
        Convenience method for subsetting initial periods of time series data
        based on a date offset

        Parameters
        ----------
        offset : string, DateOffset, dateutil.relativedelta

        Examples
        --------
        ts.last('10D') -> First 10 days

        Returns
        -------
        subset : type of caller
        """
        from pandas.tseries.frequencies import to_offset
        if not isinstance(self.index, DatetimeIndex):
            raise NotImplementedError("'first' only supports a DatetimeIndex "
                                      "index")

        if len(self.index) == 0:
            return self

        offset = to_offset(offset)
        end_date = end = self.index[0] + offset

        # Tick-like, e.g. 3 weeks
        if not offset.isAnchored() and hasattr(offset, '_inc'):
            if end_date in self.index:
                end = self.index.searchsorted(end_date, side='left')

        return self.ix[:end]

    def last(self, offset):
        """
        Convenience method for subsetting final periods of time series data
        based on a date offset

        Parameters
        ----------
        offset : string, DateOffset, dateutil.relativedelta

        Examples
        --------
        ts.last('5M') -> Last 5 months

        Returns
        -------
        subset : type of caller
        """
        from pandas.tseries.frequencies import to_offset
        if not isinstance(self.index, DatetimeIndex):
            raise NotImplementedError("'last' only supports a DatetimeIndex "
                                      "index")

        if len(self.index) == 0:
            return self

        offset = to_offset(offset)

        start_date = start = self.index[-1] - offset
        start = self.index.searchsorted(start_date, side='right')
        return self.ix[start:]

    _shared_docs['align'] = (
        """
        Align two object on their axes with the
        specified join method for each axis Index

        Parameters
        ----------
        other : DataFrame or Series
        join : {'outer', 'inner', 'left', 'right'}, default 'outer'
        axis : allowed axis of the other object, default None
            Align on index (0), columns (1), or both (None)
        level : int or level name, default None
            Broadcast across a level, matching Index values on the
            passed MultiIndex level
        copy : boolean, default True
            Always returns new objects. If copy=False and no reindexing is
            required then original objects are returned.
        fill_value : scalar, default np.NaN
            Value to use for missing values. Defaults to NaN, but can be any
            "compatible" value
        method : str, default None
        limit : int, default None
        fill_axis : %(axes_single_arg)s, default 0
            Filling axis, method and limit
        broadcast_axis : %(axes_single_arg)s, default None
            Broadcast values along this axis, if aligning two objects of
            different dimensions

            .. versionadded:: 0.17.0

        Returns
        -------
        (left, right) : (%(klass)s, type of other)
            Aligned objects
        """
    )

    @Appender(_shared_docs['align'] % _shared_doc_kwargs)
    def align(self, other, join='outer', axis=None, level=None, copy=True,
              fill_value=None, method=None, limit=None, fill_axis=0,
              broadcast_axis=None):
        from pandas import DataFrame, Series
        method = com._clean_fill_method(method)

        if broadcast_axis == 1 and self.ndim != other.ndim:
            if isinstance(self, Series):
                # this means other is a DataFrame, and we need to broadcast self
                df = DataFrame(dict((c, self) for c in other.columns),
                               **other._construct_axes_dict())
                return df._align_frame(other, join=join, axis=axis, level=level,
                                       copy=copy, fill_value=fill_value,
                                       method=method, limit=limit,
                                       fill_axis=fill_axis)
            elif isinstance(other, Series):
                # this means self is a DataFrame, and we need to broadcast other
                df = DataFrame(dict((c, other) for c in self.columns),
                               **self._construct_axes_dict())
                return self._align_frame(df, join=join, axis=axis, level=level,
                                         copy=copy, fill_value=fill_value,
                                         method=method, limit=limit,
                                         fill_axis=fill_axis)

        if axis is not None:
            axis = self._get_axis_number(axis)
        if isinstance(other, DataFrame):
            return self._align_frame(other, join=join, axis=axis, level=level,
                                     copy=copy, fill_value=fill_value,
                                     method=method, limit=limit,
                                     fill_axis=fill_axis)
        elif isinstance(other, Series):
            return self._align_series(other, join=join, axis=axis, level=level,
                                      copy=copy, fill_value=fill_value,
                                      method=method, limit=limit,
                                      fill_axis=fill_axis)
        else:  # pragma: no cover
            raise TypeError('unsupported type: %s' % type(other))

    def _align_frame(self, other, join='outer', axis=None, level=None,
                     copy=True, fill_value=np.nan, method=None, limit=None,
                     fill_axis=0):
        # defaults
        join_index, join_columns = None, None
        ilidx, iridx = None, None
        clidx, cridx = None, None

        if axis is None or axis == 0:
            if not self.index.equals(other.index):
                join_index, ilidx, iridx = \
                    self.index.join(other.index, how=join, level=level,
                                    return_indexers=True)

        if axis is None or axis == 1:
            if not self.columns.equals(other.columns):
                join_columns, clidx, cridx = \
                    self.columns.join(other.columns, how=join, level=level,
                                      return_indexers=True)

        left = self._reindex_with_indexers({0: [join_index, ilidx],
                                            1: [join_columns, clidx]},
                                           copy=copy, fill_value=fill_value,
                                           allow_dups=True)
        right = other._reindex_with_indexers({0: [join_index, iridx],
                                              1: [join_columns, cridx]},
                                             copy=copy, fill_value=fill_value,
                                             allow_dups=True)

        if method is not None:
            left = left.fillna(axis=fill_axis, method=method, limit=limit)
            right = right.fillna(axis=fill_axis, method=method, limit=limit)

        return left.__finalize__(self), right.__finalize__(other)

    def _align_series(self, other, join='outer', axis=None, level=None,
                      copy=True, fill_value=None, method=None, limit=None,
                      fill_axis=0):
        from pandas import DataFrame

        # series/series compat
        if isinstance(self, ABCSeries) and isinstance(other, ABCSeries):
            if axis:
                raise ValueError('cannot align series to a series other than '
                                 'axis 0')

            # equal
            if  self.index.equals(other.index):
                join_index, lidx, ridx = None, None, None
            else:
                join_index, lidx, ridx = self.index.join(other.index, how=join,
                                                         level=level,
                                                         return_indexers=True)

            left = self._reindex_indexer(join_index, lidx, copy)
            right = other._reindex_indexer(join_index, ridx, copy)

        else:
            # one has > 1 ndim
            fdata = self._data
            if axis == 0:
                join_index = self.index
                lidx, ridx = None, None
                if not self.index.equals(other.index):
                    join_index, lidx, ridx = \
                        self.index.join(other.index, how=join, level=level,
                                        return_indexers=True)

                if lidx is not None:
                    fdata = fdata.reindex_indexer(join_index, lidx, axis=1)

            elif axis == 1:
                join_index = self.columns
                lidx, ridx = None, None
                if not self.columns.equals(other.index):
                    join_index, lidx, ridx = \
                        self.columns.join(other.index, how=join, level=level,
                                          return_indexers=True)

                if lidx is not None:
                    fdata = fdata.reindex_indexer(join_index, lidx, axis=0)
            else:
                raise ValueError('Must specify axis=0 or 1')

            if copy and fdata is self._data:
                fdata = fdata.copy()

            left = DataFrame(fdata)

            if ridx is None:
                right = other
            else:
                right = other.reindex(join_index, level=level)

        # fill
        fill_na = notnull(fill_value) or (method is not None)
        if fill_na:
            left = left.fillna(fill_value, method=method, limit=limit, axis=fill_axis)
            right = right.fillna(fill_value, method=method, limit=limit)
        return (left.__finalize__(self), right.__finalize__(other))

    _shared_docs['where'] = ("""
        Return an object of same shape as self and whose corresponding
        entries are from self where cond is %(cond)s and otherwise are from other.

        Parameters
        ----------
        cond : boolean %(klass)s or array
        other : scalar or %(klass)s
        inplace : boolean, default False
            Whether to perform the operation in place on the data
        axis : alignment axis if needed, default None
        level : alignment level if needed, default None
        try_cast : boolean, default False
            try to cast the result back to the input type (if possible),
        raise_on_error : boolean, default True
            Whether to raise on invalid data types (e.g. trying to where on
            strings)

        Returns
        -------
        wh : same type as caller
        """)
    @Appender(_shared_docs['where'] % dict(_shared_doc_kwargs, cond="True"))
    def where(self, cond, other=np.nan, inplace=False, axis=None, level=None,
              try_cast=False, raise_on_error=True):

        if isinstance(cond, NDFrame):
            cond, _ = cond.align(self, join='right', broadcast_axis=1)
        else:
            if not hasattr(cond, 'shape'):
                raise ValueError('where requires an ndarray like object for '
                                 'its condition')
            if cond.shape != self.shape:
                raise ValueError(
                    'Array conditional must be same shape as self')
            cond = self._constructor(cond, **self._construct_axes_dict())

        if inplace:
            cond = -(cond.fillna(True).astype(bool))
        else:
            cond = cond.fillna(False).astype(bool)

        # try to align
        try_quick = True
        if hasattr(other, 'align'):

            # align with me
            if other.ndim <= self.ndim:

                _, other = self.align(other, join='left',
                                      axis=axis, level=level,
                                      fill_value=np.nan)

                # if we are NOT aligned, raise as we cannot where index
                if (axis is None and
                        not all([other._get_axis(i).equals(ax)
                                 for i, ax in enumerate(self.axes)])):
                    raise InvalidIndexError

            # slice me out of the other
            else:
                raise NotImplemented(
                    "cannot align with a higher dimensional NDFrame"
                )

        elif is_list_like(other):

            if self.ndim == 1:

                # try to set the same dtype as ourselves
                try:
                    new_other = np.array(other, dtype=self.dtype)
                except ValueError:
                    new_other = np.array(other)

                # we can end up comparing integers and m8[ns]
                # which is a numpy no no
                is_i8 = com.needs_i8_conversion(self.dtype)
                if is_i8:
                    matches = False
                else:
                    matches = (new_other == np.array(other))

                if matches is False or not matches.all():

                    # coerce other to a common dtype if we can
                    if com.needs_i8_conversion(self.dtype):
                        try:
                            other = np.array(other, dtype=self.dtype)
                        except:
                            other = np.array(other)
                    else:
                        other = np.asarray(other)
                        other = np.asarray(other, dtype=np.common_type(other, new_other))

                    # we need to use the new dtype
                    try_quick = False
                else:
                    other = new_other
            else:

                other = np.array(other)

        if isinstance(other, np.ndarray):

            if other.shape != self.shape:

                if self.ndim == 1:

                    icond = cond.values

                    # GH 2745 / GH 4192
                    # treat like a scalar
                    if len(other) == 1:
                        other = np.array(other[0])

                    # GH 3235
                    # match True cond to other
                    elif len(cond[icond]) == len(other):

                        # try to not change dtype at first (if try_quick)
                        if try_quick:

                            try:
                                new_other = _values_from_object(self).copy()
                                new_other[icond] = other
                                other = new_other
                            except:
                                try_quick = False

                        # let's create a new (if we failed at the above
                        # or not try_quick
                        if not try_quick:

                            dtype, fill_value = _maybe_promote(other.dtype)
                            new_other = np.empty(len(icond), dtype=dtype)
                            new_other.fill(fill_value)
                            com._maybe_upcast_putmask(new_other, icond, other)
                            other = new_other

                    else:
                        raise ValueError(
                            'Length of replacements must equal series length')

                else:
                    raise ValueError('other must be the same shape as self '
                                     'when an ndarray')

            # we are the same shape, so create an actual object for alignment
            else:
                other = self._constructor(other, **self._construct_axes_dict())

        if axis is None:
            axis = 0

        if self.ndim == getattr(other, 'ndim', 0):
            align = True
        else:
            align = (self._get_axis_number(axis) == 1)

        block_axis = self._get_block_manager_axis(axis)

        if inplace:
            # we may have different type blocks come out of putmask, so
            # reconstruct the block manager

            self._check_inplace_setting(other)
            new_data = self._data.putmask(mask=cond, new=other, align=align,
                                          inplace=True, axis=block_axis,
                                          transpose=self._AXIS_REVERSED)
            self._update_inplace(new_data)

        else:
            new_data = self._data.where(other=other, cond=cond, align=align,
                                        raise_on_error=raise_on_error,
                                        try_cast=try_cast, axis=block_axis,
                                        transpose=self._AXIS_REVERSED)

            return self._constructor(new_data).__finalize__(self)

    @Appender(_shared_docs['where'] % dict(_shared_doc_kwargs, cond="False"))
    def mask(self, cond, other=np.nan, inplace=False, axis=None, level=None,
             try_cast=False, raise_on_error=True):
        return self.where(~cond, other=other, inplace=inplace, axis=axis,
            level=level, try_cast=try_cast, raise_on_error=raise_on_error)

    _shared_docs['shift'] = ("""
        Shift index by desired number of periods with an optional time freq

        Parameters
        ----------
        periods : int
            Number of periods to move, can be positive or negative
        freq : DateOffset, timedelta, or time rule string, optional
            Increment to use from datetools module or time rule (e.g. 'EOM').
            See Notes.
        axis : %(axes_single_arg)s

        Notes
        -----
        If freq is specified then the index values are shifted but the data
        is not realigned. That is, use freq if you would like to extend the
        index when shifting and preserve the original data.

        Returns
        -------
        shifted : %(klass)s
    """)
    @Appender(_shared_docs['shift'] % _shared_doc_kwargs)
    def shift(self, periods=1, freq=None, axis=0):
        if periods == 0:
            return self

        block_axis = self._get_block_manager_axis(axis)
        if freq is None:
            new_data = self._data.shift(periods=periods, axis=block_axis)
        else:
            return self.tshift(periods, freq)

        return self._constructor(new_data).__finalize__(self)

    def slice_shift(self, periods=1, axis=0):
        """
        Equivalent to `shift` without copying data. The shifted data will
        not include the dropped periods and the shifted axis will be smaller
        than the original.

        Parameters
        ----------
        periods : int
            Number of periods to move, can be positive or negative

        Notes
        -----
        While the `slice_shift` is faster than `shift`, you may pay for it
        later during alignment.

        Returns
        -------
        shifted : same type as caller
        """
        if periods == 0:
            return self

        if periods > 0:
            vslicer = slice(None, -periods)
            islicer = slice(periods, None)
        else:
            vslicer = slice(-periods, None)
            islicer = slice(None, periods)

        new_obj = self._slice(vslicer, axis=axis)
        shifted_axis = self._get_axis(axis)[islicer]
        new_obj.set_axis(axis, shifted_axis)

        return new_obj.__finalize__(self)

    def tshift(self, periods=1, freq=None, axis=0):
        """
        Shift the time index, using the index's frequency if available

        Parameters
        ----------
        periods : int
            Number of periods to move, can be positive or negative
        freq : DateOffset, timedelta, or time rule string, default None
            Increment to use from datetools module or time rule (e.g. 'EOM')
        axis : int or basestring
            Corresponds to the axis that contains the Index

        Notes
        -----
        If freq is not specified then tries to use the freq or inferred_freq
        attributes of the index. If neither of those attributes exist, a
        ValueError is thrown

        Returns
        -------
        shifted : NDFrame
        """

        index = self._get_axis(axis)
        if freq is None:
            freq = getattr(index, 'freq', None)

        if freq is None:
            freq = getattr(index, 'inferred_freq', None)

        if freq is None:
            msg = 'Freq was not given and was not set in the index'
            raise ValueError(msg)

        if periods == 0:
            return self

        if isinstance(freq, string_types):
            freq = datetools.to_offset(freq)

        block_axis = self._get_block_manager_axis(axis)
        if isinstance(index, PeriodIndex):
            orig_freq = datetools.to_offset(index.freq)
            if freq == orig_freq:
                new_data = self._data.copy()
                new_data.axes[block_axis] = index.shift(periods)
            else:
                msg = ('Given freq %s does not match PeriodIndex freq %s' %
                       (freq.rule_code, orig_freq.rule_code))
                raise ValueError(msg)
        else:
            new_data = self._data.copy()
            new_data.axes[block_axis] = index.shift(periods, freq)

        return self._constructor(new_data).__finalize__(self)

    def truncate(self, before=None, after=None, axis=None, copy=True):
        """Truncates a sorted NDFrame before and/or after some particular
        dates.

        Parameters
        ----------
        before : date
            Truncate before date
        after : date
            Truncate after date
        axis : the truncation axis, defaults to the stat axis
        copy : boolean, default is True,
            return a copy of the truncated section

        Returns
        -------
        truncated : type of caller
        """

        if axis is None:
            axis = self._stat_axis_number
        axis = self._get_axis_number(axis)
        ax = self._get_axis(axis)

        # if we have a date index, convert to dates, otherwise
        # treat like a slice
        if ax.is_all_dates:
            from pandas.tseries.tools import to_datetime
            before = to_datetime(before)
            after = to_datetime(after)

        if before is not None and after is not None:
            if before > after:
                raise ValueError('Truncate: %s must be after %s' %
                                 (after, before))

        slicer = [slice(None, None)] * self._AXIS_LEN
        slicer[axis] = slice(before, after)
        result = self.ix[tuple(slicer)]

        if isinstance(ax, MultiIndex):
            setattr(result, self._get_axis_name(axis),
                    ax.truncate(before, after))

        if copy:
            result = result.copy()

        return result

    def tz_convert(self, tz, axis=0, level=None, copy=True):
        """
        Convert tz-aware axis to target time zone.

        Parameters
        ----------
        tz : string or pytz.timezone object
        axis : the axis to convert
        level : int, str, default None
            If axis ia a MultiIndex, convert a specific level. Otherwise
            must be None
        copy : boolean, default True
            Also make a copy of the underlying data

        Returns
        -------

        Raises
        ------
        TypeError
            If the axis is tz-naive.
        """
        axis = self._get_axis_number(axis)
        ax = self._get_axis(axis)

        def _tz_convert(ax, tz):
            if not hasattr(ax, 'tz_convert'):
                if len(ax) > 0:
                    ax_name = self._get_axis_name(axis)
                    raise TypeError('%s is not a valid DatetimeIndex or PeriodIndex' %
                                    ax_name)
                else:
                    ax = DatetimeIndex([],tz=tz)
            else:
                ax = ax.tz_convert(tz)
            return ax

        # if a level is given it must be a MultiIndex level or
        # equivalent to the axis name
        if isinstance(ax, MultiIndex):
            level = ax._get_level_number(level)
            new_level = _tz_convert(ax.levels[level], tz)
            ax = ax.set_levels(new_level, level=level)
        else:
            if level not in (None, 0, ax.name):
                raise ValueError("The level {0} is not valid".format(level))
            ax =  _tz_convert(ax, tz)

        result = self._constructor(self._data, copy=copy)
        result.set_axis(axis,ax)
        return result.__finalize__(self)

    @deprecate_kwarg(old_arg_name='infer_dst', new_arg_name='ambiguous',
                     mapping={True: 'infer', False: 'raise'})
    def tz_localize(self, tz, axis=0, level=None, copy=True,
                    ambiguous='raise'):
        """
        Localize tz-naive TimeSeries to target time zone

        Parameters
        ----------
        tz : string or pytz.timezone object
        axis : the axis to localize
        level : int, str, default None
            If axis ia a MultiIndex, localize a specific level. Otherwise
            must be None
        copy : boolean, default True
            Also make a copy of the underlying data
        ambiguous : 'infer', bool-ndarray, 'NaT', default 'raise'
            - 'infer' will attempt to infer fall dst-transition hours based on order
            - bool-ndarray where True signifies a DST time, False designates
              a non-DST time (note that this flag is only applicable for ambiguous times)
            - 'NaT' will return NaT where there are ambiguous times
            - 'raise' will raise an AmbiguousTimeError if there are ambiguous times
        infer_dst : boolean, default False (DEPRECATED)
            Attempt to infer fall dst-transition hours based on order

        Returns
        -------

        Raises
        ------
        TypeError
            If the TimeSeries is tz-aware and tz is not None.
        """
        axis = self._get_axis_number(axis)
        ax = self._get_axis(axis)

        def _tz_localize(ax, tz, ambiguous):
            if not hasattr(ax, 'tz_localize'):
                if len(ax) > 0:
                    ax_name = self._get_axis_name(axis)
                    raise TypeError('%s is not a valid DatetimeIndex or PeriodIndex' %
                                    ax_name)
                else:
                    ax = DatetimeIndex([],tz=tz)
            else:
                ax = ax.tz_localize(tz, ambiguous=ambiguous)
            return ax

        # if a level is given it must be a MultiIndex level or
        # equivalent to the axis name
        if isinstance(ax, MultiIndex):
            level = ax._get_level_number(level)
            new_level = _tz_localize(ax.levels[level], tz, ambiguous)
            ax = ax.set_levels(new_level, level=level)
        else:
            if level not in (None, 0, ax.name):
                raise ValueError("The level {0} is not valid".format(level))
            ax =  _tz_localize(ax, tz, ambiguous)

        result = self._constructor(self._data, copy=copy)
        result.set_axis(axis,ax)
        return result.__finalize__(self)

    #----------------------------------------------------------------------
    # Numeric Methods
    def abs(self):
        """
        Return an object with absolute value taken. Only applicable to objects
        that are all numeric

        Returns
        -------
        abs: type of caller
        """
        return np.abs(self)

    _shared_docs['describe'] = """
        Generate various summary statistics, excluding NaN values.

        Parameters
        ----------
        percentiles : array-like, optional
            The percentiles to include in the output. Should all
            be in the interval [0, 1]. By default `percentiles` is
            [.25, .5, .75], returning the 25th, 50th, and 75th percentiles.
        include, exclude : list-like, 'all', or None (default)
            Specify the form of the returned result. Either:

            - None to both (default). The result will include only numeric-typed
              columns or, if none are, only categorical columns.
            - A list of dtypes or strings to be included/excluded.
              To select all numeric types use numpy numpy.number. To select
              categorical objects use type object. See also the select_dtypes
              documentation. eg. df.describe(include=['O'])
            - If include is the string 'all', the output column-set will
              match the input one.

        Returns
        -------
        summary: %(klass)s of summary statistics

        Notes
        -----
        The output DataFrame index depends on the requested dtypes:

        For numeric dtypes, it will include: count, mean, std, min,
        max, and lower, 50, and upper percentiles.

        For object dtypes (e.g. timestamps or strings), the index
        will include the count, unique, most common, and frequency of the
        most common. Timestamps also include the first and last items.

        For mixed dtypes, the index will be the union of the corresponding
        output types. Non-applicable entries will be filled with NaN.
        Note that mixed-dtype outputs can only be returned from mixed-dtype
        inputs and appropriate use of the include/exclude arguments.

        If multiple values have the highest count, then the
        `count` and `most common` pair will be arbitrarily chosen from
        among those with the highest count.

        The include, exclude arguments are ignored for Series.

        See also
        --------
        DataFrame.select_dtypes
        """

    @Appender(_shared_docs['describe'] % _shared_doc_kwargs)
    def describe(self, percentiles=None, include=None, exclude=None ):
        if self.ndim >= 3:
            msg = "describe is not implemented on on Panel or PanelND objects."
            raise NotImplementedError(msg)

        if percentiles is not None:
            # get them all to be in [0, 1]
            self._check_percentile(percentiles)
            percentiles = np.asarray(percentiles)
        else:
            percentiles = np.array([0.25, 0.5, 0.75])

        # median should always be included
        if (percentiles != 0.5).all():  # median isn't included
            lh = percentiles[percentiles < .5]
            uh = percentiles[percentiles > .5]
            percentiles = np.hstack([lh, 0.5, uh])

        def pretty_name(x):
            x *= 100
            if x == int(x):
                return '%.0f%%' % x
            else:
                return '%.1f%%' % x

        def describe_numeric_1d(series, percentiles):
            stat_index = (['count', 'mean', 'std', 'min'] +
                  [pretty_name(x) for x in percentiles] + ['max'])
            d = ([series.count(), series.mean(), series.std(), series.min()] +
                 [series.quantile(x) for x in percentiles] + [series.max()])
            return pd.Series(d, index=stat_index, name=series.name)


        def describe_categorical_1d(data):
            names = ['count', 'unique']
            objcounts = data.value_counts()
            result = [data.count(), len(objcounts[objcounts!=0])]
            if result[1] > 0:
                top, freq = objcounts.index[0], objcounts.iloc[0]

                if data.dtype == object or com.is_categorical_dtype(data.dtype):
                    names += ['top', 'freq']
                    result += [top, freq]

                elif com.is_datetime64_dtype(data):
                    asint = data.dropna().values.view('i8')
                    names += ['top', 'freq', 'first', 'last']
                    result += [lib.Timestamp(top), freq,
                               lib.Timestamp(asint.min()),
                               lib.Timestamp(asint.max())]

            return pd.Series(result, index=names, name=data.name)

        def describe_1d(data, percentiles):
            if com.is_numeric_dtype(data):
                return describe_numeric_1d(data, percentiles)
            elif com.is_timedelta64_dtype(data):
                return describe_numeric_1d(data, percentiles)
            else:
                return describe_categorical_1d(data)

        if self.ndim == 1:
            return describe_1d(self, percentiles)
        elif (include is None) and (exclude is None):
            if len(self._get_numeric_data()._info_axis) > 0:
                # when some numerics are found, keep only numerics
                data = self.select_dtypes(include=[np.number, np.bool])
            else:
                data = self
        elif include == 'all':
            if exclude != None:
                msg = "exclude must be None when include is 'all'"
                raise ValueError(msg)
            data = self
        else:
            data = self.select_dtypes(include=include, exclude=exclude)

        ldesc = [describe_1d(s, percentiles) for _, s in data.iteritems()]
        # set a convenient order for rows
        names = []
        ldesc_indexes = sorted([x.index for x in ldesc], key=len)
        for idxnames in ldesc_indexes:
            for name in idxnames:
                if name not in names:
                    names.append(name)
        d = pd.concat(ldesc, join_axes=pd.Index([names]), axis=1)
        return d

    def _check_percentile(self, q):
        """ Validate percentiles. Used by describe and quantile """

        msg = ("percentiles should all be in the interval [0, 1]. "
               "Try {0} instead.")
        q = np.asarray(q)
        if q.ndim == 0:
            if not 0 <= q <= 1:
                raise ValueError(msg.format(q / 100.0))
        else:
            if not all(0 <= qs <= 1 for qs in q):
                raise ValueError(msg.format(q / 100.0))
        return q

    _shared_docs['pct_change'] = """
        Percent change over given number of periods.

        Parameters
        ----------
        periods : int, default 1
            Periods to shift for forming percent change
        fill_method : str, default 'pad'
            How to handle NAs before computing percent changes
        limit : int, default None
            The number of consecutive NAs to fill before stopping
        freq : DateOffset, timedelta, or offset alias string, optional
            Increment to use from time series API (e.g. 'M' or BDay())

        Returns
        -------
        chg : %(klass)s

        Notes
        -----

        By default, the percentage change is calculated along the stat
        axis: 0, or ``Index``, for ``DataFrame`` and 1, or ``minor`` for
        ``Panel``. You can change this with the ``axis`` keyword argument.
        """

    @Appender(_shared_docs['pct_change'] % _shared_doc_kwargs)
    def pct_change(self, periods=1, fill_method='pad', limit=None, freq=None,
                   **kwargs):
        # TODO: Not sure if above is correct - need someone to confirm.
        axis = self._get_axis_number(kwargs.pop('axis', self._stat_axis_name))
        if fill_method is None:
            data = self
        else:
            data = self.fillna(method=fill_method, limit=limit)

        rs = (data.div(data.shift(periods=periods, freq=freq,
                                  axis=axis, **kwargs)) - 1)
        if freq is None:
            mask = com.isnull(_values_from_object(self))
            np.putmask(rs.values, mask, np.nan)
        return rs

    def _agg_by_level(self, name, axis=0, level=0, skipna=True, **kwargs):
        grouped = self.groupby(level=level, axis=axis)
        if hasattr(grouped, name) and skipna:
            return getattr(grouped, name)(**kwargs)
        axis = self._get_axis_number(axis)
        method = getattr(type(self), name)
        applyf = lambda x: method(x, axis=axis, skipna=skipna, **kwargs)
        return grouped.aggregate(applyf)

    @classmethod
    def _add_numeric_operations(cls):
        """ add the operations to the cls; evaluate the doc strings again """

        axis_descr = "{%s}" % ', '.join([
            "{0} ({1})".format(a, i) for i, a in enumerate(cls._AXIS_ORDERS)
        ])
        name = (cls._constructor_sliced.__name__
                if cls._AXIS_LEN > 1 else 'scalar')

        _num_doc = """

%(desc)s

Parameters
----------
axis : """ + axis_descr + """
skipna : boolean, default True
    Exclude NA/null values. If an entire row/column is NA, the result
    will be NA
level : int or level name, default None
        If the axis is a MultiIndex (hierarchical), count along a
        particular level, collapsing into a """ + name + """
numeric_only : boolean, default None
    Include only float, int, boolean data. If None, will attempt to use
    everything, then use only numeric data

Returns
-------
%(outname)s : """ + name + " or " + cls.__name__ + " (if level specified)\n"

        _bool_doc = """

%(desc)s

Parameters
----------
axis : """ + axis_descr + """
skipna : boolean, default True
    Exclude NA/null values. If an entire row/column is NA, the result
    will be NA
level : int or level name, default None
        If the axis is a MultiIndex (hierarchical), count along a
        particular level, collapsing into a """ + name + """
bool_only : boolean, default None
    Include only boolean data. If None, will attempt to use everything,
    then use only boolean data

Returns
-------
%(outname)s : """ + name + " or " + cls.__name__ + " (if level specified)\n"

        _cnum_doc = """

Parameters
----------
axis : """ + axis_descr + """
skipna : boolean, default True
    Exclude NA/null values. If an entire row/column is NA, the result
    will be NA

Returns
-------
%(outname)s : """ + name + "\n"

        def _make_stat_function(name, desc, f):

            @Substitution(outname=name, desc=desc)
            @Appender(_num_doc)
            def stat_func(self, axis=None, skipna=None, level=None,
                          numeric_only=None, **kwargs):
                if skipna is None:
                    skipna = True
                if axis is None:
                    axis = self._stat_axis_number
                if level is not None:
                    return self._agg_by_level(name, axis=axis, level=level,
                                              skipna=skipna)
                return self._reduce(f, name, axis=axis,
                                    skipna=skipna, numeric_only=numeric_only)
            stat_func.__name__ = name
            return stat_func

        cls.sum = _make_stat_function(
            'sum', 'Return the sum of the values for the requested axis',
            nanops.nansum)
        cls.mean = _make_stat_function(
            'mean', 'Return the mean of the values for the requested axis',
            nanops.nanmean)
        cls.skew = _make_stat_function(
            'skew',
            'Return unbiased skew over requested axis\nNormalized by N-1',
            nanops.nanskew)
        cls.kurt = _make_stat_function(
            'kurt',
            'Return unbiased kurtosis over requested axis using Fisher''s '
            'definition of\nkurtosis (kurtosis of normal == 0.0). Normalized '
            'by N-1\n',
            nanops.nankurt)
        cls.kurtosis = cls.kurt
        cls.prod = _make_stat_function(
            'prod', 'Return the product of the values for the requested axis',
            nanops.nanprod)
        cls.product = cls.prod
        cls.median = _make_stat_function(
            'median', 'Return the median of the values for the requested axis',
            nanops.nanmedian)
        cls.max = _make_stat_function('max', """
This method returns the maximum of the values in the object. If you
want the *index* of the maximum, use ``idxmax``. This is the
equivalent of the ``numpy.ndarray`` method ``argmax``.""", nanops.nanmax)
        cls.min = _make_stat_function('min', """
This method returns the minimum of the values in the object. If you
want the *index* of the minimum, use ``idxmin``. This is the
equivalent of the ``numpy.ndarray`` method ``argmin``.""", nanops.nanmin)

        def _make_logical_function(name, desc, f):

            @Substitution(outname=name, desc=desc)
            @Appender(_bool_doc)
            def logical_func(self, axis=None, bool_only=None, skipna=None,
                             level=None, **kwargs):
                if skipna is None:
                    skipna = True
                if axis is None:
                    axis = self._stat_axis_number
                if level is not None:
                    if bool_only is not None:
                        raise NotImplementedError(
                            "Option bool_only is not implemented with option "
                            "level.")
                    return self._agg_by_level(name, axis=axis, level=level,
                                              skipna=skipna)
                return self._reduce(f, axis=axis, skipna=skipna,
                                    numeric_only=bool_only, filter_type='bool',
                                    name=name)
            logical_func.__name__ = name
            return logical_func

        cls.any = _make_logical_function(
            'any', 'Return whether any element is True over requested axis',
            nanops.nanany)
        cls.all = _make_logical_function(
            'all', 'Return whether all elements are True over requested axis',
            nanops.nanall)

        @Substitution(outname='mad',
                      desc="Return the mean absolute deviation of the values "
                           "for the requested axis")
        @Appender(_num_doc)
        def mad(self,  axis=None, skipna=None, level=None):
            if skipna is None:
                skipna = True
            if axis is None:
                axis = self._stat_axis_number
            if level is not None:
                return self._agg_by_level('mad', axis=axis, level=level,
                                          skipna=skipna)

            data = self._get_numeric_data()
            if axis == 0:
                demeaned = data - data.mean(axis=0)
            else:
                demeaned = data.sub(data.mean(axis=1), axis=0)
            return np.abs(demeaned).mean(axis=axis, skipna=skipna)
        cls.mad = mad

        def _make_stat_function_ddof(name, desc, f):

            @Substitution(outname=name, desc=desc)
            @Appender(_num_doc)
            def stat_func(self, axis=None, skipna=None, level=None, ddof=1,
                          numeric_only=None, **kwargs):
                if skipna is None:
                    skipna = True
                if axis is None:
                    axis = self._stat_axis_number
                if level is not None:
                    return self._agg_by_level(name, axis=axis, level=level,
                                              skipna=skipna, ddof=ddof)
                return self._reduce(f, name, axis=axis,
                                    numeric_only=numeric_only,
                                    skipna=skipna, ddof=ddof)
            stat_func.__name__ = name
            return stat_func

        cls.sem = _make_stat_function_ddof(
            'sem',
            "Return unbiased standard error of the mean over "
            "requested axis.\n\nNormalized by N-1 by default. "
            "This can be changed using the ddof argument",
            nanops.nansem)
        cls.var = _make_stat_function_ddof(
            'var',
            "Return unbiased variance over requested "
            "axis.\n\nNormalized by N-1 by default. "
            "This can be changed using the ddof argument",
            nanops.nanvar)
        cls.std = _make_stat_function_ddof(
            'std',
            "Return unbiased standard deviation over requested "
            "axis.\n\nNormalized by N-1 by default. "
            "This can be changed using the ddof argument",
            nanops.nanstd)

        @Substitution(outname='compounded',
                      desc="Return the compound percentage of the values for "
                           "the requested axis")
        @Appender(_num_doc)
        def compound(self, axis=None, skipna=None, level=None):
            if skipna is None:
                skipna = True
            return (1 + self).prod(axis=axis, skipna=skipna, level=level) - 1
        cls.compound = compound

        def _make_cum_function(name, accum_func, mask_a, mask_b):

            @Substitution(outname=name)
            @Appender("Return cumulative {0} over requested axis.".format(name)
                      + _cnum_doc)
            def func(self, axis=None, dtype=None, out=None, skipna=True,
                     **kwargs):
                if axis is None:
                    axis = self._stat_axis_number
                else:
                    axis = self._get_axis_number(axis)

                y = _values_from_object(self).copy()

                if skipna and issubclass(y.dtype.type,
                                         (np.datetime64, np.timedelta64)):
                    result = accum_func(y, axis)
                    mask = isnull(self)
                    np.putmask(result, mask, pd.tslib.iNaT)
                elif skipna and not issubclass(y.dtype.type, (np.integer, np.bool_)):
                    mask = isnull(self)
                    np.putmask(y, mask, mask_a)
                    result = accum_func(y, axis)
                    np.putmask(result, mask, mask_b)
                else:
                    result = accum_func(y, axis)

                d = self._construct_axes_dict()
                d['copy'] = False
                return self._constructor(result, **d).__finalize__(self)

            func.__name__ = name
            return func

        cls.cummin = _make_cum_function(
            'min', lambda y, axis: np.minimum.accumulate(y, axis),
            np.inf, np.nan)
        cls.cumsum = _make_cum_function(
            'sum', lambda y, axis: y.cumsum(axis), 0., np.nan)
        cls.cumprod = _make_cum_function(
            'prod', lambda y, axis: y.cumprod(axis), 1., np.nan)
        cls.cummax = _make_cum_function(
            'max', lambda y, axis: np.maximum.accumulate(y, axis),
            -np.inf, np.nan)

# install the indexerse
for _name, _indexer in indexing.get_indexers_list():
    NDFrame._create_indexer(_name, _indexer)
