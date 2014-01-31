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
from pandas.core.indexing import _maybe_convert_indices
from pandas.tseries.index import DatetimeIndex
from pandas.tseries.period import PeriodIndex
from pandas.core.internals import BlockManager
import pandas.core.array as pa
import pandas.core.common as com
import pandas.core.datetools as datetools
from pandas import compat, _np_version_under1p7
from pandas.compat import map, zip, lrange, string_types, isidentifier
from pandas.core.common import (isnull, notnull, is_list_like,
                                _values_from_object, _maybe_promote, ABCSeries,
                                SettingWithCopyError, SettingWithCopyWarning)
import pandas.core.nanops as nanops
from pandas.util.decorators import Appender, Substitution
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
                       'is_copy', 'str', '_subtyp', '_index', '_default_kind',
                       '_default_fill_value','__array_struct__','__array_interface__']
    _internal_names_set = set(_internal_names)
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
            dtype = np.dtype(dtype)

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

        # do not copy BlockManager unless explicitly done
        if copy and dtype is None:
            mgr = mgr.copy()
        elif dtype is not None:
            # avoid copy if we can
            if len(mgr.blocks) > 1 or mgr.blocks[0].values.dtype != dtype:
                mgr = mgr.astype(dtype)
        return mgr

    #----------------------------------------------------------------------
    # Construction

    @property
    def _constructor(self):
        raise NotImplementedError

    def __unicode__(self):
        # unicode representation based upon iterating over self
        # (since, by definition, `PandasContainers` are iterable)
        prepr = '[%s]' % ','.join(map(com.pprint_thing, self))
        return '%s(%s)' % (self.__class__.__name__, prepr)

    def _local_dir(self):
        """ add the string-like attributes from the info_axis """
        return [c for c in self._info_axis
                if isinstance(c, string_types) and isidentifier(c)]

    @property
    def _constructor_sliced(self):
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

        axes = dict([(a, kwargs.get(a)) for a in self._AXIS_ORDERS])
        return axes, kwargs

    @classmethod
    def _from_axes(cls, data, axes):
        # for construction from BlockManager
        if isinstance(data, BlockManager):
            return cls(data)
        else:
            if cls._AXIS_REVERSED:
                axes = axes[::-1]
            d = cls._construct_axes_dict_from(cls, axes, copy=False)
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

    def _get_resolvers(self):
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
        "tuple of axis dimensions"
        return tuple(len(self._get_axis(a)) for a in self._AXIS_ORDERS)

    @property
    def axes(self):
        "index(es) of the NDFrame"
        # we do it this way because if we have reversed axes, then
        # the block manager shows then reversed
        return [self._get_axis(a) for a in self._AXIS_ORDERS]

    @property
    def ndim(self):
        "Number of axes / array dimensions"
        return self._data.ndim

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
        if kwargs.get('copy') or (len(args) and args[-1]):
            new_values = new_values.copy()
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
        copy = kwargs.get('copy', True)
        inplace = kwargs.get('inplace', False)

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
            result._data = result._data.rename(f, axis=baxis, copy=copy)
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
        arr = operator.neg(_values_from_object(self))
        return self._wrap_array(arr, self.axes, copy=False)

    def __invert__(self):
        arr = operator.inv(_values_from_object(self))
        return self._wrap_array(arr, self.axes, copy=False)

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
                      DeprecationWarning)
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

    def _wrap_array(self, arr, axes, copy=False):
        d = self._construct_axes_dict_from(self, axes, copy=copy)
        return self._constructor(arr, **d).__finalize__(self)

    def __array__(self, dtype=None):
        return _values_from_object(self)

    def __array_wrap__(self, result):
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
        return self._data

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
        path_or_buf : the path (string) or buffer to put the store
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

    def save(self, path):  # TODO remove in 0.14
        "Deprecated. Use to_pickle instead"
        import warnings
        from pandas.io.pickle import to_pickle
        warnings.warn("save is deprecated, use to_pickle", FutureWarning)
        return to_pickle(self, path)

    def load(self, path):  # TODO remove in 0.14
        "Deprecated. Use read_pickle instead."
        import warnings
        from pandas.io.pickle import read_pickle
        warnings.warn("load is deprecated, use pd.read_pickle", FutureWarning)
        return read_pickle(path)

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

            setattr(cls, name, property(_indexer))

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
        except (KeyError, ValueError):
            return default

    def __getitem__(self, item):
        return self._get_item_cache(item)

    def _get_item_cache(self, item):
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

    def _box_item_values(self, key, values):
        raise NotImplementedError

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
        cacher = getattr(self, '_cacher', None)
        return cacher is not None

    def _maybe_update_cacher(self, clear=False):
        """ see if we need to update our parent cacher
            if clear, then clear our cache """
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

        # check if we are a copy
        self._check_setitem_copy(stacklevel=5, t='referant')

        if clear:
            self._clear_item_cache()

    def _clear_item_cache(self, i=None):
        if i is not None:
            self._item_cache.pop(i, None)
        else:
            self._item_cache.clear()

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

    def _check_setitem_copy(self, stacklevel=4, t='setting'):
        """ validate if we are doing a settitem on a chained copy.

        If you call this function, be sure to set the stacklevel such that the
        user will see the error *at the level of setting*"""
        if self.is_copy:

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

            if t == 'referant':
                t = ("A value is trying to be set on a copy of a slice from a "
                     "DataFrame")
            else:
                t = ("A value is trying to be set on a copy of a slice from a "
                     "DataFrame.\nTry using .loc[row_index,col_indexer] = value "
                     "instead")
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

        # check/convert indicies here
        if convert:
            axis = self._get_axis_number(axis)
            indices = _maybe_convert_indices(
                indices, len(self._get_axis(axis)))

        baxis = self._get_block_manager_axis(axis)
        if baxis == 0:
            labels = self._get_axis(axis)
            new_items = labels.take(indices)
            new_data = self._data.reindex_axis(new_items, indexer=indices,
                                               axis=baxis)
        else:
            new_data = self._data.take(indices, axis=baxis)

        result = self._constructor(new_data).__finalize__(self)

        # maybe set copy if we didn't actually change the index
        if is_copy and not result._get_axis(axis).equals(self._get_axis(axis)):
            result._set_is_copy(self)

        return result

    def xs(self, key, axis=0, level=None, copy=True, drop_level=True):
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
        copy : boolean, default True
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
        >>> s = df.xs('a', copy=False)
        >>> s['A'] = 100
        >>> df
             A  B  C
        a  100  5  2
        b    4  0  9
        c    9  7  3


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

        """
        axis = self._get_axis_number(axis)
        labels = self._get_axis(axis)
        if level is not None:
            loc, new_ax = labels.get_loc_level(key, level=level,
                                               drop_level=drop_level)

            if not copy and not isinstance(loc, slice):
                raise ValueError('Cannot retrieve view (copy=False)')

            # level = 0
            loc_is_slice = isinstance(loc, slice)
            if not loc_is_slice:
                indexer = [slice(None)] * self.ndim
                indexer[axis] = loc
                indexer = tuple(indexer)
            else:
                indexer = loc
                lev_num = labels._get_level_number(level)
                if labels.levels[lev_num].inferred_type == 'integer':
                    indexer = self.index[loc]

            # select on the correct axis
            if axis == 1 and loc_is_slice:
                indexer = slice(None), indexer
            result = self.ix[indexer]
            setattr(result, result._get_axis_name(axis), new_ax)
            return result

        if axis == 1:
            data = self[key]
            if copy:
                data = data.copy()
            return data

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
            new_values, copy = self._data.fast_2d_xs(loc, copy=copy)
            result = Series(new_values, index=self.columns,
                            name=self.index[loc])

        else:
            result = self[loc]
            result.index = new_index

        result._set_is_copy(self)
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

    def reindex_like(self, other, method=None, copy=True, limit=None):
        """ return an object with matching indicies to myself

        Parameters
        ----------
        other : Object
        method : string or None
        copy : boolean, default True
        limit : int, default None
            Maximum size gap to forward or backward fill

        Notes
        -----
        Like calling s.reindex(index=other.index, columns=other.columns,
                               method=...)

        Returns
        -------
        reindexed : same as input
        """
        d = other._construct_axes_dict(method=method)
        return self.reindex(**d)

    def drop(self, labels, axis=0, level=None, inplace=False, **kwargs):
        """
        Return new object with labels in requested axis removed

        Parameters
        ----------
        labels : single label or list-like
        axis : int or axis name
        level : int or name, default None
            For MultiIndex
        inplace : bool, default False
            If True, do operation inplace and return None.

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
                new_axis = axis.drop(labels, level=level)
            else:
                new_axis = axis.drop(labels)
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
                indexer = -lib.ismember(axis.get_level_values(level),
                                        set(labels))
            else:
                indexer = -axis.isin(labels)

            slicer = [slice(None)] * self.ndim
            slicer[self._get_axis_number(axis_name)] = indexer

            result = self.ix[tuple(slicer)]

        if inplace:
            self._update_inplace(result)
        else:
            return result

    def _update_inplace(self, result):
        "replace self internals with result."
        # NOTE: This does *not* call __finalize__ and that's an explicit
        # decision that we may revisit in the future.
        self._reset_cache()
        self._clear_item_cache()
        self._data = getattr(result,'_data',result)
        self._maybe_update_cacher()

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

    def sort_index(self, axis=0, ascending=True):
        """
        Sort object by labels (along an axis)

        Parameters
        ----------
        axis : {0, 1}
            Sort index/rows versus columns
        ascending : boolean, default True
            Sort ascending vs. descending

        Returns
        -------
        sorted_obj : type of caller
        """
        axis = self._get_axis_number(axis)
        axis_name = self._get_axis_name(axis)
        labels = self._get_axis(axis)

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
        method : {'backfill', 'bfill', 'pad', 'ffill', None}, default None
            Method to use for filling holes in reindexed DataFrame
            pad / ffill: propagate last valid observation forward to next valid
            backfill / bfill: use NEXT valid observation to fill gap
        copy : boolean, default True
            Return a new object, even if the passed indexes are the same
        level : int or name
            Broadcast across a level, matching Index values on the
            passed MultiIndex level
        fill_value : scalar, default np.NaN
            Value to use for missing values. Defaults to NaN, but can be any
            "compatible" value
        limit : int, default None
            Maximum size gap to forward or backward fill
        takeable : boolean, default False
            treat the passed as positional values

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
        method = com._clean_fill_method(kwargs.get('method'))
        level = kwargs.get('level')
        copy = kwargs.get('copy', True)
        limit = kwargs.get('limit')
        fill_value = kwargs.get('fill_value', np.nan)
        takeable = kwargs.get('takeable', False)

        self._consolidate_inplace()

        # check if we are a multi reindex
        if self._needs_reindex_multi(axes, method, level):
            try:
                return self._reindex_multi(axes, copy, fill_value)
            except:
                pass

        # if all axes that are requested to reindex are equal, then only copy
        # if indicated must have index names equal here as well as values
        if all([self._get_axis(axis).identical(ax)
                for axis, ax in axes.items() if ax is not None]):
            if copy:
                return self.copy()
            return self

        # perform the reindex on the axes
        return self._reindex_axes(axes, level, limit,
                                  method, fill_value, copy,
                                  takeable=takeable).__finalize__(self)

    def _reindex_axes(self, axes, level, limit, method, fill_value, copy,
                      takeable=False):
        """ perform the reinxed for all the axes """
        obj = self
        for a in self._AXIS_ORDERS:
            labels = axes[a]
            if labels is None:
                continue

            # convert to an index if we are not a multi-selection
            if level is None:
                labels = _ensure_index(labels)

            axis = self._get_axis_number(a)
            ax = self._get_axis(a)
            try:
                new_index, indexer = ax.reindex(
                    labels, level=level, limit=limit, method=method,
                    takeable=takeable)
            except (ValueError):

                # catch trying to reindex a non-monotonic index with a
                # specialized indexer e.g. pad, so fallback to the regular
                # indexer this will show up on reindexing a not-naturally
                # ordering series,
                # e.g.
                # Series(
                #     [1,2,3,4], index=['a','b','c','d']
                # ).reindex(['c','b','g'], method='pad')
                new_index, indexer = ax.reindex(
                    labels, level=level, limit=limit, method=None,
                    takeable=takeable)

            obj = obj._reindex_with_indexers(
                {axis: [new_index, indexer]}, method=method,
                fill_value=fill_value, limit=limit, copy=copy)

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
        index : array-like, optional
            New labels / index to conform to. Preferably an Index object to
            avoid duplicating data
        axis : %(axes_single_arg)s
        method : {'backfill', 'bfill', 'pad', 'ffill', None}, default None
            Method to use for filling holes in reindexed object.
            pad / ffill: propagate last valid observation forward to next valid
            backfill / bfill: use NEXT valid observation to fill gap
        copy : boolean, default True
            Return a new object, even if the passed indexes are the same
        level : int or name
            Broadcast across a level, matching Index values on the
            passed MultiIndex level
        limit : int, default None
            Maximum size gap to forward or backward fill

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
        method = com._clean_fill_method(method)
        new_index, indexer = axis_values.reindex(
            labels, method, level, limit=limit, copy_if_needed=True)
        return self._reindex_with_indexers(
            {axis: [new_index, indexer]}, method=method, fill_value=fill_value,
            limit=limit, copy=copy).__finalize__(self)

    def _reindex_with_indexers(self, reindexers, method=None,
                               fill_value=np.nan, limit=None, copy=False,
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

            # reindex the axis
            if method is not None:
                new_data = new_data.reindex_axis(
                    index, indexer=indexer, method=method, axis=baxis,
                    fill_value=fill_value, limit=limit, copy=copy)

            elif indexer is not None:
                # TODO: speed up on homogeneous DataFrame objects
                indexer = com._ensure_int64(indexer)
                new_data = new_data.reindex_indexer(index, indexer, axis=baxis,
                                                    fill_value=fill_value,
                                                    allow_dups=allow_dups)

            elif (baxis == 0 and index is not None and
                    index is not new_data.axes[baxis]):
                new_data = new_data.reindex_items(index, copy=copy,
                                                  fill_value=fill_value)

            elif (baxis > 0 and index is not None and
                    index is not new_data.axes[baxis]):
                new_data = new_data.copy(deep=copy)
                new_data.set_axis(baxis, index)

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
            return self.select(lambda x: matcher.search(x) is not None,
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
        if n > l:
            n = l
        elif n < -l:
            n = -l
        return self.iloc[:n]

    def tail(self, n=5):
        """
        Returns last n rows
        """
        l = len(self)
        if l == 0 or n == 0:
            return self
        if n > l:
            n = l
        elif n < -l:
            n = -l
        return self.iloc[-n:]

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
        for name in self._metadata:
            object.__setattr__(self, name, getattr(other, name, None))
        return self

    def __getattr__(self, name):
        """After regular attribute access, try looking up the name of a the
        info.

        This allows simpler access to columns for interactive use.
        """
        if name in self._internal_names_set:
            return object.__getattribute__(self, name)
        elif name in self._metadata:
            return object.__getattribute__(self, name)
        else:
            if name in self._info_axis:
                return self[name]
            raise AttributeError("'%s' object has no attribute '%s'" %
                                 (type(self).__name__, name))

    def __setattr__(self, name, value):
        """After regular attribute access, try looking up the name of the info
        This allows simpler access to columns for interactive use."""
        if name in self._internal_names_set:
            object.__setattr__(self, name, value)
        elif name in self._metadata:
            return object.__setattr__(self, name, value)
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

    def _consolidate_inplace(self):
        f = lambda: self._data.consolidate()
        self._data = self._protect_consolidate(f)

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
            if cons_data is self._data:
                cons_data = cons_data.copy()
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

    def _protect_consolidate(self, f):
        blocks_before = len(self._data.blocks)
        result = f()
        if len(self._data.blocks) != blocks_before:
            self._clear_item_cache()
        return result

    def _get_numeric_data(self):
        return self._constructor(
            self._data.get_numeric_data()).__finalize__(self)

    def _get_bool_data(self):
        return self._constructor(self._data.get_bool_data()).__finalize__(self)

    #----------------------------------------------------------------------
    # Internal Interface Methods

    def as_matrix(self, columns=None):
        """
        Convert the frame to its Numpy-array matrix representation. Columns
        are presented in sorted order unless a specific list of columns is
        provided.

        NOTE: the dtype will be a lower-common-denominator dtype (implicit
              upcasting) that is to say if the dtypes (even of numeric types)
              are mixed, the one that accommodates all will be chosen use this
              with care if you are not dealing with the blocks

              e.g. if the dtypes are float16,float32         -> float32
                                     float16,float32,float64 -> float64
                                     int32,uint8             -> int32


        Returns
        -------
        values : ndarray
            If the caller is heterogeneous and contains booleans or objects,
            the result will be of dtype=object
        """
        self._consolidate_inplace()
        if self._AXIS_REVERSED:
            return self._data.as_matrix(columns).T
        return self._data.as_matrix(columns)

    @property
    def values(self):
        "Numpy representation of NDFrame"
        return self.as_matrix()

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
        return Series(self._data.get_dtypes(),index=self._info_axis)

    @property
    def ftypes(self):
        """
        Return the ftypes (indication of sparse/dense and dtype)
        in this object.
        """
        from pandas import Series
        return Series(self._data.get_ftypes(),index=self._info_axis)

    def as_blocks(self, columns=None):
        """
        Convert the frame to a dict of dtype -> Constructor Types that each has
        a homogeneous dtype.

        are presented in sorted order unless a specific list of columns is
        provided.

        NOTE: the dtypes of the blocks WILL BE PRESERVED HERE (unlike in
              as_matrix)

        Parameters
        ----------
        columns : array-like
            Specific column order

        Returns
        -------
        values : a list of Object
        """
        self._consolidate_inplace()

        bd = dict()
        for b in self._data.blocks:
            b = b.reindex_items_from(columns or b.items)
            bd[str(b.dtype)] = self._constructor(
                BlockManager([b], [b.items, self.index])).__finalize__(self)
        return bd

    @property
    def blocks(self):
        "Internal property, property synonym for as_blocks()"
        return self.as_blocks()

    def astype(self, dtype, copy=True, raise_on_error=True):
        """
        Cast object to input numpy.dtype
        Return a copy when copy = True (be really careful with this!)

        Parameters
        ----------
        dtype : numpy.dtype or Python type
        raise_on_error : raise on invalid input

        Returns
        -------
        casted : type of caller
        """

        mgr = self._data.astype(
            dtype, copy=copy, raise_on_error=raise_on_error)
        return self._constructor(mgr).__finalize__(self)

    def copy(self, deep=True):
        """
        Make a copy of this object

        Parameters
        ----------
        deep : boolean, default True
            Make a deep copy, i.e. also copy data

        Returns
        -------
        copy : type of caller
        """
        data = self._data
        if deep:
            data = data.copy()
        return self._constructor(data).__finalize__(self)

    def convert_objects(self, convert_dates=True, convert_numeric=False,
                        convert_timedeltas=True, copy=True):
        """
        Attempt to infer better dtype for object columns

        Parameters
        ----------
        convert_dates : if True, attempt to soft convert dates, if 'coerce',
            force conversion (and non-convertibles get NaT)
        convert_numeric : if True attempt to coerce to numbers (including
            strings), non-convertibles get NaN
        convert_timedeltas : if True, attempt to soft convert timedeltas, if 'coerce',
            force conversion (and non-convertibles get NaT)
        copy : Boolean, if True, return copy, default is True

        Returns
        -------
        converted : asm as input object
        """
        return self._constructor(
            self._data.convert(convert_dates=convert_dates,
                               convert_numeric=convert_numeric,
                               convert_timedeltas=convert_timedeltas,
                               copy=copy)).__finalize__(self)

    #----------------------------------------------------------------------
    # Filling NA's

    def fillna(self, value=None, method=None, axis=0, inplace=False,
               limit=None, downcast=None):
        """
        Fill NA/NaN values using the specified method

        Parameters
        ----------
        method : {'backfill', 'bfill', 'pad', 'ffill', None}, default None
            Method to use for filling holes in reindexed Series
            pad / ffill: propagate last valid observation forward to next valid
            backfill / bfill: use NEXT valid observation to fill gap
        value : scalar, dict, or Series
            Value to use to fill holes (e.g. 0), alternately a dict/Series of
            values specifying which value to use for each index (for a Series) or
            column (for a DataFrame). (values not in the dict/Series will not be
            filled). This value cannot be a list.
        axis : {0, 1}, default 0
            * 0: fill column-by-column
            * 1: fill row-by-row
        inplace : boolean, default False
            If True, fill in place. Note: this will modify any
            other views on this object, (e.g. a no-copy slice for a column in a
            DataFrame).
        limit : int, default None
            Maximum size gap to forward or backward fill
        downcast : dict, default is None
            a dict of item->dtype of what to downcast if possible,
            or the string 'infer' which will try to downcast to an appropriate
            equal type (e.g. float64 to int64 if possible)

        See also
        --------
        reindex, asfreq

        Returns
        -------
        filled : same type as caller
        """
        if isinstance(value, (list, tuple)):
            raise TypeError('"value" parameter must be a scalar or dict, but '
                            'you passed a "{0}"'.format(type(value).__name__))
        self._consolidate_inplace()

        axis = self._get_axis_number(axis)
        method = com._clean_fill_method(method)

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

            if self.ndim == 1 and value is not None:
                if isinstance(value, (dict, com.ABCSeries)):
                    from pandas import Series
                    value = Series(value)

                new_data = self._data.fillna(value, inplace=inplace,
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
                    obj.fillna(v, inplace=True)
                return result
            else:
                new_data = self._data.fillna(value, inplace=inplace,
                                             downcast=downcast)

        if inplace:
            self._update_inplace(new_data)
        else:
            return self._constructor(new_data).__finalize__(self)

    def ffill(self, axis=0, inplace=False, limit=None, downcast=None):
        "Synonym for NDFrame.fillna(method='ffill')"
        return self.fillna(method='ffill', axis=axis, inplace=inplace,
                           limit=limit, downcast=downcast)

    def bfill(self, axis=0, inplace=False, limit=None, downcast=None):
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
                    to_rep_dict[k] = v.keys()
                    value_dict[k] = v.values()

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
                    new_data = self._data
                    for c, src in compat.iteritems(to_replace):
                        if c in value and c in self:
                            new_data = new_data.replace(src, value[c],
                                                        filter=[c],
                                                        inplace=inplace,
                                                        regex=regex)

                # {'A': NA} -> 0
                elif not com.is_list_like(value):
                    new_data = self._data
                    for k, src in compat.iteritems(to_replace):
                        if k in self:
                            new_data = new_data.replace(src, value,
                                                        filter=[k],
                                                        inplace=inplace,
                                                        regex=regex)
                else:
                    raise TypeError('Fill value must be scalar, dict, or '
                                    'Series')

            elif com.is_list_like(to_replace):  # [NA, ''] -> [0, 'missing']
                if com.is_list_like(value):
                    if len(to_replace) != len(value):
                        raise ValueError('Replacement lists must match '
                                         'in length. Expecting %d got %d ' %
                                         (len(to_replace), len(value)))

                    new_data = self._data.replace_list(to_replace, value,
                                                       inplace=inplace,
                                                       regex=regex)

                else:  # [NA, ''] -> 0
                    new_data = self._data.replace(to_replace, value,
                                                  inplace=inplace, regex=regex)
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
                            new_data = new_data.replace(to_replace, v,
                                                        filter=[k],
                                                        inplace=inplace,
                                                        regex=regex)

                elif not com.is_list_like(value):  # NA -> 0
                    new_data = self._data.replace(to_replace, value,
                                                  inplace=inplace, regex=regex)
                else:
                    msg = ('Invalid "to_replace" type: '
                           '{0!r}').format(type(to_replace).__name__)
                    raise TypeError(msg)  # pragma: no cover

        new_data = new_data.convert(copy=not inplace, convert_numeric=False)

        if inplace:
            self._update_inplace(new_data)
        else:
            return self._constructor(new_data).__finalize__(self)

    def interpolate(self, method='linear', axis=0, limit=None, inplace=False,
                    downcast='infer', **kwargs):
        """
        Interpolate values according to different methods.

        Parameters
        ----------
        method : {'linear', 'time', 'values', 'index' 'nearest', 'zero',
                  'slinear', 'quadratic', 'cubic', 'barycentric', 'krogh',
                  'polynomial', 'spline' 'piecewise_polynomial', 'pchip'}

            * 'linear': ignore the index and treat the values as equally
              spaced. default
            * 'time': interpolation works on daily and higher resolution
              data to interpolate given length of interval
            * 'index': use the actual numerical values of the index
            * 'nearest', 'zero', 'slinear', 'quadratic', 'cubic',
              'barycentric', 'polynomial' is passed to
              `scipy.interpolate.interp1d` with the order given both
              'polynomial' and 'spline' requre that you also specify and order
              (int) e.g. df.interpolate(method='polynomial', order=4)
            * 'krogh', 'piecewise_polynomial', 'spline', and 'pchip' are all
              wrappers around the scipy interpolation methods of similar
              names. See the scipy documentation for more on their behavior:
              http://docs.scipy.org/doc/scipy/reference/interpolate.html#univariate-interpolation
              http://docs.scipy.org/doc/scipy/reference/tutorial/interpolate.html

        axis : {0, 1}, default 0
            * 0: fill column-by-column
            * 1: fill row-by-row
        limit : int, default None.
            Maximum number of consecutive NaNs to fill.
        inplace : bool, default False
            Update the NDFrame in place if possible.
        downcast : optional, 'infer' or None, defaults to 'infer'
            Downcast dtypes if possible.

        Returns
        -------
        Series or DataFrame of same shape interpolated at the NaNs

        See Also
        --------
        reindex, replace, fillna

        Examples
        --------

        # Filling in NaNs:
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
        elif axis == 1:
            self = self.T
            ax = 1
        ax = self._get_axis_number(ax)

        if self.ndim == 2:
            alt_ax = 1 - ax
        else:
            alt_ax = ax

        if isinstance(self.index, MultiIndex) and method != 'linear':
            raise ValueError("Only `method=linear` interpolation is supported "
                             "on MultiIndexes.")

        if self._data.get_dtype_counts().get('object') == len(self.T):
            raise TypeError("Cannot interpolate with all NaNs.")

        # create/use the index
        if method == 'linear':
            index = np.arange(len(self._get_axis(alt_ax)))  # prior default
        else:
            index = self._get_axis(alt_ax)

        if pd.isnull(index).any():
            raise NotImplementedError("Interpolation with NaNs in the index "
                                      "has not been implemented. Try filling "
                                      "those NaNs before interpolating.")
        new_data = self._data.interpolate(method=method,
                                          axis=ax,
                                          index=index,
                                          values=self,
                                          limit=limit,
                                          inplace=inplace,
                                          downcast=downcast,
                                          **kwargs)

        if inplace:
            if axis == 1:
                self._update_inplace(new_data)
                self = self.T
            else:
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
        """
        return isnull(self).__finalize__(self)

    def notnull(self):
        """Return a boolean same-sized object indicating if the values are
        not null
        """
        return notnull(self).__finalize__(self)

    def clip(self, lower=None, upper=None, out=None):
        """
        Trim values at input threshold(s)

        Parameters
        ----------
        lower : float, default None
        upper : float, default None

        Returns
        -------
        clipped : Series
        """
        if out is not None:  # pragma: no cover
            raise Exception('out argument is not supported yet')

        # GH 2747 (arguments were reversed)
        if lower is not None and upper is not None:
            lower, upper = min(lower, upper), max(lower, upper)

        result = self
        if lower is not None:
            result = result.clip_lower(lower)
        if upper is not None:
            result = result.clip_upper(upper)

        return result

    def clip_upper(self, threshold):
        """
        Return copy of input with values above given value truncated

        See also
        --------
        clip

        Returns
        -------
        clipped : same type as input
        """
        if isnull(threshold):
            raise ValueError("Cannot use an NA value as a clip threshold")

        return self.where((self <= threshold) | isnull(self), threshold)

    def clip_lower(self, threshold):
        """
        Return copy of the input with values below given value truncated

        See also
        --------
        clip

        Returns
        -------
        clipped : same type as input
        """
        if isnull(threshold):
            raise ValueError("Cannot use an NA value as a clip threshold")

        return self.where((self >= threshold) | isnull(self), threshold)

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
            Sort group keys. Get better performance by turning this off
        group_keys : boolean, default True
            When calling apply, add group keys to index to identify pieces
        squeeze : boolean, default False
            reduce the dimensionaility of the return type if possible,
            otherwise return a consistent type

        Examples
        --------
        # DataFrame result
        >>> data.groupby(func, axis=0).mean()

        # DataFrame result
        >>> data.groupby(['col1', 'col2'])['col3'].mean()

        # DataFrame with hierarchical index
        >>> data.groupby(['col1', 'col2']).mean()

        Returns
        -------
        GroupBy object

        """

        from pandas.core.groupby import groupby
        axis = self._get_axis_number(axis)
        return groupby(self, by, axis=axis, level=level, as_index=as_index,
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
            raise NotImplementedError

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
            raise NotImplementedError

        if len(self.index) == 0:
            return self

        offset = to_offset(offset)

        start_date = start = self.index[-1] - offset
        start = self.index.searchsorted(start_date, side='right')
        return self.ix[start:]

    def align(self, other, join='outer', axis=None, level=None, copy=True,
              fill_value=None, method=None, limit=None, fill_axis=0):
        """
        Align two object on their axes with the
        specified join method for each axis Index

        Parameters
        ----------
        other : DataFrame or Series
        join : {'outer', 'inner', 'left', 'right'}, default 'outer'
        axis : allowed axis of the other object, default None
            Align on index (0), columns (1), or both (None)
        level : int or name
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
        fill_axis : {0, 1}, default 0
            Filling axis, method and limit

        Returns
        -------
        (left, right) : (type of input, type of other)
            Aligned objects
        """
        from pandas import DataFrame, Series
        method = com._clean_fill_method(method)

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

        left = self._reindex_with_indexers({0: [join_index,   ilidx],
                                            1: [join_columns, clidx]},
                                           copy=copy, fill_value=fill_value,
                                           allow_dups=True)
        right = other._reindex_with_indexers({0: [join_index,   iridx],
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

            join_index, lidx, ridx = self.index.join(other.index, how=join,
                                                     level=level,
                                                     return_indexers=True)

            left_result = self._reindex_indexer(join_index, lidx, copy)
            right_result = other._reindex_indexer(join_index, ridx, copy)

        else:

            # one has > 1 ndim
            fdata = self._data
            if axis == 0:
                join_index = self.index
                lidx, ridx = None, None
                if not self.index.equals(other.index):
                    join_index, lidx, ridx = self.index.join(
                        other.index, how=join, return_indexers=True)

                if lidx is not None:
                    fdata = fdata.reindex_indexer(join_index, lidx, axis=1)
            elif axis == 1:
                join_index = self.columns
                lidx, ridx = None, None
                if not self.columns.equals(other.index):
                    join_index, lidx, ridx = \
                        self.columns.join(other.index, how=join,
                                          return_indexers=True)

                if lidx is not None:
                    fdata = fdata.reindex_indexer(join_index, lidx, axis=0)
            else:
                raise ValueError('Must specify axis=0 or 1')

            if copy and fdata is self._data:
                fdata = fdata.copy()

            left_result = DataFrame(fdata)
            right_result = other if ridx is None else other.reindex(join_index)

        # fill
        fill_na = notnull(fill_value) or (method is not None)
        if fill_na:
            return (left_result.fillna(fill_value, method=method, limit=limit,
                                       axis=fill_axis),
                    right_result.fillna(fill_value, method=method,
                                        limit=limit))
        else:
            return (left_result.__finalize__(self),
                    right_result.__finalize__(other))

    def where(self, cond, other=np.nan, inplace=False, axis=None, level=None,
              try_cast=False, raise_on_error=True):
        """
        Return an object of same shape as self and whose corresponding
        entries are from self where cond is True and otherwise are from other.

        Parameters
        ----------
        cond : boolean NDFrame or array
        other : scalar or NDFrame
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
        """
        if isinstance(cond, NDFrame):
            cond = cond.reindex(**self._construct_axes_dict())
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
                new_other = np.array(other, dtype=self.dtype)
                if not (new_other == np.array(other)).all():
                    other = np.array(other)

                    # we can't use our existing dtype
                    # because of incompatibilities
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

        if inplace:
            # we may have different type blocks come out of putmask, so
            # reconstruct the block manager
            new_data = self._data.putmask(cond, other, align=axis is None,
                                          inplace=True)
            self._update_inplace(new_data)

        else:
            new_data = self._data.where(other, cond, align=axis is None,
                                        raise_on_error=raise_on_error,
                                        try_cast=try_cast)

            return self._constructor(new_data).__finalize__(self)

    def mask(self, cond):
        """
        Returns copy whose values are replaced with nan if the
        inverted condition is True

        Parameters
        ----------
        cond : boolean NDFrame or array

        Returns
        -------
        wh: same as input
        """
        return self.where(~cond, np.nan)

    def shift(self, periods=1, freq=None, axis=0, **kwds):
        """
        Shift index by desired number of periods with an optional time freq

        Parameters
        ----------
        periods : int
            Number of periods to move, can be positive or negative
        freq : DateOffset, timedelta, or time rule string, optional
            Increment to use from datetools module or time rule (e.g. 'EOM')

        Notes
        -----
        If freq is specified then the index values are shifted but the data
        if not realigned

        Returns
        -------
        shifted : same type as caller
        """
        if periods == 0:
            return self

        if freq is None and not len(kwds):
            block_axis = self._get_block_manager_axis(axis)
            indexer = com._shift_indexer(len(self), periods)
            new_data = self._data.shift(indexer, periods, axis=block_axis)
        else:
            return self.tshift(periods, freq, **kwds)

        return self._constructor(new_data).__finalize__(self)

    def tshift(self, periods=1, freq=None, axis=0, **kwds):
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
        from pandas.core.datetools import _resolve_offset

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

        offset = _resolve_offset(freq, kwds)

        if isinstance(offset, string_types):
            offset = datetools.to_offset(offset)

        block_axis = self._get_block_manager_axis(axis)
        if isinstance(index, PeriodIndex):
            orig_offset = datetools.to_offset(index.freq)
            if offset == orig_offset:
                new_data = self._data.copy()
                new_data.axes[block_axis] = index.shift(periods)
            else:
                msg = ('Given freq %s does not match PeriodIndex freq %s' %
                       (offset.rule_code, orig_offset.rule_code))
                raise ValueError(msg)
        else:
            new_data = self._data.copy()
            new_data.axes[block_axis] = index.shift(periods, offset)

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

    def tz_convert(self, tz, axis=0, copy=True):
        """
        Convert TimeSeries to target time zone. If it is time zone naive, it
        will be localized to the passed time zone.

        Parameters
        ----------
        tz : string or pytz.timezone object
        copy : boolean, default True
            Also make a copy of the underlying data

        Returns
        -------
        """
        axis = self._get_axis_number(axis)
        ax = self._get_axis(axis)

        if not hasattr(ax, 'tz_convert'):
            ax_name = self._get_axis_name(axis)
            raise TypeError('%s is not a valid DatetimeIndex or PeriodIndex' %
                            ax_name)

        new_data = self._data
        if copy:
            new_data = new_data.copy()

        new_obj = self._constructor(new_data)
        new_ax = ax.tz_convert(tz)

        if axis == 0:
            new_obj._set_axis(1, new_ax)
        elif axis == 1:
            new_obj._set_axis(0, new_ax)
            self._clear_item_cache()

        return new_obj.__finalize__(self)

    def tz_localize(self, tz, axis=0, copy=True, infer_dst=False):
        """
        Localize tz-naive TimeSeries to target time zone

        Parameters
        ----------
        tz : string or pytz.timezone object
        copy : boolean, default True
            Also make a copy of the underlying data
        infer_dst : boolean, default False
            Attempt to infer fall dst-transition times based on order

        Returns
        -------
        """
        axis = self._get_axis_number(axis)
        ax = self._get_axis(axis)

        if not hasattr(ax, 'tz_localize'):
            ax_name = self._get_axis_name(axis)
            raise TypeError('%s is not a valid DatetimeIndex or PeriodIndex' %
                            ax_name)

        new_data = self._data
        if copy:
            new_data = new_data.copy()

        new_obj = self._constructor(new_data)
        new_ax = ax.tz_localize(tz, infer_dst=infer_dst)

        if axis == 0:
            new_obj._set_axis(1, new_ax)
        elif axis == 1:
            new_obj._set_axis(0, new_ax)
            self._clear_item_cache()

        return new_obj.__finalize__(self)

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

        # suprimo numpy 1.6 hacking
        # for timedeltas
        if _np_version_under1p7:

            def _convert_timedeltas(x):
                if x.dtype.kind == 'm':
                    return np.abs(x.view('i8')).astype(x.dtype)
                return np.abs(x)

            if self.ndim == 1:
                return _convert_timedeltas(self)
            elif self.ndim == 2:
                return  self.apply(_convert_timedeltas)

        return np.abs(self)

    def pct_change(self, periods=1, fill_method='pad', limit=None, freq=None,
                   **kwds):
        """
        Percent change over given number of periods

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
        chg : same type as caller
        """
        # TODO: Not sure if above is correct - need someone to confirm.
        if fill_method is None:
            data = self
        else:
            data = self.fillna(method=fill_method, limit=limit)
        rs = data / data.shift(periods=periods, freq=freq, **kwds) - 1
        if freq is None:
            mask = com.isnull(_values_from_object(self))
            np.putmask(rs.values, mask, np.nan)
        return rs

    def _agg_by_level(self, name, axis=0, level=0, skipna=True, **kwds):
        grouped = self.groupby(level=level, axis=axis)
        if hasattr(grouped, name) and skipna:
            return getattr(grouped, name)(**kwds)
        axis = self._get_axis_number(axis)
        method = getattr(type(self), name)
        applyf = lambda x: method(x, axis=axis, skipna=skipna, **kwds)
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
level : int, default None
        If the axis is a MultiIndex (hierarchical), count along a
        particular level, collapsing into a """ + name + """
numeric_only : boolean, default None
    Include only float, int, boolean data. If None, will attempt to use
    everything, then use only numeric data

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
                return self._reduce(f, axis=axis,
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
            'Return unbiased kurtosis over requested axis\nNormalized by N-1',
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

        @Substitution(outname='mad',
                      desc="Return the mean absolute deviation of the values "
                           "for the requested axis")
        @Appender(_num_doc)
        def mad(self,  axis=None, skipna=None, level=None, **kwargs):
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

        @Substitution(outname='variance',
                      desc="Return unbiased variance over requested "
                           "axis\nNormalized by N-1")
        @Appender(_num_doc)
        def var(self, axis=None, skipna=None, level=None, ddof=1, **kwargs):
            if skipna is None:
                skipna = True
            if axis is None:
                axis = self._stat_axis_number
            if level is not None:
                return self._agg_by_level('var', axis=axis, level=level,
                                          skipna=skipna, ddof=ddof)

            return self._reduce(nanops.nanvar, axis=axis, skipna=skipna,
                                ddof=ddof)
        cls.var = var

        @Substitution(outname='stdev',
                      desc="Return unbiased standard deviation over requested "
                           "axis\nNormalized by N-1")
        @Appender(_num_doc)
        def std(self, axis=None, skipna=None, level=None, ddof=1, **kwargs):
            if skipna is None:
                skipna = True
            if axis is None:
                axis = self._stat_axis_number
            if level is not None:
                return self._agg_by_level('std', axis=axis, level=level,
                                          skipna=skipna, ddof=ddof)
            result = self.var(axis=axis, skipna=skipna, ddof=ddof)
            if getattr(result, 'ndim', 0) > 0:
                return result.apply(np.sqrt)
            return np.sqrt(result)
        cls.std = std

        @Substitution(outname='compounded',
                      desc="Return the compound percentage of the values for "
                           "the requested axis")
        @Appender(_num_doc)
        def compound(self, axis=None, skipna=None, level=None, **kwargs):
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
                if not issubclass(y.dtype.type, (np.integer, np.bool_)):
                    mask = isnull(self)
                    if skipna:
                        np.putmask(y, mask, mask_a)
                    result = accum_func(y, axis)
                    if skipna:
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
