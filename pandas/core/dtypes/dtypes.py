""" define extension dtypes """
import re
import warnings

import numpy as np
import pytz

from pandas._libs.interval import Interval
from pandas._libs.tslibs import NaT, Period, Timestamp, timezones

from pandas.core.dtypes.generic import ABCCategoricalIndex, ABCIndexClass

from pandas import compat

from .base import ExtensionDtype, _DtypeOpsMixin
from .inference import is_list_like


def register_extension_dtype(cls):
    """Class decorator to register an ExtensionType with pandas.

    .. versionadded:: 0.24.0

    This enables operations like ``.astype(name)`` for the name
    of the ExtensionDtype.

    Examples
    --------
    >>> from pandas.api.extensions import register_extension_dtype
    >>> from pandas.api.extensions import ExtensionDtype
    >>> @register_extension_dtype
    ... class MyExtensionDtype(ExtensionDtype):
    ...     pass
    """
    registry.register(cls)
    return cls


class Registry(object):
    """
    Registry for dtype inference

    The registry allows one to map a string repr of a extension
    dtype to an extension dtype. The string alias can be used in several
    places, including

    * Series and Index constructors
    * :meth:`pandas.array`
    * :meth:`pandas.Series.astype`

    Multiple extension types can be registered.
    These are tried in order.
    """
    def __init__(self):
        self.dtypes = []

    def register(self, dtype):
        """
        Parameters
        ----------
        dtype : ExtensionDtype
        """
        if not issubclass(dtype, (PandasExtensionDtype, ExtensionDtype)):
            raise ValueError("can only register pandas extension dtypes")

        self.dtypes.append(dtype)

    def find(self, dtype):
        """
        Parameters
        ----------
        dtype : PandasExtensionDtype or string

        Returns
        -------
        return the first matching dtype, otherwise return None
        """
        if not isinstance(dtype, compat.string_types):
            dtype_type = dtype
            if not isinstance(dtype, type):
                dtype_type = type(dtype)
            if issubclass(dtype_type, ExtensionDtype):
                return dtype

            return None

        for dtype_type in self.dtypes:
            try:
                return dtype_type.construct_from_string(dtype)
            except TypeError:
                pass

        return None


registry = Registry()


class PandasExtensionDtype(_DtypeOpsMixin):
    """
    A np.dtype duck-typed class, suitable for holding a custom dtype.

    THIS IS NOT A REAL NUMPY DTYPE
    """
    type = None
    subdtype = None
    kind = None
    str = None
    num = 100
    shape = tuple()
    itemsize = 8
    base = None
    isbuiltin = 0
    isnative = 0
    _cache = {}

    def __unicode__(self):
        return self.name

    def __str__(self):
        """
        Return a string representation for a particular Object

        Invoked by str(df) in both py2/py3.
        Yields Bytestring in Py2, Unicode String in py3.
        """

        if compat.PY3:
            return self.__unicode__()
        return self.__bytes__()

    def __bytes__(self):
        """
        Return a string representation for a particular object.

        Invoked by bytes(obj) in py3 only.
        Yields a bytestring in both py2/py3.
        """
        from pandas.core.config import get_option

        encoding = get_option("display.encoding")
        return self.__unicode__().encode(encoding, 'replace')

    def __repr__(self):
        """
        Return a string representation for a particular object.

        Yields Bytestring in Py2, Unicode String in py3.
        """
        return str(self)

    def __hash__(self):
        raise NotImplementedError("sub-classes should implement an __hash__ "
                                  "method")

    def __getstate__(self):
        # pickle support; we don't want to pickle the cache
        return {k: getattr(self, k, None) for k in self._metadata}

    @classmethod
    def reset_cache(cls):
        """ clear the cache """
        cls._cache = {}


class CategoricalDtypeType(type):
    """
    the type of CategoricalDtype, this metaclass determines subclass ability
    """
    pass


@register_extension_dtype
class CategoricalDtype(PandasExtensionDtype, ExtensionDtype):
    """
    Type for categorical data with the categories and orderedness

    .. versionchanged:: 0.21.0

    Parameters
    ----------
    categories : sequence, optional
        Must be unique, and must not contain any nulls.
    ordered : bool, default False

    Attributes
    ----------
    categories
    ordered

    Methods
    -------
    None

    See Also
    --------
    pandas.Categorical

    Notes
    -----
    This class is useful for specifying the type of a ``Categorical``
    independent of the values. See :ref:`categorical.categoricaldtype`
    for more.

    Examples
    --------
    >>> t = pd.CategoricalDtype(categories=['b', 'a'], ordered=True)
    >>> pd.Series(['a', 'b', 'a', 'c'], dtype=t)
    0      a
    1      b
    2      a
    3    NaN
    dtype: category
    Categories (2, object): [b < a]
    """
    # TODO: Document public vs. private API
    name = 'category'
    type = CategoricalDtypeType
    kind = 'O'
    str = '|O08'
    base = np.dtype('O')
    _metadata = ('categories', 'ordered')
    _cache = {}

    def __init__(self, categories=None, ordered=None):
        self._finalize(categories, ordered, fastpath=False)

    @classmethod
    def _from_fastpath(cls, categories=None, ordered=None):
        self = cls.__new__(cls)
        self._finalize(categories, ordered, fastpath=True)
        return self

    @classmethod
    def _from_categorical_dtype(cls, dtype, categories=None, ordered=None):
        if categories is ordered is None:
            return dtype
        if categories is None:
            categories = dtype.categories
        if ordered is None:
            ordered = dtype.ordered
        return cls(categories, ordered)

    @classmethod
    def _from_values_or_dtype(cls, values=None, categories=None, ordered=None,
                              dtype=None):
        """
        Construct dtype from the input parameters used in :class:`Categorical`.

        This constructor method specifically does not do the factorization
        step, if that is needed to find the categories. This constructor may
        therefore return ``CategoricalDtype(categories=None, ordered=None)``,
        which may not be useful. Additional steps may therefore have to be
        taken to create the final dtype.

        The return dtype is specified from the inputs in this prioritized
        order:
        1. if dtype is a CategoricalDtype, return dtype
        2. if dtype is the string 'category', create a CategoricalDtype from
           the supplied categories and ordered parameters, and return that.
        3. if values is a categorical, use value.dtype, but override it with
           categories and ordered if either/both of those are not None.
        4. if dtype is None and values is not a categorical, construct the
           dtype from categories and ordered, even if either of those is None.

        Parameters
        ----------
        values : list-like, optional
            The list-like must be 1-dimensional.
        categories : list-like, optional
            Categories for the CategoricalDtype.
        ordered : bool, optional
            Designating if the categories are ordered.
        dtype : CategoricalDtype or the string "category", optional
            If ``CategoricalDtype``, cannot be used together with
            `categories` or `ordered`.

        Returns
        -------
        CategoricalDtype

        Examples
        --------
        >>> CategoricalDtype._from_values_or_dtype()
        CategoricalDtype(categories=None, ordered=None)
        >>> CategoricalDtype._from_values_or_dtype(categories=['a', 'b'],
        ...                                        ordered=True)
        CategoricalDtype(categories=['a', 'b'], ordered=True)
        >>> dtype1 = CategoricalDtype(['a', 'b'], ordered=True)
        >>> dtype2 = CategoricalDtype(['x', 'y'], ordered=False)
        >>> c = Categorical([0, 1], dtype=dtype1, fastpath=True)
        >>> CategoricalDtype._from_values_or_dtype(c, ['x', 'y'], ordered=True,
        ...                                        dtype=dtype2)
        ValueError: Cannot specify `categories` or `ordered` together with
        `dtype`.

        The supplied dtype takes precedence over values' dtype:

        >>> CategoricalDtype._from_values_or_dtype(c, dtype=dtype2)
        CategoricalDtype(['x', 'y'], ordered=False)
        """
        from pandas.core.dtypes.common import is_categorical

        if dtype is not None:
            # The dtype argument takes precedence over values.dtype (if any)
            if isinstance(dtype, compat.string_types):
                if dtype == 'category':
                    dtype = CategoricalDtype(categories, ordered)
                else:
                    msg = "Unknown dtype {dtype!r}"
                    raise ValueError(msg.format(dtype=dtype))
            elif categories is not None or ordered is not None:
                raise ValueError("Cannot specify `categories` or `ordered` "
                                 "together with `dtype`.")
        elif is_categorical(values):
            # If no "dtype" was passed, use the one from "values", but honor
            # the "ordered" and "categories" arguments
            dtype = values.dtype._from_categorical_dtype(values.dtype,
                                                         categories, ordered)
        else:
            # If dtype=None and values is not categorical, create a new dtype.
            # Note: This could potentially have categories=None and
            # ordered=None.
            dtype = CategoricalDtype(categories, ordered)

        return dtype

    def _finalize(self, categories, ordered, fastpath=False):

        if ordered is not None:
            self.validate_ordered(ordered)

        if categories is not None:
            categories = self.validate_categories(categories,
                                                  fastpath=fastpath)

        self._categories = categories
        self._ordered = ordered

    def __setstate__(self, state):
        self._categories = state.pop('categories', None)
        self._ordered = state.pop('ordered', False)

    def __hash__(self):
        # _hash_categories returns a uint64, so use the negative
        # space for when we have unknown categories to avoid a conflict
        if self.categories is None:
            if self.ordered:
                return -1
            else:
                return -2
        # We *do* want to include the real self.ordered here
        return int(self._hash_categories(self.categories, self.ordered))

    def __eq__(self, other):
        """
        Rules for CDT equality:
        1) Any CDT is equal to the string 'category'
        2) Any CDT is equal to itself
        3) Any CDT is equal to a CDT with categories=None regardless of ordered
        4) A CDT with ordered=True is only equal to another CDT with
           ordered=True and identical categories in the same order
        5) A CDT with ordered={False, None} is only equal to another CDT with
           ordered={False, None} and identical categories, but same order is
           not required. There is no distinction between False/None.
        6) Any other comparison returns False
        """
        if isinstance(other, compat.string_types):
            return other == self.name
        elif other is self:
            return True
        elif not (hasattr(other, 'ordered') and hasattr(other, 'categories')):
            return False
        elif self.categories is None or other.categories is None:
            # We're forced into a suboptimal corner thanks to math and
            # backwards compatibility. We require that `CDT(...) == 'category'`
            # for all CDTs **including** `CDT(None, ...)`. Therefore, *all*
            # CDT(., .) = CDT(None, False) and *all*
            # CDT(., .) = CDT(None, True).
            return True
        elif self.ordered or other.ordered:
            # At least one has ordered=True; equal if both have ordered=True
            # and the same values for categories in the same order.
            return ((self.ordered == other.ordered) and
                    self.categories.equals(other.categories))
        else:
            # Neither has ordered=True; equal if both have the same categories,
            # but same order is not necessary.  There is no distinction between
            # ordered=False and ordered=None: CDT(., False) and CDT(., None)
            # will be equal if they have the same categories.
            return hash(self) == hash(other)

    def __repr__(self):
        tpl = u'CategoricalDtype(categories={}ordered={})'
        if self.categories is None:
            data = u"None, "
        else:
            data = self.categories._format_data(name=self.__class__.__name__)
        return tpl.format(data, self.ordered)

    @staticmethod
    def _hash_categories(categories, ordered=True):
        from pandas.core.util.hashing import (
            hash_array, _combine_hash_arrays, hash_tuples
        )
        from pandas.core.dtypes.common import is_datetime64tz_dtype, _NS_DTYPE

        if len(categories) and isinstance(categories[0], tuple):
            # assumes if any individual category is a tuple, then all our. ATM
            # I don't really want to support just some of the categories being
            # tuples.
            categories = list(categories)  # breaks if a np.array of categories
            cat_array = hash_tuples(categories)
        else:
            if categories.dtype == 'O':
                types = [type(x) for x in categories]
                if not len(set(types)) == 1:
                    # TODO: hash_array doesn't handle mixed types. It casts
                    # everything to a str first, which means we treat
                    # {'1', '2'} the same as {'1', 2}
                    # find a better solution
                    hashed = hash((tuple(categories), ordered))
                    return hashed

            if is_datetime64tz_dtype(categories.dtype):
                # Avoid future warning.
                categories = categories.astype(_NS_DTYPE)

            cat_array = hash_array(np.asarray(categories), categorize=False)
        if ordered:
            cat_array = np.vstack([
                cat_array, np.arange(len(cat_array), dtype=cat_array.dtype)
            ])
        else:
            cat_array = [cat_array]
        hashed = _combine_hash_arrays(iter(cat_array),
                                      num_items=len(cat_array))
        return np.bitwise_xor.reduce(hashed)

    @classmethod
    def construct_array_type(cls):
        """
        Return the array type associated with this dtype

        Returns
        -------
        type
        """
        from pandas import Categorical
        return Categorical

    @classmethod
    def construct_from_string(cls, string):
        """
        attempt to construct this type from a string, raise a TypeError if
        it's not possible """
        try:
            if string == 'category':
                return cls()
            else:
                raise TypeError("cannot construct a CategoricalDtype")
        except AttributeError:
            pass

    @staticmethod
    def validate_ordered(ordered):
        """
        Validates that we have a valid ordered parameter. If
        it is not a boolean, a TypeError will be raised.

        Parameters
        ----------
        ordered : object
            The parameter to be verified.

        Raises
        ------
        TypeError
            If 'ordered' is not a boolean.
        """
        from pandas.core.dtypes.common import is_bool
        if not is_bool(ordered):
            raise TypeError("'ordered' must either be 'True' or 'False'")

    @staticmethod
    def validate_categories(categories, fastpath=False):
        """
        Validates that we have good categories

        Parameters
        ----------
        categories : array-like
        fastpath : bool
            Whether to skip nan and uniqueness checks

        Returns
        -------
        categories : Index
        """
        from pandas import Index

        if not fastpath and not is_list_like(categories):
            msg = "Parameter 'categories' must be list-like, was {!r}"
            raise TypeError(msg.format(categories))
        elif not isinstance(categories, ABCIndexClass):
            categories = Index(categories, tupleize_cols=False)

        if not fastpath:

            if categories.hasnans:
                raise ValueError('Categorial categories cannot be null')

            if not categories.is_unique:
                raise ValueError('Categorical categories must be unique')

        if isinstance(categories, ABCCategoricalIndex):
            categories = categories.categories

        return categories

    def update_dtype(self, dtype):
        """
        Returns a CategoricalDtype with categories and ordered taken from dtype
        if specified, otherwise falling back to self if unspecified

        Parameters
        ----------
        dtype : CategoricalDtype

        Returns
        -------
        new_dtype : CategoricalDtype
        """
        if isinstance(dtype, compat.string_types) and dtype == 'category':
            # dtype='category' should not change anything
            return self
        elif not self.is_dtype(dtype):
            msg = ('a CategoricalDtype must be passed to perform an update, '
                   'got {dtype!r}').format(dtype=dtype)
            raise ValueError(msg)
        elif dtype.categories is not None and dtype.ordered is self.ordered:
            return dtype

        # dtype is CDT: keep current categories/ordered if None
        new_categories = dtype.categories
        if new_categories is None:
            new_categories = self.categories

        new_ordered = dtype.ordered
        if new_ordered is None:
            new_ordered = self.ordered

        return CategoricalDtype(new_categories, new_ordered)

    @property
    def categories(self):
        """
        An ``Index`` containing the unique categories allowed.
        """
        return self._categories

    @property
    def ordered(self):
        """
        Whether the categories have an ordered relationship.
        """
        return self._ordered

    @property
    def _is_boolean(self):
        from pandas.core.dtypes.common import is_bool_dtype

        return is_bool_dtype(self.categories)


@register_extension_dtype
class DatetimeTZDtype(PandasExtensionDtype, ExtensionDtype):

    """
    A np.dtype duck-typed class, suitable for holding a custom datetime with tz
    dtype.

    THIS IS NOT A REAL NUMPY DTYPE, but essentially a sub-class of
    np.datetime64[ns]
    """
    type = Timestamp
    kind = 'M'
    str = '|M8[ns]'
    num = 101
    base = np.dtype('M8[ns]')
    na_value = NaT
    _metadata = ('unit', 'tz')
    _match = re.compile(r"(datetime64|M8)\[(?P<unit>.+), (?P<tz>.+)\]")
    _cache = {}

    def __init__(self, unit="ns", tz=None):
        """
        An ExtensionDtype for timezone-aware datetime data.

        Parameters
        ----------
        unit : str, default "ns"
            The precision of the datetime data. Currently limited
            to ``"ns"``.
        tz : str, int, or datetime.tzinfo
            The timezone.

        Raises
        ------
        pytz.UnknownTimeZoneError
            When the requested timezone cannot be found.

        Examples
        --------
        >>> pd.core.dtypes.dtypes.DatetimeTZDtype(tz='UTC')
        datetime64[ns, UTC]

        >>> pd.core.dtypes.dtypes.DatetimeTZDtype(tz='dateutil/US/Central')
        datetime64[ns, tzfile('/usr/share/zoneinfo/US/Central')]
        """
        if isinstance(unit, DatetimeTZDtype):
            unit, tz = unit.unit, unit.tz

        if unit != 'ns':
            if isinstance(unit, compat.string_types) and tz is None:
                # maybe a string like datetime64[ns, tz], which we support for
                # now.
                result = type(self).construct_from_string(unit)
                unit = result.unit
                tz = result.tz
                msg = (
                    "Passing a dtype alias like 'datetime64[ns, {tz}]' "
                    "to DatetimeTZDtype is deprecated. Use "
                    "'DatetimeTZDtype.construct_from_string()' instead."
                )
                warnings.warn(msg.format(tz=tz), FutureWarning, stacklevel=2)
            else:
                raise ValueError("DatetimeTZDtype only supports ns units")

        if tz:
            tz = timezones.maybe_get_tz(tz)
        elif tz is not None:
            raise pytz.UnknownTimeZoneError(tz)
        elif tz is None:
            raise TypeError("A 'tz' is required.")

        self._unit = unit
        self._tz = tz

    @property
    def unit(self):
        """The precision of the datetime data."""
        return self._unit

    @property
    def tz(self):
        """The timezone."""
        return self._tz

    @classmethod
    def construct_array_type(cls):
        """
        Return the array type associated with this dtype

        Returns
        -------
        type
        """
        from pandas.core.arrays import DatetimeArray
        return DatetimeArray

    @classmethod
    def construct_from_string(cls, string):
        """
        Construct a DatetimeTZDtype from a string.

        Parameters
        ----------
        string : str
            The string alias for this DatetimeTZDtype.
            Should be formatted like ``datetime64[ns, <tz>]``,
            where ``<tz>`` is the timezone name.

        Examples
        --------
        >>> DatetimeTZDtype.construct_from_string('datetime64[ns, UTC]')
        datetime64[ns, UTC]
        """
        if isinstance(string, compat.string_types):
            msg = "Could not construct DatetimeTZDtype from '{}'"
            try:
                match = cls._match.match(string)
                if match:
                    d = match.groupdict()
                    return cls(unit=d['unit'], tz=d['tz'])
            except Exception:
                # TODO(py3): Change this pass to `raise TypeError(msg) from e`
                pass
            raise TypeError(msg.format(string))

        raise TypeError("Could not construct DatetimeTZDtype")

    def __unicode__(self):
        return "datetime64[{unit}, {tz}]".format(unit=self.unit, tz=self.tz)

    @property
    def name(self):
        """A string representation of the dtype."""
        return str(self)

    def __hash__(self):
        # make myself hashable
        # TODO: update this.
        return hash(str(self))

    def __eq__(self, other):
        if isinstance(other, compat.string_types):
            return other == self.name

        return (isinstance(other, DatetimeTZDtype) and
                self.unit == other.unit and
                str(self.tz) == str(other.tz))

    def __setstate__(self, state):
        # for pickle compat.
        self._tz = state['tz']
        self._unit = state['unit']


@register_extension_dtype
class PeriodDtype(ExtensionDtype, PandasExtensionDtype):
    """
    A Period duck-typed class, suitable for holding a period with freq dtype.

    THIS IS NOT A REAL NUMPY DTYPE, but essentially a sub-class of np.int64.
    """
    type = Period
    kind = 'O'
    str = '|O08'
    base = np.dtype('O')
    num = 102
    _metadata = ('freq',)
    _match = re.compile(r"(P|p)eriod\[(?P<freq>.+)\]")
    _cache = {}

    def __new__(cls, freq=None):
        """
        Parameters
        ----------
        freq : frequency
        """

        if isinstance(freq, PeriodDtype):
            return freq

        elif freq is None:
            # empty constructor for pickle compat
            return object.__new__(cls)

        from pandas.tseries.offsets import DateOffset
        if not isinstance(freq, DateOffset):
            freq = cls._parse_dtype_strict(freq)

        try:
            return cls._cache[freq.freqstr]
        except KeyError:
            u = object.__new__(cls)
            u.freq = freq
            cls._cache[freq.freqstr] = u
            return u

    @classmethod
    def _parse_dtype_strict(cls, freq):
        if isinstance(freq, compat.string_types):
            if freq.startswith('period[') or freq.startswith('Period['):
                m = cls._match.search(freq)
                if m is not None:
                    freq = m.group('freq')
            from pandas.tseries.frequencies import to_offset
            freq = to_offset(freq)
            if freq is not None:
                return freq

        raise ValueError("could not construct PeriodDtype")

    @classmethod
    def construct_from_string(cls, string):
        """
        Strict construction from a string, raise a TypeError if not
        possible
        """
        from pandas.tseries.offsets import DateOffset

        if (isinstance(string, compat.string_types) and
            (string.startswith('period[') or
             string.startswith('Period[')) or
                isinstance(string, DateOffset)):
            # do not parse string like U as period[U]
            # avoid tuple to be regarded as freq
            try:
                return cls(freq=string)
            except ValueError:
                pass
        raise TypeError("could not construct PeriodDtype")

    def __unicode__(self):
        return compat.text_type(self.name)

    @property
    def name(self):
        return str("period[{freq}]".format(freq=self.freq.freqstr))

    @property
    def na_value(self):
        return NaT

    def __hash__(self):
        # make myself hashable
        return hash(str(self))

    def __eq__(self, other):
        if isinstance(other, compat.string_types):
            return other == self.name or other == self.name.title()

        return isinstance(other, PeriodDtype) and self.freq == other.freq

    @classmethod
    def is_dtype(cls, dtype):
        """
        Return a boolean if we if the passed type is an actual dtype that we
        can match (via string or type)
        """

        if isinstance(dtype, compat.string_types):
            # PeriodDtype can be instantiated from freq string like "U",
            # but doesn't regard freq str like "U" as dtype.
            if dtype.startswith('period[') or dtype.startswith('Period['):
                try:
                    if cls._parse_dtype_strict(dtype) is not None:
                        return True
                    else:
                        return False
                except ValueError:
                    return False
            else:
                return False
        return super(PeriodDtype, cls).is_dtype(dtype)

    @classmethod
    def construct_array_type(cls):
        from pandas.core.arrays import PeriodArray

        return PeriodArray


@register_extension_dtype
class IntervalDtype(PandasExtensionDtype, ExtensionDtype):
    """
    A Interval duck-typed class, suitable for holding an interval

    THIS IS NOT A REAL NUMPY DTYPE
    """
    name = 'interval'
    kind = None
    str = '|O08'
    base = np.dtype('O')
    num = 103
    _metadata = ('subtype',)
    _match = re.compile(r"(I|i)nterval\[(?P<subtype>.+)\]")
    _cache = {}

    def __new__(cls, subtype=None):
        """
        Parameters
        ----------
        subtype : the dtype of the Interval
        """
        from pandas.core.dtypes.common import (
            is_categorical_dtype, is_string_dtype, pandas_dtype)

        if isinstance(subtype, IntervalDtype):
            return subtype
        elif subtype is None:
            # we are called as an empty constructor
            # generally for pickle compat
            u = object.__new__(cls)
            u.subtype = None
            return u
        elif (isinstance(subtype, compat.string_types) and
              subtype.lower() == 'interval'):
            subtype = None
        else:
            if isinstance(subtype, compat.string_types):
                m = cls._match.search(subtype)
                if m is not None:
                    subtype = m.group('subtype')

            try:
                subtype = pandas_dtype(subtype)
            except TypeError:
                raise TypeError("could not construct IntervalDtype")

        if is_categorical_dtype(subtype) or is_string_dtype(subtype):
            # GH 19016
            msg = ('category, object, and string subtypes are not supported '
                   'for IntervalDtype')
            raise TypeError(msg)

        try:
            return cls._cache[str(subtype)]
        except KeyError:
            u = object.__new__(cls)
            u.subtype = subtype
            cls._cache[str(subtype)] = u
            return u

    @classmethod
    def construct_array_type(cls):
        """
        Return the array type associated with this dtype

        Returns
        -------
        type
        """
        from pandas.core.arrays import IntervalArray
        return IntervalArray

    @classmethod
    def construct_from_string(cls, string):
        """
        attempt to construct this type from a string, raise a TypeError
        if its not possible
        """
        if (isinstance(string, compat.string_types) and
            (string.startswith('interval') or
             string.startswith('Interval'))):
            return cls(string)

        msg = "a string needs to be passed, got type {typ}"
        raise TypeError(msg.format(typ=type(string)))

    @property
    def type(self):
        return Interval

    def __unicode__(self):
        if self.subtype is None:
            return "interval"
        return "interval[{subtype}]".format(subtype=self.subtype)

    def __hash__(self):
        # make myself hashable
        return hash(str(self))

    def __eq__(self, other):
        if isinstance(other, compat.string_types):
            return other.lower() in (self.name.lower(), str(self).lower())
        elif not isinstance(other, IntervalDtype):
            return False
        elif self.subtype is None or other.subtype is None:
            # None should match any subtype
            return True
        else:
            from pandas.core.dtypes.common import is_dtype_equal
            return is_dtype_equal(self.subtype, other.subtype)

    @classmethod
    def is_dtype(cls, dtype):
        """
        Return a boolean if we if the passed type is an actual dtype that we
        can match (via string or type)
        """

        if isinstance(dtype, compat.string_types):
            if dtype.lower().startswith('interval'):
                try:
                    if cls.construct_from_string(dtype) is not None:
                        return True
                    else:
                        return False
                except ValueError:
                    return False
            else:
                return False
        return super(IntervalDtype, cls).is_dtype(dtype)
