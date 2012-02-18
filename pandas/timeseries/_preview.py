"""
This file contains local copies of functions to be released in numpy 1.4.0.
Once we drop support for numpy 1.3.x we can eliminate this file.
"""

__all__ = ['genfromtxt']

import itertools
import warnings
from operator import itemgetter
from __builtin__ import bool, int, long, float, complex, object, unicode, str

import numpy as np
import numpy.core.numeric as nx

#####--------------------------------------------------------------------------
#---- numpy.lib._iotools functions ---
#####--------------------------------------------------------------------------

def _is_string_like(obj):
    """
    Check whether obj behaves like a string.
    """
    try:
        obj + ''
    except (TypeError, ValueError):
        return False
    return True


def _to_filehandle(fname, flag='r', return_opened=False):
    """
    Returns the filehandle corresponding to a string or a file.
    If the string ends in '.gz', the file is automatically unzipped.

    Parameters
    ----------
    fname : string, filehandle
        Name of the file whose filehandle must be returned.
    flag : string, optional
        Flag indicating the status of the file ('r' for read, 'w' for write).
    return_opened : boolean, optional
        Whether to return the opening status of the file.
    """
    if _is_string_like(fname):
        if fname.endswith('.gz'):
            import gzip
            fhd = gzip.open(fname, flag)
        elif fname.endswith('.bz2'):
            import bz2
            fhd = bz2.BZ2File(fname)
        else:
            fhd = file(fname, flag)
        opened = True
    elif hasattr(fname, 'seek'):
        fhd = fname
        opened = False
    else:
        raise ValueError('fname must be a string or file handle')
    if return_opened:
        return fhd, opened
    return fhd


def has_nested_fields(ndtype):
    """
    Returns whether one or several fields of a dtype are nested.

    Parameters
    ----------
    ndtype : dtype
        Data-type of a structured array.

    Raises
    ------
    AttributeError : If `ndtype` does not have a `names` attribute.

    Examples
    --------
    >>> dt = np.dtype([('name', 'S4'), ('x', float), ('y', float)])
    >>> np.lib._iotools.has_nested_fields(dt)
    False

    """
    for name in ndtype.names or ():
        if ndtype[name].names:
            return True
    return False


def flatten_dtype(ndtype, flatten_base=False):
    """
    Unpack a structured data-type by collapsing nested fields and/or fields
    with a shape.

    Note that the field names are lost.

    Parameters
    ----------
    ndtype : dtype
        The datatype to collapse
    flatten_base : {False, True}, optional
        Whether to transform a field with a shape into several fields or not.

    Examples
    --------
    >>> dt = np.dtype([('name', 'S4'), ('x', float), ('y', float),
    ...                ('block', int, (2, 3))])
    >>> np.lib._iotools.flatten_dtype(dt)
    [dtype('|S4'), dtype('float64'), dtype('float64'), dtype('int32')]
    >>> np.lib._iotools.flatten_dtype(dt, flatten_base=True)
    [dtype('|S4'), dtype('float64'), dtype('float64'), dtype('int32'),
     dtype('int32'), dtype('int32'), dtype('int32'), dtype('int32'),
     dtype('int32')]

    """
    names = ndtype.names
    if names is None:
        if flatten_base:
            return [ndtype.base] * int(np.prod(ndtype.shape))
        return [ndtype.base]
    else:
        types = []
        for field in names:
            (typ, _) = ndtype.fields[field]
            flat_dt = flatten_dtype(typ, flatten_base)
            types.extend(flat_dt)
        return types






class LineSplitter:
    """
    Object to split a string at a given delimiter or at given places.

    Parameters
    ----------
    delimiter : str, int, or sequence of ints, optional
        If a string, character used to delimit consecutive fields.
        If an integer or a sequence of integers, width(s) of each field.
    comment : str, optional
        Character used to mark the beginning of a comment. Default is '#'.
    autostrip : bool, optional
        Whether to strip each individual field. Default is True.

    """

    def autostrip(self, method):
        """
        Wrapper to strip each member of the output of `method`.

        Parameters
        ----------
        method : function
            Function that takes a single argument and returns a sequence of
            strings.

        Returns
        -------
        wrapped : function
            The result of wrapping `method`. `wrapped` takes a single input
            argument and returns a list of strings that are stripped of
            white-space.

        """
        return lambda input: [_.strip() for _ in method(input)]
    #
    def __init__(self, delimiter=None, comments='#', autostrip=True):
        self.comments = comments
        # Delimiter is a character
        if (delimiter is None) or _is_string_like(delimiter):
            delimiter = delimiter or None
            _handyman = self._delimited_splitter
        # Delimiter is a list of field widths
        elif hasattr(delimiter, '__iter__'):
            _handyman = self._variablewidth_splitter
            idx = np.cumsum([0] + list(delimiter))
            delimiter = [slice(i, j) for (i, j) in zip(idx[:-1], idx[1:])]
        # Delimiter is a single integer
        elif int(delimiter):
            (_handyman, delimiter) = (self._fixedwidth_splitter, int(delimiter))
        else:
            (_handyman, delimiter) = (self._delimited_splitter, None)
        self.delimiter = delimiter
        if autostrip:
            self._handyman = self.autostrip(_handyman)
        else:
            self._handyman = _handyman
    #
    def _delimited_splitter(self, line):
        line = line.split(self.comments)[0].strip()
        if not line:
            return []
        return line.split(self.delimiter)
    #
    def _fixedwidth_splitter(self, line):
        line = line.split(self.comments)[0]
        if not line:
            return []
        fixed = self.delimiter
        slices = [slice(i, i + fixed) for i in range(len(line))[::fixed]]
        return [line[s] for s in slices]
    #
    def _variablewidth_splitter(self, line):
        line = line.split(self.comments)[0]
        if not line:
            return []
        slices = self.delimiter
        return [line[s] for s in slices]
    #
    def __call__(self, line):
        return self._handyman(line)



class NameValidator:
    """
    Object to validate a list of strings to use as field names.

    The strings are stripped of any non alphanumeric character, and spaces
    are replaced by '_'. During instantiation, the user can define a list of
    names to exclude, as well as a list of invalid characters. Names in the
    exclusion list are appended a '_' character.

    Once an instance has been created, it can be called with a list of names,
    and a list of valid names will be created.
    The `__call__` method accepts an optional keyword "default" that sets
    the default name in case of ambiguity. By default this is 'f', so
    that names will default to `f0`, `f1`, etc.

    Parameters
    ----------
    excludelist : sequence, optional
        A list of names to exclude. This list is appended to the default list
        ['return', 'file', 'print']. Excluded names are appended an underscore:
        for example, `file` becomes `file_` if supplied.
    deletechars : str, optional
        A string combining invalid characters that must be deleted from the
        names.
    casesensitive : {True, False, 'upper', 'lower'}, optional
        * If True, field names are case-sensitive.
        * If False or 'upper', field names are converted to upper case.
        * If 'lower', field names are converted to lower case.

        The default value is True.

    Notes
    -----
    Calling an instance of `NameValidator` is the same as calling its method
    `validate`.

    Examples
    --------
    >>> validator = np.lib._iotools.NameValidator()
    >>> validator(['file', 'field2', 'with space', 'CaSe'])
    ['file_', 'field2', 'with_space', 'CaSe']

    >>> validator = np.lib._iotools.NameValidator(excludelist=['excl'],
                                                  deletechars='q',
                                                  case_sensitive='False')
    >>> validator(['excl', 'field2', 'no_q', 'with space', 'CaSe'])
    ['excl_', 'field2', 'no_', 'with_space', 'case']

    """
    #
    defaultexcludelist = ['return', 'file', 'print']
    defaultdeletechars = set("""~!@#$%^&*()-=+~\|]}[{';: /?.>,<""")
    #
    def __init__(self, excludelist=None, deletechars=None, case_sensitive=None):
        # Process the exclusion list ..
        if excludelist is None:
            excludelist = []
        excludelist.extend(self.defaultexcludelist)
        self.excludelist = excludelist
        # Process the list of characters to delete
        if deletechars is None:
            delete = self.defaultdeletechars
        else:
            delete = set(deletechars)
        delete.add('"')
        self.deletechars = delete
        # Process the case option .....
        if (case_sensitive is None) or (case_sensitive is True):
            self.case_converter = lambda x: x
        elif (case_sensitive is False) or ('u' in case_sensitive):
            self.case_converter = lambda x: x.upper()
        elif 'l' in case_sensitive:
            self.case_converter = lambda x: x.lower()
        else:
            self.case_converter = lambda x: x

    def validate(self, names, defaultfmt="f%i", nbfields=None):
        """
        Validate a list of strings to use as field names for a structured array.

        Parameters
        ----------
        names : sequence of str
            Strings to be validated.
        defaultfmt : str, optional
            Default format string, used if validating a given string reduces its
            length to zero.
        nboutput : integer, optional
            Final number of validated names, used to expand or shrink the initial
            list of names.

        Returns
        -------
        validatednames : list of str
            The list of validated field names.

        Notes
        -----
        A `NameValidator` instance can be called directly, which is the same as
        calling `validate`. For examples, see `NameValidator`.

        """
        # Initial checks ..............
        if (names is None):
            if (nbfields is None):
                return None
            names = []
        if isinstance(names, basestring):
            names = [names, ]
        if nbfields is not None:
            nbnames = len(names)
            if (nbnames < nbfields):
                names = list(names) + [''] * (nbfields - nbnames)
            elif (nbnames > nbfields):
                names = names[:nbfields]
        # Set some shortcuts ...........
        deletechars = self.deletechars
        excludelist = self.excludelist
        case_converter = self.case_converter
        # Initializes some variables ...
        validatednames = []
        seen = dict()
        nbempty = 0
        #
        for item in names:
            item = case_converter(item)
            item = item.strip().replace(' ', '_')
            item = ''.join([c for c in item if c not in deletechars])
            if item == '':
                item = defaultfmt % nbempty
                while item in names:
                    nbempty += 1
                    item = defaultfmt % nbempty
                nbempty += 1
            elif item in excludelist:
                item += '_'
            cnt = seen.get(item, 0)
            if cnt > 0:
                validatednames.append(item + '_%d' % cnt)
            else:
                validatednames.append(item)
            seen[item] = cnt + 1
        return tuple(validatednames)
    #
    def __call__(self, names, defaultfmt="f%i", nbfields=None):
        return self.validate(names, defaultfmt=defaultfmt, nbfields=nbfields)



def str2bool(value):
    """
    Tries to transform a string supposed to represent a boolean to a boolean.

    Parameters
    ----------
    value : str
        The string that is transformed to a boolean.

    Returns
    -------
    boolval : bool
        The boolean representation of `value`.

    Raises
    ------
    ValueError
        If the string is not 'True' or 'False' (case independent)

    Examples
    --------
    >>> np.lib._iotools.str2bool('TRUE')
    True
    >>> np.lib._iotools.str2bool('false')
    False

    """
    value = value.upper()
    if value == 'TRUE':
        return True
    elif value == 'FALSE':
        return False
    else:
        raise ValueError("Invalid boolean")


class ConverterError(Exception):
    """
    Exception raised when an error occurs in a converter for string values.

    """
    pass

class ConverterLockError(ConverterError):
    """
    Exception raised when an attempt is made to upgrade a locked converter.

    """
    pass

class ConversionWarning(UserWarning):
    """
    Warning issued when a string converter has a problem.

    Notes
    -----
    In `genfromtxt` a `ConversionWarning` is issued if raising exceptions
    is explicitly suppressed with the "invalid_raise" keyword.

    """
    pass



class StringConverter:
    """
    Factory class for function transforming a string into another object (int,
    float).

    After initialization, an instance can be called to transform a string
    into another object. If the string is recognized as representing a missing
    value, a default value is returned.

    Attributes
    ----------
    func : function
        Function used for the conversion.
    default : any
        Default value to return when the input corresponds to a missing value.
    type : type
        Type of the output.
    _status : int
        Integer representing the order of the conversion.
    _mapper : sequence of tuples
        Sequence of tuples (dtype, function, default value) to evaluate in
        order.
    _locked : bool
        Holds `locked` parameter.

    Parameters
    ----------
    dtype_or_func : {None, dtype, function}, optional
        If a `dtype`, specifies the input data type, used to define a basic
        function and a default value for missing data. For example, when
        `dtype` is float, the `func` attribute is set to `float` and the
        default value to `np.nan`.
        If a function, this function is used to convert a string to another
        object. In this case, it is recommended to give an associated default
        value as input.
    default : any, optional
        Value to return by default, that is, when the string to be converted
        is flagged as missing. If not given, `StringConverter` tries to supply
        a reasonable default value.
    missing_values : sequence of str, optional
        Sequence of strings indicating a missing value.
    locked : bool, optional
        Whether the StringConverter should be locked to prevent automatic
        upgrade or not. Default is False.

    """
    #
    _mapper = [(nx.bool_, str2bool, False),
               (nx.integer, int, -1),
               (nx.floating, float, nx.nan),
               (complex, complex, nx.nan + 0j),
               (nx.string_, str, '???')]
    (_defaulttype, _defaultfunc, _defaultfill) = zip(*_mapper)
    #
    @classmethod
    def _getsubdtype(cls, val):
        """Returns the type of the dtype of the input variable."""
        return np.array(val).dtype.type
    #
    @classmethod
    def upgrade_mapper(cls, func, default=None):
        """
    Upgrade the mapper of a StringConverter by adding a new function and its
    corresponding default.

    The input function (or sequence of functions) and its associated default
    value (if any) is inserted in penultimate position of the mapper.
    The corresponding type is estimated from the dtype of the default value.

    Parameters
    ----------
    func : var
        Function, or sequence of functions

    Examples
    --------
    >>> import dateutil.parser
    >>> import datetime
    >>> dateparser = datetustil.parser.parse
    >>> defaultdate = datetime.date(2000, 1, 1)
    >>> StringConverter.upgrade_mapper(dateparser, default=defaultdate)
        """
        # Func is a single functions
        if hasattr(func, '__call__'):
            cls._mapper.insert(-1, (cls._getsubdtype(default), func, default))
            return
        elif hasattr(func, '__iter__'):
            if isinstance(func[0], (tuple, list)):
                for _ in func:
                    cls._mapper.insert(-1, _)
                return
            if default is None:
                default = [None] * len(func)
            else:
                default = list(default)
                default.append([None] * (len(func) - len(default)))
            for (fct, dft) in zip(func, default):
                cls._mapper.insert(-1, (cls._getsubdtype(dft), fct, dft))
    #
    def __init__(self, dtype_or_func=None, default=None, missing_values=None,
                 locked=False):
        # Defines a lock for upgrade
        self._locked = bool(locked)
        # No input dtype: minimal initialization
        if dtype_or_func is None:
            self.func = str2bool
            self._status = 0
            self.default = default or False
            ttype = np.bool
        else:
            # Is the input a np.dtype ?
            try:
                self.func = None
                ttype = np.dtype(dtype_or_func).type
            except TypeError:
                # dtype_or_func must be a function, then
                if not hasattr(dtype_or_func, '__call__'):
                    errmsg = "The input argument `dtype` is neither a function"\
                             " or a dtype (got '%s' instead)"
                    raise TypeError(errmsg % type(dtype_or_func))
                # Set the function
                self.func = dtype_or_func
                # If we don't have a default, try to guess it or set it to None
                if default is None:
                    try:
                        default = self.func('0')
                    except ValueError:
                        default = None
                ttype = self._getsubdtype(default)
            # Set the status according to the dtype
            _status = -1
            for (i, (deftype, func, default_def)) in enumerate(self._mapper):
                if np.issubdtype(ttype, deftype):
                    _status = i
                    if default is None:
                        self.default = default_def
                    else:
                        self.default = default
                    break
            if _status == -1:
                # We never found a match in the _mapper...
                _status = 0
                self.default = default
            self._status = _status
            # If the input was a dtype, set the function to the last we saw
            if self.func is None:
                self.func = func
            # If the status is 1 (int), change the function to smthg more robust
            if self.func == self._mapper[1][1]:
                self.func = lambda x : int(float(x))
        # Store the list of strings corresponding to missing values.
        if missing_values is None:
            self.missing_values = set([''])
        else:
            if isinstance(missing_values, basestring):
                missing_values = missing_values.split(",")
            self.missing_values = set(list(missing_values) + [''])
        #
        self._callingfunction = self._strict_call
        self.type = ttype
        self._checked = False
        self._initial_default = default
    #
    def _loose_call(self, value):
        try:
            return self.func(value)
        except ValueError:
            return self.default
    #
    def _strict_call(self, value):
        try:
            return self.func(value)
        except ValueError:
            if value.strip() in self.missing_values:
                if not self._status:
                    self._checked = False
                return self.default
            raise ValueError("Cannot convert string '%s'" % value)
    #
    def __call__(self, value):
        return self._callingfunction(value)
    #
    def upgrade(self, value):
        """
        Try to find the best converter for a given string, and return the result.

        The supplied string `value` is converted by testing different
        converters in order. First the `func` method of the `StringConverter`
        instance is tried, if this fails other available converters are tried.
        The order in which these other converters are tried is determined by the
        `_status` attribute of the instance.

        Parameters
        ----------
        value : str
            The string to convert.

        Returns
        -------
        out : any
            The result of converting `value` with the appropriate converter.

        """
        self._checked = True
        try:
            self._strict_call(value)
        except ValueError:
            # Raise an exception if we locked the converter...
            if self._locked:
                errmsg = "Converter is locked and cannot be upgraded"
                raise ConverterLockError(errmsg)
            _statusmax = len(self._mapper)
            # Complains if we try to upgrade by the maximum
            _status = self._status
            if _status == _statusmax:
                errmsg = "Could not find a valid conversion function"
                raise ConverterError(errmsg)
            elif _status < _statusmax - 1:
                _status += 1
            (self.type, self.func, default) = self._mapper[_status]
            self._status = _status
            if self._initial_default is not None:
                self.default = self._initial_default
            else:
                self.default = default
            self.upgrade(value)

    def iterupgrade(self, value):
        self._checked = True
        if not hasattr(value, '__iter__'):
            value = (value,)
        _strict_call = self._strict_call
        try:
            map(_strict_call, value)
        except ValueError:
            # Raise an exception if we locked the converter...
            if self._locked:
                errmsg = "Converter is locked and cannot be upgraded"
                raise ConverterLockError(errmsg)
            _statusmax = len(self._mapper)
            # Complains if we try to upgrade by the maximum
            _status = self._status
            if _status == _statusmax:
                raise ConverterError("Could not find a valid conversion function")
            elif _status < _statusmax - 1:
                _status += 1
            (self.type, self.func, default) = self._mapper[_status]
            if self._initial_default is not None:
                self.default = self._initial_default
            else:
                self.default = default
            self._status = _status
            self.iterupgrade(value)

    def update(self, func, default=None, missing_values='', locked=False):
        """
        Set StringConverter attributes directly.

        Parameters
        ----------
        func : function
            Conversion function.
        default : any, optional
            Value to return by default, that is, when the string to be converted
            is flagged as missing. If not given, `StringConverter` tries to supply
            a reasonable default value.
        missing_values : sequence of str, optional
            Sequence of strings indicating a missing value.
        locked : bool, optional
            Whether the StringConverter should be locked to prevent automatic
            upgrade or not. Default is False.

        Notes
        -----
        `update` takes the same parameters as the constructor of `StringConverter`,
        except that `func` does not accept a `dtype` whereas `dtype_or_func` in
        the constructor does.

        """
        self.func = func
        self._locked = locked
        # Don't reset the default to None if we can avoid it
        if default is not None:
            self.default = default
            self.type = self._getsubdtype(default)
        else:
            try:
                tester = func('1')
            except (TypeError, ValueError):
                tester = None
            self.type = self._getsubdtype(tester)
        # Add the missing values to the existing set
        if missing_values is not None:
            if _is_string_like(missing_values):
                self.missing_values.add(missing_values)
            elif hasattr(missing_values, '__iter__'):
                for val in missing_values:
                    self.missing_values.add(val)
        else:
            self.missing_values = []



def easy_dtype(ndtype, names=None, defaultfmt="f%i", **validationargs):
    """
    Convenience function to create a `np.dtype` object.

    The function processes the input `dtype` and matches it with the given
    names.

    Parameters
    ----------
    ndtype : var
        Definition of the dtype. Can be any string or dictionary
        recognized by the `np.dtype` function, or a sequence of types.
    names : str or sequence, optional
        Sequence of strings to use as field names for a structured dtype.
        For convenience, `names` can be a string of a comma-separated list of
        names.
    defaultfmt : str, optional
        Format string used to define missing names, such as ``"f%i"``
        (default) or ``"fields_%02i"``.
    validationargs : optional
        A series of optional arguments used to initialize a `NameValidator`.

    Examples
    --------
    >>> np.lib._iotools.easy_dtype(float)
    dtype('float64')
    >>> np.lib._iotools.easy_dtype("i4, f8")
    dtype([('f0', '<i4'), ('f1', '<f8')])
    >>> np.lib._iotools.easy_dtype("i4, f8", defaultfmt="field_%03i")
    dtype([('field_000', '<i4'), ('field_001', '<f8')])

    >>> np.lib._iotools.easy_dtype((int, float, float), names="a,b,c")
    dtype([('a', '<i8'), ('b', '<f8'), ('c', '<f8')])
    >>> np.lib._iotools.easy_dtype(float, names="a,b,c")
    dtype([('a', '<f8'), ('b', '<f8'), ('c', '<f8')])

    """
    try:
        ndtype = np.dtype(ndtype)
    except TypeError:
        validate = NameValidator(**validationargs)
        nbfields = len(ndtype)
        if names is None:
            names = [''] * len(ndtype)
        elif isinstance(names, basestring):
            names = names.split(",")
        names = validate(names, nbfields=nbfields, defaultfmt=defaultfmt)
        ndtype = np.dtype(dict(formats=ndtype, names=names))
    else:
        nbtypes = len(ndtype)
        # Explicit names
        if names is not None:
            validate = NameValidator(**validationargs)
            if isinstance(names, basestring):
                names = names.split(",")
            # Simple dtype: repeat to match the nb of names
            if nbtypes == 0:
                formats = tuple([ndtype.type] * len(names))
                names = validate(names, defaultfmt=defaultfmt)
                ndtype = np.dtype(zip(names, formats))
            # Structured dtype: just validate the names as needed
            else:
                ndtype.names = validate(names, nbfields=nbtypes,
                                        defaultfmt=defaultfmt)
        # No implicit names
        elif (nbtypes > 0):
            validate = NameValidator(**validationargs)
            # Default initial names : should we change the format ?
            if (ndtype.names == tuple("f%i" % i for i in range(nbtypes))) and \
               (defaultfmt != "f%i"):
                ndtype.names = validate([''] * nbtypes, defaultfmt=defaultfmt)
            # Explicit initial names : just validate
            else:
                ndtype.names = validate(ndtype.names, defaultfmt=defaultfmt)
    return ndtype


#####--------------------------------------------------------------------------
#---- numpy.lib.io functions ---
#####--------------------------------------------------------------------------

_string_like = _is_string_like

def genfromtxt(fname, dtype=float, comments='#', delimiter=None,
               skiprows=0, skip_header=0, skip_footer=0, converters=None,
               missing='', missing_values=None, filling_values=None,
               usecols=None, names=None, excludelist=None, deletechars=None,
               autostrip=False, case_sensitive=True, defaultfmt="f%i",
               unpack=None, usemask=False, loose=True, invalid_raise=True):
    """
    Load data from a text file, with missing values handled as specified.

    Each line past the first `skiprows` lines is split at the `delimiter`
    character, and characters following the `comments` character are discarded.

    Parameters
    ----------
    fname : file or str
        File or filename to read.  If the filename extension is `.gz` or
        `.bz2`, the file is first decompressed.
    dtype : dtype, optional
        Data type of the resulting array.
        If None, the dtypes will be determined by the contents of each
        column, individually.
    comments : str, optional
        The character used to indicate the start of a comment.
        All the characters occurring on a line after a comment are discarded
    delimiter : str, int, or sequence, optional
        The string used to separate values.  By default, any consecutive
        whitespaces act as delimiter.  An integer or sequence of integers
        can also be provided as width(s) of each field.
    skip_header : int, optional
        The numbers of lines to skip at the beginning of the file.
    skip_footer : int, optional
        The numbers of lines to skip at the end of the file
    converters : variable or None, optional
        The set of functions that convert the data of a column to a value.
        The converters can also be used to provide a default value
        for missing data: ``converters = {3: lambda s: float(s or 0)}``.
    missing_values : variable or None, optional
        The set of strings corresponding to missing data.
    filling_values : variable or None, optional
        The set of values to be used as default when the data are missing.
    usecols : sequence or None, optional
        Which columns to read, with 0 being the first.  For example,
        ``usecols = (1, 4, 5)`` will extract the 2nd, 5th and 6th columns.
    names : {None, True, str, sequence}, optional
        If `names` is True, the field names are read from the first valid line
        after the first `skiprows` lines.
        If `names` is a sequence or a single-string of comma-separated names,
        the names will be used to define the field names in a structured dtype.
        If `names` is None, the names of the dtype fields will be used, if any.
    excludelist : sequence, optional
        A list of names to exclude. This list is appended to the default list
        ['return','file','print']. Excluded names are appended an underscore:
        for example, `file` would become `file_`.
    deletechars : str, optional
        A string combining invalid characters that must be deleted from the
        names.
    defaultfmt : str, optional
        A format used to define default field names, such as "f%i" or "f_%02i".
    autostrip : bool, optional
        Whether to automatically strip white spaces from the variables.
    case_sensitive : {True, False, 'upper', 'lower'}, optional
        If True, field names are case sensitive.
        If False or 'upper', field names are converted to upper case.
        If 'lower', field names are converted to lower case.
    unpack : bool, optional
        If True, the returned array is transposed, so that arguments may be
        unpacked using ``x, y, z = loadtxt(...)``
    usemask : bool, optional
        If True, return a masked array.
        If False, return a regular array.
    invalid_raise : bool, optional
        If True, an exception is raised if an inconsistency is detected in the
        number of columns.
        If False, a warning is emitted and the offending lines are skipped.

    Returns
    -------
    out : ndarray
        Data read from the text file. If `usemask` is True, this is a
        masked array.

    See Also
    --------
    numpy.loadtxt : equivalent function when no data is missing.

    Notes
    -----
    * When spaces are used as delimiters, or when no delimiter has been given
      as input, there should not be any missing data between two fields.
    * When the variables are named (either by a flexible dtype or with `names`,
      there must not be any header in the file (else a ValueError
      exception is raised).
    * Individual values are not stripped of spaces by default.
      When using a custom converter, make sure the function does remove spaces.

    Examples
    ---------
    >>> from StringIO import StringIO
    >>> import numpy as np

    Comma delimited file with mixed dtype

    >>> s = StringIO("1,1.3,abcde")
    >>> data = np.genfromtxt(s, dtype=[('myint','i8'),('myfloat','f8'),
        ('mystring','S5')], delimiter=",")
    >>> data
    array((1, 1.3, 'abcde'),
          dtype=[('myint', '<i8'), ('myfloat', '<f8'), ('mystring', '|S5')])

    Using dtype = None

    >>> s.seek(0) # needed for StringIO example only
    >>> data = np.genfromtxt(s, dtype=None,
        names = ['myint','myfloat','mystring'], delimiter=",")
    >>> data
    array((1, 1.3, 'abcde'),
          dtype=[('myint', '<i8'), ('myfloat', '<f8'), ('mystring', '|S5')])

    Specifying dtype and names

    >>> s.seek(0)
    >>> data = np.genfromtxt(s, dtype="i8,f8,S5",
        names=['myint','myfloat','mystring'], delimiter=",")
    >>> data
    array((1, 1.3, 'abcde'),
          dtype=[('myint', '<i8'), ('myfloat', '<f8'), ('mystring', '|S5')])

    An example with fixed-width columns

    >>> s = StringIO("11.3abcde")
    >>> data = np.genfromtxt(s, dtype=None, names=['intvar','fltvar','strvar'],
            delimiter=[1,3,5])
    >>> data
    array((1, 1.3, 'abcde'),
          dtype=[('intvar', '<i8'), ('fltvar', '<f8'), ('strvar', '|S5')])

    """
    #
    if usemask:
        from numpy.ma import MaskedArray, make_mask_descr
    # Check the input dictionary of converters
    user_converters = converters or {}
    if not isinstance(user_converters, dict):
        errmsg = "The input argument 'converter' should be a valid dictionary "\
            "(got '%s' instead)"
        raise TypeError(errmsg % type(user_converters))

    # Initialize the filehandle, the LineSplitter and the NameValidator
    if isinstance(fname, basestring):
        fhd = np.lib._datasource.open(fname)
    elif not hasattr(fname, 'read'):
        raise TypeError("The input should be a string or a filehandle. "\
                        "(got %s instead)" % type(fname))
    else:
        fhd = fname
    split_line = LineSplitter(delimiter=delimiter, comments=comments,
                              autostrip=autostrip)._handyman
    validate_names = NameValidator(excludelist=excludelist,
                                   deletechars=deletechars,
                                   case_sensitive=case_sensitive)

    # Get the first valid lines after the first skiprows ones ..
    if skiprows:
        warnings.warn("The use of `skiprows` is deprecated.\n"\
                      "Please use `skip_header` instead.",
                      DeprecationWarning)
        skip_header = skiprows
    # Skip the first `skip_header` rows
    for i in xrange(skip_header):
        fhd.readline()
    # Keep on until we find the first valid values
    first_values = None
    while not first_values:
        first_line = fhd.readline()
        if first_line == '':
            raise IOError('End-of-file reached before encountering data.')
        if names is True:
            if comments in first_line:
                first_line = ''.join(first_line.split(comments)[1])
        first_values = split_line(first_line)
    # Should we take the first values as names ?
    if names is True:
        fval = first_values[0].strip()
        if fval in comments:
            del first_values[0]

    # Check the columns to use
    if usecols is not None:
        try:
            usecols = [_.strip() for _ in usecols.split(",")]
        except AttributeError:
            try:
                usecols = list(usecols)
            except TypeError:
                usecols = [usecols, ]
    nbcols = len(usecols or first_values)

    # Check the names and overwrite the dtype.names if needed
    if names is True:
        names = validate_names([_.strip() for _ in first_values])
        first_line = ''
    elif _is_string_like(names):
        names = validate_names([_.strip() for _ in names.split(',')])
    elif names:
        names = validate_names(names)
    # Get the dtype
    if dtype is not None:
        dtype = easy_dtype(dtype, defaultfmt=defaultfmt, names=names)
        names = dtype.names
    # Make sure the names is a list (for 2.5)
    if names is not None:
        names = list(names)


    if usecols:
        for (i, current) in enumerate(usecols):
            # if usecols is a list of names, convert to a list of indices
            if _is_string_like(current):
                usecols[i] = names.index(current)
            elif current < 0:
                usecols[i] = current + len(first_values)
        # If the dtype is not None, make sure we update it
        if (dtype is not None) and (len(dtype) > nbcols):
            descr = dtype.descr
            dtype = np.dtype([descr[_] for _ in usecols])
            names = list(dtype.names)
        # If `names` is not None, update the names
        elif (names is not None) and (len(names) > nbcols):
            names = [names[_] for _ in usecols]


    # Process the missing values ...............................
    # Rename missing_values for convenience
    user_missing_values = missing_values or ()

    # Define the list of missing_values (one column: one list)
    missing_values = [list(['']) for _ in range(nbcols)]

    # We have a dictionary: process it field by field
    if isinstance(user_missing_values, dict):
        # Loop on the items
        for (key, val) in user_missing_values.items():
            # Is the key a string ?
            if _is_string_like(key):
                try:
                    # Transform it into an integer
                    key = names.index(key)
                except ValueError:
                    # We couldn't find it: the name must have been dropped, then
                    continue
            # Redefine the key as needed if it's a column number
            if usecols:
                try:
                    key = usecols.index(key)
                except ValueError:
                    pass
            # Transform the value as a list of string
            if isinstance(val, (list, tuple)):
                val = [str(_) for _ in val]
            else:
                val = [str(val), ]
            # Add the value(s) to the current list of missing
            if key is None:
                # None acts as default
                for miss in missing_values:
                    miss.extend(val)
            else:
                missing_values[key].extend(val)
    # We have a sequence : each item matches a column
    elif isinstance(user_missing_values, (list, tuple)):
        for (value, entry) in zip(user_missing_values, missing_values):
            value = str(value)
            if value not in entry:
                entry.append(value)
    # We have a string : apply it to all entries
    elif isinstance(user_missing_values, basestring):
        user_value = user_missing_values.split(",")
        for entry in missing_values:
            entry.extend(user_value)
    # We have something else: apply it to all entries
    else:
        for entry in missing_values:
            entry.extend([str(user_missing_values)])

    # Process the deprecated `missing`
    if missing != '':
        warnings.warn("The use of `missing` is deprecated.\n"\
                      "Please use `missing_values` instead.",
                      DeprecationWarning)
        values = [str(_) for _ in missing.split(",")]
        for entry in missing_values:
            entry.extend(values)

    # Process the filling_values ...............................
    # Rename the input for convenience
    user_filling_values = filling_values or []
    # Define the default
    filling_values = [None] * nbcols
    # We have a dictionary : update each entry individually
    if isinstance(user_filling_values, dict):
        for (key, val) in user_filling_values.items():
            if _is_string_like(key):
                try:
                    # Transform it into an integer
                    key = names.index(key)
                except ValueError:
                    # We couldn't find it: the name must have been dropped, then
                    continue
            # Redefine the key if it's a column number and usecols is defined
            if usecols:
                try:
                    key = usecols.index(key)
                except ValueError:
                    pass
            # Add the value to the list
            filling_values[key] = val
    # We have a sequence : update on a one-to-one basis
    elif isinstance(user_filling_values, (list, tuple)):
        n = len(user_filling_values)
        if (n <= nbcols):
            filling_values[:n] = user_filling_values
        else:
            filling_values = user_filling_values[:nbcols]
    # We have something else : use it for all entries
    else:
        filling_values = [user_filling_values] * nbcols

    # Initialize the converters ................................
    if dtype is None:
        # Note: we can't use a [...]*nbcols, as we would have 3 times the same
        # ... converter, instead of 3 different converters.
        converters = [StringConverter(None, missing_values=miss, default=fill)
                      for (miss, fill) in zip(missing_values, filling_values)]
    else:
        dtype_flat = flatten_dtype(dtype, flatten_base=True)
        # Initialize the converters
        if len(dtype_flat) > 1:
            # Flexible type : get a converter from each dtype
            zipit = zip(dtype_flat, missing_values, filling_values)
            converters = [StringConverter(dt, locked=True,
                                          missing_values=miss, default=fill)
                           for (dt, miss, fill) in zipit]
        else:
            # Set to a default converter (but w/ different missing values)
            zipit = zip(missing_values, filling_values)
            converters = [StringConverter(dtype, locked=True,
                                          missing_values=miss, default=fill)
                          for (miss, fill) in zipit]
    # Update the converters to use the user-defined ones
    uc_update = []
    for (i, conv) in user_converters.items():
        # If the converter is specified by column names, use the index instead
        if _is_string_like(i):
            try:
                i = names.index(i)
            except ValueError:
                continue
        elif usecols:
            try:
                i = usecols.index(i)
            except ValueError:
                # Unused converter specified
                continue
        converters[i].update(conv, locked=True,
                             default=filling_values[i],
                             missing_values=missing_values[i],)
        uc_update.append((i, conv))
    # Make sure we have the corrected keys in user_converters...
    user_converters.update(uc_update)

    miss_chars = [_.missing_values for _ in converters]


    # Initialize the output lists ...
    # ... rows
    rows = []
    append_to_rows = rows.append
    # ... masks
    if usemask:
        masks = []
        append_to_masks = masks.append
    # ... invalid
    invalid = []
    append_to_invalid = invalid.append

    # Parse each line
    for (i, line) in enumerate(itertools.chain([first_line, ], fhd)):
        values = split_line(line)
        nbvalues = len(values)
        # Skip an empty line
        if nbvalues == 0:
            continue
        # Select only the columns we need
        if usecols:
            try:
                values = [values[_] for _ in usecols]
            except IndexError:
                append_to_invalid((i, nbvalues))
                continue
        elif nbvalues != nbcols:
            append_to_invalid((i, nbvalues))
            continue
        # Store the values
        append_to_rows(tuple(values))
        if usemask:
            append_to_masks(tuple([v.strip() in m
                                   for (v, m) in zip(values, missing_values)]))

    # Strip the last skip_footer data
    if skip_footer > 0:
        rows = rows[:-skip_footer]
        if usemask:
            masks = masks[:-skip_footer]

    # Upgrade the converters (if needed)
    if dtype is None:
        for (i, converter) in enumerate(converters):
            current_column = map(itemgetter(i), rows)
            try:
                converter.iterupgrade(current_column)
            except ConverterLockError:
                errmsg = "Converter #%i is locked and cannot be upgraded: " % i
                current_column = itertools.imap(itemgetter(i), rows)
                for (j, value) in enumerate(current_column):
                    try:
                        converter.upgrade(value)
                    except (ConverterError, ValueError):
                        errmsg += "(occurred line #%i for value '%s')"
                        errmsg %= (j + 1 + skip_header, value)
                        raise ConverterError(errmsg)

    # Check that we don't have invalid values
    if len(invalid) > 0:
        nbrows = len(rows)
        # Construct the error message
        template = "    Line #%%i (got %%i columns instead of %i)" % nbcols
        if skip_footer > 0:
            nbrows -= skip_footer
            errmsg = [template % (i + skip_header + 1, nb)
                      for (i, nb) in invalid if i < nbrows]
        else:
            errmsg = [template % (i + skip_header + 1, nb)
                      for (i, nb) in invalid]
        if len(errmsg):
            errmsg.insert(0, "Some errors were detected !")
            errmsg = "\n".join(errmsg)
            # Raise an exception ?
            if invalid_raise:
                raise ValueError(errmsg)
            # Issue a warning ?
            else:
                warnings.warn(errmsg, ConversionWarning)

    # Convert each value according to the converter:
    # We want to modify the list in place to avoid creating a new one...
#    if loose:
#        conversionfuncs = [conv._loose_call for conv in converters]
#    else:
#        conversionfuncs = [conv._strict_call for conv in converters]
#    for (i, vals) in enumerate(rows):
#        rows[i] = tuple([convert(val)
#                         for (convert, val) in zip(conversionfuncs, vals)])
    if loose:
        rows = zip(*(map(converter._loose_call, map(itemgetter(i), rows))
                     for (i, converter) in enumerate(converters)))
    else:
        rows = zip(*(map(converter._strict_call, map(itemgetter(i), rows))
                     for (i, converter) in enumerate(converters)))
    # Reset the dtype
    data = rows
    if dtype is None:
        # Get the dtypes from the types of the converters
        column_types = [conv.type for conv in converters]
        # Find the columns with strings...
        strcolidx = [i for (i, v) in enumerate(column_types)
                     if v in (type('S'), np.string_)]
        # ... and take the largest number of chars.
        for i in strcolidx:
            column_types[i] = "|S%i" % max(len(row[i]) for row in data)
        #
        if names is None:
            # If the dtype is uniform, don't define names, else use ''
            base = set([c.type for c in converters if c._checked])
            if len(base) == 1:
                (ddtype, mdtype) = (list(base)[0], np.bool)
            else:
                ddtype = [(defaultfmt % i, dt)
                          for (i, dt) in enumerate(column_types)]
                if usemask:
                    mdtype = [(defaultfmt % i, np.bool)
                              for (i, dt) in enumerate(column_types)]
        else:
            ddtype = zip(names, column_types)
            mdtype = zip(names, [np.bool] * len(column_types))
        output = np.array(data, dtype=ddtype)
        if usemask:
            outputmask = np.array(masks, dtype=mdtype)
    else:
        # Overwrite the initial dtype names if needed
        if names and dtype.names:
            dtype.names = names
        # Case 1. We have a structured type
        if len(dtype_flat) > 1:
            # Nested dtype, eg  [('a', int), ('b', [('b0', int), ('b1', 'f4')])]
            # First, create the array using a flattened dtype:
            # [('a', int), ('b1', int), ('b2', float)]
            # Then, view the array using the specified dtype.
            if 'O' in (_.char for _ in dtype_flat):
                if has_nested_fields(dtype):
                    errmsg = "Nested fields involving objects "\
                             "are not supported..."
                    raise NotImplementedError(errmsg)
                else:
                    output = np.array(data, dtype=dtype)
            else:
                rows = np.array(data, dtype=[('', _) for _ in dtype_flat])
                output = rows.view(dtype)
            # Now, process the rowmasks the same way
            if usemask:
                rowmasks = np.array(masks,
                                    dtype=np.dtype([('', np.bool)
                                    for t in dtype_flat]))
                # Construct the new dtype
                mdtype = make_mask_descr(dtype)
                outputmask = rowmasks.view(mdtype)
        # Case #2. We have a basic dtype
        else:
            # We used some user-defined converters
            if user_converters:
                ishomogeneous = True
                descr = []
                for (i, ttype) in enumerate([conv.type for conv in converters]):
                    # Keep the dtype of the current converter
                    if i in user_converters:
                        ishomogeneous &= (ttype == dtype.type)
                        if ttype == np.string_:
                            ttype = "|S%i" % max(len(row[i]) for row in data)
                        descr.append(('', ttype))
                    else:
                        descr.append(('', dtype))
                # So we changed the dtype ?
                if not ishomogeneous:
                    # We have more than one field
                    if len(descr) > 1:
                        dtype = np.dtype(descr)
                    # We have only one field: drop the name if not needed.
                    else:
                        dtype = np.dtype(ttype)
            #
            output = np.array(data, dtype)
            if usemask:
                if dtype.names:
                    mdtype = [(_, np.bool) for _ in dtype.names]
                else:
                    mdtype = np.bool
                outputmask = np.array(masks, dtype=mdtype)
    # Try to take care of the missing data we missed
    names = output.dtype.names
    if usemask and names:
        for (name, conv) in zip(names or (), converters):
            missing_values = [conv(_) for _ in conv.missing_values if _ != '']
            for mval in missing_values:
                outputmask[name] |= (output[name] == mval)
    # Construct the final array
    if usemask:
        output = output.view(MaskedArray)
        output._mask = outputmask
    if unpack:
        return output.squeeze().T
    return output.squeeze()
