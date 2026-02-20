"""Utilities for handling different array flavors in PyTables.

Variables
=========

`__docformat`__
    The format of documentation strings in this module.
`internal_flavor`
    The flavor used internally by PyTables.
`all_flavors`
    List of all flavors available to PyTables.
`alias_map`
    Maps old flavor names to the most similar current flavor.
`description_map`
    Maps flavors to short descriptions of their supported objects.
`identifier_map`
    Maps flavors to functions that can identify their objects.

    The function associated with a given flavor will return a true
    value if the object passed to it can be identified as being of
    that flavor.

    See the `flavor_of()` function for a friendlier interface to
    flavor identification.

`converter_map`
    Maps (source, destination) flavor pairs to converter functions.

    Converter functions get an array of the source flavor and return
    an array of the destination flavor.

    See the `array_of_flavor()` and `flavor_to_flavor()` functions for
    friendlier interfaces to flavor conversion.

"""

from __future__ import annotations

import warnings
from typing import Any, Literal
from collections.abc import Callable, Sequence

import numpy as np
import numpy.typing as npt

from .exceptions import FlavorError, FlavorWarning

# @compatibility with numpy < 1.25
try:
    from numpy.exceptions import VisibleDeprecationWarning
except ImportError:
    from numpy import VisibleDeprecationWarning

__docformat__ = "reStructuredText"
"""The format of documentation strings in this module."""

FlavorType = Literal["numpy", "python"]

internal_flavor: FlavorType = "numpy"
"""The flavor used internally by PyTables."""

# This is very slightly slower than a set for a small number of values
# in terms of (infrequent) lookup time, but allows `flavor_of()`
# (which may be called much more frequently) to check for flavors in
# order, beginning with the most common one.
all_flavors: list[FlavorType] = []  # filled as flavors are registered
"""List of all flavors available to PyTables."""

alias_map: dict[str, FlavorType] = {}  # filled as flavors are registered
"""Maps old flavor names to the most similar current flavor."""

description_map: dict[FlavorType, str] = {}  # filled as flavors are registered
"""Maps flavors to short descriptions of their supported objects."""

# filled as flavors are registered
identifier_map: dict[FlavorType, Callable[[Any], FlavorType]] = {}
"""Maps flavors to functions that can identify their objects.

The function associated with a given flavor will return a true value
if the object passed to it can be identified as being of that flavor.

See the `flavor_of()` function for a friendlier interface to flavor
identification.
"""

# filled as flavors are registered
converter_map: dict[tuple[FlavorType, FlavorType], Callable[[Any], Any]] = {}
"""Maps (source, destination) flavor pairs to converter functions.

Converter functions get an array of the source flavor and return an
array of the destination flavor.

See the `array_of_flavor()` and `flavor_to_flavor()` functions for
friendlier interfaces to flavor conversion.
"""


def check_flavor(flavor: FlavorType) -> None:
    """Raise a ``FlavorError`` if the `flavor` is not valid."""
    if flavor not in all_flavors:
        available_flavs = ", ".join(flav for flav in all_flavors)
        raise FlavorError(
            "flavor ``%s`` is unsupported or unavailable; "
            "available flavors in this system are: %s"
            % (flavor, available_flavs)
        )


def array_of_flavor2(
    array: npt.ArrayLike, src_flavor: FlavorType, dst_flavor: FlavorType
) -> Any | list[Any] | np.ndarray:
    """Get a version of the given `array` in a different flavor.

    The input `array` must be of the given `src_flavor`, and the
    returned array will be of the indicated `dst_flavor`.  Both
    flavors may be the same, but it is not guaranteed that the
    returned array will be the same object as the input one in this
    case.

    If the conversion is not supported, a ``FlavorError`` is raised.

    """
    convkey = (src_flavor, dst_flavor)
    if convkey not in converter_map:
        raise FlavorError(
            "conversion from flavor ``%s`` to flavor ``%s`` "
            "is unsupported or unavailable in this system"
            % (src_flavor, dst_flavor)
        )

    convfunc = converter_map[convkey]
    return convfunc(array)


def flavor_to_flavor(
    array: npt.ArrayLike, src_flavor: FlavorType, dst_flavor: FlavorType
) -> Any | list[Any] | np.ndarray | npt.ArrayLike:
    """Get a version of the given `array` in a different flavor.

    The input `array` must be of the given `src_flavor`, and the
    returned array will be of the indicated `dst_flavor` (see below
    for an exception to this).  Both flavors may be the same, but it
    is not guaranteed that the returned array will be the same object
    as the input one in this case.

    If the conversion is not supported, a `FlavorWarning` is issued
    and the input `array` is returned as is.

    """
    try:
        return array_of_flavor2(array, src_flavor, dst_flavor)
    except FlavorError as fe:
        warnings.warn(
            "%s; returning an object of the ``%s`` flavor instead"
            % (fe.args[0], src_flavor),
            FlavorWarning,
        )
        return array


def internal_to_flavor(
    array: npt.ArrayLike, dst_flavor: FlavorType
) -> Any | list[Any] | np.ndarray:
    """Get a version of the given `array` in a different `dst_flavor`.

    The input `array` must be of the internal flavor, and the returned
    array will be of the given `dst_flavor`.  See `flavor_to_flavor()`
    for more information.

    """
    return flavor_to_flavor(array, internal_flavor, dst_flavor)


def array_as_internal(
    array: npt.ArrayLike, src_flavor: FlavorType
) -> np.ndarray:
    """Get a version of the given `array` in the internal flavor.

    The input `array` must be of the given `src_flavor`, and the
    returned array will be of the internal flavor.

    If the conversion is not supported, a ``FlavorError`` is raised.

    """
    return array_of_flavor2(array, src_flavor, internal_flavor)


def flavor_of(array: npt.ArrayLike) -> FlavorType:
    """Identify the flavor of a given `array`.

    If the `array` can not be matched with any flavor, a ``TypeError``
    is raised.

    """
    for flavor in all_flavors:
        if identifier_map[flavor](array):
            return flavor
    type_name = type(array).__name__
    supported_descs = "; ".join(description_map[fl] for fl in all_flavors)
    raise TypeError(
        "objects of type ``%s`` are not supported in this context, sorry; "
        "supported objects are: %s" % (type_name, supported_descs)
    )


def array_of_flavor(
    array: npt.ArrayLike, dst_flavor: FlavorType
) -> Any | list[Any] | np.ndarray:
    """Get a version of the given `array` in a different `dst_flavor`.

    The flavor of the input `array` is guessed, and the returned array
    will be of the given `dst_flavor`.

    If the conversion is not supported, a ``FlavorError`` is raised.

    """
    return array_of_flavor2(array, flavor_of(array), dst_flavor)


def restrict_flavors(keep: Sequence[FlavorType] = ("python",)) -> None:
    """Disable all flavors except those in keep.

    Providing an empty keep sequence implies disabling all flavors (but the
    internal one).  If the sequence is not specified, only optional flavors are
    disabled.

    .. important:: Once you disable a flavor, it can not be enabled again.

    """
    remove = set(all_flavors) - set(keep) - {internal_flavor}
    for flavor in remove:
        _disable_flavor(flavor)


# Flavor registration
#
# The order in which flavors appear in `all_flavors` determines the
# order in which they will be tested for by `flavor_of()`, so place
# most frequent flavors first.
all_flavors.append("numpy")  # this is the internal flavor

all_flavors.append("python")  # this is always supported


def _register_aliases() -> None:
    """Register aliases of *available* flavors."""
    for flavor in all_flavors:
        aliases = eval("_%s_aliases" % flavor)
        for alias in aliases:
            alias_map[alias] = flavor


def _register_descriptions() -> None:
    """Register descriptions of *available* flavors."""
    for flavor in all_flavors:
        description_map[flavor] = eval("_%s_desc" % flavor)


def _register_identifiers() -> None:
    """Register identifier functions of *available* flavors."""
    for flavor in all_flavors:
        identifier_map[flavor] = eval("_is_%s" % flavor)


def _register_converters() -> None:
    """Register converter functions between *available* flavors."""

    def identity(array: Any) -> Any:
        return array

    for src_flavor in all_flavors:
        for dst_flavor in all_flavors:
            # Converters with the same source and destination flavor
            # are used when available, since they may perform some
            # optimizations on the resulting array (e.g. making it
            # contiguous).  Otherwise, an identity function is used.
            convfunc = None
            try:
                convfunc = eval(f"_conv_{src_flavor}_to_{dst_flavor}")
            except NameError:
                if src_flavor == dst_flavor:
                    convfunc = identity
            if convfunc:
                converter_map[(src_flavor, dst_flavor)] = convfunc


def _register_all() -> None:
    """Register all *available* flavors."""
    _register_aliases()
    _register_descriptions()
    _register_identifiers()
    _register_converters()


def _deregister_aliases(flavor: FlavorType) -> None:
    """Deregister aliases of a given `flavor` (no checks)."""
    rm_aliases = []
    for an_alias, a_flavor in alias_map.items():
        if a_flavor == flavor:
            rm_aliases.append(an_alias)
    for an_alias in rm_aliases:
        del alias_map[an_alias]


def _deregister_description(flavor: FlavorType) -> None:
    """Deregister description of a given `flavor` (no checks)."""
    del description_map[flavor]


def _deregister_identifier(flavor: FlavorType) -> None:
    """Deregister identifier function of a given `flavor` (no checks)."""
    del identifier_map[flavor]


def _deregister_converters(flavor: FlavorType) -> None:
    """Deregister converter functions of a given `flavor` (no checks)."""
    rm_flavor_pairs = []
    for flavor_pair in converter_map:
        if flavor in flavor_pair:
            rm_flavor_pairs.append(flavor_pair)
    for flavor_pair in rm_flavor_pairs:
        del converter_map[flavor_pair]


def _disable_flavor(flavor: FlavorType) -> None:
    """Completely disable the given `flavor` (no checks)."""
    _deregister_aliases(flavor)
    _deregister_description(flavor)
    _deregister_identifier(flavor)
    _deregister_converters(flavor)
    all_flavors.remove(flavor)


# Implementation of flavors
_python_aliases: list[str] = [
    "List",
    "Tuple",
    "Int",
    "Float",
    "String",
    "VLString",
    "Object",
]
_python_desc = "homogeneous list or tuple, " "integer, float, complex or bytes"


def _is_python(array: Any) -> bool:
    return isinstance(array, (tuple, list, int, float, complex, bytes))


_numpy_aliases: list[str] = []
_numpy_desc = "NumPy array, record or scalar"


if np.lib.NumpyVersion(np.__version__) >= np.lib.NumpyVersion("1.19.0"):

    def toarray(array: npt.ArrayLike, *args, **kwargs) -> np.ndarray:
        """Convert the input to a numpy array if needed."""
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            try:
                array = np.array(array, *args, **kwargs)
            except VisibleDeprecationWarning:
                raise ValueError(
                    "cannot guess the desired dtype from the input"
                )

        return array

else:
    toarray = np.array


def _is_numpy(array: Any) -> bool:
    return isinstance(array, (np.ndarray, np.generic))


def _numpy_contiguous(
    convfunc: Callable[[npt.ArrayLike], np.ndarray]
) -> Callable[[npt.ArrayLike], np.ndarray]:
    """Decorate `convfunc` to return a *contiguous* NumPy array.

    Note: When arrays are 0-strided, the copy is avoided.  This allows
    to use `array` to still carry info about the dtype and shape.
    """

    def conv_to_numpy(array: npt.ArrayLike) -> np.ndarray:
        nparr = convfunc(array)
        if (
            hasattr(nparr, "flags")
            and not nparr.flags.contiguous
            and sum(nparr.strides) != 0
        ):
            nparr = nparr.copy()  # copying the array makes it contiguous
        return nparr

    conv_to_numpy.__name__ = convfunc.__name__
    conv_to_numpy.__doc__ = convfunc.__doc__
    return conv_to_numpy


@_numpy_contiguous
def _conv_numpy_to_numpy(array: npt.ArrayLike) -> np.ndarray:
    # Passes contiguous arrays through and converts scalars into
    # scalar arrays.
    nparr = np.asarray(array)
    if nparr.dtype.kind == "U":
        # from Python 3 loads of common strings are disguised as Unicode
        try:
            # try to convert to basic 'S' type
            return nparr.astype("S")
        except UnicodeEncodeError:
            pass
            # pass on true Unicode arrays downstream in case it can be
            # handled in the future
    return nparr


@_numpy_contiguous
def _conv_python_to_numpy(array: npt.ArrayLike) -> np.ndarray:
    nparr = toarray(array)
    if nparr.dtype.kind == "U":
        # from Python 3 loads of common strings are disguised as Unicode
        try:
            # try to convert to basic 'S' type
            return nparr.astype("S")
        except UnicodeEncodeError:
            pass
            # pass on true Unicode arrays downstream in case it can be
            # handled in the future
    return nparr


def _conv_numpy_to_python(array: np.ndarray) -> Any | list[Any]:
    if array.shape != ():
        # Lists are the default for returning multidimensional objects
        array = array.tolist()
    else:
        # 0-dim or scalar case
        array = array.item()
    return array


# Now register everything related with *available* flavors.
_register_all()


def _test() -> None:
    """Run ``doctest`` on this module."""
    import doctest

    doctest.testmod()


if __name__ == "__main__":
    _test()
