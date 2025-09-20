"""Functionality related with filters in a PyTables file."""

from __future__ import annotations

import warnings
from typing import Any, Literal, TYPE_CHECKING

import numpy as np
from packaging.version import Version

from . import utilsextension
from .exceptions import FiltersWarning

if TYPE_CHECKING:
    from .leaf import Leaf

blosc_version = Version(utilsextension.which_lib_version("blosc")[1])
blosc2_version = Version(utilsextension.which_lib_version("blosc2")[1])
blosc_compcode_to_compname = utilsextension.blosc_compcode_to_compname_
blosc2_compcode_to_compname = utilsextension.blosc2_compcode_to_compname_


__docformat__ = "reStructuredText"
"""The format of documentation strings in this module."""

all_complibs = ["zlib", "lzo", "bzip2", "blosc", "blosc2"]
all_complibs += [
    f"blosc:{cname}" for cname in utilsextension.blosc_compressor_list()
]
all_complibs += [
    f"blosc2:{cname}" for cname in utilsextension.blosc2_compressor_list()
]


"""List of all compression libraries."""

foreign_complibs = ["szip"]
"""List of known but unsupported compression libraries."""

default_complib = "zlib"
"""The default compression library."""


_shuffle_flag = 0x1
_fletcher32_flag = 0x2
_rounding_flag = 0x4
_bitshuffle_flag = 0x8


class Filters:
    """Container for filter properties.

    This class is meant to serve as a container that keeps information about
    the filter properties associated with the chunked leaves, that is Table,
    CArray, EArray and VLArray.

    Instances of this class can be directly compared for equality.

    Parameters
    ----------
    complevel : int
        Specifies a compression level for data. The allowed
        range is 0-9. A value of 0 (the default) disables
        compression.
    complib : str
        Specifies the compression library to be used. Right now, 'zlib' (the
        default), 'lzo', 'bzip2', 'blosc' and 'blosc2' are supported.
        Additional compressors for Blosc like 'blosc:blosclz' ('blosclz' is
        the default in case the additional compressor is not specified),
        'blosc:lz4', 'blosc:lz4hc', 'blosc:zlib' and 'blosc:zstd' are
        supported too.
        Also, additional compressors for Blosc2 like 'blosc2:blosclz'
        ('blosclz' is the default in case the additional compressor is not
        specified), 'blosc2:lz4', 'blosc2:lz4hc', 'blosc2:zlib' and
        'blosc2:zstd' are supported too.
        Specifying a compression library which is not available
        in the system issues a FiltersWarning and sets the library to the
        default one.
    shuffle : bool
        Whether to use the *Shuffle* filter in the HDF5 library.
        This is normally used to improve the compression ratio.
        A false value disables shuffling and a true one enables
        it. The default value depends on whether compression is
        enabled or not; if compression is enabled, shuffling defaults
        to be enabled, else shuffling is disabled. Shuffling can only
        be used when compression is enabled.
    bitshuffle : bool
        Whether to use the *BitShuffle* filter in the Blosc/Blosc2
        libraries. This is normally used to improve the compression
        ratio. A false value disables bitshuffling and a true one
        enables it. The default value is disabled.
    fletcher32 : bool
        Whether to use the *Fletcher32* filter in the HDF5 library.
        This is used to add a checksum on each data chunk. A false
        value (the default) disables the checksum.
    least_significant_digit : int
        If specified, data will be truncated (quantized). In conjunction
        with enabling compression, this produces 'lossy', but
        significantly more efficient compression. For example, if
        *least_significant_digit=1*, data will be quantized using
        ``around(scale*data)/scale``, where ``scale = 2**bits``, and
        bits is determined so that a precision of 0.1 is retained (in
        this case bits=4). Default is *None*, or no quantization.

        .. note::

            quantization is only applied if some form of compression is
            enabled

    Examples
    --------
    This is a small example on using the Filters class::

        import numpy as np
        import tables as tb

        fileh = tb.open_file('test5.h5', mode='w')
        atom = Float32Atom()
        filters = Filters(complevel=1, complib='blosc', fletcher32=True)
        arr = fileh.create_earray(fileh.root, 'earray', atom, (0,2),
                                 "A growable array", filters=filters)

        # Append several rows in only one call
        arr.append(np.array([[1., 2.],
                             [2., 3.],
                             [3., 4.]], dtype=np.float32))

        # Print information on that enlargeable array
        print("Result Array:")
        print(repr(arr))
        fileh.close()

    This enforces the use of the Blosc library, a compression level of 1 and a
    Fletcher32 checksum filter as well. See the output of this example::

        Result Array:
        /earray (EArray(3, 2), fletcher32, shuffle, blosc(1)) 'A growable ...
        type = float32
        shape = (3, 2)
        itemsize = 4
        nrows = 3
        extdim = 0
        flavor = 'numpy'
        byteorder = 'little'

    .. rubric:: Filters attributes

    .. attribute:: fletcher32

        Whether the *Fletcher32* filter is active or not.

    .. attribute:: complevel

        The compression level (0 disables compression).

    .. attribute:: complib

        The compression filter used (irrelevant when compression is not
        enabled).

    .. attribute:: shuffle

        Whether the *Shuffle* filter is active or not.

    .. attribute:: bitshuffle

        Whether the *BitShuffle* filter is active or not (Blosc/Blosc2 only).

    """

    @property
    def shuffle_bitshuffle(self) -> Literal[0, 1, 2]:
        """Encode NoShuffle (0), Shuffle (1) and BitShuffle (2) filters."""
        if self.shuffle and self.bitshuffle:
            raise ValueError(
                "Shuffle and BitShuffle cannot be active at the same time"
            )
        if not (self.shuffle or self.bitshuffle):
            return 0
        if self.shuffle:
            return 1
        if self.bitshuffle:
            return 2

    @classmethod
    def _from_leaf(cls, leaf: Leaf) -> Filters:
        # Get a dictionary with all the filters
        parent = leaf._v_parent
        filters_dict = utilsextension.get_filters(
            parent._v_objectid, leaf._v_name
        )
        if filters_dict is None:
            filters_dict = {}  # not chunked

        # Keyword arguments are all off
        kwargs = {
            "complevel": 0,
            "shuffle": False,
            "bitshuffle": False,
            "fletcher32": False,
            "least_significant_digit": None,
            "_new": False,
        }
        for name, values in filters_dict.items():
            if name == "deflate":
                name = "zlib"
            if name in all_complibs:
                kwargs["complib"] = name
                if name in ("blosc", "blosc2"):
                    kwargs["complevel"] = values[4]
                    if values[5] == 1:
                        # Shuffle filter is internal to blosc/blosc2
                        kwargs["shuffle"] = True
                    elif values[5] == 2:
                        # Shuffle filter is internal to blosc/blosc2
                        kwargs["bitshuffle"] = True
                    # From Blosc 1.3 on, parameter 6 is used for the compressor
                    if len(values) > 6:
                        if name == "blosc":
                            cname = blosc_compcode_to_compname(values[6])
                            kwargs["complib"] = "blosc:%s" % cname
                        else:
                            cname = blosc2_compcode_to_compname(values[6])
                            kwargs["complib"] = "blosc2:%s" % cname
                else:
                    kwargs["complevel"] = values[0]
            elif name in foreign_complibs:
                kwargs["complib"] = name
                kwargs["complevel"] = 1  # any nonzero value will do
            elif name in ["shuffle", "fletcher32"]:
                kwargs[name] = True
        return cls(**kwargs)

    @classmethod
    def _unpack(cls, packed: int) -> Filters:
        """Create a new `Filters` object from a packed version.

        >>> Filters._unpack(0)
        Filters(complevel=0, shuffle=False, bitshuffle=False, \
fletcher32=False, least_significant_digit=None)
        >>> Filters._unpack(0x101)
        Filters(complevel=1, complib='zlib', shuffle=False, \
bitshuffle=False, fletcher32=False, least_significant_digit=None)
        >>> Filters._unpack(0x30109)
        Filters(complevel=9, complib='zlib', shuffle=True, \
bitshuffle=False, fletcher32=True, least_significant_digit=None)
        >>> Filters._unpack(0x3010A)
        Traceback (most recent call last):
          ...
        ValueError: compression level must be between 0 and 9
        >>> Filters._unpack(0x1)
        Traceback (most recent call last):
          ...
        ValueError: invalid compression library id: 0

        """
        kwargs = {"_new": False}

        # Byte 0: compression level.
        kwargs["complevel"] = complevel = packed & 0xFF
        packed >>= 8

        # Byte 1: compression library id (0 for none).
        if complevel > 0:
            complib_id = int(packed & 0xFF)
            if not (0 < complib_id <= len(all_complibs)):
                raise ValueError(
                    f"invalid compression library id: {complib_id}"
                )
            kwargs["complib"] = all_complibs[complib_id - 1]
        packed >>= 8

        # Byte 2: parameterless filters.
        kwargs["shuffle"] = packed & _shuffle_flag
        kwargs["bitshuffle"] = packed & _bitshuffle_flag
        kwargs["fletcher32"] = packed & _fletcher32_flag
        has_rounding = packed & _rounding_flag
        packed >>= 8

        # Byte 3: least significant digit.
        if has_rounding:
            kwargs["least_significant_digit"] = np.int8(packed & 0xFF)
        else:
            kwargs["least_significant_digit"] = None

        return cls(**kwargs)

    def _pack(self) -> np.int64:
        """Pack the `Filters` object into a 64-bit NumPy integer."""
        packed = np.int64(0)

        # Byte 3: least significant digit.
        if self.least_significant_digit is not None:
            # assert isinstance(self.least_significant_digit, np.int8)
            packed |= self.least_significant_digit
        packed <<= 8

        # Byte 2: parameterless filters.
        if self.shuffle:
            packed |= _shuffle_flag
        if self.bitshuffle:
            packed |= _bitshuffle_flag
        if self.fletcher32:
            packed |= _fletcher32_flag
        if self.least_significant_digit:
            packed |= _rounding_flag
        packed <<= 8

        # Byte 1: compression library id (0 for none).
        if self.complevel > 0:
            packed |= all_complibs.index(self.complib) + 1
        packed <<= 8

        # Byte 0: compression level.
        packed |= self.complevel

        return packed

    def __init__(
        self,
        complevel: int = 0,
        complib: Literal[
            "zlib", "lzo", "bzip2", "blosc", "blosc2"
        ] = default_complib,
        shuffle: bool = True,
        bitshuffle: bool = False,
        fletcher32: bool = False,
        least_significant_digit: int | None = None,
        _new: bool = True,
    ) -> None:

        if not (0 <= complevel <= 9):
            raise ValueError("compression level must be between 0 and 9")

        if _new and complevel > 0:
            # These checks are not performed when loading filters from disk.
            if complib not in all_complibs:
                raise ValueError(
                    "compression library ``%s`` is not supported; "
                    "it must be one of: %s"
                    % (complib, ", ".join(all_complibs))
                )
            if utilsextension.which_lib_version(complib) is None:
                warnings.warn(
                    "compression library ``%s`` is not available; "
                    "using ``%s`` instead" % (complib, default_complib),
                    FiltersWarning,
                )
                complib = default_complib  # always available

        complevel = int(complevel)
        complib = str(complib)
        shuffle = bool(shuffle)
        bitshuffle = bool(bitshuffle)
        fletcher32 = bool(fletcher32)
        if least_significant_digit is not None:
            least_significant_digit = np.int8(least_significant_digit)

        if complevel == 0:
            # Override some inputs when compression is not enabled.
            complib = None  # make it clear there is no compression
            shuffle = False  # shuffling and not compressing makes no sense
            least_significant_digit = None
        elif complib not in all_complibs:
            # Do not try to use a meaningful level for unsupported libs.
            complevel = -1

        self.complevel = complevel
        """The compression level (0 disables compression)."""

        self.complib = complib
        """The compression filter used (irrelevant when compression is
        not enabled).
        """

        self.shuffle = shuffle
        """Whether the *Shuffle* filter is active or not."""

        self.bitshuffle = bitshuffle
        """Whether the *BitShuffle* filter is active or not."""

        if (
            self.complib
            and self.bitshuffle
            and not self.complib.startswith("blosc")
        ):
            raise ValueError("BitShuffle can only be used inside Blosc/Blosc2")

        if self.shuffle and self.bitshuffle:
            # BitShuffle has priority in case both are specified
            self.shuffle = False

        self.fletcher32 = fletcher32
        """Whether the *Fletcher32* filter is active or not."""

        self.least_significant_digit = least_significant_digit
        """The least significant digit to which data shall be truncated."""

    def __repr__(self) -> str:
        args = []
        if self.complevel >= 0:  # meaningful compression level
            args.append(f"complevel={self.complevel}")
        if self.complevel != 0:  # compression enabled (-1 or > 0)
            args.append(f"complib={self.complib!r}")
        args.append(f"shuffle={self.shuffle}")
        args.append(f"bitshuffle={self.bitshuffle}")
        args.append(f"fletcher32={self.fletcher32}")
        args.append(f"least_significant_digit={self.least_significant_digit}")
        return f'{self.__class__.__name__}({", ".join(args)})'

    def __str__(self) -> str:
        return repr(self)

    def __eq__(self, other: Any) -> bool:
        if not isinstance(other, self.__class__):
            return False
        for attr in self.__dict__:
            if getattr(self, attr) != getattr(other, attr):
                return False
        return True

    # XXX: API incompatible change for PyTables 3 line
    # Overriding __eq__ blocks inheritance of __hash__ in 3.x
    # def __hash__(self):
    #    return hash((self.__class__, self.complevel, self.complib,
    #                 self.shuffle, self.bitshuffle, self.fletcher32))

    def copy(self, **override) -> Filters:
        """Get a copy of the filters, possibly overriding some arguments.

        Constructor arguments to be overridden must be passed as keyword
        arguments.

        Using this method is recommended over replacing the attributes of an
        instance, since instances of this class may become immutable in the
        future::

            >>> filters1 = Filters()
            >>> filters2 = filters1.copy()
            >>> filters1 == filters2
            True
            >>> filters1 is filters2
            False
            >>> filters3 = filters1.copy(complevel=1) #doctest: +ELLIPSIS
            Traceback (most recent call last):
            ...
            ValueError: compression library ``None`` is not supported...
            >>> filters3 = filters1.copy(complevel=1, complib='zlib')
            >>> print(filters1)
            Filters(complevel=0, shuffle=False, bitshuffle=False, \
fletcher32=False, least_significant_digit=None)
            >>> print(filters3)
            Filters(complevel=1, complib='zlib', shuffle=False, \
bitshuffle=False, fletcher32=False, least_significant_digit=None)
            >>> filters1.copy(foobar=42) #doctest: +ELLIPSIS
            Traceback (most recent call last):
            ...
            TypeError: ...__init__() got an unexpected keyword argument ...

        """
        newargs = self.__dict__.copy()
        newargs.update(override)
        return self.__class__(**newargs)


def _test() -> None:
    """Run ``doctest`` on this module."""
    import doctest

    doctest.testmod()


if __name__ == "__main__":
    _test()
