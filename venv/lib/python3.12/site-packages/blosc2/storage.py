#######################################################################
# Copyright (c) 2019-present, Blosc Development Team <blosc@blosc.org>
# All rights reserved.
#
# This source code is licensed under a BSD-style license (found in the
# LICENSE file in the root directory of this source tree)
#######################################################################

import contextlib
import warnings
from dataclasses import asdict, dataclass, field, fields

import blosc2


def default_nthreads():
    return blosc2.nthreads


def default_filters():
    return [
        blosc2.Filter.NOFILTER,
        blosc2.Filter.NOFILTER,
        blosc2.Filter.NOFILTER,
        blosc2.Filter.NOFILTER,
        blosc2.Filter.NOFILTER,
        blosc2.Filter.SHUFFLE,
    ]


def default_filters_meta():
    return [0] * 6


@dataclass
class CParams:
    """Dataclass for hosting the different compression parameters.

    Parameters
    ----------
    codec: :class:`Codec` or int
        The compressor code. Default is :py:obj:`Codec.ZSTD <Codec>`.
    codec_meta: int
        The metadata for the compressor code. Default is 0.
    clevel: int
        The compression level from 0 (no compression) to 9
        (maximum compression). Default is 1.
    use_dict: bool
        Whether to use dictionaries when compressing
        (only for :py:obj:`blosc2.Codec.ZSTD <Codec>`). Default is `False`.
    typesize: int
        The data type size, ranging from 1 to 255. Default is 8.
    nthreads: int
        The number of threads to use internally. By default, the
        value of :py:obj:`blosc2.nthreads` is used. If not set with
        :func:`blosc2.set_nthreads`, blosc2 computes a good guess for it.
    blocksize: int
        The requested size of the compressed blocks. If set to 0 (the default)
        blosc2 will choose the size automatically.
    splitmode: :class:`SplitMode`
        The split mode for the blocks.
        The default value is :py:obj:`SplitMode.AUTO_SPLIT <SplitMode>`.
    filters: :class:`Filter` or int list or None
        The sequence of filters. Default: [:py:obj:`Filter.NOFILTER <Filter>`,
        :py:obj:`Filter.NOFILTER <Filter>`, :py:obj:`Filter.NOFILTER <Filter>`, :py:obj:`Filter.NOFILTER <Filter>`,
        :py:obj:`Filter.NOFILTER <Filter>`, :py:obj:`Filter.SHUFFLE <Filter>`].
    filters_meta: list
        The metadata for filters. Default: `[0, 0, 0, 0, 0, 0]`.
    tuner: :class:`Tuner`
        The tuner to use. Default: :py:obj:`Tuner.STUNE <Tuner>`.
    """

    codec: blosc2.Codec | int = blosc2.Codec.ZSTD
    codec_meta: int = 0
    clevel: int = 5
    use_dict: bool = False
    typesize: int = 8
    nthreads: int = field(default_factory=default_nthreads)
    blocksize: int = 0
    splitmode: blosc2.SplitMode = blosc2.SplitMode.AUTO_SPLIT
    filters: list[blosc2.Filter | int] = field(default_factory=default_filters)
    filters_meta: list[int] = field(default_factory=default_filters_meta)
    tuner: blosc2.Tuner = blosc2.Tuner.STUNE

    def __post_init__(self):
        # C2Array sends metadata (like codec, filters, splitmode and tuner) as ints
        if not isinstance(self.codec, blosc2.Codec):
            with contextlib.suppress(ValueError):
                # User-defined codecs may have no entries in Codec
                self.codec = blosc2.Codec(self.codec)
        if not isinstance(self.splitmode, blosc2.SplitMode):
            with contextlib.suppress(ValueError):
                self.splitmode = blosc2.SplitMode(self.splitmode)
        if not isinstance(self.tuner, blosc2.Tuner):
            with contextlib.suppress(ValueError):
                self.tuner = blosc2.Tuner(self.tuner)

        if len(self.filters) > 6:
            raise ValueError("Number of filters exceeds 6")
        if len(self.filters) < len(self.filters_meta):
            self.filters_meta = self.filters_meta[: len(self.filters)]
            # There is no need to raise a warning here
            # warnings.warn("Changed `filters_meta` length to match `filters` length")
        if len(self.filters) > len(self.filters_meta):
            raise ValueError("Number of filters cannot exceed number of filters meta")

        for i, filter_i in enumerate(self.filters):
            if not isinstance(filter_i, blosc2.Filter):
                with contextlib.suppress(ValueError):
                    # User-defined filters may have no entries in Filter
                    self.filters[i] = blosc2.Filter(filter_i)
            if self.filters_meta[i] == 0 and self.filters[i] == blosc2.Filter.BYTEDELTA:
                self.filters_meta[i] = self.typesize


@dataclass
class DParams:
    """Dataclass for hosting the different decompression parameters.

    Parameters
    ----------
    nthreads: int
        The number of threads to use internally. By default, the
        value of :py:obj:`blosc2.nthreads` is used. If not set with
        :func:`blosc2.set_nthreads`, blosc2 computes a good guess for it.
    """

    nthreads: int = field(default_factory=default_nthreads)


@dataclass
class Storage:
    """Dataclass for hosting the different storage parameters.

    Parameters
    ----------
    contiguous: bool
        Indicates whether the chunks are stored contiguously.
        Default is True when :paramref:`urlpath` is not None;
        False otherwise.
    urlpath: str or pathlib.Path, optional
        If the storage is persistent, the name of the file (when
        `contiguous = True`) or the directory (if `contiguous = False`).
        If the storage is in-memory, then this field is `None`.
    mode: str, optional
        Persistence mode: 'r' means read only (must exist);
        'a' means read/write (create if it doesn't exist);
        'w' means create (overwrite if it exists). Default is 'a'.
    mmap_mode: str, optional
        If set, the file will be memory-mapped instead of using the default
        I/O functions and the `mode` argument will be ignored. The memory-mapping
        modes are similar to those used by the
        `numpy.memmap <https://numpy.org/doc/stable/reference/generated/numpy.memmap.html>`_
        function, but it is possible to extend the file:

        .. list-table::
            :widths: 10 90
            :header-rows: 1

            * - mode
              - description
            * - 'r'
              - Open an existing file for reading only.
            * - 'r+'
              - Open an existing file for reading and writing. Use this mode if you want
                to append data to an existing schunk file.
            * - 'w+'
              - Create or overwrite an existing file for reading and writing. Use this
                mode if you want to create a new schunk.
            * - 'c'
              - Open an existing file in copy-on-write mode: all changes affect the data
                in memory but changes are not saved to disk. The file on disk is
                read-only. On Windows, the size of the mapping cannot change.

        Only contiguous storage can be memory-mapped. Hence, `urlpath` must point to a
        file (and not a directory).

        .. note::
            Memory-mapped files are opened once, and their contents remain in (virtual)
            memory for the lifetime of the schunk. Using memory-mapped I/O can be faster
            than the default I/O functions, depending on the use case. While
            reading performance is generally better, writing performance may be
            slower in some cases on certain systems. Memory-mapped files
            can be especially beneficial when operating with network file systems
            (like NFS).

            This is currently a beta feature (especially for write operations) and we
            recommend trying it out and reporting any issues you may encounter.

    initial_mapping_size: int, optional
        The initial size of the mapping for the memory-mapped file when writes are
        allowed (r+ w+, or c mode). Once a file is memory-mapped and extended beyond the
        initial mapping size, the file must be remapped, which may be expensive. This
        parameter allows decoupling the mapping size from the actual file size to
        reserve memory early for future writes and avoid remappings. The memory is only
        reserved virtually and does not occupy physical memory unless actual writes
        occur. Since the virtual address space is large enough, it is ok to be generous
        with this parameter (with special consideration on Windows, see note below).
        For best performance, set this to the maximum expected size of the compressed
        data (see example in :obj:`SChunk.__init__ <blosc2.schunk.SChunk.__init__>`).
        The size is in bytes.

        Default: 1 GiB.

        .. note::
            On Windows, the size of the mapping is directly coupled to the file size.
            When the schunk is destroyed, the file size will be truncated to the
            actual size of the schunk.

    meta: dict or None
        A dictionary with different metalayers.  Each entry represents a metalayer:

            key: bytes or str
                The name of the metalayer.
            value: object
                The metalayer object that will be serialized using msgpack.
    """

    contiguous: bool = None
    urlpath: str = None
    mode: str = "a"
    mmap_mode: str = None
    initial_mapping_size: int = None
    meta: dict = None

    def __post_init__(self):
        if self.contiguous is None:
            self.contiguous = self.urlpath is not None
        # Check for None values
        for f in fields(self):
            if getattr(self, f.name) is None and f.name not in [
                "urlpath",
                "mmap_mode",
                "initial_mapping_size",
                "meta",
            ]:
                setattr(self, f.name, getattr(Storage(), f.name))
                warnings.warn(f"`{f.name}` field value changed from `None` to `{getattr(self, f.name)}`")


# Defaults for compression params
cparams_dflts = asdict(CParams())
"""
Compression params defaults.
"""

# Defaults for decompression params
dparams_dflts = asdict(DParams())
"""
Decompression params defaults.
"""
# Default for storage
storage_dflts = asdict(Storage())
"""
Storage params defaults. This is meant only for :ref:`SChunk <SChunk>` or :ref:`NDArray <NDArray>`.
"""
