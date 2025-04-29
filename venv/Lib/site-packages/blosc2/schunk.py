#######################################################################
# Copyright (c) 2019-present, Blosc Development Team <blosc@blosc.org>
# All rights reserved.
#
# This source code is licensed under a BSD-style license (found in the
# LICENSE file in the root directory of this source tree)
#######################################################################
from __future__ import annotations

import os
import pathlib
from collections import namedtuple
from collections.abc import Iterator, Mapping, MutableMapping
from dataclasses import asdict
from typing import Any, NamedTuple

import numpy as np
from msgpack import packb, unpackb

import blosc2
from blosc2 import SpecialValue, blosc2_ext


class vlmeta(MutableMapping, blosc2_ext.vlmeta):
    """
    Class providing access to user metadata on an :ref:`SChunk`.
    It is available via the `.vlmeta` property of an :ref:`SChunk`.
    """

    def __init__(self, schunk, urlpath, mode, mmap_mode, initial_mapping_size):
        self.urlpath = urlpath
        self.mode = mode
        self.mmap_mode = mmap_mode
        self.initial_mapping_size = initial_mapping_size
        super().__init__(schunk)

    def __setitem__(self, name, content):
        blosc2_ext.check_access_mode(self.urlpath, self.mode)
        # If name is a slice, assume that content is a dictionary and copy all the items
        if isinstance(name, slice):
            if name.start is None and name.stop is None:
                for k, v in content.items():
                    self.set_vlmeta(k, v)
                return
            raise NotImplementedError("Slicing is not supported, unless [:]")
        cparams = {"typesize": 1}
        content = packb(
            content,
            default=blosc2_ext.encode_tuple,
            strict_types=True,
            use_bin_type=True,
        )
        super().set_vlmeta(name, content, **cparams)

    def __getitem__(self, name):
        if isinstance(name, slice):
            if name.start is None and name.stop is None:
                # Return all the vlmetalayers
                return self.getall()
            raise NotImplementedError("Slicing is not supported, unless [:]")
        return unpackb(super().get_vlmeta(name), list_hook=blosc2_ext.decode_tuple)

    def __delitem__(self, name):
        blosc2_ext.check_access_mode(self.urlpath, self.mode)
        super().del_vlmeta(name)

    def __len__(self):
        return super().nvlmetalayers()

    def __iter__(self):
        keys = super().get_names()
        yield from keys

    def getall(self):
        """
        Return all the variable length metalayers as a dictionary

        """
        return super().to_dict()

    def __repr__(self):
        return repr(self.getall())

    def __str__(self):
        return str(self.getall())


class Meta(Mapping):
    """
    Class providing access to fixed-length metadata on an :ref:`SChunk`.
    It is available via the `.meta` property of an :ref:`SChunk`.
    """

    def get(self, key: str, default: Any = None) -> Any:
        """Return the value for `key` if `key` is in the dictionary, else return `default`.
        If `default` is not given, it defaults to ``None``."""
        return self.get(key, default)

    def __init__(self, schunk):
        self.schunk = schunk

    def __contains__(self, key: str) -> bool:
        """Check if the `key` metalayer exists or not."""
        return blosc2_ext.meta__contains__(self.schunk, key)

    def __delitem__(self, key: str) -> None:
        raise NotImplementedError("Cannot remove a metalayer")

    def __setitem__(self, key: str, value: bytes) -> None:
        """Update the `key` metalayer with `value`.

        Parameters
        ----------
        key: str
            The name of the metalayer to update.
        value: bytes
            The buffer containing the new content for the metalayer.

            ..warning: Note that the *length* of the metalayer cannot change,
            otherwise an exception will be raised.
        """
        value = packb(value, default=blosc2_ext.encode_tuple, strict_types=True, use_bin_type=True)
        blosc2_ext.meta__setitem__(self.schunk, key, value)

    def __getitem__(self, item: str | slice) -> bytes | dict[str, bytes]:
        """Return the specified metalayer.

        Parameters
        ----------
        item: str or slice
            The name of the metalayer to return.  If a slice is passed,
            and start and stop are None ([:]), all the metalayers are returned;
            else, a NotImplementedError is raised.

        Returns
        -------
        bytes or dict
            The buffer containing the metalayer information. If a slice is passed,
            a dictionary with all the metalayers is returned.
        """
        if isinstance(item, slice):
            if item.start is None and item.stop is None:
                return self.getall()
            raise NotImplementedError("Slicing is not supported, unless [:]")
        if self.__contains__(item):
            return unpackb(
                blosc2_ext.meta__getitem__(self.schunk, item),
                list_hook=blosc2_ext.decode_tuple,
            )
        else:
            raise KeyError(f"{item} not found")

    def keys(self) -> list[str]:
        """Return the metalayers keys."""
        return blosc2_ext.meta_keys(self.schunk)

    def values(self):
        raise NotImplementedError("Values can not be accessed")

    def items(self):
        raise NotImplementedError("Items can not be accessed")

    def __iter__(self) -> Iterator[str]:
        """Iter over the keys of the metalayers."""
        return iter(self.keys())

    def __len__(self) -> int:
        """Return the number of metalayers."""
        return blosc2_ext.meta__len__(self.schunk)

    def getall(self):
        """
        Return all the variable length metalayers as a dictionary

        """
        return {key: self[key] for key in self.keys()}

    def __repr__(self):
        return repr(self.getall())

    def __str__(self):
        return str(self.getall())


class SChunk(blosc2_ext.SChunk):
    def __init__(  # noqa: C901
        self,
        chunksize: int | None = None,
        data: object = None,
        **kwargs: dict | blosc2.CParams | blosc2.Storage | blosc2.DParams,
    ) -> None:
        """Create a new super-chunk, or open an existing one.

        Parameters
        ----------
        chunksize: int, optional
            The size, in bytes, of the chunks in the super-chunk. If not provided,
            it is set automatically to a reasonable value.

        data: bytes-like object, optional
            The data to be split into different chunks of size :paramref:`chunksize`.
            If None, the Schunk instance will be empty initially.

        kwargs: dict, optional
            Storage parameters. The default values are in :class:`blosc2.Storage`.
            Supported keyword arguments:
                storage: :class:`blosc2.Storage` or dict
                    All the storage parameters that you want to use as
                    a :class:`blosc2.Storage` or dict instance.
                cparams: :class:`blosc2.CParams` or dict
                    All the compression parameters that you want to use as
                    a :class:`blosc2.CParams` or dict instance.
                dparams: :class:`blosc2.DParams` or dict
                    All the decompression parameters that you want to use as
                    a :class:`blosc2.DParams` or dict instance.
                others: Any
                    If `storage` is not passed, all the parameters of a :class:`blosc2.Storage`
                    can be passed as keyword arguments.

        Examples
        --------
        >>> import blosc2
        >>> import numpy as np
        >>> import os.path
        >>> import shutil
        >>> import tempfile
        >>> cparams = blosc2.CParams()
        >>> dparams = blosc2.DParams()
        >>> storage = blosc2.Storage(contiguous=True)
        >>> schunk = blosc2.SChunk(cparams=cparams, dparams=dparams, storage=storage)

        In the following, we will write and read a super-chunk to and from disk
        via memory-mapped files.

        >>> a = np.arange(3, dtype=np.int64)
        >>> chunksize = a.size * a.itemsize
        >>> n_chunks = 2
        >>> tmpdirname = tempfile.mkdtemp()
        >>> urlpath = os.path.join(tmpdirname, 'schunk.b2frame')

        Optional: we intend to write 2 chunks of 24 bytes each, and we expect
        the compressed size to be smaller than the original size. Therefore, we
        generously set the initial size of the mapping to 48 bytes
        effectively avoiding remappings.

        >>> initial_mapping_size = chunksize * n_chunks
        >>> schunk_mmap = blosc2.SChunk(
        ...     chunksize=chunksize,
        ...     mmap_mode="w+",
        ...     initial_mapping_size=initial_mapping_size,
        ...     urlpath=urlpath,
        ... )
        >>> schunk_mmap.append_data(a)
        1
        >>> schunk_mmap.append_data(a * 2)
        2

        Optional: explicitly close the file and free the mapping.

        >>> del schunk_mmap

        Reading the data back again via memory-mapped files:

        >>> schunk_mmap = blosc2.open(urlpath, mmap_mode="r")
        >>> np.frombuffer(schunk_mmap.decompress_chunk(0), dtype=np.int64).tolist()
        [0, 1, 2]
        >>> np.frombuffer(schunk_mmap.decompress_chunk(1), dtype=np.int64).tolist()
        [0, 2, 4]
        >>> shutil.rmtree(tmpdirname)
        """
        # Check only allowed kwarg are passed
        allowed_kwargs = [
            "urlpath",
            "contiguous",
            "cparams",
            "dparams",
            "_schunk",
            "meta",
            "mode",
            "mmap_mode",
            "initial_mapping_size",
            "_is_view",
            "storage",
        ]
        for kwarg in kwargs:
            if kwarg not in allowed_kwargs:
                raise ValueError(f"{kwarg} is not supported as keyword argument")
        if kwargs.get("storage") is not None:
            if any(key in list(blosc2.Storage.__annotations__) for key in kwargs):
                raise AttributeError(
                    "Cannot pass both `storage` and other kwargs already included in Storage"
                )
            storage = kwargs.get("storage")
            if isinstance(storage, blosc2.Storage):
                kwargs = {**kwargs, **asdict(storage)}
            else:
                kwargs = {**kwargs, **storage}

        if isinstance(kwargs.get("cparams"), blosc2.CParams):
            kwargs["cparams"] = asdict(kwargs.get("cparams"))

        if isinstance(kwargs.get("dparams"), blosc2.DParams):
            kwargs["dparams"] = asdict(kwargs.get("dparams"))

        urlpath = kwargs.get("urlpath")
        if "contiguous" not in kwargs:
            # Make contiguous true for disk, else sparse (for in-memory performance)
            kwargs["contiguous"] = urlpath is not None

        # This a private param to get an SChunk from a blosc2_schunk*
        sc = kwargs.pop("_schunk", None)

        # If not passed, set a sensible typesize
        itemsize = data.itemsize if data is not None and hasattr(data, "itemsize") else 1
        if "cparams" in kwargs:
            if "typesize" not in kwargs["cparams"]:
                cparams = kwargs.pop("cparams").copy()
                cparams["typesize"] = itemsize
                kwargs["cparams"] = cparams
        else:
            kwargs["cparams"] = {"typesize": itemsize}

        # chunksize handling
        if chunksize is None:
            chunksize = 2**24
            if data is not None:
                if hasattr(data, "itemsize"):
                    chunksize = data.size * data.itemsize
                    # Make that a multiple of typesize
                    chunksize = chunksize // data.itemsize * data.itemsize
                else:
                    chunksize = len(data)
            # Use a cap of 256 MB (modern boxes should all have this RAM available)
            if chunksize > 2**28:
                chunksize = 2**28

        super().__init__(_schunk=sc, chunksize=chunksize, data=data, **kwargs)
        self._vlmeta = vlmeta(
            super().c_schunk, self.urlpath, self.mode, self.mmap_mode, self.initial_mapping_size
        )
        self._cparams = super().get_cparams()
        self._dparams = super().get_dparams()

    @property
    def cparams(self) -> blosc2.CParams:
        """
        :class:`blosc2.CParams` instance with the compression parameters.
        """
        return self._cparams

    @cparams.setter
    def cparams(self, value: blosc2.CParams) -> None:
        super().update_cparams(value)
        self._cparams = super().get_cparams()

    @property
    def dparams(self) -> blosc2.DParams:
        """
        :class:`blosc2.DParams` instance with the decompression parameters.
        """
        return self._dparams

    @dparams.setter
    def dparams(self, value: blosc2.DParams) -> None:
        super().update_dparams(value)
        self._dparams = super().get_dparams()

    @property
    def meta(self) -> Meta:
        """
        Access to the fixed-length metadata of the `SChunk`.
        """
        return Meta(self)

    @property
    def vlmeta(self) -> vlmeta:
        """
        Access to the variable-length metadata of the `SChunk`.
        """
        return self._vlmeta

    @property
    def chunkshape(self) -> int:
        """
        Number of elements per chunk.
        """
        return self.chunksize // self.typesize

    @property
    def chunksize(self) -> int:
        """
        Number of bytes in each chunk.
        """
        return super().chunksize

    @property
    def blocksize(self) -> int:
        """The block size (in bytes)."""
        return super().blocksize

    @property
    def nchunks(self) -> int:
        """The number of chunks."""
        return super().nchunks

    @property
    def cratio(self) -> float:
        """
        Compression ratio.
        """
        if self.cbytes == 0:
            return 0.0
        return self.nbytes / self.cbytes

    @property
    def nbytes(self) -> int:
        """
        Amount of uncompressed data bytes.
        """
        return super().nbytes

    @property
    def cbytes(self) -> int:
        """
        Amount of compressed data bytes (data size + chunk headers size).
        """
        return super().cbytes

    @property
    def typesize(self) -> int:
        """
        Type size of the `SChunk`.
        """
        return super().typesize

    @property
    def urlpath(self) -> str:
        """
        Path where the `SChunk` is stored.
        """
        return super().urlpath

    @property
    def contiguous(self) -> bool:
        """
        Whether the `SChunk` is stored contiguously or sparsely.
        """
        return super().contiguous

    def append_data(self, data: object) -> int:
        """Append a data buffer to the SChunk.

        The data buffer must be of size `chunksize` specified in
        :func:`SChunk.__init__ <blosc2.schunk.SChunk.__init__>`.

        Parameters
        ----------
        data: bytes-like object
            The data to be compressed and added as a chunk.

        Returns
        -------
        out: int
            The number of chunks in the SChunk.

        Raises
        ------
        RunTimeError
            If the :paramref:`data` could not be appended.

        Examples
        --------
        >>> import blosc2
        >>> import numpy as np
        >>> schunk = blosc2.SChunk(chunksize=200*1000*4)
        >>> data = np.arange(200 * 1000, dtype='int32')
        >>> schunk.append_data(data)
        1
        """
        blosc2_ext.check_access_mode(self.urlpath, self.mode)
        return super().append_data(data)

    def fill_special(
        self,
        nitems: int,
        special_value: blosc2.SpecialValue,
        value: bytes | int | float | bool | None = None,
    ) -> int:
        """Fill the SChunk with a special value. The SChunk must be empty.

        Parameters
        ----------
        nitems: int
            The number of items to fill with the special value.
        special_value: SpecialValue
            The special value to be used for filling the SChunk.
        value: bytes, int, float, bool (optional)
            The value to fill the SChunk. This parameter is only supported if
            :paramref:`special_value` is ``blosc2.SpecialValue.VALUE``.

        Returns
        -------
        out: int
            The number of chunks in the SChunk.

        Raises
        ------
        RunTimeError
            If the SChunk could not be filled with the special value.

        Examples
        --------
        >>> import blosc2
        >>> import numpy as np
        >>> import time
        >>> nitems = 100_000_000
        >>> dtype = np.dtype(np.float64)
        >>> # Measure the time to create SChunk from a NumPy array
        >>> t0 = time.time()
        >>> data = np.full(nitems, np.pi, dtype)
        >>> cparams = blosc2.CParams(typesize=dtype.itemsize)
        >>> schunk = blosc2.SChunk(data=data, cparams=cparams)
        >>> t = (time.time() - t0) * 1000.
        >>> f"Time creating a schunk with a numpy array: {t:10.3f} ms"
        Time creating a schunk with a numpy array:    710.273 ms
        >>> # Measure the time to create SChunk using fill_special
        >>> t0 = time.time()
        >>> cparams = blosc2.CParams(typesize=dtype.itemsize)
        >>> schunk = blosc2.SChunk(cparams=cparams)
        >>> schunk.fill_special(nitems, blosc2.SpecialValue.VALUE, np.pi)
        >>> t = (time.time() - t0) * 1000.
        >>> f"Time passing directly the value to `fill_special`: {t:10.3f} ms"
        Time passing directly the value to `fill_special`:      2.109 ms
        """
        if not isinstance(special_value, SpecialValue) or special_value == SpecialValue.NOT_SPECIAL:
            raise TypeError("special_value must be a SpecialValue instance other than NOT_SPECIAL")
        if special_value == SpecialValue.VALUE and value is None:
            raise ValueError("value cannot be None when special_value is VALUE")

        nchunks = super().fill_special(nitems, special_value.value, value)
        if nchunks < 0:
            raise RuntimeError("Unable to fill with special values")
        return nchunks

    def decompress_chunk(self, nchunk: int, dst: object = None) -> str | bytes:
        """Decompress the chunk given by its index :paramref:`nchunk`.

        Parameters
        ----------
        nchunk: int
            The index of the chunk that will be decompressed.
        dst: NumPy object or bytearray
            The destination NumPy object or bytearray to fill, the length
            of which must be greater than 0. The user must ensure
            that it has enough capacity to host the decompressed
            chunk. Default is None, meaning that a new bytes object
            is created, filled and returned.

        Returns
        -------
        out: str or bytes
            The decompressed chunk as a Python str or bytes object if
            :paramref:`dst` is `None`. Otherwise, it returns `None` because the
            result will already be in :paramref:`dst`.

        Raises
        ------
        RunTimeError
            If a problem is detected.

        Examples
        --------
        >>> import blosc2
        >>> cparams = blosc2.CParams(typesize=1)
        >>> schunk = blosc2.SChunk(cparams=cparams)
        >>> buffer = b"wermqeoir23"
        >>> schunk.append_data(buffer)
        1
        >>> schunk.decompress_chunk(0)
        b'wermqeoir23'
        >>> # Construct a mutable bytearray object
        >>> bytes_obj = bytearray(len(buffer))
        >>> schunk.decompress_chunk(0, dst=bytes_obj)
        >>> bytes_obj == buffer
        True
        """
        return super().decompress_chunk(nchunk, dst)

    def get_chunk(self, nchunk: int) -> bytes:
        """Return the compressed chunk that is in the SChunk.

        Parameters
        ----------
        nchunk: int
            The index of the chunk that will be returned.

        Returns
        -------
        out: bytes object
            The compressed chunk.

        Raises
        ------
        RunTimeError
            If a problem is detected.

        Examples
        --------
        >>> import blosc2
        >>> import numpy as np
        >>> # Create an SChunk with 3 chunks
        >>> nchunks = 3
        >>> data = np.arange(200 * 1000 * nchunks, dtype=np.int32)
        >>> cparams = blosc2.CParams(typesize=4)
        >>> schunk = blosc2.SChunk(data=data, cparams=cparams)
        >>> # Retrieve the first chunk (index 0)
        >>> chunk = schunk.get_chunk(0)
        >>> # Check the type and length of the compressed chunk
        >>> type(chunk)
        <class 'bytes'>
        >>> len(chunk)
        10552
        """
        return super().get_chunk(nchunk)

    def delete_chunk(self, nchunk: int) -> int:
        """Delete the specified chunk from the SChunk.

        Parameters
        ----------
        nchunk: int
            The index of the chunk that will be removed.

        Returns
        -------
        out: int
            The number of chunks in the SChunk.

        Raises
        ------
        RunTimeError
            If a problem is detected.

        Examples
        --------
        >>> import blosc2
        >>> import numpy as np
        >>> # Create an SChunk with 3 chunks
        >>> nchunks = 3
        >>> data = np.arange(200 * 1000 * nchunks, dtype=np.int32)
        >>> cparams = blosc2.CParams(typesize=4)
        >>> schunk = blosc2.SChunk(chunksize=200 * 1000 * 4, data=data, cparams=cparams)
        >>> # Check the number of chunks before deletion
        >>> schunk.nchunks
        3
        >>>  # Delete the second chunk (index 1)
        >>> schunk.delete_chunk(1)
        >>>  # Check the number of chunks after deletion
        >>> schunk.nchunks
        2
        """
        blosc2_ext.check_access_mode(self.urlpath, self.mode)
        return super().delete_chunk(nchunk)

    def insert_chunk(self, nchunk: int, chunk: bytes) -> int:
        """Insert an already compressed chunk into the SChunk.

        Parameters
        ----------
        nchunk: int
            The index at which the chunk will be inserted.
        chunk: bytes object
            The compressed chunk.

        Returns
        -------
        out: int
            The number of chunks in the SChunk.

        Raises
        ------
        RunTimeError
            If a problem is detected.

        Examples
        --------
        >>> import blosc2
        >>> import numpy as np
        >>> # Create an SChunk with 2 chunks
        >>> data = np.arange(400 * 1000, dtype=np.int32)
        >>> cparams = blosc2.CParams(typesize=4)
        >>> schunk = blosc2.SChunk(chunksize=200*1000*4, data=data, cparams=cparams)
        >>> # Get a compressed chunk from the SChunk
        >>> chunk = schunk.get_chunk(0)
        >>> # Insert a chunk in the second position (index 1)"
        >>> schunk.insert_chunk(1, chunk)
        >>> # Verify the total number of chunks after insertion
        >>> schunk.nchunks
        3
        """
        blosc2_ext.check_access_mode(self.urlpath, self.mode)
        return super().insert_chunk(nchunk, chunk)

    def insert_data(self, nchunk: int, data: object, copy: bool) -> int:
        """Insert the data in the specified position in the SChunk.

        Parameters
        ----------
        nchunk: int
            The index at which the chunk will be inserted.
        data: bytes object
            The data that will be compressed and inserted as a chunk.
        copy: bool
            Whether to make an internal copy of the chunk to insert it or not.

        Returns
        -------
        out: int
            The number of chunks in the SChunk.

        Raises
        ------
        RunTimeError
            If a problem is detected.

        Examples
        --------
        >>> import blosc2
        >>> import numpy as np
        >>> # Create an SChunk with 2 chunks
        >>> data = np.arange(400 * 1000, dtype=np.int32)
        >>> cparams = blosc2.CParams(typesize=4)
        >>> schunk = blosc2.SChunk(chunksize=200*1000*4, data=data, cparams=cparams)
        >>> # Create a new array to insert into the second chunk of the SChunk
        >>> new_data = np.arange(200 * 1000, dtype=np.int32)
        >>> # Insert the new data at position 1, compressing it
        >>> schunk.insert_data(1, new_data, copy=True)
        >>> # Verify the total number of chunks after insertion
        >>> schunk.nchunks
        3
        """
        blosc2_ext.check_access_mode(self.urlpath, self.mode)
        return super().insert_data(nchunk, data, copy)

    def update_chunk(self, nchunk: int, chunk: bytes) -> int:
        """Update an existing chunk in the SChunk.

        Parameters
        ----------
        nchunk: int
            The index of the chunk to be updated.
        chunk: bytes object
            The new compressed chunk that will replace the old chunk's content.

        Returns
        -------
        out: int
            The number of chunks in the SChunk.

        Raises
        ------
        RunTimeError
            If a problem is detected.

        Examples
        --------
        >>> import blosc2
        >>> import numpy as np
        >>> nchunks = 5
        >>> chunk_size = 200 * 1000 * 4
        >>> data = np.arange(nchunks * chunk_size // 4, dtype=np.int32)
        >>> cparams = blosc2.CParams(typesize=4)
        >>> schunk = blosc2.SChunk(chunksize=chunk_size, data=data, cparams=cparams)
        >>> f"Initial number of chunks: {schunk.nchunks}"
        Initial number of chunks: 5
        >>> c_index = 1
        >>> new_data = np.full(chunk_size // 4, fill_value=c_index, dtype=np.int32).tobytes()
        >>> compressed_data = blosc2.compress2(new_data, typesize=4)
        >>> # Update the 2nd chunk (index 1) with new data
        >>> nchunks = schunk.update_chunk(c_index, compressed_data)
        >>> f"Number of chunks after update: {nchunks}"
        Number of chunks after update: 5
        """
        blosc2_ext.check_access_mode(self.urlpath, self.mode)
        return super().update_chunk(nchunk, chunk)

    def update_data(self, nchunk: int, data: object, copy: bool) -> int:
        """Update the chunk in the specified position with the given data.

        Parameters
        ----------
        nchunk: int
            The index of the chunk to be updated.
        data: bytes object
            The data to be compressed and will replace the old chunk.
        copy: bool
             Whether to make an internal copy of the chunk before updating it.

        Returns
        -------
        out: int
            The number of chunks in the SChunk.

        Raises
        ------
        RunTimeError
            If a problem is detected.

        Examples
        --------
        >>> import blosc2
        >>> import numpy as np
        >>> nchunks = 4
        >>> chunk_size = 200 * 1000 * 4
        >>> data = np.arange(nchunks * chunk_size // 4, dtype=np.int32)
        >>> cparams = blosc2.CParams(typesize=4)
        >>> schunk = blosc2.SChunk(chunksize=chunk_size, data=data, cparams=cparams)
        >>> f"Initial number of chunks: {schunk.nchunks}"
        Initial number of chunks: 4
        >>> c_index = 1 # Update the 2nd chunk (index 1)
        >>> new_data = np.full(chunk_size // 4, fill_value=c_index, dtype=np.int32).tobytes()
        >>> nchunks = schunk.update_data(c_index, new_data, copy=True)
        >>> f"Number of chunks after update: {schunk.nchunks}"
        Number of chunks after update: 4
        """
        blosc2_ext.check_access_mode(self.urlpath, self.mode)
        return super().update_data(nchunk, data, copy)

    def get_slice(self, start: int = 0, stop: int | None = None, out: object = None) -> str | bytes | None:
        """Get a slice from :paramref:`start` to :paramref:`stop`.

        Parameters
        ----------
        start: int
            The starting index of the slice. Default is 0.
        stop: int
            The ending index of the slice (exclusive).
            Default is until the SChunk ends.
        out: bytes-like object or bytearray
            The target object (supporting the
            `Buffer Protocol <https://docs.python.org/3/c-api/buffer.html>`_) to fill.
            Verify that the buffer has enough space for the decompressed data.
            If `None` is provided, a new bytes object will be created, filled,
            and returned.

        Returns
        -------
        out: str or bytes or None
            The decompressed slice a Python str or bytes object if
            :paramref:`out` is `None`. Otherwise, it returns `None` since the result
            will already be in :paramref:`out`.

        Raises
        ------
        ValueError
            If the size to get is negative.
            If there is not enough space in :paramref:`out`.
            If :paramref:`start` is greater or equal to the number of items in the SChunk.
        RunTimeError
            If a problem is detected.

        See Also
        --------
        :func:`__getitem__`

        Examples
        --------
        >>> import blosc2
        >>> import numpy as np
        >>> nchunks = 4
        >>> chunk_size = 200 * 1000 * 4
        >>> data = np.arange(nchunks * chunk_size // 4, dtype=np.int32)
        >>> cparams = blosc2.CParams(typesize=4)
        >>> schunk = blosc2.SChunk(data=data, cparams=cparams)
        >>> # Define the slice parameters
        >>> start_index = 200 * 1000
        >>> stop_index = 2 * 200 * 1000
        >>> # Prepare an output buffer
        >>> slice_size = stop_index - start_index
        >>> out_buffer = bytearray(slice_size * 4)  # Ensure the buffer is large enough
        >>> result = schunk.get_slice(start=start_index, stop=stop_index, out=out_buffer)
        >>> # Convert bytearray to NumPy array for easier inspection
        >>> slice_array = np.frombuffer(out_buffer, dtype=np.int32)
        >>> f"Slice data: {slice_array[:10]} ..."  # Print the first 10 elements
        Slice data: [200000 200001 200002 200003 200004 200005 200006 200007 200008 200009] ...
        """
        return super().get_slice(start, stop, out)

    def __len__(self) -> int:
        """
        Return the number of items in the SChunk.
        """
        return self.nbytes // self.typesize

    def __getitem__(self, item: int | slice) -> str | bytes:
        """Get a slice from the SChunk.

        Parameters
        ----------
        item: int or slice
            The index or slice for the data. Note that the step parameter is not honored.

        Returns
        -------
        out: str or bytes
            The decompressed slice as a Python str or bytes object.

        Raises
        ------
        ValueError
            If the size to get is negative.
            If :paramref:`item`.start is greater than or equal to the number of
            items in the SChunk.
        RunTimeError
            If a problem is detected.
        IndexError
            If `step` is not 1.

        See Also
        --------
        :func:`get_slice`

        Examples
        --------
        >>> import blosc2
        >>> import numpy as np
        >>> nchunks = 4
        >>> chunk_size = 200 * 1000 * 4
        >>> data = np.arange(nchunks * chunk_size // 4, dtype=np.int32)
        >>> cparams = blosc2.CParams(typesize=4)
        >>> schunk = blosc2.SChunk(chunksize=chunk_size, data=data, cparams=cparams)
        >>> # Use __getitem__ to retrieve the same slice of data from the SChunk
        >>> res = schunk[150:155]
        >>> f"Slice data: {np.frombuffer(res, dtype=np.int32)}"
        Slice data: [150 151 152 153 154]
        """
        if isinstance(item, int):
            if item == -1:
                return self.get_slice(item)
            return self.get_slice(item, item + 1)
        if item.step is not None and item.step != 1:
            raise IndexError("`step` must be 1")
        return self.get_slice(item.start, item.stop)

    def __setitem__(self, key: int | slice, value: object) -> None:
        """Set slice to :paramref:`value`.

        Parameters
        ----------
        key: int or slice
            The index of the slice to update. Note that step parameter is not honored.
        value: bytes-like object
            An object supporting the
            `Buffer Protocol <https://docs.python.org/3/c-api/buffer.html>`_ used to
            fill the slice.

        Returns
        -------
        out: None

        Raises
        ------
        ValueError
            If the object cannot be modified.
            If the size to get is negative.
            If there is not enough space in :paramref:`value` to update the slice.
            If :paramref:`start` is greater than the number of items in the SChunk.
        RunTimeError
            If a problem is detected.
        IndexError
            If `step` is not 1.

        Notes
        -----
        This method can also be used to append new data if :paramref:`key`.stop
        is greater than the number of items in the SChunk.

        Examples
        --------
        >>> import blosc2
        >>> import numpy as np
        >>> nchunks = 4
        >>> chunk_size = 200 * 1000 * 4
        >>> data = np.arange(nchunks * chunk_size // 4, dtype=np.int32)
        >>> cparams = blosc2.CParams(typesize=4)
        >>> schunk = blosc2.SChunk(data=data, cparams=cparams)
        >>> # Create a new array of values to update the slice (values from 1000 to 1999 multiplied by 2)
        >>> start_ = 1000
        >>> stop = 2000
        >>> new_values = np.arange(start_, stop, dtype=np.int32) * 2
        >>> schunk[start_:stop] = new_values
        >>> # Retrieve the updated slice using the slicing syntax
        >>> retrieved_slice = np.frombuffer(schunk[start_:stop], dtype=np.int32)
        >>> f"First 10 values of the updated slice: {retrieved_slice[:10]}"
        >>> f"Last 10 values of the updated slice: {retrieved_slice[-10:]}"
        First 10 values of the updated slice: [2000 2002 2004 2006 2008 2010 2012 2014 2016 2018]
        Last 10 values of the updated slice: [3980 3982 3984 3986 3988 3990 3992 3994 3996 3998]
        """
        if key.step is not None and key.step != 1:
            raise IndexError("`step` must be 1")
        blosc2_ext.check_access_mode(self.urlpath, self.mode)
        return super().set_slice(start=key.start, stop=key.stop, value=value)

    def to_cframe(self) -> bytes:
        """Get a bytes object containing the serialized :ref:`SChunk` instance.

        Returns
        -------
        out: bytes
            The buffer containing the serialized :ref:`SChunk` instance.

        See Also
        --------
        :func:`~blosc2.schunk_from_cframe`

        Examples
        --------
        >>> import blosc2
        >>> import numpy as np
        >>> nchunks = 4
        >>> chunk_size = 200 * 1000 * 4
        >>> data = np.arange(nchunks * chunk_size // 4, dtype=np.int32)
        >>> cparams = blosc2.CParams(typesize=4)
        >>> schunk = blosc2.SChunk(data=data, cparams=cparams)
        >>> # Serialize the SChunk instance to a bytes object
        >>> serialized_schunk = schunk.to_cframe()
        >>> f"Serialized SChunk length: {len(serialized_schunk)} bytes"
        Serialized SChunk length: 14129 bytes
        >>> # Create a new SChunk from the serialized data
        >>> deserialized_schunk = blosc2.schunk_from_cframe(serialized_schunk)
        >>> start = 500
        >>> stop = 505
        >>> sl_bytes = deserialized_schunk[start:stop]
        >>> sl = np.frombuffer(sl_bytes, dtype=np.int32)
        >>> res = data[start:stop]
        >>> f"Original slice: {res}"
        Original slice: [500 501 502 503 504]
        >>> f"Deserialized slice: {sl}"
        Deserialized slice: [500 501 502 503 504]
        """
        return super().to_cframe()

    def iterchunks(self, dtype: np.dtype) -> Iterator[np.ndarray]:
        """
        Iterate over the :paramref:`self` chunks of the SChunk.

        Parameters
        ----------
        dtype: np.dtype
            The data type to use for the decompressed chunks.

        Yields
        ------
        chunk: NumPy ndarray
           The decompressed chunk.

        Examples
        --------
        >>> import blosc2
        >>> import numpy as np
        >>> # Create sample data and an SChunk
        >>> data = np.arange(400 * 1000, dtype=np.int32)
        >>> cparams = blosc2.CParams(typesize=4)
        >>> schunk = blosc2.SChunk(data=data, cparams=cparams)
        >>> # Iterate over chunks using the iterchunks method
        >>> for chunk in schunk.iterchunks(dtype=np.int32):
        >>>     f"Chunk shape: {chunk.shape} "
        >>>     f"First 5 elements of chunk: {chunk[:5]}"
        Chunk shape: (400000,)
        First 5 elements of chunk: [0 1 2 3 4]
        """
        out = np.empty(self.chunkshape, dtype)
        for i in range(0, len(self), self.chunkshape):
            self.get_slice(i, i + self.chunkshape, out)
            yield out

    def iterchunks_info(
        self,
    ) -> Iterator[
        NamedTuple(
            "info",
            nchunk=int,
            cratio=float,
            special=blosc2.SpecialValue,
            repeated_value=bytes | None,
            lazychunk=bytes,
        )
    ]:
        """
        Iterate over the chunks of the SChunk, providing info on index and special values.

        Yields
        ------
        info: namedtuple
            A namedtuple with the following fields:

                nchunk: int
                    The index of the chunk.
                cratio: float
                    The compression ratio of the chunk.
                special: :class:`~blosc2.SpecialValue`
                    The special value enum of the chunk; if 0, the chunk is not special.
                repeated_value: bytes or None
                    The repeated value for the chunk; if not SpecialValue.VALUE, it is None.
                lazychunk: bytes
                    A buffer with the complete lazy chunk.

        Examples
        --------
        >>> import blosc2
        >>> import numpy as np
        >>> # Create sample data and an SChunk
        >>> data = np.arange(400 * 1000, dtype=np.int32)
        >>> cparams = blosc2.CParams(typesize=4)
        >>> schunk = blosc2.SChunk(data=data, cparams=cparams)
        >>> # Iterate over chunks and print detailed information
        >>> for chunk_info in schunk.iterchunks_info():
        >>>     f"Chunk index: {chunk_info.nchunk}"
        >>>     f"Compression ratio: {chunk_info.cratio:.2f}"
        >>>     f"Special value: {chunk_info.special.name}"
        >>>     f"Repeated value: {chunk_info.repeated_value[:10] if chunk_info.repeated_value else None}"
        Chunk index: 0
        Compression ratio: 223.56
        Special value: NOT_SPECIAL
        Repeated value: None
        """
        ChunkInfo = namedtuple("ChunkInfo", ["nchunk", "cratio", "special", "repeated_value", "lazychunk"])
        for nchunk in range(self.nchunks):
            lazychunk = self.get_lazychunk(nchunk)
            # Blosc2 flags are encoded at the end of the header
            # (see https://github.com/Blosc/c-blosc2/blob/main/README_CHUNK_FORMAT.rst)
            is_special = (lazychunk[31] & 0x70) >> 4
            special = SpecialValue(is_special)
            # The special value is encoded at the end of the header
            repeated_value = lazychunk[32:] if special == SpecialValue.VALUE else None
            # Compression ratio (nbytes and cbytes are little-endian)
            cratio = (
                np.frombuffer(lazychunk[4:8], dtype="<i4")[0]
                / np.frombuffer(lazychunk[12:16], dtype="<i4")[0]
            )
            yield ChunkInfo(nchunk, cratio, special, repeated_value, lazychunk)

    def postfilter(self, input_dtype: np.dtype, output_dtype: np.dtype = None) -> None:
        """Decorator to set a function as a postfilter.

        The postfilter function will be executed each time after decompressing
        blocks of data. It will receive three parameters:

        * the input `ndarray` to be read from
        * the output `ndarray` to be filled out
        * the offset inside the `SChunk` instance where the corresponding block begins (see example below).

        Parameters
        ----------
        input_dtype: np.dtype
            Data type of the input that will receive the postfilter function.
        output_dtype: np.dtype
            Data type of the output that will receive and fill the postfilter function.
            If None (default) it will be set to :paramref:`input_dtype`.

        Returns
        -------
        out: None

        Notes
        -----
        * `nthreads` must be 1 when decompressing.

        * The :paramref:`input_dtype` itemsize must be the same as the
          :paramref:`output_dtype` itemsize.

        See Also
        --------
        :meth:`remove_postfilter`
        :meth:`prefilter`

        Examples
        --------
        .. code-block:: python

            # Create SChunk
            input_dtype = np.dtype(np.int64)
            cparams = blosc2.CParams(typesize=input_dtype.itemsize)
            dparams = blosc2.DParams(nthreads=1)
            schunk = blosc2.SChunk(
                chunksize=20_000 * input_dtype.itemsize, cparams=cparams, dparams=dparams
            )


            # Create postfilter and associate it to the schunk
            @schunk.postfilter(input_dtype)
            def postfilter(input, output, offset):
                output[:] = offset + np.arange(input.size)
        """

        def initialize(func):
            super(SChunk, self)._set_postfilter(func, input_dtype, output_dtype)

            def exec_func(*args):
                func(*args)

            return exec_func

        return initialize

    def remove_postfilter(self, func_name: str, _new_ctx: bool = True) -> None:
        """Remove the postfilter from the `SChunk` instance.

        Parameters
        ----------
        func_name: str
            The name of the postfilter function to remove.

        Returns
        -------
        out: None

        Examples
        --------
        >>> import blosc2
        >>> import numpy as np
        >>> dtype = np.dtype(np.int32)
        >>> cparams = blosc2.CParams(typesize=dtype.itemsize)
        >>> dparams = blosc2.DParams(nthreads=1)
        >>> data = np.arange(500, dtype=np.int32)
        >>> schunk = blosc2.SChunk(data=data, cparams=cparams, dparams=dparams)
        >>> # Define the postfilter function
        >>> @schunk.postfilter(dtype)
        >>> def postfilter(input, output, offset):
        >>>     output[:] = input + offset + np.arange(input.size)
        >>> out = np.empty(data.size, dtype=dtype)
        >>> schunk.get_slice(out=out)
        >>> f"Data slice with postfilter applied (first 8 elements): {out[:8]}"
        Data slice with postfilter applied (first 8 elements): [ 0  2  4  6  8 10 12 14]
        >>> schunk.remove_postfilter('postfilter')
        >>> retrieved_data = np.empty(data.size, dtype=dtype)
        >>> schunk.get_slice(out=retrieved_data)
        >>> f"Original data (first 8 elements): {data[:8]}"
        Original data (first 8 elements): [0 1 2 3 4 5 6 7]
        """
        return super().remove_postfilter(func_name)

    def filler(self, inputs_tuple: tuple[tuple], schunk_dtype: np.dtype, nelem: int | None = None) -> None:
        """Decorator to set a filler function.

        This function will fill :paramref:`self` according to :paramref:`nelem`.
        It will receive three parameters: a tuple with the inputs as `ndarrays`
        from which to read, the `ndarray` to fill :paramref:`self` and the
        offset inside the `SChunk` instance where the corresponding block
        begins (see example below).

        Parameters
        ----------
        inputs_tuple: tuple of tuples
            Tuple containing a tuple for each argument that the function will receive, along with their
            corresponding np.dtype.
            Supported operand types are :ref:`SChunk`, `ndarray` and
            Python scalars.
        schunk_dtype: np.dtype
            The data type to use to fill :paramref:`self`.
        nelem: int
            Number of elements to append to :paramref:`self`. If None (default) it
            will be the number of elements from the operands.

        Returns
        -------
        out: None

        Notes
        -----
        * Compression `nthreads` must be 1 when using this.
        * This does not need to be removed from the created `SChunk` instance.

        See Also
        --------
        :meth:`prefilter`

        Examples
        --------
        .. code-block:: python

            # Set the compression and decompression parameters
            schunk_dtype = np.dtype(np.float64)
            cparams = blosc2.CParams(typesize=schunk_dtype.itemsize, nthreads=1)
            # Create empty SChunk
            schunk = blosc2.SChunk(chunksize=20_000 * schunk_dtype.itemsize, cparams=cparams)

            # Create operands
            op_dtype = np.dtype(np.int32)
            data = np.full(20_000 * 3, 12, dtype=op_dtype)
            schunk_op = blosc2.SChunk(chunksize=20_000 * op_dtype.itemsize, data=data)


            # Create filler
            @schunk.filler(((schunk_op, op_dtype), (np.e, np.float32)), schunk_dtype)
            def filler(inputs_tuple, output, offset):
                output[:] = inputs_tuple[0] - inputs_tuple[1]

        """

        def initialize(func):
            if self.nbytes != 0:
                raise ValueError("Cannot apply a filler to a non empty SChunk")
            nelem_ = blosc2_ext.nelem_from_inputs(inputs_tuple, nelem)
            super(SChunk, self)._set_filler(func, id(inputs_tuple), schunk_dtype)
            chunksize = self.chunksize
            written_nbytes = 0
            nbytes = nelem_ * self.typesize
            while written_nbytes < nbytes:
                chunk = np.zeros(chunksize // self.typesize, dtype=schunk_dtype)
                self.append_data(chunk)
                written_nbytes += chunksize
                if (nbytes - written_nbytes) < self.chunksize:
                    chunksize = nbytes - written_nbytes
            self.remove_prefilter(func.__name__)

            def exec_func(*args):
                func(*args)

            return exec_func

        return initialize

    def prefilter(self, input_dtype: np.dtype, output_dtype: np.dtype = None) -> None:
        """Decorator to set a function as a prefilter.

        This function will be executed each time before compressing the data.
        It will receive three parameters:

        * The actual data as a `ndarray` from which to read,
        * The `ndarray` to be filled,
        * The offset inside the `SChunk` instance where the corresponding block begins (see example below).

        Parameters
        ----------
        input_dtype: np.dtype
            Data type of the input that will be processed the prefilter function.
        output_dtype: np.dtype, optional
            Data type of the output that will be filled by the prefilter function.
            If None (default), it will be the same as :paramref:`input_dtype`.

        Returns
        -------
        out: None

        Notes
        -----
        * `nthreads` must be 1 when compressing.

        * The :paramref:`input_dtype` itemsize must be the same as the
          :paramref:`output_dtype` itemsize.

        See Also
        --------
        :meth:`remove_prefilter`
        :meth:`postfilter`
        :meth:`filler`

        Examples
        --------
        .. code-block:: python

            # Set the compression and decompression parameters
            input_dtype = np.dtype(np.int32)
            output_dtype = np.dtype(np.float32)
            cparams = blosc2.CParams(typesize=output_dtype.itemsize, nthreads=1)
            # Create schunk
            schunk = blosc2.SChunk(chunksize=200 * 1000 * input_dtype.itemsize, cparams=cparams)


            # Set prefilter with decorator
            @schunk.prefilter(input_dtype, output_dtype)
            def prefilter(input, output, offset):
                output[:] = input - np.pi
        """

        def initialize(func):
            super(SChunk, self)._set_prefilter(func, input_dtype, output_dtype)

            def exec_func(*args):
                func(*args)

            return exec_func

        return initialize

    def remove_prefilter(self, func_name: str, _new_ctx: bool = True) -> None:
        """Remove the prefilter from the `SChunk` instance.

        Parameters
        ----------
        func_name: str
            Name of the prefilter function.

        Returns
        -------
        out: None

        Examples
        --------
        >>> import blosc2
        >>> import numpy as np
        >>> dtype = np.dtype(np.int32)
        >>> cparams = blosc2.CParams(typesize=dtype.itemsize, nthreads=1)
        >>> data = np.arange(1000, dtype=np.int32)
        >>> output_dtype = np.float32
        >>> schunk = blosc2.SChunk(cparams=cparams)
        >>> # Define the prefilter function
        >>> @schunk.prefilter(dtype, output_dtype)
        >>> def prefilter(input, output, offset):
        >>>     output[:] = input - np.pi
        >>> schunk[:1000] = data
        >>> # Retrieve and convert compressed data with the prefilter to a NumPy array.
        >>> compressed_array_with_filter = np.frombuffer(schunk.get_slice(), dtype=output_dtype)
        >>> f"Compressed data with prefilter applied (first 8 elements): {compressed_array_with_filter[:8]}"
        Compressed data with prefilter applied (first 8 elements): [-3.1415927  -2.1415927  -1.1415926  -0.14159265  0.8584073   1.8584074
         2.8584073   3.8584073 ]
        >>> schunk.remove_prefilter('prefilter')
        >>> schunk[:1000] = data
        >>> compressed_array_without_filter = np.frombuffer(schunk.get_slice(), dtype=dtype)
        >>> f"Compressed data without prefilter (first 8 elements): {compressed_array_without_filter[:8]}"
        Compressed data without prefilter (first 8 elements): [0. 1. 2. 3. 4. 5. 6. 7.]
        """
        return super().remove_prefilter(func_name)

    def __dealloc__(self):
        super().__dealloc__()


def open(
    urlpath: str | pathlib.Path | blosc2.URLPath, mode: str = "a", offset: int = 0, **kwargs: dict
) -> blosc2.SChunk | blosc2.NDArray | blosc2.C2Array | blosc2.LazyArray | blosc2.Proxy:
    """Open a persistent :ref:`SChunk`, :ref:`NDArray`, a remote :ref:`C2Array`
    or a :ref:`Proxy`

    See the `Notes` section for more info on opening `Proxy` objects.

    Parameters
    ----------
    urlpath: str | pathlib.Path | :ref:`URLPath`
        The path where the :ref:`SChunk` (or :ref:`NDArray`)
        is stored. If it is a remote array, a :ref:`URLPath` must be passed.
    mode: str, optional
        The open mode.
    offset: int, optional
        An offset in the file where super-chunk or array data is located
        (e.g. in a file containing several such objects).
    kwargs: dict, optional
        mmap_mode: The memory mapping mode.
        initial_mapping_size: The initial size of the memory mapping.
        cparams: dict
            A dictionary with the compression parameters, which are the same that can be
            used in the :func:`~blosc2.compress2` function.
            Typesize and blocksize cannot be changed.
        dparams: dict
            A dictionary with the decompression parameters, which are the same that can
            be used in the :func:`~blosc2.decompress2` function.

    Returns
    -------
    out: :ref:`SChunk`, :ref:`NDArray` or :ref:`C2Array`
        The SChunk or NDArray (if there is a "b2nd" metalayer")
        or the C2Array if :paramref:`urlpath` is a :ref:`blosc2.URLPath <URLPath>` instance.

    Notes
    -----
    * This is just a 'logical' open, so there is no `close()` counterpart because
      currently, there is no need for it.

    * If :paramref:`urlpath` is a :ref:`URLPath` instance, :paramref:`mode`
      must be 'r', :paramref:`offset` must be 0, and kwargs cannot be passed.

    * If the original object saved in :paramref:`urlpath` is a :ref:`Proxy`,
      this function will only return a :ref:`Proxy` if its source is a local
      :ref:`SChunk`, :ref:`NDArray` or a remote :ref:`C2Array`. Otherwise,
      it will return the Python-Blosc2 container used to cache the data which
      can be a :ref:`SChunk` or a :ref:`NDArray` and may not have all the data
      initialized (e.g. if the user has not accessed to it yet).

    * When opening a :ref:`LazyExpr` keep in mind the note above regarding operands.

    Examples
    --------
    >>> import blosc2
    >>> import numpy as np
    >>> import os
    >>> import tempfile
    >>> tmpdirname = tempfile.mkdtemp()
    >>> urlpath = os.path.join(tmpdirname, 'b2frame')
    >>> storage = blosc2.Storage(contiguous=True, urlpath=urlpath, mode="w")
    >>> nelem = 20 * 1000
    >>> nchunks = 5
    >>> chunksize = nelem * 4 // nchunks
    >>> data = np.arange(nelem, dtype="int32")
    >>> # Create SChunk and append data
    >>> schunk = blosc2.SChunk(chunksize=chunksize, data=data.tobytes(), storage=storage)
    >>> # Open SChunk
    >>> sc_open = blosc2.open(urlpath=urlpath)
    >>> for i in range(nchunks):
    ...     dest = np.empty(nelem // nchunks, dtype=data.dtype)
    ...     schunk.decompress_chunk(i, dest)
    ...     dest1 = np.empty(nelem // nchunks, dtype=data.dtype)
    ...     sc_open.decompress_chunk(i, dest1)
    ...     np.array_equal(dest, dest1)
    True
    True
    True
    True
    True

    To open the same schunk memory-mapped, we simply need to pass the `mmap_mode` parameter:

    >>> sc_open_mmap = blosc2.open(urlpath=urlpath, mmap_mode="r")
    >>> sc_open.nchunks == sc_open_mmap.nchunks
    True
    >>> all(sc_open.decompress_chunk(i, dest1) == sc_open_mmap.decompress_chunk(i, dest1) for i in range(nchunks))
    True
    """
    if isinstance(urlpath, blosc2.URLPath):
        if mode != "r" or offset != 0 or kwargs != {}:
            raise NotImplementedError(
                "Cannot open a C2Array with mode != 'r', or offset != 0 or some kwargs"
            )
        return blosc2.C2Array(urlpath.path, urlbase=urlpath.urlbase, auth_token=urlpath.auth_token)

    if isinstance(urlpath, pathlib.PurePath):
        urlpath = str(urlpath)
    if not os.path.exists(urlpath):
        raise FileNotFoundError(f"No such file or directory: {urlpath}")

    res = blosc2_ext.open(urlpath, mode, offset, **kwargs)

    meta = getattr(res, "schunk", res).meta
    if "proxy-source" in meta:
        proxy_src = meta["proxy-source"]
        if proxy_src["local_abspath"] is not None:
            src = blosc2.open(proxy_src["local_abspath"])
            return blosc2.Proxy(src, _cache=res)
        elif proxy_src["urlpath"] is not None:
            src = blosc2.C2Array(proxy_src["urlpath"][0], proxy_src["urlpath"][1], proxy_src["urlpath"][2])
            return blosc2.Proxy(src, _cache=res)
        elif not proxy_src["caterva2_env"]:
            raise RuntimeError("Could not find the source when opening a Proxy")

    if isinstance(res, blosc2.NDArray) and "LazyArray" in res.schunk.meta:
        return blosc2._open_lazyarray(res)
    else:
        return res
