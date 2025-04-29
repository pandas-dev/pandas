#######################################################################
# Copyright (c) 2019-present, Blosc Development Team <blosc@blosc.org>
# All rights reserved.
#
# This source code is licensed under a BSD-style license (found in the
# LICENSE file in the root directory of this source tree)
#######################################################################
from abc import ABC, abstractmethod

import numpy as np

import blosc2


class ProxyNDSource(ABC):
    """
    Base interface for NDim sources in :ref:`Proxy`.
    """

    @property
    @abstractmethod
    def shape(self) -> tuple:
        """
        The shape of the source.
        """
        pass

    @property
    @abstractmethod
    def chunks(self) -> tuple:
        """
        The chunk shape of the source.
        """
        pass

    @property
    @abstractmethod
    def blocks(self) -> tuple:
        """
        The block shape of the source.
        """
        pass

    @property
    @abstractmethod
    def dtype(self) -> np.dtype:
        """
        The dtype of the source.
        """
        pass

    @property
    def cparams(self) -> blosc2.CParams:
        """
        The compression parameters of the source.

        This property is optional and can be overridden if the source has a
        different compression configuration.
        """
        return blosc2.CParams(typesize=self.dtype.itemsize)

    @abstractmethod
    def get_chunk(self, nchunk: int) -> bytes:
        """
        Return the compressed chunk in :paramref:`self`.

        Parameters
        ----------
        nchunk: int
            The unidimensional index of the chunk to retrieve.

        Returns
        -------
        out: bytes object
            The compressed chunk.
        """
        pass

    async def aget_chunk(self, nchunk: int) -> bytes:
        """
        Return the compressed chunk in :paramref:`self` asynchronously.

        Parameters
        ----------
        nchunk: int
            The index of the chunk to retrieve.

        Returns
        -------
        out: bytes object
            The compressed chunk.

        Notes
        -----
        This method is optional, and only available if the source has an async
        `aget_chunk` method.
        """
        raise NotImplementedError(
            "aget_chunk is only available if the source has an async aget_chunk method"
        )


class ProxySource(ABC):
    """
    Base interface for sources of :ref:`Proxy` that are not NDim objects.
    """

    @property
    @abstractmethod
    def nbytes(self) -> int:
        """
        The total number of bytes in the source.
        """
        pass

    @property
    @abstractmethod
    def chunksize(self) -> tuple:
        """
        The chunksize of the source.
        """
        pass

    @property
    @abstractmethod
    def typesize(self) -> int:
        """
        The typesize of the source.
        """
        pass

    @property
    def cparams(self) -> blosc2.CParams:
        """
        The compression parameters of the source.

        This property is optional and can be overridden if the source has a
        different compression configuration.
        """
        return blosc2.CParams(typesize=self.typesize)

    @abstractmethod
    def get_chunk(self, nchunk: int) -> bytes:
        """
        Return the compressed chunk in :paramref:`self`.

        Parameters
        ----------
        nchunk: int
            The index of the chunk to retrieve.

        Returns
        -------
        out: bytes object
            The compressed chunk.
        """
        pass

    async def aget_chunk(self, nchunk: int) -> bytes:
        """
        Return the compressed chunk in :paramref:`self` asynchronously.

        Parameters
        ----------
        nchunk: int
            The index of the chunk to retrieve.

        Returns
        -------
        out: bytes object
            The compressed chunk.

        Notes
        -----
        This method is optional and only available if the source has an async
        `aget_chunk` method.
        """
        raise NotImplementedError(
            "aget_chunk is only available if the source has an async aget_chunk method"
        )


class Proxy(blosc2.Operand):
    """Proxy (with cache support) for an object following the :ref:`ProxySource` interface.

    This can be used to cache chunks of a regular data container which follows the
    :ref:`ProxySource` or :ref:`ProxyNDSource` interfaces.
    """

    def __init__(
        self, src: ProxySource or ProxyNDSource, urlpath: str | None = None, mode="a", **kwargs: dict
    ):
        """
        Create a new :ref:`Proxy` to serve as a cache to save accessed chunks locally.

        Parameters
        ----------
        src: :ref:`ProxySource` or :ref:`ProxyNDSource`
            The original container.
        urlpath: str, optional
            The urlpath where to save the container that will work as a cache.
        mode: str, optional
            "a" means read/write (create if it doesn't exist); "w" means create
            (overwrite if it exists). Default is "a".
        kwargs: dict, optional
            Keyword arguments supported:

                vlmeta: dict or None
                    A dictionary with different variable length metalayers.  One entry per metalayer:
                        key: bytes or str
                            The name of the metalayer.
                        value: object
                            The metalayer object that will be serialized using msgpack.

        """
        self.src = src
        self.urlpath = urlpath
        if kwargs is None:
            kwargs = {}
        self._cache = kwargs.pop("_cache", None)

        if self._cache is None:
            meta_val = {
                "local_abspath": None,
                "urlpath": None,
                "caterva2_env": kwargs.pop("caterva2_env", False),
            }
            container = getattr(self.src, "schunk", self.src)
            if hasattr(container, "urlpath"):
                meta_val["local_abspath"] = container.urlpath
            elif isinstance(self.src, blosc2.C2Array):
                meta_val["urlpath"] = (self.src.path, self.src.urlbase, self.src.auth_token)
            meta = {"proxy-source": meta_val}
            if hasattr(self.src, "shape"):
                self._cache = blosc2.empty(
                    self.src.shape,
                    self.src.dtype,
                    chunks=self.src.chunks,
                    blocks=self.src.blocks,
                    cparams=self.src.cparams,
                    urlpath=urlpath,
                    mode=mode,
                    meta=meta,
                )
            else:
                self._cache = blosc2.SChunk(
                    chunksize=self.src.chunksize,
                    cparams=self.src.cparams,
                    urlpath=urlpath,
                    mode=mode,
                    meta=meta,
                )
                self._cache.fill_special(self.src.nbytes // self.src.typesize, blosc2.SpecialValue.UNINIT)
        self._schunk_cache = getattr(self._cache, "schunk", self._cache)
        vlmeta = kwargs.get("vlmeta")
        if vlmeta:
            for key in vlmeta:
                self._schunk_cache.vlmeta[key] = vlmeta[key]

    def fetch(self, item: slice | list[slice] | None = None) -> blosc2.NDArray | blosc2.schunk.SChunk:
        """
        Get the container used as cache with the requested data updated.

        Parameters
        ----------
        item: slice or list of slices, optional
            If not None, only the chunks that intersect with the slices
            in items will be retrieved if they have not been already.

        Returns
        -------
        out: :ref:`NDArray` or :ref:`SChunk`
            The local container used to cache the already requested data.

        Examples
        --------
        >>> import numpy as np
        >>> import blosc2
        >>> data = np.arange(20).reshape(10, 2)
        >>> ndarray = blosc2.asarray(data)
        >>> proxy = blosc2.Proxy(ndarray)
        >>> slice_data = proxy.fetch((slice(0, 3), slice(0, 2)))
        >>> slice_data[:3, :2]
        [[0 1]
        [2 3]
        [4 5]]
        """
        if item is None:
            # Full realization
            for info in self._schunk_cache.iterchunks_info():
                if info.special != blosc2.SpecialValue.NOT_SPECIAL:
                    chunk = self.src.get_chunk(info.nchunk)
                    self._schunk_cache.update_chunk(info.nchunk, chunk)
        else:
            # Get only a slice
            nchunks = blosc2.get_slice_nchunks(self._cache, item)
            for info in self._schunk_cache.iterchunks_info():
                if info.nchunk in nchunks and info.special != blosc2.SpecialValue.NOT_SPECIAL:
                    chunk = self.src.get_chunk(info.nchunk)
                    self._schunk_cache.update_chunk(info.nchunk, chunk)

        return self._cache

    async def afetch(self, item: slice | list[slice] | None = None) -> blosc2.NDArray | blosc2.schunk.SChunk:
        """
        Retrieve the cache container with the requested data updated asynchronously.

        Parameters
        ----------
        item: slice or list of slices, optional
            If provided, only the chunks intersecting with the specified slices
            will be retrieved if they have not been already.

        Returns
        -------
        out: :ref:`NDArray` or :ref:`SChunk`
            The local container used to cache the already requested data.

        Notes
        -----
        This method is only available if the :ref:`ProxySource` or :ref:`ProxyNDSource`
        have an async `aget_chunk` method.

        Examples
        --------
        >>> import numpy as np
        >>> import blosc2
        >>> import asyncio
        >>> from blosc2 import ProxyNDSource
        >>> class MyProxySource(ProxyNDSource):
        >>>     def __init__(self, data):
        >>>         # If the next source is multidimensional, it must have the attributes:
        >>>         self.data = data
        >>>         f"Data shape: {self.shape}, Chunks: {self.chunks}"
        >>>         f"Blocks: {self.blocks}, Dtype: {self.dtype}"
        >>>     @property
        >>>     def shape(self):
        >>>         return self.data.shape
        >>>     @property
        >>>     def chunks(self):
        >>>         return self.data.chunks
        >>>     @property
        >>>     def blocks(self):
        >>>         return self.data.blocks
        >>>     @property
        >>>     def dtype(self):
        >>>         return self.data.dtype
        >>>     # This method must be present
        >>>     def get_chunk(self, nchunk):
        >>>         return self.data.get_chunk(nchunk)
        >>>     # This method is optional
        >>>     async def aget_chunk(self, nchunk):
        >>>         await asyncio.sleep(0.1) # Simulate an asynchronous operation
        >>>         return self.data.get_chunk(nchunk)
        >>> data = np.arange(20).reshape(4, 5)
        >>> chunks = [2, 5]
        >>> blocks = [1, 5]
        >>> data = blosc2.asarray(data, chunks=chunks, blocks=blocks)
        >>> source = MyProxySource(data)
        >>> proxy = blosc2.Proxy(source)
        >>> async def fetch_data():
        >>>     # Fetch a slice of the data from the proxy asynchronously
        >>>     slice_data = await proxy.afetch(slice(0, 2))
        >>>     # Note that only data fetched is shown, the rest is uninitialized
        >>>     slice_data[:]
        >>> asyncio.run(fetch_data())
        >>> # Using getitem to get a slice of the data
        >>> result = proxy[1:2, 1:3]
        >>> f"Proxy getitem: {result}"
        Data shape: (4, 5), Chunks: (2, 5)
        Blocks: (1, 5), Dtype: int64
        [[0 1 2 3 4]
        [5 6 7 8 9]
        [0 0 0 0 0]
        [0 0 0 0 0]]
        Proxy getitem: [[6 7]]
        """
        if not callable(getattr(self.src, "aget_chunk", None)):
            raise NotImplementedError("afetch is only available if the source has an aget_chunk method")
        if item is None:
            # Full realization
            for info in self._schunk_cache.iterchunks_info():
                if info.special != blosc2.SpecialValue.NOT_SPECIAL:
                    chunk = await self.src.aget_chunk(info.nchunk)
                    self._schunk_cache.update_chunk(info.nchunk, chunk)
        else:
            # Get only a slice
            nchunks = blosc2.get_slice_nchunks(self._cache, item)
            for info in self._schunk_cache.iterchunks_info():
                if info.nchunk in nchunks and info.special != blosc2.SpecialValue.NOT_SPECIAL:
                    chunk = await self.src.aget_chunk(info.nchunk)
                    self._schunk_cache.update_chunk(info.nchunk, chunk)

        return self._cache

    def __getitem__(self, item: slice | list[slice]) -> np.ndarray:
        """
        Get a slice as a numpy.ndarray using the :ref:`Proxy`.

        Parameters
        ----------
        item: slice or list of slices
            The slice of the desired data.

        Returns
        -------
        out: numpy.ndarray
            An array with the data slice.

        Examples
        --------
        >>> import numpy as np
        >>> import blosc2
        >>> data = np.arange(25).reshape(5, 5)
        >>> ndarray = blosc2.asarray(data)
        >>> proxy = blosc2.Proxy(ndarray)
        >>> proxy[0:3, 0:3]
        [[ 0  1  2]
        [ 5  6  7]
        [10 11 12]
        [20 21 22]]
        >>> proxy[2:5, 2:5]
        [[12 13 14]
        [17 18 19]
        [22 23 24]]
        """
        # Populate the cache
        self.fetch(item)
        return self._cache[item]

    @property
    def dtype(self) -> np.dtype:
        """The dtype of :paramref:`self` or None if the data is unidimensional"""
        return self._cache.dtype if isinstance(self._cache, blosc2.NDArray) else None

    @property
    def shape(self) -> tuple[int]:
        """The shape of :paramref:`self`"""
        return self._cache.shape if isinstance(self._cache, blosc2.NDArray) else len(self._cache)

    @property
    def schunk(self) -> blosc2.schunk.SChunk:
        """The :ref:`SChunk` of the cache"""
        return self._schunk_cache

    @property
    def cparams(self) -> blosc2.CParams:
        """The compression parameters of the cache"""
        return self._cache.cparams

    @property
    def info(self) -> str:
        """The info of the cache"""
        if isinstance(self._cache, blosc2.NDArray):
            return self._cache.info
        raise NotImplementedError("info is only available if the source is a NDArray")

    def __str__(self):
        return f"Proxy({self.src}, urlpath={self.urlpath})"

    @property
    def vlmeta(self) -> blosc2.schunk.vlmeta:
        """
        Get the vlmeta of the cache.

        See Also
        --------
        :py:attr:`blosc2.schunk.SChunk.vlmeta`
        """
        return self._schunk_cache.vlmeta

    @property
    def fields(self) -> dict:
        """
        Dictionary with the fields of :paramref:`self`.

        Returns
        -------
        fields: dict
            A dictionary with the fields of the :ref:`Proxy`.

        See Also
        --------
        :ref:`NDField`

        Examples
        --------
        >>> import numpy as np
        >>> import blosc2
        >>> data = np.ones(16, dtype=[('field1', 'i4'), ('field2', 'f4')]).reshape(4, 4)
        >>> ndarray = blosc2.asarray(data)
        >>> proxy = blosc2.Proxy(ndarray)
        >>>  # Get a dictionary of fields from the proxy, where each field can be accessed individually
        >>> fields_dict = proxy.fields
        >>> for field_name, field_proxy in fields_dict.items():
        >>>     print(f"Field name: {field_name}, Field data: {field_proxy}")
        Field name: field1, Field data: <blosc2.proxy.ProxyNDField object at 0x114472d20>
        Field name: field2, Field data: <blosc2.proxy.ProxyNDField object at 0x10e215be0>
        >>> fields_dict['field2'][:]
        [[1. 1. 1. 1.]
         [1. 1. 1. 1.]
         [1. 1. 1. 1.]
         [1. 1. 1. 1.]]
        """
        _fields = getattr(self._cache, "fields", None)
        if _fields is None:
            return None
        return {key: ProxyNDField(self, key) for key in _fields}


class ProxyNDField(blosc2.Operand):
    def __init__(self, proxy: Proxy, field: str):
        self.proxy = proxy
        self.field = field
        self.dtype = proxy.dtype[field]
        self.shape = proxy.shape

    def __getitem__(self, item: slice | list[slice]) -> np.ndarray:
        """
        Get a slice as a numpy.ndarray using the `field` in `proxy`.

        Parameters
        ----------
        item: slice or list of slices
            The slice of the desired data.

        Returns
        -------
        out: numpy.ndarray
            An array with the data slice.
        """
        # Get the data and return the corresponding field
        nparr = self.proxy[item]
        return nparr[self.field]


class SimpleProxy(blosc2.Operand):
    """
    Simple proxy for a NumPy array (or similar) that can be used with the Blosc2 compute engine.

    This only supports the __getitem__ method. No caching is performed.

    Examples
    --------
    >>> import numpy as np
    >>> import blosc2
    >>> a = np.arange(20, dtype=np.float32).reshape(4, 5)
    >>> proxy = blosc2.SimpleProxy(a)
    >>> proxy[1:3, 2:4]
    [[ 7.  8.]
     [12. 13.]]
    """

    def __init__(self, src, chunks: tuple | None = None, blocks: tuple | None = None):
        if not hasattr(src, "shape") or not hasattr(src, "dtype"):
            # If the source is not a NumPy array, convert it to one
            src = np.asarray(src)
        self.src = src
        self.dtype = src.dtype
        self.shape = src.shape
        # Compute reasonable values for chunks and blocks
        cparams = blosc2.CParams(clevel=0)
        self.chunks, self.blocks = blosc2.compute_chunks_blocks(
            self.shape, chunks, blocks, self.dtype, **{"cparams": cparams}
        )

    def __getitem__(self, item: slice | list[slice]) -> np.ndarray:
        """
        Get a slice as a numpy.ndarray using the :ref:`ProxyNumPy`.

        Parameters
        ----------
        item

        Returns
        -------
        out: numpy.ndarray
            An array with the data slice.
        """
        return self.src[item]


def jit(func=None, *, out=None, **kwargs):  # noqa: C901
    """
    Prepare a function so that it can be used with the Blosc2 compute engine.

    The inputs of the function can be any combination of NumPy/NDArray arrays
    and scalars.  The function will be called with the NumPy arrays replaced by
    :ref:`SimpleProxy` objects, whereas NDArray objects will be used as is.

    The returned value will be a NumPy array if all arguments are NumPy arrays
    or if not kwargs are provided. Else, the return value will be a NDArray
    created using the provided kwargs.

    Parameters
    ----------
    func: callable
        The function to be prepared for the Blosc2 compute engine.
    out: np.ndarray, NDArray, optional
        The output array where the result will be stored.
    **kwargs: dict, optional
        Additional keyword arguments supported by the :func:`empty` constructor.

    Returns
    -------
    wrapper

    Notes
    -----
    * Although many NumPy functions are supported, some may not be implemented yet.
      If you find a function that is not supported, please open an issue.
    * `out` and `kwargs` parameters are not supported for all expressions
      (e.g. when using a reduction as the last function).  In this case, you can
      still use the `out` parameter of the reduction function for some custom
      control over the output.

    Examples
    --------
    >>> import numpy as np
    >>> import blosc2
    >>> @blosc2.jit
    >>> def compute_expression(a, b, c):
    >>>     return np.sum(((a ** 3 + np.sin(a * 2)) > 2 * c) & (b > 0), axis=1)
    >>> a = np.arange(20, dtype=np.float32).reshape(4, 5)
    >>> b = np.arange(20).reshape(4, 5)
    >>> c = np.arange(5)
    >>> compute_expression(a, b, c)
    [5 5 5 5]
    """

    def decorator(func):  # noqa: C901
        def wrapper(*args, **func_kwargs):
            # Get some kwargs in decorator for SimpleProxy constructor
            proxy_kwargs = {"chunks": kwargs.get("chunks"), "blocks": kwargs.get("blocks")}

            # Wrap the arguments in SimpleProxy objects if they are not NDArrays
            new_args = []
            for arg in args:
                if issubclass(type(arg), blosc2.Operand):
                    new_args.append(arg)
                else:
                    new_args.append(SimpleProxy(arg, **proxy_kwargs))
            # The same for the keyword arguments
            for key, value in func_kwargs.items():
                if issubclass(type(value), blosc2.Operand):
                    continue
                func_kwargs[key] = SimpleProxy(value, **proxy_kwargs)

            # Call function with the new arguments
            retval = func(*new_args, **func_kwargs)

            # Treat return value
            # If it is a numpy array, return it as is
            if isinstance(retval, np.ndarray):
                if kwargs and any(kwargs[key] is not None for key in kwargs):
                    # But if kwargs are provided, return a NDArray instead
                    return blosc2.asarray(retval, **kwargs)
                return retval

            # In some instances, the return value is not a LazyExpr
            # (e.g. using a reduction as the last function, and using an `out` param)
            if not isinstance(retval, blosc2.LazyExpr):
                return retval

            # If the return value is a LazyExpr, compute it
            if out is not None:
                return retval.compute(out=out, **kwargs)
            if kwargs and any(kwargs[key] is not None for key in kwargs):
                return retval.compute(**kwargs)
            # If no kwargs are provided, return a numpy array
            return retval[()]

        return wrapper

    if func is None:
        return decorator
    else:
        return decorator(func)
