import asyncio
import contextlib

import aiohttp
import aiohttp.client_exceptions
import botocore.response
from botocore.response import (
    ReadTimeoutError,
    ResponseStreamingError,
)

from aiobotocore import parsers


class AioReadTimeoutError(ReadTimeoutError, asyncio.TimeoutError):
    pass


class AioStreamingBodyBase(
    botocore.response.StreamingBody, contextlib.AbstractAsyncContextManager
):
    """Shared async base for aiobotocore response bodies.

    Subclasses ``botocore.response.StreamingBody`` to provide the
    backend-agnostic async API (iteration, line splitting, content-length
    validation, async context manager) and to block the inherited sync API
    surface. Concrete subclasses implement the backend-specific I/O
    (``read``, ``readinto``, ``readable``, ``close``):
    :class:`AioStreamingBody` for aiohttp and :class:`AioHttpxStreamingBody`
    for httpx.

    Access the underlying raw response object via the ``raw_stream`` property.
    """

    _DEFAULT_CHUNK_SIZE = 1024

    @property
    def raw_stream(self):
        """Access the underlying raw HTTP response object."""
        return self._raw_stream

    async def __aenter__(self):
        return self

    async def __aexit__(self, exc_type, exc_val, exc_tb):
        await self.aclose()

    async def readlines(self):
        return [line async for line in self.iter_lines()]

    def __aiter__(self):
        """Return an async iterator to yield 1k chunks from the raw stream."""
        return self.iter_chunks(self._DEFAULT_CHUNK_SIZE)

    async def __anext__(self):
        """Return the next 1k chunk from the raw stream."""
        current_chunk = await self.read(self._DEFAULT_CHUNK_SIZE)
        if current_chunk:
            return current_chunk
        raise StopAsyncIteration

    anext = __anext__

    async def iter_lines(self, chunk_size=_DEFAULT_CHUNK_SIZE, keepends=False):
        """Return an async iterator to yield lines from the raw stream."""
        pending = b''
        async for chunk in self.iter_chunks(chunk_size):
            lines = (pending + chunk).splitlines(True)
            for line in lines[:-1]:
                yield line.splitlines(keepends)[0]
            pending = lines[-1]
        if pending:
            yield pending.splitlines(keepends)[0]

    async def iter_chunks(self, chunk_size=_DEFAULT_CHUNK_SIZE):
        """Return an async iterator to yield chunks of chunk_size bytes
        from the raw stream.
        """
        while True:
            current_chunk = await self.read(chunk_size)
            if current_chunk == b"":
                break
            yield current_chunk

    def tell(self):
        return self._amount_read

    # Block sync methods inherited from botocore.StreamingBody / IOBase
    # that would otherwise call urllib3-shaped APIs on an async stream.
    def __iter__(self):
        raise TypeError(
            f"{type(self).__name__} is async; use 'async for' or __aiter__"
        )

    def __next__(self):
        raise TypeError(f"{type(self).__name__} is async; use __anext__")

    next = __next__

    def __enter__(self):
        raise TypeError(f"{type(self).__name__} is async; use 'async with'")

    def __exit__(self, exc_type, exc_val, exc_tb):
        raise TypeError(f"{type(self).__name__} is async; use 'async with'")

    def set_socket_timeout(self, timeout):
        raise NotImplementedError(
            "set_socket_timeout is not supported for async streaming bodies; "
            "configure timeouts on the aiobotocore Config or HTTP session "
            "instead."
        )


class AioStreamingBody(AioStreamingBodyBase):
    """Async response body for the aiohttp backend."""

    def readable(self):
        return not self._raw_stream.content.at_eof()

    async def read(self, amt=None):
        """Read at most amt bytes from the stream.

        If the amt argument is omitted, read all data.
        """
        try:
            chunk = await self._raw_stream.content.read(
                amt if amt is not None else -1
            )
        except asyncio.TimeoutError as e:
            raise AioReadTimeoutError(
                endpoint_url=self._raw_stream.url, error=e
            )
        except aiohttp.client_exceptions.ClientConnectionError as e:
            raise ResponseStreamingError(error=e)

        self._amount_read += len(chunk)
        if amt is None or (not chunk and amt > 0):
            self._verify_content_length()
        return chunk

    async def readinto(self, b: bytearray):
        """Read bytes into a pre-allocated, writable bytes-like object b,
        and return the number of bytes read.
        """
        try:
            chunk = await self._raw_stream.content.read(len(b))
            amount_read = len(chunk)
            b[:amount_read] = chunk
        except asyncio.TimeoutError as e:
            raise AioReadTimeoutError(
                endpoint_url=self._raw_stream.url, error=e
            )
        except aiohttp.client_exceptions.ClientConnectionError as e:
            raise ResponseStreamingError(error=e)

        self._amount_read += amount_read
        if amount_read == 0 and len(b) > 0:
            self._verify_content_length()
        return amount_read

    def close(self):
        """Close the underlying response stream synchronously."""
        self._raw_stream.close()

    async def aclose(self):
        """Close the underlying response stream asynchronously."""
        self.close()


class AioHttpxStreamingBody(AioStreamingBodyBase):
    """Async response body for the httpx backend.

    httpx.Response does not support ``read(amt)`` natively — it yields
    chunks via ``aiter_bytes()``. This class implements an internal buffer
    so ``read(amt)``, ``readinto()``, ``iter_chunks(size)`` and
    ``iter_lines()`` all work with arbitrary sizes. The backend-agnostic
    API (iteration, line splitting, tell, content-length validation) is
    inherited from ``AioStreamingBodyBase``.
    """

    def __init__(self, raw_stream, content_length=None):
        super().__init__(raw_stream, content_length)
        self._buffer = b''
        self._stream_iter = None
        self._stream_exhausted = False

    def _ensure_stream(self):
        if self._stream_iter is None:
            # aiter_bytes decodes; counts must match wire Content-Length/checksums
            self._stream_iter = self._raw_stream.aiter_raw().__aiter__()

    async def _fill_buffer(self, min_bytes):
        """Fill internal buffer until it has at least min_bytes or the
        stream is exhausted.
        """
        self._ensure_stream()
        while len(self._buffer) < min_bytes and not self._stream_exhausted:
            try:
                chunk = await self._stream_iter.__anext__()
                self._buffer += chunk
            except StopAsyncIteration:
                self._stream_exhausted = True

    async def read(self, amt=None):
        """Read at most amt bytes from the stream.

        If the amt argument is omitted, read all data.
        """
        self._ensure_stream()

        # a negative amt means read-all to file-like callers, as it does to aiohttp
        if amt is None or amt < 0:
            chunks = [self._buffer]
            self._buffer = b''
            async for chunk in self._stream_iter:
                chunks.append(chunk)
            self._stream_exhausted = True
            result = b''.join(chunks)
        elif amt == 0:
            return b''
        else:
            await self._fill_buffer(amt)
            result = self._buffer[:amt]
            self._buffer = self._buffer[amt:]

        self._amount_read += len(result)
        if amt is None or (not result and amt > 0):
            self._verify_content_length()
        return result

    async def readinto(self, b: bytearray):
        if len(b) == 0:
            return 0

        await self._fill_buffer(len(b))
        chunk = self._buffer[: len(b)]
        amount_read = len(chunk)
        b[:amount_read] = chunk
        self._buffer = self._buffer[amount_read:]

        self._amount_read += amount_read
        if amount_read == 0 and len(b) > 0:
            self._verify_content_length()
        return amount_read

    def readable(self):
        return bool(self._buffer) or not self._stream_exhausted

    async def close(self):
        """Close the underlying httpx response (async-only — httpx has no
        synchronous close).
        """
        await self._raw_stream.aclose()

    aclose = close


# Backwards-compatibility aliases. External code imports these names from
# aiobotocore.response; keeping them here avoids a breaking rename.
# ``StreamingBody`` stays bound to the aiohttp concrete class so that
# ``isinstance(body, StreamingBody)`` continues to distinguish the aiohttp
# backend (sync ``close()``) from the httpx backend (async ``aclose()``).
StreamingBody = AioStreamingBody
HttpxStreamingBody = AioHttpxStreamingBody


async def get_response(operation_model, http_response):
    protocol = operation_model.service_model.resolved_protocol
    response_dict = {
        'headers': http_response.headers,
        'status_code': http_response.status_code,
    }
    # TODO: Unfortunately, we have to have error logic here.
    # If it looks like an error, in the streaming response case we
    # need to actually grab the contents.
    if response_dict['status_code'] >= 300:
        response_dict['body'] = await http_response.content
    elif operation_model.has_streaming_output:
        response_dict['body'] = AioStreamingBody(
            http_response.raw, response_dict['headers'].get('content-length')
        )
    else:
        response_dict['body'] = await http_response.content

    parser = parsers.create_parser(protocol)
    parsed = await parser.parse(response_dict, operation_model.output_shape)
    return http_response, parsed
