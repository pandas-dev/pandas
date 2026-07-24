import asyncio
import collections
import concurrent.futures
import contextlib
import ctypes
import logging
import os
import threading
import weakref
from io import BytesIO

from google.api_core.exceptions import NotFound
from google.cloud.storage.asyncio.async_appendable_object_writer import (
    _DEFAULT_FLUSH_INTERVAL_BYTES,
    AsyncAppendableObjectWriter,
)
from google.cloud.storage.asyncio.async_multi_range_downloader import (
    AsyncMultiRangeDownloader,
)

MRD_MAX_RANGES = 1000  # MRD supports up to 1000 ranges per request
DEFAULT_CONCURRENCY = int(os.environ.get("DEFAULT_GCSFS_CONCURRENCY", "1"))
MAX_PREFETCH_SIZE = 256 * 1024 * 1024
logger = logging.getLogger("gcsfs")


try:
    PyBytes_FromStringAndSize = ctypes.pythonapi.PyBytes_FromStringAndSize
    PyBytes_FromStringAndSize.argtypes = (ctypes.c_void_p, ctypes.c_ssize_t)
    PyBytes_FromStringAndSize.restype = ctypes.py_object

    PyBytes_AsString = ctypes.pythonapi.PyBytes_AsString
    PyBytes_AsString.argtypes = (ctypes.py_object,)
    PyBytes_AsString.restype = ctypes.c_void_p
    HAS_CPYTHON_API = True
except Exception:
    PyBytes_FromStringAndSize = None
    PyBytes_AsString = None
    HAS_CPYTHON_API = False


async def init_mrd(grpc_client, bucket_name, object_name, generation=None):
    """
    Creates the AsyncMultiRangeDownloader using an existing client.
    Wraps Google API errors into standard Python exceptions.
    """
    try:
        return await AsyncMultiRangeDownloader.create_mrd(
            grpc_client, bucket_name, object_name, generation
        )
    except NotFound:
        # We wrap the error here to match standard Python error handling
        # and avoid leaking Google API exceptions to users.
        raise FileNotFoundError(f"{bucket_name}/{object_name}")


async def download_range(offset, length, mrd):
    """
    Downloads a byte range from the file asynchronously.
    """
    # If length = 0, mrd returns till end of file, so handle that case here
    if length == 0:
        return b""
    buffer = BytesIO()
    await mrd.download_ranges([(offset, length, buffer)])
    data = buffer.getvalue()
    bytes_downloaded = len(data)

    if length != bytes_downloaded:
        logger.warning(
            f"Short read detected for {mrd.bucket_name}/{mrd.object_name}! "
            f"Requested {length} bytes but downloaded {bytes_downloaded} bytes."
        )

    logger.debug(
        f"Requested {length} bytes from offset {offset}, downloaded {bytes_downloaded} "
        f"bytes from mrd path: {mrd.bucket_name}/{mrd.object_name}"
    )
    return data


async def download_ranges(ranges, mrd):
    """
    Downloads multiple byte ranges from the file asynchronously in a single batch.

    Args:
        ranges: List of (offset, length) tuples to download. Max 1000 ranges allowed.
        mrd: AsyncMultiRangeDownloader instance

    Returns:
        List of bytes objects, one for each range
    """
    # Prepare tasks: Filter out empty ranges and create buffers immediately
    # Structure: (original_index, offset, length, buffer)
    # Calling MRD with length=0 returns till end of file. We handle zero-length
    # ranges by returning b"" without calling MRD. So only create tasks for length > 0

    if len(ranges) > MRD_MAX_RANGES:
        raise ValueError("Invalid input - number of ranges cannot be more than 1000")

    tasks = [
        (i, off, length, BytesIO())
        for i, (off, length) in enumerate(ranges)
        if length > 0
    ]

    # Execute Download
    if tasks:
        # The MRD expects list of (offset, length, buffer)
        # We extract these from our task list
        await mrd.download_ranges([(off, length, buf) for _, off, length, buf in tasks])

    # Map results back to their original positions
    results = [b""] * len(ranges)
    for i, _, _, buffer in tasks:
        results[i] = buffer.getvalue()

    # Log stats
    total_requested = sum(r[1] for r in ranges)
    total_downloaded = sum(len(r) for r in results)

    if total_requested != total_downloaded:
        logger.warning(
            f"Short read detected for {mrd.bucket_name}/{mrd.object_name}! "
            f"Requested {total_requested} bytes but downloaded {total_downloaded} bytes."
        )

    if logger.isEnabledFor(logging.DEBUG):
        requested_ranges_to_log = [(r[0], r[1]) for r in ranges]
        logger.debug(
            f"mrd path: {mrd.bucket_name}/{mrd.object_name} | "
            f"Requested {len(ranges)} ranges: {requested_ranges_to_log} | "
            f"total bytes requested: {total_requested} | "
            f"total bytes downloaded: {total_downloaded}"
        )

    return results


async def init_aaow(
    grpc_client, bucket_name, object_name, generation=None, flush_interval_bytes=None
):
    """
    Creates and opens the AsyncAppendableObjectWriter.
    """
    writer_options = {}
    # Only pass flush_interval_bytes if the user explicitly provided a
    # non-default flush interval.
    if flush_interval_bytes and flush_interval_bytes != _DEFAULT_FLUSH_INTERVAL_BYTES:
        writer_options["FLUSH_INTERVAL_BYTES"] = flush_interval_bytes
    writer = AsyncAppendableObjectWriter(
        client=grpc_client,
        bucket_name=bucket_name,
        object_name=object_name,
        generation=generation,
        writer_options=writer_options,
    )
    await writer.open()
    return writer


async def close_mrd(mrd):
    """
    Closes the AsyncMultiRangeDownloader gracefully.
    Logs a warning if closing fails, instead of raising an exception.
    """
    if mrd:
        try:
            await mrd.close()
        except Exception as e:
            logger.warning(
                f"Error closing AsyncMultiRangeDownloader for {mrd.bucket_name}/{mrd.object_name}: {e}"
            )


async def close_aaow(aaow, finalize_on_close=False):
    """
    Closes the AsyncAppendableObjectWriter gracefully.
    Logs a warning if closing fails, instead of raising an exception.
    """
    if aaow:
        try:
            await aaow.close(finalize_on_close=finalize_on_close)
        except Exception as e:
            logger.warning(
                f"Error closing AsyncAppendableObjectWriter for {aaow.bucket_name}/{aaow.object_name}: {e}"
            )


class PartialView:
    """A bounded memory writer providing robust overfill/underfill constraint validations."""

    def __init__(self, parent, start_offset, expected_size):
        self.parent = parent
        self.start_offset = start_offset
        self.expected_size = expected_size
        self.current_offset = 0
        self._view_lock = threading.Lock()

    def write(self, data):
        """
        Schedules a write operation to memory mapping.
        """
        if not isinstance(data, bytes):
            raise ValueError(f"Expected bytes, but got {type(data)}")

        size = len(data)
        with self._view_lock:
            if self.current_offset + size > self.expected_size:
                error_msg = (
                    f"Attempted to write {size} bytes "
                    f"at offset {self.current_offset}. "
                    f"Max capacity is {self.expected_size} bytes."
                )
                raise BufferError(error_msg)

            abs_offset = self.start_offset + self.current_offset
            self.current_offset += size

        return self.parent._submit_write(abs_offset, data, size)

    def close(self):
        """
        Validates boundaries enforcing complete local payload consistency.
        """
        if self.current_offset < self.expected_size:
            error_msg = (
                f"Expected {self.expected_size} bytes, "
                f"but only received {self.current_offset} bytes. "
                f"Buffer contains uninitialized data."
            )
            raise BufferError(error_msg)


class DirectMemmoveBuffer:
    """
    A buffer-like object that writes data directly to memory asynchronously.

    This class provides an interface that queues `ctypes.memmove` operations
    to a thread pool executor. It provides synchronous backpressure: if `max_pending`
    operations are currently writing, the `write()` call will safely block the
    calling thread (e.g., an asyncio loop) until capacity frees up.

    Memory allocation is natively deferred. If the payload precisely aligns
    with expected bounds sequentially, it gracefully overrides manual memmoves
    using true Zero-Copy payload replacement safely under the hood.

    Note: This class is now strictly Thread-Safe
    """

    THRESHOLD_BYTES_FOR_SCHEDULING = 128 * 1024

    def __init__(self, expected_size, executor, max_pending=5):
        """
        Initializes the DirectMemmoveBuffer.

        Args:
            expected_size (int): The total amount of bytes expected to populate memory.
            executor (concurrent.futures.Executor): The thread pool executor to run the
                memmove operations. The lifecycle of this executor is managed by the caller.
            max_pending (int, optional): The maximum number of pending write operations
                allowed in the queue. Defaults to 5.
        """
        self.expected_size = expected_size
        self.executor = executor

        # Volatile state variables. Must only be amended while holding self._lock.
        self._pending_count = 0
        self._error = None
        self._total_bytes_written = 0
        self._stop_accepting_writes = False
        self._is_closed = False

        # Track allocated (start, end) intervals to prevent overlapping views.
        self._allocated_intervals = []

        # PyBytes Native Pointers & Allocation tracking natively handled
        self._result_bytes = None
        self._start_address = None

        # Primitives:
        # 1. semaphore: Provides backpressure by limiting the number of active tasks.
        # 2. _lock: Protects mutations to the volatile state variables above.
        # 3. _done_event: Signals when the queue of active background tasks reaches zero.
        self.semaphore = threading.Semaphore(max_pending)
        self._lock = threading.Lock()
        self._done_event = threading.Event()
        self._done_event.set()

    def get_view(self, offset, size):
        """Constructs secure mapped offset references correctly handling constraint layouts."""
        if offset < 0 or offset + size > self.expected_size:
            raise ValueError("Invalid view requested: exceeds physical boundaries!")

        start = offset
        end = offset + size

        with self._lock:
            if self._stop_accepting_writes or self._is_closed:
                raise ValueError("Cannot get view on a closed/closing buffer.")

            # Enforce Write-Once memory semantics: prevent overlapping views
            for a_start, a_end in self._allocated_intervals:
                if max(start, a_start) < min(end, a_end):
                    raise ValueError(
                        f"Overlapping view requested: [{start}, {end}) "
                        f"overlaps with already allocated view [{a_start}, {a_end})"
                    )

            self._allocated_intervals.append((start, end))

        return PartialView(self, offset, size)

    def _decrement_pending(self):
        """Helper to cleanly release concurrency primitives after a task finishes."""
        self.semaphore.release()
        with self._lock:
            self._pending_count -= 1
            if self._pending_count == 0:
                self._done_event.set()

    def _submit_write(self, dest_offset, data_bytes, size):
        if size == 0:
            with self._lock:
                if self._stop_accepting_writes or self._is_closed:
                    raise ValueError("I/O operation on closed buffer.")
                if self._error:
                    raise self._error

            fut = concurrent.futures.Future()
            fut.set_result(None)
            return fut

        self.semaphore.acquire()

        try:
            with self._lock:
                if self._stop_accepting_writes or self._is_closed:
                    raise ValueError("I/O operation on closed buffer.")

                if self._error:
                    raise self._error

                if self._result_bytes is None:
                    if dest_offset == 0 and size == self.expected_size:
                        # fastpath: return buffer directly
                        self._result_bytes = data_bytes
                        self.semaphore.release()  # Release because we skip the executor
                        fut = concurrent.futures.Future()
                        fut.set_result(None)
                        self._total_bytes_written += size
                        return fut
                    if HAS_CPYTHON_API:
                        self._result_bytes = PyBytes_FromStringAndSize(
                            None, self.expected_size
                        )
                        self._start_address = PyBytes_AsString(self._result_bytes)
                    else:
                        self._result_bytes = bytearray(self.expected_size)
                        self._start_address = (
                            -1
                        )  # Dummy value to pass the defensive check below

                # Defensive programming: gracefully catch internal overwrite attempts
                if self._start_address is None:
                    raise BufferError(
                        "Attempted to execute standard write over a Zero-Copied payload."
                    )

                if self._pending_count == 0:
                    self._done_event.clear()
                self._pending_count += 1

        except BaseException:
            self.semaphore.release()
            raise

        if size <= self.THRESHOLD_BYTES_FOR_SCHEDULING:
            # Fast path, no need to send it to executor
            try:
                self._do_memmove(dest_offset, data_bytes, size)
            except BaseException:
                # The exception is already captured in self._error by _do_memmove
                pass

            fut = concurrent.futures.Future()
            local_err = self._error
            if local_err:
                fut.set_exception(local_err)
            else:
                fut.set_result(None)
            return fut
        else:
            try:
                # Slow path, schedule it on executor.
                return self.executor.submit(
                    self._do_memmove, dest_offset, data_bytes, size
                )
            except BaseException as e:
                with self._lock:
                    self._error = e
                self._decrement_pending()
                raise e

    def _do_memmove(self, dest_offset, data_bytes, size):
        try:
            with self._lock:
                if self._error:
                    return

            # Isolate pointer math to CPython only.
            # PyPy uses memory-safe native slice assignment.
            if HAS_CPYTHON_API:
                dest = self._start_address + dest_offset
                ctypes.memmove(dest, data_bytes, size)
            else:
                memoryview(self._result_bytes)[
                    dest_offset : dest_offset + size
                ] = data_bytes

            with self._lock:
                self._total_bytes_written += size

        except BaseException as e:
            with self._lock:
                if self._error is None:
                    self._error = e
            raise
        finally:
            self._decrement_pending()

    def get_value(self):
        with self._lock:
            if self._error:
                raise self._error
            if not self._is_closed:
                raise RuntimeError("Buffer is still not closed yet!")
            if self._result_bytes is None and self.expected_size == 0:
                return b""
            if self._total_bytes_written < self.expected_size:
                raise BufferError(
                    f"Buffer incomplete: Expected {self.expected_size} bytes but "
                    f"only populated {self._total_bytes_written}. Returning this "
                    f"payload would leak uninitialized memory."
                )

            if not isinstance(self._result_bytes, bytes):
                return bytes(self._result_bytes)

            return self._result_bytes

    def close(self):
        """
        Locks the buffer preventing further incoming writes, waits for all pending
        write operations to complete, and checks for errors.
        """
        with self._lock:
            self._stop_accepting_writes = True

        self._done_event.wait()
        with self._lock:
            self._is_closed = True
            if self._error:
                raise self._error


async def _close_mrds(mrds, raise_exception=False):
    """Close a list of MRDs asynchronously."""
    if not mrds:
        return
    results = await asyncio.gather(
        *(mrd.close() for mrd in mrds), return_exceptions=True
    )
    for r in results:
        if isinstance(r, Exception):
            if raise_exception:
                raise r
            logger.warning("Error closing MRD: %s", r)


class MRDPool:
    """Manages a pool of AsyncMultiRangeDownloader objects with on-demand scaling.

    When constructed by `MRDPoolCache`, the instance acts as a pool over a shared
    MRD queue and donates its MRDs back to that queue on close.
    """

    def __init__(
        self,
        gcsfs,
        bucket_name,
        object_name,
        generation,
        finalized,
        pool_size,
        cache=None,
    ):
        self.gcsfs = gcsfs
        self.bucket_name = bucket_name
        self.object_name = object_name
        self.generation = generation
        self._cache = cache
        self._key = (bucket_name, object_name, generation)
        self.pool_size = pool_size
        self._free_mrds = asyncio.Queue(maxsize=pool_size)
        self._active_count = 0
        self._lock = asyncio.Lock()
        self.details = None
        self.persisted_size = None
        self.finalized = finalized
        self._initialized = False
        self._closed = False

        self._all_mrds = []
        self._rr_index = 0
        # Maps each checked-out AsyncMultiRangeDownloader to its number of active
        # get_mrd() holders. An MRD is only requeued into _free_mrds (or closed,
        # when the pool is closing) by whichever holder releases it LAST, so an
        # MRD still being driven by a round-robin sharer is never closed/requeued
        # out from under it.
        self._inflight = {}

    def _mark_inflight(self, mrd):
        """Record one more holder of `mrd`. Called under self._lock while the MRD
        is handed to exactly one get_mrd() caller."""
        self._inflight[mrd] = self._inflight.get(mrd, 0) + 1

    def _release_inflight(self, mrd):
        """Drop one holder of `mrd`; return True iff this was the LAST holder
        (so the caller must now requeue or close it).

        Both helpers only do synchronous dict mutations with no `await`, so they
        are atomic under asyncio even though get_mrd's finally runs WITHOUT
        self._lock."""
        count = self._inflight.get(mrd, 0) - 1
        if count > 0:
            self._inflight[mrd] = count
            return False
        self._inflight.pop(mrd, None)
        return True

    async def _create_mrd(self):
        await self.gcsfs._get_grpc_client()
        mrd = await init_mrd(
            self.gcsfs.grpc_client, self.bucket_name, self.object_name, self.generation
        )
        return mrd

    async def _get_or_create_mrd(self):
        """Gets an MRD from the cache or creates a new one."""
        mrd = None
        if self._cache is not None:
            mrd = self._cache.get_idle_mrd(self._key)
        if mrd is None:
            mrd = await self._create_mrd()
        self._all_mrds.append(mrd)
        return mrd

    async def initialize(self):
        """Initializes the MRDPool by creating the first downloader instance."""
        async with self._lock:
            if self._closed:
                raise RuntimeError("Cannot initialize a closed MRDPool.")

            if not self._initialized and self._active_count == 0:
                if self.finalized:
                    mrd = await self._get_or_create_mrd()
                else:
                    # Always create a new MRD for unfinalized objects to get the up-to-date persisted_size
                    mrd = await self._create_mrd()
                    self._all_mrds.append(mrd)
                self.persisted_size = mrd.persisted_size
                self._free_mrds.put_nowait(mrd)
                self._active_count += 1

            self._initialized = True

    @contextlib.asynccontextmanager
    async def get_mrd(self):
        """
        Dynamically provisions MRDs using an async context manager.

        If a downloader is available in the pool, it is yielded immediately. If the
        pool is empty but hasn't reached `pool_size`, a new downloader is spawned
        on demand or fetched from the cache. Automatically returns the downloader
        to the free queue upon exit.

        Yields:
            AsyncMultiRangeDownloader: An active downloader ready for requests.

        Raises:
            Exception: Bubbles up any exceptions encountered during MRD creation.
        """
        mrd = None

        async with self._lock:
            if self._closed:
                raise RuntimeError("MRDPool is closed.")

            if self._free_mrds.empty():
                if self._active_count < self.pool_size:
                    self._active_count += 1
                    try:
                        mrd = await self._get_or_create_mrd()
                    except BaseException as e:
                        self._active_count -= 1
                        raise e
                elif self._all_mrds:
                    # Pool is full and the queue is empty: share a busy MRD in
                    # round-robin fashion. The MRD now has multiple holders;
                    # refcounting ensures it is requeued/closed only once the
                    # LAST holder is done with it.
                    mrd = self._all_mrds[self._rr_index]
                    self._rr_index = (self._rr_index + 1) % len(self._all_mrds)

            if mrd is None:
                # If the queue was non-empty, this gets an MRD immediately without blocking.
                # If the queue was empty (pool is full and sharing is disabled), this blocks
                # until a holder returns an MRD.
                # NOTE: the lock is intentionally held across this await -- get_mrd's finally
                # returns MRDs via put_nowait WITHOUT the lock, so a waiter blocked
                # here is still unblocked by a concurrent release (no deadlock).
                mrd = await self._free_mrds.get()

            self._mark_inflight(mrd)

        try:
            yield mrd
        finally:
            # Intentionally lock-free (see note above). Only the holder that
            # releases the MRD last requeues or closes it, so a round-robin
            # sharer is never torn down by a peer or by close().
            if self._release_inflight(mrd):
                if self._closed:
                    await close_mrd(mrd)
                else:
                    self._free_mrds.put_nowait(mrd)

    async def close(self):
        """
        Cleanly shut down all MRDs.

        Iterates through all instantiated downloaders and releases them back to
        the cache if available, otherwise closes them.

        In-flight MRDs are not touched here; the last get_mrd() holder closes them on return once _closed is set.
        """
        async with self._lock:
            if self._closed:
                return
            self._closed = True

            free_mrds = []
            while not self._free_mrds.empty():
                free_mrds.append(self._free_mrds.get_nowait())

            try:
                if self._cache is not None:
                    await self._cache.release(self._key, free_mrds)
                else:
                    await _close_mrds(free_mrds, raise_exception=True)
            finally:
                self._all_mrds.clear()


def _drain_queue(q):
    if q is None:
        return []
    items = list(q)
    q.clear()
    return items


class MRDPoolCache:
    """Filesystem-level cache of MRD pools.

    Keyed by (bucket, object, generation). Idle pools are kept in an LRU cache
    and evicted when exceeding `max_idle_pools`.

    Lifecycle:
    1. `get()` returns an `MRDPool`.
    2. When the pool is closed, it returns its MRDs to this cache via `release()`.
    3. When a key's refcount hits zero, it becomes eligible for LRU eviction.
    """

    def __init__(self, gcsfs, max_idle_pools: int = 16, max_queue_size: int = 8):
        """
        Initializes the MRDPoolCache.

        Args:
            gcsfs (ExtendedGcsFileSystem): The filesystem instance.
            max_idle_pools (int, optional): Maximum number of idle pools to retain. Defaults to 16.
            max_queue_size (int, optional): Maximum number of idle MRDs per key. Defaults to 8.
        """
        self._gcsfs = weakref.ref(gcsfs)
        self._max_idle_pools = max_idle_pools
        self._max_queue_size = max_queue_size
        self._mrd_queues = {}
        self._refcounts = {}
        self._evictable_keys = collections.OrderedDict()
        self._closed = False

    def get_idle_mrd(self, key):
        """Gets an MRD from the queue for the given key."""
        if self._closed:
            return None
        queue = self._mrd_queues.get(key)
        if queue:
            return queue.popleft()
        return None

    def _incref(self, key):
        """Mark `key` as in use: ensure its queue exists, bump refcount,
        and remove the key from the evictable set so it can't be LRU'd out
        while a caller still holds the pool.
        """
        if key not in self._mrd_queues:
            self._mrd_queues[key] = collections.deque()
        self._refcounts[key] = self._refcounts.get(key, 0) + 1
        self._evictable_keys.pop(key, None)

    def _decref(self, key):
        """Release one reference on `key`. When the last reference goes,
        mark the key evictable and run LRU eviction. Returns MRDs whose
        keys were evicted and must be closed by the caller.
        """
        refcount = self._refcounts.get(key, 0) - 1
        if refcount > 0:
            self._refcounts[key] = refcount
            return []

        self._refcounts.pop(key, None)
        if self._closed:
            return []

        self._evictable_keys[key] = None
        mrds_to_close = []
        while len(self._evictable_keys) > self._max_idle_pools:
            evict_key, _ = self._evictable_keys.popitem(last=False)
            mrds_to_close.extend(_drain_queue(self._mrd_queues.pop(evict_key, None)))
        return mrds_to_close

    async def get(self, bucket_name, object_name, generation, pool_size):
        """
        Gets an MRDPool for the specified object.

        Args:
            bucket_name (str): Name of the bucket.
            object_name (str): Name of the object.
            generation (int): Object generation.
            pool_size (int): Requested pool size.

        Returns:
            MRDPool: An initialized MRDPool instance.
        """
        if self._closed:
            raise RuntimeError("MRDPoolCache is closed.")
        fs = self._gcsfs()
        if fs is None:
            raise RuntimeError("ExtendedGcsFileSystem has been garbage collected.")

        info = await fs._info(f"{bucket_name}/{object_name}", generation=generation)
        if generation is None:
            generation = info.get("generation")
        key = (bucket_name, object_name, generation)
        finalized = info.get("timeFinalized") is not None

        self._incref(key)
        mrd_pool = MRDPool(
            fs,
            bucket_name,
            object_name,
            generation,
            finalized,
            pool_size,
            cache=self,
        )
        if info is not None:
            mrd_pool.details = info

        try:
            await mrd_pool.initialize()
        except BaseException:
            # Init failed. `mrd_pool.close()` donates any partial MRDs back
            # via release() and drops the refcount we just took. If that was
            # the last reference, purge the key entirely.
            await mrd_pool.close()
            mrds_to_close = []
            if key not in self._refcounts:
                self._evictable_keys.pop(key, None)
                mrds_to_close = _drain_queue(self._mrd_queues.pop(key, None))
            await _close_mrds(mrds_to_close, raise_exception=False)
            raise

        return mrd_pool

    async def release(self, key, mrds):
        """
        Releases MRDs back to the cache or closes them if necessary.

        Args:
            key (tuple): Cache key (bucket, object, generation).
            mrds (list): List of MRDs to release.
        """
        mrds_to_close = []
        mrd_queue = self._mrd_queues.get(key)
        if mrd_queue is not None:
            for mrd in mrds:
                if len(mrd_queue) < self._max_queue_size:
                    mrd_queue.append(mrd)
                else:
                    mrds_to_close.append(mrd)
        else:
            mrds_to_close.extend(mrds)

        mrds_to_close.extend(self._decref(key))
        await _close_mrds(mrds_to_close, raise_exception=False)

    async def close(self):
        """
        Closes the cache and all pooled MRDs.
        """
        if self._closed:
            return
        mrds_to_close = []
        for q in self._mrd_queues.values():
            mrds_to_close.extend(_drain_queue(q))
        self._mrd_queues.clear()
        self._refcounts.clear()
        self._evictable_keys.clear()
        self._closed = True
        await _close_mrds(mrds_to_close, raise_exception=True)
