import asyncio
import contextlib
import logging
import os
import uuid
import weakref
from concurrent.futures import ThreadPoolExecutor
from enum import Enum
from glob import has_magic

import aiohttp
import fsspec
from fsspec import asyn
from fsspec.callbacks import NoOpCallback
from google.api_core import exceptions as api_exceptions
from google.api_core.client_info import ClientInfo
from google.api_core.client_options import ClientOptions
from google.auth.credentials import AnonymousCredentials
from google.cloud import storage_control_v2
from google.cloud.storage.asyncio.async_appendable_object_writer import (
    AsyncAppendableObjectWriter,
)
from google.cloud.storage.asyncio.async_grpc_client import AsyncGrpcClient
from google.cloud.storage.asyncio.async_multi_range_downloader import (
    AsyncMultiRangeDownloader,
)

from gcsfs import __version__ as version
from gcsfs import zb_hns_utils
from gcsfs._dircache import HnsDirCacheUpdater
from gcsfs.concurrency import split_range
from gcsfs.core import GCSFile, GCSFileSystem
from gcsfs.retry import DEFAULT_RETRY_CONFIG, get_storage_control_retry_config
from gcsfs.zb_hns_utils import DirectMemmoveBuffer, MRDPool
from gcsfs.zonal_file import ZonalFile

logger = logging.getLogger("gcsfs")

USER_AGENT = "python-gcsfs"
STORAGE_CONTROL_RPC_TIMEOUT = 30.0


class BucketType(Enum):
    ZONAL_HIERARCHICAL = "ZONAL_HIERARCHICAL"
    HIERARCHICAL = "HIERARCHICAL"
    NON_HIERARCHICAL = "NON_HIERARCHICAL"
    UNKNOWN = "UNKNOWN"


gcs_file_types = {
    BucketType.ZONAL_HIERARCHICAL: ZonalFile,
    BucketType.NON_HIERARCHICAL: GCSFile,
    BucketType.HIERARCHICAL: GCSFile,
    BucketType.UNKNOWN: GCSFile,
}


@contextlib.asynccontextmanager
async def _get_mrd_from_pool_or_mrd(mrd_or_pool):
    """
    Helper function to yield an AsyncMultiRangeDownloader
    whether a single instance or an MRDPool is provided.
    """
    if isinstance(mrd_or_pool, MRDPool):
        async with mrd_or_pool.get_mrd() as m:
            yield m
    elif isinstance(mrd_or_pool, AsyncMultiRangeDownloader):
        yield mrd_or_pool
    else:
        raise TypeError(
            f"Expected MRDPool or AsyncMultiRangeDownloader, got {type(mrd_or_pool)}"
        )


async def _get_mrd_size(mrd_or_pool):
    """Helper to extract the persisted_size from either a pool or a single MRD."""
    if mrd_or_pool is None:
        return None
    async with _get_mrd_from_pool_or_mrd(mrd_or_pool) as m:
        return m.persisted_size


class ExtendedGcsFileSystem(HnsDirCacheUpdater, GCSFileSystem):
    """
    This class will be used when GCSFS_EXPERIMENTAL_ZB_HNS_SUPPORT env variable is set to true.
    ExtendedGcsFileSystem is a subclass of GCSFileSystem that adds new logic for bucket types
    including zonal and hierarchical. For buckets without special properties, it forwards requests
    to the parent class GCSFileSystem for default processing.
    """

    def __init__(
        self,
        *args,
        finalize_on_close=False,
        mrd_pool_cache_size=16,
        max_mrd_pool_cache_queue_size=8,
        **kwargs,
    ):
        """
        Parameters
        ----------
        finalize_on_close : bool, default False
            By default, files in zonal buckets are left unfinalized to allow appends.
        mrd_pool_cache_size : int, default 16
            Maximum number of idle pools to retain in the cache.
        max_mrd_pool_cache_queue_size : int, default 8
            Maximum number of idle MRDs per key in the cache.
        **kwargs : dict
            Additional arguments passed to GCSFileSystem.
            Supports retry configuration overrides for Storage Control API:
            - retry_timeout: Total time to spend retrying (seconds).
            - retry_initial: Initial delay between retries (seconds).
            - retry_maximum: Maximum delay between retries (seconds).
            - retry_multiplier: Multiplier for delay between retries.
            These map to `google.api_core.retry.AsyncRetry` arguments (without 'retry_' prefix).
        """
        valid_keys = DEFAULT_RETRY_CONFIG.keys()
        self.retry_config = {
            k[6:]: v
            for k, v in kwargs.items()
            if k.startswith("retry_") and k[6:] in valid_keys and v is not None
        }
        super().__init__(*args, **kwargs)
        # By default, files in zonal buckets are left unfinalized to allow appends.
        self.finalize_on_close = finalize_on_close
        self._grpc_client = None
        self._storage_control_client = None
        # Adds user-passed credentials to ExtendedGcsFileSystem to pass to gRPC/Storage Control clients.
        # We unwrap the nested credentials here because self.credentials is a GCSFS wrapper,
        # but the clients expect the underlying google.auth credentials object.
        self.credential = self.credentials.credentials
        # When token="anon", self.credentials.credentials is None. This is
        # often used for testing with emulators. However, the gRPC and storage
        # control clients require a credentials object for initialization.
        # We explicitly use AnonymousCredentials() to allow unauthenticated access.
        if self.credentials.token == "anon":
            self.credential = AnonymousCredentials()
        self._storage_layout_cache = {}
        self._memmove_executor = ThreadPoolExecutor(
            max_workers=kwargs.get("memmove_max_workers", 8)
        )
        weakref.finalize(self, self._memmove_executor.shutdown)
        self._mrd_pool_cache = zb_hns_utils.MRDPoolCache(
            self,
            max_idle_pools=mrd_pool_cache_size,
            max_queue_size=max_mrd_pool_cache_queue_size,
        )
        weakref.finalize(
            self,
            self._finalize_mrd_pool_cache,
            self.loop,
            self._mrd_pool_cache,
        )

    async def _get_threshold_for_disk_reads(self, bucket):
        if await self._is_zonal_bucket(bucket):
            return (
                5 * 1024 * 1024
            )  # Thanks to our in house, zero copy DirectMemmoveBuffer
        return await super()._get_threshold_for_disk_reads(bucket)

    @staticmethod
    def _finalize_mrd_pool_cache(loop, cache):
        """Tear down the MRDPoolCache when ExtendedGcsFileSystem is garbage collected."""
        if cache is None or getattr(cache, "_closed", False):
            return

        try:
            current_loop = asyncio.get_running_loop()
        except RuntimeError:
            current_loop = None

        if loop and loop.is_running():
            asyncio.run_coroutine_threadsafe(cache.close(), loop)
        elif current_loop is not None and current_loop.is_running():
            asyncio.run_coroutine_threadsafe(cache.close(), current_loop)
        elif asyn.loop[0] is not None and asyn.loop[0].is_running():
            try:
                asyn.sync(asyn.loop[0], cache.close, timeout=5.0)
            except fsspec.FSTimeoutError:
                pass

    @property
    def _user_project(self):
        """Value used for billing - enabling "requestor pays" access"""
        if self.requester_pays:
            return (
                self.requester_pays
                if isinstance(self.requester_pays, str)
                else self.project
            )
        return None

    def _get_retry_config(self, **kwargs):
        return get_storage_control_retry_config(self.retry_config, **kwargs)

    @property
    def grpc_client(self):
        if self.asynchronous and self._grpc_client is None:
            raise RuntimeError(
                "Please await _get_grpc_client() before accessing grpc_client"
            )
        if self._grpc_client is None:
            self._grpc_client = asyn.sync(self.loop, self._get_grpc_client)
        return self._grpc_client

    async def _get_grpc_client(self):
        if self._grpc_client is None:
            client_options = ClientOptions(quota_project_id=self._user_project)
            if self._location:
                # client_options expects only the host:port, without any protocol or path components.
                endpoint = self._location.split("://")[-1].split("/")[0]
                client_options.api_endpoint = endpoint
            self._grpc_client = AsyncGrpcClient(
                credentials=self.credential,
                client_info=ClientInfo(user_agent=f"{USER_AGENT}/{version}"),
                client_options=client_options,
            )
        return self._grpc_client

    async def _get_control_plane_client(self):
        if self._storage_control_client is None:

            # Initialize the storage control plane client for bucket
            # metadata operations
            transport_cls = (
                storage_control_v2.StorageControlAsyncClient.get_transport_class(
                    "grpc_asyncio"
                )
            )
            channel_kwargs = {
                "credentials": self.credential,
                "options": [("grpc.primary_user_agent", f"{USER_AGENT}/{version}")],
                "quota_project_id": self._user_project,
            }
            if self._location:
                # Extract host:port safely (strips protocol and trailing URL paths if any).
                endpoint = self._location.split("://")[-1].split("/")[0]
                channel_kwargs["host"] = endpoint

            channel = transport_cls.create_channel(**channel_kwargs)

            transport = transport_cls(channel=channel)
            self._storage_control_client = storage_control_v2.StorageControlAsyncClient(
                transport=transport
            )
        return self._storage_control_client

    async def _close_resources(self):
        """
        Close gRPC clients, channels, and other resources.

        Order matters: pooled MRDs ride on the gRPC channel, so the MRD pool
        cache must be drained BEFORE the gRPC transport is closed. The storage
        control client owns a separate channel and is independent.
        """
        if self._mrd_pool_cache is not None:
            try:
                await self._mrd_pool_cache.close()
            except Exception as e:
                logger.warning(f"Failed to close MRDPoolCache: {e}")
        if self._storage_control_client is not None:
            try:
                await self._storage_control_client.transport.close()
            except Exception as e:
                logger.warning(f"Failed to close storage_control_client: {e}")
            self._storage_control_client = None
        if self._grpc_client is not None:
            try:
                await self._grpc_client.grpc_client.transport.close()
            except Exception as e:
                logger.warning(f"Failed to close grpc_client: {e}")
            self._grpc_client = None

    async def _lookup_bucket_type(self, bucket):
        if bucket in self._storage_layout_cache:
            return self._storage_layout_cache[bucket]
        bucket_type = await self._get_bucket_type(bucket)
        # Don't cache UNKNOWN type.
        # This ensures that subsequent operations will retry the lookup,
        # allowing it to recover when the transient error resolves.
        if bucket_type == BucketType.UNKNOWN:
            return bucket_type
        self._storage_layout_cache[bucket] = bucket_type
        return self._storage_layout_cache[bucket]

    _sync_lookup_bucket_type = asyn.sync_wrapper(_lookup_bucket_type)

    async def _get_bucket_type(self, bucket):
        try:
            client = await self._get_control_plane_client()
            bucket_name_value = f"projects/_/buckets/{bucket}/storageLayout"
            logger.debug(f"get_storage_layout request for name: {bucket_name_value}")
            response = await client.get_storage_layout(
                name=bucket_name_value,
                retry=self._get_retry_config(),
                timeout=STORAGE_CONTROL_RPC_TIMEOUT,
            )

            if response.location_type == "zone":
                return BucketType.ZONAL_HIERARCHICAL
            if (
                response.hierarchical_namespace
                and response.hierarchical_namespace.enabled
            ):
                return BucketType.HIERARCHICAL
            return BucketType.NON_HIERARCHICAL
        except api_exceptions.NotFound:
            logger.warning(
                f"Error: Bucket {bucket} not found or you lack permissions for "
                f"storage layout api used to detect bucket type. Falling back to GCSFileSystem."
            )
            return BucketType.UNKNOWN
        except Exception as e:
            logger.warning(
                f"Could not determine bucket type for bucket name {bucket}: {e}, falling back to GCSFileSystem"
            )
            # Default to UNKNOWN in case bucket type is not obtained
            return BucketType.UNKNOWN

    def _open(
        self,
        path,
        mode="rb",
        block_size=None,
        cache_options=None,
        acl=None,
        consistency=None,
        metadata=None,
        autocommit=True,
        fixed_key_metadata=None,
        generation=None,
        **kwargs,
    ):
        """
        Open a file.
        """
        bucket, _, _ = self.split_path(path)
        bucket_type = self._sync_lookup_bucket_type(bucket)

        return gcs_file_types[bucket_type](
            self,
            path,
            mode,
            block_size=block_size or self.default_block_size,
            cache_options=cache_options,
            consistency=consistency or self.consistency,
            metadata=metadata,
            acl=acl,
            autocommit=autocommit,
            fixed_key_metadata=fixed_key_metadata,
            generation=generation,
            finalize_on_close=kwargs.pop("finalize_on_close", self.finalize_on_close),
            **kwargs,
        )

    # Replacement method for _process_limits to support new params (offset and length) for MRD.
    async def _process_limits_to_offset_and_length(
        self, path, start, end, file_size=None
    ):
        """
        Calculates the read offset and length from start and end parameters.

        Args:
            path (str): The path to the file.
            start (int | None): The starting byte position.
            end (int | None): The ending byte position.
            file_size (int | None): The total size of the file. If None, it will be fetched via _info().

        Returns:
            tuple: A tuple containing (offset, length).
        """
        size = file_size

        async def _get_size():
            nonlocal size
            if size is None:
                size = (await self._info(path))["size"]
            return size

        if start is None:
            offset = 0
        elif start < 0:
            offset = max(0, await _get_size() + start)
        else:
            offset = start

        if end is None:
            effective_end = await _get_size()
        elif end < 0:
            effective_end = await _get_size() + end
        else:
            effective_end = end

        # If the requested end is before/ same as the start, return empty.
        if effective_end <= offset:
            return offset, 0
        else:
            length = effective_end - offset  # Normal case
            s = await _get_size()
            if effective_end > s:
                length = max(0, s - offset)  # Clamp and ensure non-negative

        return offset, length

    sync_process_limits_to_offset_and_length = asyn.sync_wrapper(
        _process_limits_to_offset_and_length
    )

    async def _is_zonal_bucket(self, bucket):
        bucket_type = await self._lookup_bucket_type(bucket)
        return bucket_type == BucketType.ZONAL_HIERARCHICAL

    async def _fetch_range_split(
        self,
        path,
        start,
        chunk_lengths,
        concurrency,
        mrd=None,
        size=None,
        **kwargs,
    ):
        """
        Reading multiple adjacent ranges concurrently.

        Delegates concurrent fetching of individual chunks directly to `_cat_file`.
        """
        file_size = size or await _get_mrd_size(mrd)
        if file_size is None:
            logger.warning(
                f"AsyncMultiRangeDownloader (MRD) for {path} has no 'persisted_size'. "
                "Falling back to _info() to get the file size."
            )
            file_size = (await self._info(path))["size"]

        start_offset = start if start is not None else 0
        if start_offset >= file_size or start_offset + sum(chunk_lengths) > file_size:
            raise RuntimeError("Request not satisfiable.")

        pool_created_here = False
        bucket, object_name, generation = self.split_path(path)

        if mrd is None:
            # If no mrd is provided, we create one with pool size equal to passed concurrency.
            pool_size = min(len(chunk_lengths), concurrency)
            mrd = await self._mrd_pool_cache.get(
                bucket, object_name, generation, pool_size=pool_size
            )
            pool_created_here = True

        tasks = []
        try:
            current_offset = start_offset

            cat_kwargs = kwargs.copy()

            for length in chunk_lengths:
                end_offset = current_offset + length
                tasks.append(
                    asyncio.create_task(
                        self._cat_file(
                            path,
                            start=current_offset,
                            end=end_offset,
                            mrd=mrd,
                            # Distribute the concurrency budget proportionally.
                            # Since these outer tasks are already concurrent, this is typically 1.
                            # However, if a large chunk dominates the total size, it receives
                            # higher concurrency to prevent it from becoming a bottleneck.
                            concurrency=max(
                                1, length * concurrency // sum(chunk_lengths)
                            ),
                            **cat_kwargs,
                        )
                    )
                )
                current_offset = end_offset

            results = await asyncio.gather(*tasks, return_exceptions=True)

            # Bubble up any exceptions encountered during concurrent fetching
            for res in results:
                if isinstance(res, Exception):
                    raise res

            return results
        except BaseException:
            for t in tasks:
                if not t.done():
                    t.cancel()
            await asyncio.gather(*tasks, return_exceptions=True)
            raise
        finally:
            if pool_created_here:
                await mrd.close()

    async def _concurrent_mrd_fetch(self, offset, length, concurrency, mrd_or_pool):
        """Helper to handle concurrent chunk downloads cleanly."""
        ranges = split_range(length, concurrency, self.MIN_CHUNK_SIZE_FOR_CONCURRENCY)

        tasks = []
        views = []
        has_error = False

        # The master buffer manages its own allocation under the hood
        master_buffer = DirectMemmoveBuffer(length, self._memmove_executor)

        async def _download(o, s, view, mrd_or_pool):
            async with _get_mrd_from_pool_or_mrd(mrd_or_pool) as m_client:
                if logger.isEnabledFor(logging.DEBUG):
                    logger.debug(
                        f"mrd path: {m_client.object_name} | "
                        f"Requested range: [({o}, {s})]"
                    )
                await m_client.download_ranges([(o, s, view)])

        for relative_offset, actual_size in ranges:
            part_offset = offset + relative_offset

            # Give each task a restricted view of the master buffer
            view = master_buffer.get_view(part_offset - offset, actual_size)
            views.append(view)

            tasks.append(
                asyncio.create_task(
                    _download(part_offset, actual_size, view, mrd_or_pool)
                )
            )

        try:
            results = await asyncio.gather(*tasks, return_exceptions=True)
            for res in results:
                if isinstance(res, Exception):
                    has_error = True
                    raise res
            for view in views:
                view.close()
        except BaseException:
            has_error = True
            for t in tasks:
                if not t.done():
                    t.cancel()
            await asyncio.gather(*tasks, return_exceptions=True)
            raise
        finally:
            try:
                master_buffer.close()
            except Exception:
                # If we are already handling a network/download exception,
                # ignore the exception from buffer (which is just a symptom of the drop).
                # If there's no download error, this means the buffer logic
                # itself failed, so we must surface the error.
                if not has_error:
                    raise

        return master_buffer.get_value()

    async def _cat_file(
        self,
        path,
        start=None,
        end=None,
        mrd=None,
        concurrency=zb_hns_utils.DEFAULT_CONCURRENCY,
        **kwargs,
    ):
        """Fetch a file's contents as bytes, with an optimized path for Zonal buckets.

        This method overrides the parent `_cat_file` to read objects in Zonal buckets using gRPC.

        Args:
            path (str): The full GCS path to the file (e.g., "bucket/object").
            start (int, optional): The starting byte position to read from.
            end (int, optional): The ending byte position to read to.
            mrd (AsyncMultiRangeDownloader, MRDPool, optional): An existing multi-range
                downloader instance or a pool of MRD. If not provided, a new one will be created for Zonal buckets.
            concurrency (int, optional): The max number of concurrent request to fetch the data.

        Returns:
            bytes: The content of the file or file range.
        """
        pool_created_here = False

        # A new MRDPool is required when read is done directly by the
        # GCSFilesystem class without creating a GCSFile object first.
        if mrd is None:
            bucket, object_name, generation = self.split_path(path)
            if not await self._is_zonal_bucket(bucket):
                # Fall back to default implementation if not a zonal bucket
                return await super()._cat_file(
                    path, start=start, end=end, concurrency=concurrency, **kwargs
                )

            # Instantiate an MRDPool locally for this call
            mrd = await self._mrd_pool_cache.get(
                bucket, object_name, generation, pool_size=concurrency
            )
            pool_created_here = True

        try:
            file_size = await _get_mrd_size(mrd)
            if file_size is None:
                logger.warning(
                    f"AsyncMultiRangeDownloader (MRD) for {path} has no 'persisted_size'. "
                    "Falling back to _info() to get the file size. "
                    "This may result in incorrect behavior for unfinalized objects."
                )
                file_size = (await self._info(path))["size"]

            offset, length = await self._process_limits_to_offset_and_length(
                path, start, end, file_size
            )

            if length == 0:
                return b""

            return await self._concurrent_mrd_fetch(
                offset,
                length,
                concurrency,
                mrd,
            )

        finally:
            # If we created a temporary pool specifically for this _cat_file call, clean it up
            if pool_created_here:
                await mrd.close()

    async def _is_bucket_hns_enabled(self, bucket):
        """Checks if a bucket has Hierarchical Namespace enabled."""
        try:
            bucket_type = await self._lookup_bucket_type(bucket)
        except Exception as e:
            logger.warning(
                f"Could not determine if bucket '{bucket}' is HNS-enabled, falling back to default non-HNS: {e}",
                stack_info=True,
            )
            return False

        return bucket_type in [BucketType.ZONAL_HIERARCHICAL, BucketType.HIERARCHICAL]

    async def _mv(self, path1, path2, **kwargs):
        """
        Move a file or directory. Overrides the parent `_mv` to provide an
        optimized, atomic implementation for renaming folders and moving files
        in HNS-enabled buckets. Falls back to the parent's object-level
        copy-and-delete implementation for non-HNS buckets.
        """
        if path1 == path2:
            logger.debug(
                "%s mv: The paths are the same, so no files/directories were moved.",
                self,
            )
            return

        if (
            isinstance(path1, list)
            or isinstance(path2, list)
            or (isinstance(path1, str) and has_magic(path1))
        ):
            return await super()._mv(path1, path2, **kwargs)

        bucket1, key1, _ = self.split_path(path1)
        bucket2, key2, _ = self.split_path(path2)

        is_hns = await self._is_bucket_hns_enabled(bucket1)

        if not is_hns:
            logger.debug(
                f"Not an HNS bucket. Falling back to object-level mv for '{path1}' to '{path2}'."
            )
            return await super()._mv(path1, path2, **kwargs)

        try:
            info1 = await self._info(path1)
            is_folder = info1.get("type") == "directory"

            # We only use HNS rename if the source is a folder and the move is
            # within the same bucket.
            if is_folder and bucket1 == bucket2 and key1:
                logger.debug(
                    f"Using HNS-aware folder rename for '{path1}' to '{path2}'."
                )
                source_folder_name = f"projects/_/buckets/{bucket1}/folders/{key1}"
                destination_folder_id = key2 or key1.rstrip("/").split("/")[-1]

                request = storage_control_v2.RenameFolderRequest(
                    name=source_folder_name,
                    destination_folder_id=destination_folder_id,
                    request_id=str(uuid.uuid4()),
                )

                logger.debug(f"rename_folder request: {request}")
                client = await self._get_control_plane_client()
                operation = await client.rename_folder(
                    request=request,
                    retry=self._get_retry_config(),
                    timeout=STORAGE_CONTROL_RPC_TIMEOUT,
                )
                await operation.result()
                self._update_dircache_after_rename(path1, path2)

                logger.debug(
                    "Successfully renamed folder from '%s' to '%s'", path1, path2
                )
                return
            elif not is_folder:
                await self._mv_file(path1, path2)
                return
        except Exception as e:
            if isinstance(e, FileNotFoundError):
                # If the source doesn't exist, fail fast.
                raise
            if isinstance(e, api_exceptions.NotFound):
                raise FileNotFoundError(
                    f"Source '{path1}' not found for move operation."
                ) from e
            if isinstance(e, api_exceptions.Conflict):
                # This occurs if the destination folder already exists.
                # Raise FileExistsError for fsspec compatibility.
                raise FileExistsError(
                    f"HNS rename failed due to conflict for '{path1}' to '{path2}'"
                ) from e
            if isinstance(e, api_exceptions.FailedPrecondition):
                raise OSError(f"HNS rename failed: {e}") from e

            logger.warning(f"Could not perform HNS-aware mv: {e}")

        logger.debug(f"Falling back to object-level mv for '{path1}' to '{path2}'.")
        return await super()._mv(path1, path2, **kwargs)

    mv = asyn.sync_wrapper(_mv)

    async def _list_objects(self, path, prefix="", versions=False, **kwargs):
        try:
            return await super()._list_objects(
                path, prefix=prefix, versions=versions, **kwargs
            )
        except FileNotFoundError:
            bucket, key, _ = self.split_path(path)
            if key and await self._is_bucket_hns_enabled(bucket):
                try:
                    await self._get_directory_info(path, bucket, key, None)
                    return []
                except (FileNotFoundError, Exception):
                    pass
            raise

    async def _mkdir(
        self,
        path,
        create_parents=False,
        enable_hierarchical_namespace=False,
        placement=None,
        location=None,
        **kwargs,
    ):
        """
        Create a directory or bucket.

        If the path refers to a bucket (no object key), a new bucket is created.
        If the path refers to a directory (includes object key), a directory is created.

        Parameters
        ----------
        path : str
            Path to create.
        create_parents : bool
            If True, create parent directories if they do not exist.
            If the path includes a bucket that does not exist, the bucket will also be created.
        enable_hierarchical_namespace : bool
            If True, and a bucket is being created, the bucket will have Hierarchical
            Namespace (HNS) enabled.
        placement : str, optional
            If set to a zone (e.g. "us-central1-a"), a Zonal bucket is created.
            Zonal buckets are HNS-enabled by default.
            When creating a Zonal bucket, `location` must be passed as a
            region (e.g. "us-central1"). If `location` is not specified, it defaults
            to `self.default_location`. The zone specified in `placement` must belong
            to the region specified in `location`.
        location : str, optional
            Location where buckets are created, like 'US' or 'EUROPE-WEST3'.
            If not provided, defaults to `self.default_location`.
        **kwargs : dict
            Additional arguments passed to the bucket creation API.

        Notes
        -----
        - For HNS-enabled buckets (including Zonal buckets), this method creates a
          native folder object.
        - If `create_parents` is False and a parent directory does not exist in an
          HNS/Zonal bucket, a FileNotFoundError is raised.
        - For non-HNS buckets, this falls back to the parent implementation. Since
          standard GCS has no true directories, `mkdir` on a path with a key is
          typically a no-op unless `create_parents=True` triggers bucket creation.
        """
        path = self._strip_protocol(path)
        bucket, key, _ = self.split_path(path)

        # Determine if we are requesting creation of a Zonal or HNS bucket
        should_create_zonal_bucket = placement is not None
        should_create_hns_bucket = (
            enable_hierarchical_namespace or should_create_zonal_bucket
        )

        # Prepare arguments for bucket creation
        bucket_kwargs = kwargs.copy()
        if location:
            bucket_kwargs["location"] = location
        if should_create_zonal_bucket:
            bucket_kwargs["customPlacementConfig"] = {"dataLocations": [placement]}
            bucket_kwargs["storageClass"] = "RAPID"

        if should_create_hns_bucket:
            bucket_kwargs["hierarchicalNamespace"] = {"enabled": True}
            # HNS buckets require uniform bucket-level access.
            bucket_kwargs["iamConfiguration"] = {
                "uniformBucketLevelAccess": {"enabled": True}
            }
            # When uniformBucketLevelAccess is enabled, ACLs cannot be used.
            # We must explicitly set them to None to prevent the parent
            # method from using default ACLs.
            bucket_kwargs["acl"] = None
            bucket_kwargs["default_acl"] = None

        # Case 1: Path is just a bucket
        if not key:
            return await super()._mkdir(
                path, create_parents=create_parents, **bucket_kwargs
            )

        # Case 2: Path is a folder
        is_hns_bucket = False

        # If creating parents and HNS/Zonal requested, ensure bucket exists with correct config
        if create_parents and should_create_hns_bucket:
            if not await self._exists(bucket):
                await super()._mkdir(bucket, create_parents=True, **bucket_kwargs)
                is_hns_bucket = True

        if not is_hns_bucket:
            is_hns_bucket = await self._is_bucket_hns_enabled(bucket)

        if is_hns_bucket:
            return await self._create_hns_folder(path, bucket, key, create_parents)

        return await super()._mkdir(
            path, create_parents=create_parents, **bucket_kwargs
        )

    async def _create_hns_folder(
        self, path, bucket, key, create_parents, exist_ok=True
    ):
        logger.debug(f"Using HNS-aware mkdir for '{path}'.")

        # Preemptively check if the path already exists to ensure fsspec compatibility.
        # This check is required because GCS natively allows a file and a folder
        # with the exact same name to co-exist at the same level. Consequently, the GCS
        # server's create_folder API will succeed silently and not raise a Conflict
        # if a file already exists at the target path. To enforce standard POSIX/fsspec
        # filesystem semantics (forbidding file/directory name collisions), we need to
        # explicitly verify client-side that no file occupies the target path.
        try:
            info = await self._info(path)
            if info["type"] != "directory":
                raise FileExistsError(f"A file already exists at the path: {path}")
            if not exist_ok:
                raise FileExistsError(f"Directory already exists: {path}")
            return
        except FileNotFoundError:
            pass

        parent = f"projects/_/buckets/{bucket}"
        folder_id = key.rstrip("/") + "/"
        request = storage_control_v2.CreateFolderRequest(
            parent=parent,
            folder_id=folder_id,
            recursive=create_parents,
            request_id=str(uuid.uuid4()),
        )
        try:
            logger.debug(f"create_folder request: {request}")
            client = await self._get_control_plane_client()
            await client.create_folder(
                request=request,
                retry=self._get_retry_config(),
                timeout=STORAGE_CONTROL_RPC_TIMEOUT,
            )
            # Instead of invalidating the parent cache, update it to add the new entry.
            self._cache_add_entry(
                self._parent(path),
                self._directory_cache_entry(path, key.rstrip("/")),
            )
        except api_exceptions.Conflict as e:
            logger.debug(f"Conflict detected for path: {path}: {e}")
            # Under race conditions, folder might have been created concurrently
            if not exist_ok:
                raise FileExistsError(f"Directory already exists: {path}") from e
        except api_exceptions.FailedPrecondition as e:
            # This error can occur if create_parents=False and the parent dir doesn't exist.
            # Translate it to FileNotFoundError for fsspec compatibility.
            raise FileNotFoundError(
                f"mkdir for '{path}' failed due to a precondition error: {e}"
            ) from e

    mkdir = asyn.sync_wrapper(_mkdir)

    async def _makedirs(self, path, exist_ok=False):
        """Recursively make directories.

        For HNS-enabled buckets, this natively creates the folder objects recursively.
        For standard GCS buckets, since folders are simulated, this is a no-op.

        Parameters
        ----------
        path : str
            Leaf directory name.
        exist_ok : bool (False)
            If False, will error if the target directory already exists.
            If True, will succeed silently if it exists as a directory, but
            will still error if it exists as a file.

        Raises
        ------
        FileNotFoundError
            If the target bucket does not exist.
        FileExistsError
            If `exist_ok` is False and the directory already exists, or if a
            file already exists at the target path regardless of `exist_ok`.
        """
        path = self._strip_protocol(path)
        bucket, key, _ = self.split_path(path)

        if bucket in ["", "/"]:
            raise ValueError("Cannot create root bucket")

        # First, check if the bucket exists
        if not await self._exists(bucket):
            raise FileNotFoundError(f"Bucket does not exist: {bucket}")

        if not key:
            # Bucket-only case: makedirs is not for bucket creation.
            # Since the bucket exists, raise FileExistsError if not exist_ok.
            if not exist_ok:
                raise FileExistsError(f"Bucket already exists: {bucket}")
            return

        # Check if the bucket is HNS-enabled
        is_hns = await self._is_bucket_hns_enabled(bucket)
        if not is_hns:
            # For non-HNS buckets, directory creation is a no-op
            return

        # Recursively create the native folders in the HNS bucket
        await self._create_hns_folder(
            path, bucket, key, create_parents=True, exist_ok=exist_ok
        )

    makedirs = asyn.sync_wrapper(_makedirs)

    async def _get_directory_info(self, path, bucket, key, generation):
        """
        Override to use Storage Control API's get_folder for HNS buckets.
        For HNS, we avoid calling _ls (_list_objects) entirely.
        """
        is_hns = await self._is_bucket_hns_enabled(bucket)

        # If bucket is HNS, use get folder metadata api to determine a directory
        if is_hns:
            try:
                # folder_id is the path relative to the bucket
                folder_id = key.rstrip("/")
                folder_resource_name = (
                    f"projects/_/buckets/{bucket}/folders/{folder_id}"
                )

                request = storage_control_v2.GetFolderRequest(
                    name=folder_resource_name, request_id=str(uuid.uuid4())
                )

                # Verify existence using get_folder API
                client = await self._get_control_plane_client()
                response = await client.get_folder(
                    request=request,
                    retry=self._get_retry_config(),
                    timeout=STORAGE_CONTROL_RPC_TIMEOUT,
                )

                # If successful, return directory metadata
                return {
                    "bucket": bucket,
                    "name": path,
                    "size": 0,
                    "storageClass": "DIRECTORY",
                    "type": "directory",
                    "ctime": response.create_time,
                    "mtime": response.update_time,
                    "metageneration": response.metageneration,
                }
            except api_exceptions.NotFound:
                # If get_folder fails, the folder does not exist.
                raise FileNotFoundError(path)
            except Exception as e:
                # Log unexpected errors
                logger.error(f"Error fetching folder metadata for {path}: {e}")
                raise e

        # Fallback to standard GCS behavior for non-HNS buckets
        return await super()._get_directory_info(path, bucket, key, generation)

    async def _rmdir(self, path):
        """
        Deletes an empty directory. Overrides the parent `_rmdir` to delete
        empty directories in HNS-enabled buckets.
        """
        path = self._strip_protocol(path)
        bucket, key, _ = self.split_path(path)

        # The parent _rmdir is only for deleting buckets. If key is empty,
        # given path is a bucket, we can fall back.
        if not key:
            return await super()._rmdir(path)

        is_hns = await self._is_bucket_hns_enabled(bucket)
        if not is_hns:
            return await super()._rmdir(path)

        # In HNS buckets, a placeholder object (e.g., 'a/b/c/') might exist,
        # which would cause rmdir('a/b/c') to fail because the directory is not empty.
        # To handle this, we first attempt to delete the placeholder object.
        # If it doesn't exist, _rm_file will raise a FileNotFoundError, which we can
        # safely ignore and proceed with the directory deletion.
        #
        # Note: This may delete the placeholder even if the directory contains
        # other files and the final `delete_folder` call fails. This side
        # effect is acceptable because placeholder objects are used to simulate
        # folders and are not strictly necessary in HNS-enabled buckets, which
        # have native folder entities.
        try:
            placeholder_path = f"{path.rstrip('/')}/"

            await self._rm_file(placeholder_path)
            logger.debug(
                f"Removed placeholder object '{placeholder_path}' before rmdir."
            )
        except FileNotFoundError:
            # This is expected if no placeholder object exists and can be safely ignored.
            pass

        try:
            logger.debug(f"Using HNS-aware rmdir for '{path}'.")
            folder_name = f"projects/_/buckets/{bucket}/folders/{key.rstrip('/')}"
            request = storage_control_v2.DeleteFolderRequest(
                name=folder_name,
                request_id=str(uuid.uuid4()),
            )

            logger.debug(f"delete_folder request: {request}")
            client = await self._get_control_plane_client()
            await client.delete_folder(
                request=request,
                retry=self._get_retry_config(),
                timeout=STORAGE_CONTROL_RPC_TIMEOUT,
            )

            # Remove the directory from the cache and from its parent's listing.
            self.dircache.pop(path, None)
            # Remove the deleted directory entry from the parent's listing.
            self._cache_drop_entries(self._parent(path), {path})
            return
        except api_exceptions.NotFound as e:
            # This can happen if the directory does not exist, or if the path
            # points to a file object instead of a directory.
            raise FileNotFoundError(f"rmdir failed for path: {path}: {e}") from e
        except api_exceptions.FailedPrecondition as e:
            # This can happen if the directory is not empty.
            raise OSError(
                f"Pre condition failed for rmdir for path: {path}: {e}"
            ) from e
        except Exception as e:
            logger.error(f"HNS rmdir: Failed to delete folder '{path}': {e}")
            raise

    rmdir = asyn.sync_wrapper(_rmdir)

    # TODO: This method is only added to be used in rm method, can be deleted once
    # rm method is integrated with recursive API
    async def _expand_path_with_details(
        self, path, recursive=False, maxdepth=None, detail=False, assume_literal=False
    ):
        """
        Expand path with details, similar to `_expand_path` but returning full details.

        This method is used in `_rm` to get a list of files and directories to delete.
        It handles special characters in filenames by avoiding re-globbing concrete paths
        found during recursive expansion.

        Parameters
        ----------
        path : str or list
            The path or list of paths to expand.
        recursive : bool, default False
            Whether to recursively expand directories.
        maxdepth : int, optional
            Maximum depth to traverse for expansion.
        detail : bool, default False
            If True, returns a list of dictionaries with file details;
            otherwise, returns a list of path strings.
        assume_literal : bool, default False
            If True, treats the paths as literal even if they contain magic characters.
            This is used in recursive calls to prevent re-globbing concrete paths.
        """
        if maxdepth is not None and maxdepth < 1:
            raise ValueError("maxdepth must be at least 1")

        if isinstance(path, str):
            return await self._expand_path_with_details(
                [path],
                recursive,
                maxdepth,
                detail=detail,
                assume_literal=assume_literal,
            )
        else:
            out = {} if detail else set()
            path = [self._strip_protocol(p) for p in path]

            for p in path:
                if not assume_literal and has_magic(p):
                    bit = await self._glob(p, maxdepth=maxdepth, detail=detail)
                    if detail:
                        out.update(bit)
                        bit_paths = list(bit.keys())
                    else:
                        bit_set = set(bit)
                        out |= bit_set
                        bit_paths = list(bit_set)

                    if recursive:
                        if maxdepth is not None and maxdepth <= 1:
                            continue
                        rec = await self._expand_path_with_details(
                            bit_paths,
                            recursive=recursive,
                            maxdepth=maxdepth - 1 if maxdepth is not None else None,
                            detail=detail,
                            assume_literal=True,
                        )
                        if detail:
                            for info in rec:
                                out[info["name"]] = info
                        else:
                            out |= set(rec)
                    continue
                elif recursive:
                    rec = await self._find(
                        p, maxdepth=maxdepth, withdirs=True, detail=detail
                    )
                    if detail:
                        out.update(rec)
                    else:
                        out |= set(rec)

                if p not in out:
                    if detail:
                        try:
                            info = await self._info(p)
                            out[p] = info
                        except (FileNotFoundError, OSError):
                            pass
                    elif recursive is False or (await self._exists(p)):
                        out.add(p)

            if detail:
                out = list(out.values())
            else:
                out = sorted(out)

        if not out:
            raise FileNotFoundError(path)
        return out

    async def _rm(self, path, recursive=False, maxdepth=None, batchsize=20):
        """
        Deletes files and directories.

        This method overrides the parent `_rm` to correctly handle directory
        deletion in HNS-enabled buckets. Input paths are grouped by bucket and
        each bucket is routed independently: HNS buckets use HNS-specific
        deletion, while non-HNS buckets fall back to the parent implementation.
        This means a single call may mix paths from HNS and flat buckets, and
        the per-bucket groups are deleted concurrently.

        For HNS buckets, it first expands the path to get a list of all files
        and directories, then categorizes them. It deletes files in batches
        and then deletes directories individually.

        Args:
            path (str or list): The path(s) to delete.
            recursive (bool): If True, deletes directories and their contents.
            maxdepth (int, optional): The maximum depth to traverse for deletion.
            batchsize (int): The number of files to delete in a single batch request.
        """
        # Normalize to a list so single-path and multi-path inputs share a
        # single code path.
        if isinstance(path, str):
            path = [path]

        if not path:
            return []

        # Group paths by bucket so each bucket is routed independently to the
        # correct deletion strategy (HNS vs. flat), while same-bucket paths are
        # still batched together.
        grouped = {}
        for p in path:
            bucket, _, _ = self.split_path(p)
            grouped.setdefault(bucket, []).append(p)

        # Buckets are independent of one another, so delete each group
        # concurrently.
        bucket_results = await asyn._run_coros_in_chunks(
            [
                self._rm_bucket_paths(
                    bucket,
                    bucket_paths,
                    recursive=recursive,
                    maxdepth=maxdepth,
                    batchsize=batchsize,
                )
                for bucket, bucket_paths in grouped.items()
            ],
            return_exceptions=True,
        )

        results = []
        succeeded = False
        for res in bucket_results:
            if isinstance(res, FileNotFoundError):
                # A bucket group matched nothing; tolerated unless every group
                # fails.
                continue
            if isinstance(res, Exception):
                raise res
            succeeded = True
            if res:
                results.extend(res)

        if not succeeded:
            raise FileNotFoundError(path)

        return results

    async def _rm_bucket_paths(
        self, bucket, path, recursive=False, maxdepth=None, batchsize=20
    ):
        """Helper method to handle the rm operation for paths within a single bucket."""
        if not await self._is_bucket_hns_enabled(bucket):
            return await super()._rm(
                path, recursive=recursive, maxdepth=maxdepth, batchsize=batchsize
            )

        paths = await self._expand_path_with_details(
            path, recursive=recursive, maxdepth=maxdepth, detail=True
        )

        # Separate files and directories based on their type. Directories are
        # sorted in reverse so they are deleted from the deepest first.
        files = list({p["name"] for p in paths if p["type"] == "file"})
        dirs = sorted(
            {p["name"] for p in paths if p["type"] == "directory"},
            reverse=True,
        )

        return await self._perform_rm(files, dirs, path, batchsize=batchsize)

    async def _perform_rm(self, files, dirs, path, batchsize):
        """
        Helper method to perform the deletion of files and directories.

        Args:
            files (list[str]): A list of file paths to delete.
            dirs (list[str]): A list of directory paths to delete, sorted from
                              deepest to shallowest.
            path (str): The original path for the rm operation, for error reporting.
            batchsize (int): The number of files to delete in a single batch request.

        Returns:
            list: A list of exceptions that occurred during deletion.
        """
        # If no files or directories were found to delete, raise FileNotFoundError.
        if not files and not dirs:
            raise FileNotFoundError(path)

        exs = await self._delete_files(files, batchsize)
        # For directories, we must delete them from deepest to shallowest
        # to avoid race conditions where a parent is deleted before its child.
        # We group directories by depth and delete them level by level.
        dirs_by_depth = {}
        for d in dirs:
            depth = d.count("/")
            dirs_by_depth.setdefault(depth, []).append(d)

        for depth in sorted(dirs_by_depth.keys(), reverse=True):
            level_dirs = dirs_by_depth[depth]
            results = await asyn._run_coros_in_chunks(
                [self._rmdir(d) for d in level_dirs],
                batch_size=batchsize,
                return_exceptions=True,
            )
            for res in results:
                if isinstance(res, Exception):
                    exs.append(res)

        errors = [
            ex
            for ex in exs
            if isinstance(ex, Exception)
            and not isinstance(ex, (FileNotFoundError, api_exceptions.NotFound))
            and "No such object" not in str(ex)
        ]

        if errors:
            raise errors[0]

        # Filter out non-critical "not found" errors from the final list.
        # A successful rm should return an empty list.
        return [
            e
            for e in exs
            if e is not None
            and not isinstance(e, (FileNotFoundError, api_exceptions.NotFound))
            and "No such object" not in str(e)
        ]

    rm = asyn.sync_wrapper(_rm)

    async def _find(
        self,
        path,
        withdirs=False,
        detail=False,
        prefix="",
        versions=False,
        maxdepth=None,
        **kwargs,
    ):
        """
        HNS-aware find. Overrides the parent to correctly list empty folders in HNS buckets.

        For HNS buckets this method uses a hybrid approach for fetching files and directories:
        1. It fetches all files in a single recursive API call (like the parent).
        2. It concurrently fetches all folder objects (including empty ones)
           using a recursive walk with the Storage Control API.
        3. It merges these results for a complete listing.

        For buckets with flat structure, it falls back to the parent implementation.
        """
        path = self._strip_protocol(path)
        if maxdepth is not None and maxdepth < 1:
            raise ValueError("maxdepth must be at least 1")
        bucket, _, _ = self.split_path(path)

        is_hns = await self._is_bucket_hns_enabled(bucket)
        if not is_hns:
            # Use parent implementation if bucket is not HNS.
            return await super()._find(
                path,
                withdirs=withdirs,
                detail=detail,
                prefix=prefix,
                versions=versions,
                maxdepth=maxdepth,
                **kwargs,
            )

        # Hybrid approach for HNS enabled buckets
        # 1. Fetch all files from super find() method by passing withdirs as False.
        # We pass maxdepth as None here to ensure we fetch all files for caching,
        # and then filter by maxdepth at the end of this method.
        files_task = asyncio.create_task(
            super()._find(
                path,
                withdirs=False,  # Fetch files only
                detail=True,  # Get full details for merging and populating cache
                prefix=prefix,
                versions=versions,
                maxdepth=None,
                update_cache=False,  # Defer caching until merging files and folders
                **kwargs,
            )
        )

        # 2. Fetch all folders recursively. This is necessary to find all folders,
        # especially empty ones.
        folders_task = asyncio.create_task(
            self._get_all_folders(path, bucket, prefix=prefix)
        )
        # 3. Run tasks concurrently and merge results.
        files_result, folders_result = await asyncio.gather(files_task, folders_task)

        # Always update the cache with both files and folders for consistency.
        cacheable_objects = list(files_result.values()) + folders_result
        self._get_dirs_and_update_cache(path, cacheable_objects, prefix=prefix)

        if not withdirs:
            # If not including directories, the final output should only contain files.
            all_objects = list(files_result.values())
        else:
            all_objects = cacheable_objects

        all_objects.sort(key=lambda o: o["name"])

        # Final filtering and formatting. `all_objects` now contains a complete
        # list of all files and folders.
        if maxdepth:
            depth = path.rstrip("/").count("/") + maxdepth
            all_objects = [o for o in all_objects if o["name"].count("/") <= depth]

        if detail:
            if versions:
                return {
                    (
                        f"{o['name']}#{o['generation']}"
                        if "generation" in o
                        else o["name"]
                    ): o
                    for o in all_objects
                }
            return {o["name"]: o for o in all_objects}

        if versions:
            return [
                f"{o['name']}#{o['generation']}" if "generation" in o else o["name"]
                for o in all_objects
            ]
        return [o["name"] for o in all_objects]

    async def _get_all_folders(self, path, bucket, prefix=""):
        """
        Recursively fetches all folder objects under a given path using the
        Storage Control API.
        """
        _, base_path, _ = self.split_path(path)
        base_path = "" if not base_path else base_path.rstrip("/") + "/"
        full_prefix = f"{base_path}{prefix}"

        # To find folders matching a partial prefix, we need to list their parent.
        if full_prefix and not full_prefix.endswith("/"):
            partition = full_prefix.rpartition("/")
            start_dir = partition[0] + partition[1]
        else:
            start_dir = full_prefix

        folders = []
        client = await self._get_control_plane_client()
        parent = f"projects/_/buckets/{bucket}"

        request = storage_control_v2.ListFoldersRequest(
            parent=parent, prefix=start_dir, request_id=str(uuid.uuid4())
        )
        logger.debug(f"list_folders request: {request}")

        try:
            async for folder in await client.list_folders(
                request=request,
                retry=self._get_retry_config(),
                timeout=STORAGE_CONTROL_RPC_TIMEOUT,
            ):
                entry = self._create_folder_entry(bucket, folder)
                _, key, _ = self.split_path(entry["name"])
                # Key from _create_folder_entry does not have a trailing slash
                key_with_slash = key + "/"

                if key.startswith(full_prefix) or key_with_slash.startswith(
                    full_prefix
                ):
                    folders.append(entry)
        except api_exceptions.NotFound:
            # If the start_dir itself doesn't exist, we just return empty
            pass

        return folders

    def _create_folder_entry(self, bucket, folder):
        """Helper to create a dictionary representing a folder entry."""
        path = f"{bucket}/{folder.name.split('/folders/')[1]}".rstrip("/")
        _, key, _ = self.split_path(path)
        return {
            "Key": key,
            "Size": 0,
            "name": path,
            "size": 0,
            "type": "directory",
            "storageClass": "DIRECTORY",
            "ctime": folder.create_time,
            "mtime": folder.update_time,
            "metageneration": folder.metageneration,
        }

    async def _put_file(
        self,
        lpath,
        rpath,
        metadata=None,
        consistency=None,
        content_type=None,
        chunksize=50 * 2**20,
        callback=None,
        fixed_key_metadata=None,
        mode="overwrite",
        **kwargs,
    ):
        """Upload a local file.

        This method is optimized for Zonal buckets, using gRPC for uploads.
        In zonal buckets, file is left *unfinalized* by default unless
        `finalize_on_close` is set to True.
        For non-Zonal buckets, it delegates to the parent class's implementation.

        Parameters
        ----------
        lpath: str
            Path to the local file to be uploaded.
        rpath: str
            Path on GCS to upload the file to.
        metadata: dict, optional
            Unsupported for Zonal buckets and will be ignored.
        consistency: str, optional
            Unsupported for Zonal buckets and will be ignored.
        content_type: str, optional
            Unsupported for Zonal buckets and will be ignored except for the default.
        chunksize: int, optional
            The size of chunks to upload data in.
        callback: fsspec.callbacks.Callback, optional
            Callback to monitor the upload progress.
        fixed_key_metadata: dict, optional
            Unsupported for Zonal buckets and will be ignored.
        mode: str, optional
            The write mode, either 'overwrite' or 'create'.
        """
        bucket, key, generation = self.split_path(rpath)
        if not await self._is_zonal_bucket(bucket):
            return await super()._put_file(
                lpath,
                rpath,
                metadata=metadata,
                consistency=consistency,
                content_type=content_type,
                chunksize=chunksize,
                callback=callback,
                fixed_key_metadata=fixed_key_metadata,
                mode=mode,
                **kwargs,
            )

        if os.path.isdir(lpath):
            return

        if generation:
            raise ValueError("Cannot write to specific object generation")

        if (
            metadata
            or fixed_key_metadata
            or consistency
            or (content_type and content_type != "application/octet-stream")
        ):
            logger.warning(
                "Zonal buckets do not support content_type, metadata, "
                "fixed_key_metadata or consistency during upload. "
                "These parameters will be ignored."
            )
        await self._get_grpc_client()
        # Works for both 'overwrite' and 'create' modes
        writer = await zb_hns_utils.init_aaow(self.grpc_client, bucket, key)

        try:
            with open(lpath, "rb") as f:
                await writer.append_from_file(f, block_size=chunksize)
        finally:
            finalize_on_close = kwargs.get("finalize_on_close", self.finalize_on_close)
            await zb_hns_utils.close_aaow(writer, finalize_on_close=finalize_on_close)

        await self._write_file_cache_update(rpath)

    async def _pipe_file(
        self,
        path,
        data,
        metadata=None,
        consistency=None,
        content_type="application/octet-stream",
        fixed_key_metadata=None,
        chunksize=50 * 2**20,
        mode="overwrite",
        **kwargs,
    ):
        """Upload bytes to a file.

        This method is optimized for Zonal buckets, using gRPC for uploads.
        In zonal buckets, file is left *unfinalized* by default unless
        `finalize_on_close` is set to True.
        For non-Zonal buckets, it delegates to the parent class's implementation.

        Parameters
        ----------
        path: str
            Path to the file to be written.
        data: bytes
            The content to write to the file.
        metadata: dict, optional
            Unsupported for Zonal buckets and will be ignored.
        consistency: str, optional
            Unsupported for Zonal buckets and will be ignored.
        content_type: str, optional
            Unsupported for Zonal buckets and will be ignored, except for the default.
        fixed_key_metadata: dict, optional
            Unsupported for Zonal buckets and will be ignored.
        chunksize: int, optional
            The size of chunks to upload data in.
        mode: str, optional
            The write mode, either 'overwrite' or 'create'.
        """
        bucket, key, generation = self.split_path(path)
        if not await self._is_zonal_bucket(bucket):
            return await super()._pipe_file(
                path,
                data,
                metadata=metadata,
                consistency=consistency,
                content_type=content_type,
                fixed_key_metadata=fixed_key_metadata,
                chunksize=chunksize,
                mode=mode,
            )

        if (
            metadata
            or fixed_key_metadata
            or (content_type and content_type != "application/octet-stream")
        ):
            logger.warning(
                "Zonal buckets do not support content_type, metadata or "
                "fixed_key_metadata during upload. These parameters will be ignored."
            )
        await self._get_grpc_client()
        # Works for both 'overwrite' and 'create' modes
        writer = await zb_hns_utils.init_aaow(self.grpc_client, bucket, key)
        try:
            with memoryview(data) as data_view:
                for i in range(0, len(data_view), chunksize):
                    await writer.append(data_view[i : i + chunksize])
        finally:
            finalize_on_close = kwargs.get("finalize_on_close", self.finalize_on_close)
            await zb_hns_utils.close_aaow(writer, finalize_on_close=finalize_on_close)

        await self._write_file_cache_update(path)

    async def _get_file_request(
        self, rpath, lpath, *args, headers=None, callback=None, **kwargs
    ):
        """
        Downloads a file from GCS to a local path.

        For Zonal buckets, it uses gRPC client for optimized downloads.
        For Standard buckets, it delegates to the parent class implementation.

        Parameters
        ----------
        rpath: str
            Path on GCS to download the file from.
        lpath: str
            Path to the local file to be downloaded.
        callback: fsspec.callbacks.Callback, optional
            Callback to monitor the download progress.
        **kwargs:
            For Zonal buckets, `chunksize` bytes (int) can be provided to control
            the download chunk size (default is 128KB).
        """
        bucket, key, path_generation = self.split_path(rpath)

        if not await self._is_zonal_bucket(bucket):
            return await super()._get_file_request(
                rpath, lpath, *args, headers=headers, callback=callback, **kwargs
            )

        if os.path.isdir(lpath):
            return

        # Coalesce generation: URL parsed vs kwargs
        generation = path_generation or kwargs.get("generation")
        callback = callback or NoOpCallback()

        mrd_pool = await self._mrd_pool_cache.get(bucket, key, generation, pool_size=1)
        try:
            async with mrd_pool.get_mrd() as mrd:
                size = mrd.persisted_size
                if size is None:
                    logger.warning(
                        f"AsyncMultiRangeDownloader (MRD) for {rpath} has no 'persisted_size'. "
                        "Falling back to _info() to get the file size. "
                        "This may result in incorrect behavior for unfinalized objects."
                    )
                    size = (await self._info(rpath, **kwargs)).get("size", 0)

                callback.set_size(size)

                lparent = os.path.dirname(lpath) or os.curdir
                os.makedirs(lparent, exist_ok=True)

                chunksize = kwargs.get("chunksize", 4096 * 32)  # 128KB default
                offset = 0

                with open(lpath, "wb") as f2:
                    while True:
                        if offset >= size:
                            break

                        data = await zb_hns_utils.download_range(
                            offset=offset, length=chunksize, mrd=mrd
                        )
                        if not data:
                            break

                        f2.write(data)
                        offset += len(data)
                        callback.relative_update(len(data))

                if offset != size:
                    raise aiohttp.ClientError(
                        f"Expected {size} bytes, but only received {offset} bytes"
                    )
        except Exception as e:
            # Clean up the corrupted file before raising error
            if os.path.exists(lpath):
                os.remove(lpath)
            raise e
        finally:
            await mrd_pool.close()

    async def _get_file_concurrent(
        self,
        rpath,
        lpath,
        concurrency,
        chunk_size,
        max_prefetch_size,
        headers=None,
        callback=None,
        fetcher_fn=None,
        **kwargs,
    ):
        bucket, key, path_generation = self.split_path(rpath)

        # Delegate standard buckets to the original implementation
        if not await self._is_zonal_bucket(bucket):
            return await super()._get_file_concurrent(
                rpath,
                lpath,
                concurrency,
                chunk_size,
                max_prefetch_size,
                headers=headers,
                callback=callback,
                fetcher_fn=fetcher_fn,
                **kwargs,
            )

        generation = path_generation or kwargs.get("generation")

        # Initialize the MRDPool once for this concurrent operation
        mrd_pool = await self._mrd_pool_cache.get(
            bucket, key, generation, pool_size=concurrency
        )

        # Define a custom fetcher that passes the pool to _cat_file
        async def custom_fetcher(start, size, split_factor=1):
            return await self._cat_file(
                rpath,
                start=start,
                end=start + size,
                mrd=mrd_pool,  # Inject the shared pool here
                concurrency=split_factor,
                headers=headers,
                **kwargs,
            )

        try:
            return await super()._get_file_concurrent(
                rpath,
                lpath,
                concurrency,
                chunk_size,
                max_prefetch_size,
                headers=headers,
                callback=callback,
                fetcher_fn=custom_fetcher,  # Pass our custom fetcher up the chain
                **kwargs,
            )
        finally:
            # Ensure the pool is closed when the download completes or fails
            await mrd_pool.close()

    async def _do_list_objects(
        self,
        path,
        max_results=None,
        delimiter="/",
        prefix="",
        versions=False,
        **kwargs,
    ):
        """
        Lists objects in a bucket.

        For HNS-enabled buckets, it sets `includeFoldersAsPrefixes` to True
        when the delimiter is '/'.
        """
        bucket, _, _ = self.split_path(path)
        if await self._is_bucket_hns_enabled(bucket) and delimiter == "/":
            kwargs["includeFoldersAsPrefixes"] = "true"

        return await super()._do_list_objects(
            path,
            max_results=max_results,
            delimiter=delimiter,
            prefix=prefix,
            versions=versions,
            **kwargs,
        )

    async def _cp_file(self, path1, path2, acl=None, **kwargs):
        """Duplicate remote file.

        For Standard GCS buckets, falls back to the parent class's implementation

        Zonal Bucket Support:
        Server-side copy is currently NOT supported for Zonal buckets because
        the `RewriteObject` API is unavailable for them.

        The following scenarios will raise a `NotImplementedError`:
        * Intra-zonal: Copying within the same Zonal bucket.
        * Inter-zonal: Copying between two different Zonal buckets.
        * Mixed: Copying between a Zonal bucket and a Standard bucket.

        """
        b1, _, _ = self.split_path(path1)
        b2, _, _ = self.split_path(path2)

        is_zonal_source, is_zonal_dest = await asyncio.gather(
            self._is_zonal_bucket(b1), self._is_zonal_bucket(b2)
        )

        # 1. Standard -> Standard (Delegate to core implementation)
        if not is_zonal_source and not is_zonal_dest:
            return await super()._cp_file(path1, path2, acl=acl, **kwargs)

        # 2. Zonal Scenarios (Currently Unsupported)
        raise NotImplementedError(
            "Server-side copy involving Zonal buckets is not supported. "
            "Zonal objects do not support rewrite."
        )

    async def _merge(self, path, paths, acl=None):
        """Concatenate objects within a single bucket.

        For Standard GCS buckets, falls back to the parent class's implementation

        Zonal Bucket Support:
        Server-side compose is currently NOT supported for Zonal buckets.
        """
        bucket, _, _ = self.split_path(path)

        if await self._is_zonal_bucket(bucket):
            raise NotImplementedError(
                "Server-side compose/merge is not supported for Zonal buckets."
            )

        return await super()._merge(path, paths, acl=acl)

    merge = asyn.sync_wrapper(_merge)


async def upload_chunk(fs, location, data, offset, size, content_type):
    """
    Uploads a chunk of data using AsyncAppendableObjectWriter for zonal buckets.
    Finalizes the upload when the total uploaded data size reaches the specified size.
    Delegates to core upload_chunk implementation for Non-Zonal buckets.
    """
    # If `location` is an HTTP resumable-upload URL (string), delegate to core upload_chunk
    # for Standard buckets.
    if isinstance(location, (str, bytes)):
        from gcsfs.core import upload_chunk as core_upload_chunk

        return await core_upload_chunk(fs, location, data, offset, size, content_type)

    if not isinstance(location, AsyncAppendableObjectWriter):
        raise TypeError(
            "upload_chunk for Zonal buckets expects an AsyncAppendableObjectWriter instance."
        )

    if offset or size or content_type:
        logger.warning(
            "Zonal buckets do not support offset, or content_type during upload. These parameters will be ignored."
        )

    if not location._is_stream_open:
        raise ValueError("Writer is closed. Please initiate a new upload.")

    try:
        await location.append(data)
    except Exception as e:
        logger.error(
            f"Error uploading chunk at offset {location.offset}: {e}. Closing stream."
        )
        # Don't finalize the upload on error
        await zb_hns_utils.close_aaow(location, finalize_on_close=False)
        raise

    if (location.offset or 0) >= size:
        logger.debug("Uploaded data is equal or greater than size. Finalizing upload.")
        await zb_hns_utils.close_aaow(location, finalize_on_close=True)


async def initiate_upload(
    fs,
    bucket,
    key,
    content_type="application/octet-stream",
    metadata=None,
    fixed_key_metadata=None,
    mode="overwrite",
    kms_key_name=None,
):
    """
    Initiates an upload for Zonal buckets by creating an AsyncAppendableObjectWriter.
    Delegates to core initiate_upload implementation for Non-Zonal buckets.

    Parameters
    ----------
    fs: GCSFileSystem
        The GCS filesystem instance.
    bucket: str
        The target bucket name.
    key: str
        The target object key.
    content_type: str, optional
        Unsupported for Zonal buckets and will be ignored, except for the default.
    metadata: dict, optional
        Unsupported for Zonal buckets and will be ignored.
    fixed_key_metadata: dict, optional
        Unsupported for Zonal buckets and will be ignored.
    mode: str, optional
        The write mode, either 'overwrite' or 'create'.
    kms_key_name: str, optional
        Unsupported for Zonal buckets and will be ignored.
    """
    if not await fs._is_zonal_bucket(bucket):
        from gcsfs.core import initiate_upload as core_initiate_upload

        return await core_initiate_upload(
            fs,
            bucket,
            key,
            content_type,
            metadata,
            fixed_key_metadata,
            mode,
            kms_key_name,
        )

    if (
        metadata
        or fixed_key_metadata
        or kms_key_name
        or (content_type and content_type != "application/octet-stream")
    ):
        logger.warning(
            "Zonal buckets do not support content_type, metadata, fixed_key_metadata, "
            "or kms_key_name during upload. These parameters will be ignored."
        )

    await fs._get_grpc_client()
    # If generation is not passed to init_aaow, it creates a new object and overwrites if object already exists.
    # Hence it works for both 'overwrite' and 'create' modes.
    return await zb_hns_utils.init_aaow(fs.grpc_client, bucket, key)


async def simple_upload(
    fs,
    bucket,
    key,
    datain,
    metadatain=None,
    consistency=None,
    content_type="application/octet-stream",
    fixed_key_metadata=None,
    mode="overwrite",
    kms_key_name=None,
    **kwargs,
):
    """
    Performs a simple, single-request upload to Zonal bucket using gRPC.
    In zonal buckets, file is left *unfinalized* by default unless
    `finalize_on_close` is set to True.
    Delegates to core simple_upload implementation for Non-Zonal buckets.

    Parameters
    ----------
    fs: GCSFileSystem
        The GCS filesystem instance.
    bucket: str
        The target bucket name.
    key: str
        The target object key.
    datain: bytes
        The data to be uploaded.
    metadatain: dict, optional
        Unsupported for Zonal buckets and will be ignored.
    consistency: str, optional
        Unsupported for Zonal buckets and will be ignored.
    content_type: str, optional
        Unsupported for Zonal buckets and will be ignored, except for the default.
    fixed_key_metadata: dict, optional
        Unsupported for Zonal buckets and will be ignored.
    mode: str, optional
        The write mode, either 'overwrite' or 'create'.
    kms_key_name: str, optional
        Unsupported for Zonal buckets and will be ignored.
    """
    if not await fs._is_zonal_bucket(bucket):
        from gcsfs.core import simple_upload as core_simple_upload

        return await core_simple_upload(
            fs,
            bucket,
            key,
            datain,
            metadatain,
            consistency,
            content_type,
            fixed_key_metadata,
            mode,
            kms_key_name,
        )

    if (
        metadatain
        or fixed_key_metadata
        or kms_key_name
        or consistency
        or (content_type and content_type != "application/octet-stream")
    ):
        logger.warning(
            "Zonal buckets do not support content_type, metadatain, fixed_key_metadata, "
            "consistency or kms_key_name during upload. These parameters will be ignored."
        )
    await fs._get_grpc_client()
    # If generation is not passed to init_aaow, it creates a new object and overwrites if object already exists.
    # Hence it works for both 'overwrite' and 'create' modes.
    writer = await zb_hns_utils.init_aaow(fs.grpc_client, bucket, key)
    try:
        await writer.append(datain)
    finally:
        default_finalize = getattr(fs, "finalize_on_close", False)
        finalize_on_close = kwargs.get("finalize_on_close", default_finalize)
        await zb_hns_utils.close_aaow(writer, finalize_on_close=finalize_on_close)
