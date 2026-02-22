import asyncio
import logging
import os
from enum import Enum
from functools import partial
from glob import has_magic
from io import BytesIO

from fsspec import asyn
from fsspec.callbacks import NoOpCallback
from google.api_core import exceptions as api_exceptions
from google.api_core import gapic_v1
from google.api_core.client_info import ClientInfo
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
from gcsfs.core import GCSFile, GCSFileSystem
from gcsfs.zonal_file import ZonalFile

logger = logging.getLogger("gcsfs")

USER_AGENT = "python-gcsfs"


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


class ExtendedGcsFileSystem(GCSFileSystem):
    """
    This class will be used when GCSFS_EXPERIMENTAL_ZB_HNS_SUPPORT env variable is set to true.
    ExtendedGcsFileSystem is a subclass of GCSFileSystem that adds new logic for bucket types
    including zonal and hierarchical. For buckets without special properties, it forwards requests
    to the parent class GCSFileSystem for default processing.
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
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
            self._grpc_client = AsyncGrpcClient(
                credentials=self.credential,
                client_info=ClientInfo(user_agent=f"{USER_AGENT}/{version}"),
            )
        return self._grpc_client

    async def _get_control_plane_client(self):
        if self._storage_control_client is None:

            # Initialize the storage control plane client for bucket
            # metadata operations
            client_info = gapic_v1.client_info.ClientInfo(
                user_agent=f"{USER_AGENT}/{version}"
            )
            # The HNS RenameFolder operation began failing with an "input/output error"
            # after an authentication library change caused it to send a
            # `quota_project_id` from application default credentials. The
            # RenameFolder API rejects requests with this parameter.
            #
            # This workaround explicitly removes the `quota_project_id` to prevent
            # the API from rejecting the request. A long-term fix is in progress
            # in the GCS backend to relax this restriction.
            #
            # TODO: Remove this workaround once the GCS backend fix is deployed.
            creds = self.credential
            if hasattr(creds, "with_quota_project"):
                creds = creds.with_quota_project(None)

            self._storage_control_client = storage_control_v2.StorageControlAsyncClient(
                credentials=creds, client_info=client_info
            )
        return self._storage_control_client

    async def _lookup_bucket_type(self, bucket):
        if bucket in self._storage_layout_cache:
            return self._storage_layout_cache[bucket]
        bucket_type = await self._get_bucket_type(bucket)
        # Dont cache UNKNOWN type
        if bucket_type == BucketType.UNKNOWN:
            return bucket_type
        self._storage_layout_cache[bucket] = bucket_type
        return self._storage_layout_cache[bucket]

    _sync_lookup_bucket_type = asyn.sync_wrapper(_lookup_bucket_type)

    async def _get_bucket_type(self, bucket):
        try:
            await self._get_control_plane_client()
            bucket_name_value = f"projects/_/buckets/{bucket}/storageLayout"
            logger.debug(f"get_storage_layout request for name: {bucket_name_value}")
            response = await self._storage_control_client.get_storage_layout(
                name=bucket_name_value
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
            logger.warning(f"Error: Bucket {bucket} not found or you lack permissions.")
            return BucketType.UNKNOWN
        except Exception as e:
            logger.error(
                f"Could not determine bucket type for bucket name {bucket}: {e}"
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

        Raises:
            ValueError: If the calculated range is invalid.
        """
        size = file_size

        if start is None:
            offset = 0
        elif start < 0:
            size = (await self._info(path))["size"] if size is None else size
            offset = size + start
        else:
            offset = start

        if end is None:
            size = (await self._info(path))["size"] if size is None else size
            effective_end = size
        elif end < 0:
            size = (await self._info(path))["size"] if size is None else size
            effective_end = size + end
        else:
            effective_end = end

        if offset < 0:
            raise ValueError(f"Calculated start offset ({offset}) cannot be negative.")
        if effective_end < offset:
            raise ValueError(
                f"Calculated end position ({effective_end}) cannot be before start offset ({offset})."
            )
        elif effective_end == offset:
            length = 0  # Handle zero-length slice
        else:
            length = effective_end - offset  # Normal case
            size = (await self._info(path))["size"] if size is None else size
            if effective_end > size:
                length = max(0, size - offset)  # Clamp and ensure non-negative

        return offset, length

    sync_process_limits_to_offset_and_length = asyn.sync_wrapper(
        _process_limits_to_offset_and_length
    )

    async def _is_zonal_bucket(self, bucket):
        bucket_type = await self._lookup_bucket_type(bucket)
        return bucket_type == BucketType.ZONAL_HIERARCHICAL

    async def _fetch_range_split(
        self, path, start=None, chunk_lengths=None, mrd=None, size=None, **kwargs
    ):
        """
        Reading multiple reads in one large stream.

        Optimized for Zonal Buckets:
        Leverages AsyncMultiRangeDownloader.download_ranges() to fetch all requested
        'chunk_lengths' (chunks) concurrently in a single batch request, significantly
        improving performance for ReadAheadV2.
        """

        bucket, object_name, generation = self.split_path(path)
        mrd_created = False
        try:
            if mrd is None:
                # Check before creating MRD
                if not await self._is_zonal_bucket(bucket):
                    raise RuntimeError(
                        "Internal error, this method is only supported for zonal buckets!"
                    )

                await self._get_grpc_client()
                mrd = await AsyncMultiRangeDownloader.create_mrd(
                    self.grpc_client, bucket, object_name, generation
                )
                mrd_created = True

            file_size = size or mrd.persisted_size
            if file_size is None:
                logger.warning(
                    "AsyncMultiRangeDownloader (MRD) has no 'persisted_size'. "
                    "Falling back to _info() to get the file size."
                )
                file_size = (await self._info(path))["size"]

            if chunk_lengths:
                start_offset = start if start is not None else 0
                current_offset = start_offset

                if start_offset >= file_size or (
                    chunk_lengths is not None
                    and start_offset + sum(chunk_lengths) > file_size
                ):
                    raise RuntimeError("Request not satisfiable.")

                buffers = []  # To hold the results in order
                read_ranges = []  # To pass to MRD

                for length in chunk_lengths:
                    buf = BytesIO()
                    buffers.append(buf)

                    if length > 0:
                        read_ranges.append((current_offset, length, buf))

                    current_offset += length

                if read_ranges:
                    await mrd.download_ranges(read_ranges)

                return [b.getvalue() for b in buffers]
            else:
                end = kwargs.get("end")
                offset, length = await self._process_limits_to_offset_and_length(
                    path, start, end, file_size
                )

                data = await zb_hns_utils.download_range(
                    offset=offset, length=length, mrd=mrd
                )
                return [data]
        finally:
            if mrd_created:
                await mrd.close()

    async def _cat_file(self, path, start=None, end=None, mrd=None, **kwargs):
        """Fetch a file's contents as bytes, with an optimized path for Zonal buckets.

        This method overrides the parent `_cat_file` to read objects in Zonal buckets using gRPC.

        Args:
            path (str): The full GCS path to the file (e.g., "bucket/object").
            start (int, optional): The starting byte position to read from.
            end (int, optional): The ending byte position to read to.
            mrd (AsyncMultiRangeDownloader, optional): An existing multi-range
                downloader instance. If not provided, a new one will be created for Zonal buckets.

        Returns:
            bytes: The content of the file or file range.
        """
        try:
            mrd_created = False

            # A new MRD is required when read is done directly by the
            # GCSFilesystem class without creating a GCSFile object first.
            if mrd is None:
                bucket, object_name, generation = self.split_path(path)
                # Fall back to default implementation if not a zonal bucket
                if not await self._is_zonal_bucket(bucket):
                    return await super()._cat_file(path, start=start, end=end, **kwargs)

                await self._get_grpc_client()
                mrd = await AsyncMultiRangeDownloader.create_mrd(
                    self.grpc_client, bucket, object_name, generation
                )
                mrd_created = True

            file_size = mrd.persisted_size
            if file_size is None:
                logger.warning(
                    "AsyncMultiRangeDownloader (MRD) exists but has no 'persisted_size'. "
                    "Falling back to _info() to get the file size. "
                    "This may result in incorrect behavior for unfinalized objects."
                )
            offset, length = await self._process_limits_to_offset_and_length(
                path, start, end, file_size
            )

            return await zb_hns_utils.download_range(
                offset=offset, length=length, mrd=mrd
            )
        finally:
            # Explicit cleanup if we created the MRD
            if mrd_created:
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

    def _update_dircache_after_rename(self, path1, path2):
        """
        Performs a targeted update of the directory cache after a successful
        folder rename operation.

        This involves three main steps:
        1. Removing the source folder and all its descendants from the cache.
        2. Removing the source folder's entry from its parent's listing.
        3. Adding the new destination folder's entry to its parent's listing.

        Args:
            path1 (str): The source path that was renamed.
            path2 (str): The destination path.
        """
        # 1. Find and remove all descendant paths of the source from the cache.
        source_prefix = f"{path1.rstrip('/')}/"
        for key in list(self.dircache):
            if key.startswith(source_prefix):
                self.dircache.pop(key, None)

        # 2. Remove the old source entry from its parent's listing.
        self.dircache.pop(path1, None)
        parent1 = self._parent(path1)
        if parent1 in self.dircache:
            self.dircache[parent1] = [
                e for e in self.dircache[parent1] if e.get("name") != path1
            ]

        # 3. Invalidate the destination path and update its parent's cache.
        self.dircache.pop(path2, None)
        parent2 = self._parent(path2)
        if parent2 in self.dircache:
            _, key2, _ = self.split_path(path2)
            new_entry = {
                "Key": key2,
                "Size": 0,
                "name": path2,
                "size": 0,
                "type": "directory",
                "storageClass": "DIRECTORY",
            }
            self.dircache[parent2].append(new_entry)

    async def _mv(self, path1, path2, **kwargs):
        """
        Move a file or directory. Overrides the parent `_mv` to provide an
        optimized, atomic implementation for renaming folders in HNS-enabled
        buckets. Falls back to the parent's object-level copy-and-delete
        implementation for files or for non-HNS buckets.
        """
        if path1 == path2:
            logger.debug(
                "%s mv: The paths are the same, so no files/directories were moved.",
                self,
            )
            return

        bucket1, key1, _ = self.split_path(path1)
        bucket2, key2, _ = self.split_path(path2)

        is_hns = await self._is_bucket_hns_enabled(bucket1)

        if not is_hns:
            logger.debug(
                f"Not an HNS bucket. Falling back to object-level mv for '{path1}' to '{path2}'."
            )
            return await self.loop.run_in_executor(
                None, partial(super().mv, path1, path2, **kwargs)
            )

        try:
            info1 = await self._info(path1)
            is_folder = info1.get("type") == "directory"

            # We only use HNS rename if the source is a folder and the move is
            # within the same bucket.
            if is_folder and bucket1 == bucket2 and key1 and key2:
                logger.info(
                    f"Using HNS-aware folder rename for '{path1}' to '{path2}'."
                )
                source_folder_name = f"projects/_/buckets/{bucket1}/folders/{key1}"
                destination_folder_id = key2

                request = storage_control_v2.RenameFolderRequest(
                    name=source_folder_name,
                    destination_folder_id=destination_folder_id,
                )

                logger.debug(f"rename_folder request: {request}")
                await self._storage_control_client.rename_folder(request=request)
                self._update_dircache_after_rename(path1, path2)

                logger.info(
                    "Successfully renamed folder from '%s' to '%s'", path1, path2
                )
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
        # TODO: Check feasibility to call async copy and rm methods instead of sync mv method
        return await self.loop.run_in_executor(
            None, partial(super().mv, path1, path2, **kwargs)
        )

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
        self, path, create_parents=False, enable_hierarchical_namespace=False, **kwargs
    ):
        """
        If the path does not contain an object key, a new bucket is created.
        If `enable_hierarchical_namespace` is True, the bucket will have Hierarchical Namespace enabled.

        For HNS-enabled buckets, this method creates a folder object. If
        `create_parents` is True, any missing parent folders are also created.

        If bucket doesn't exist, enable_hierarchical_namespace and create_parents are set to True
        and the path includes a key then HNS-enabled bucket will be created
        and also the folders within that bucket.

        If `create_parents` is False and a parent does not exist, a
        FileNotFoundError is raised.

        For non-HNS buckets, it falls back to the parent implementation which
        may involve creating a bucket or doing nothing (as GCS has no true empty directories).
        """
        path = self._strip_protocol(path)
        if enable_hierarchical_namespace:
            kwargs["hierarchicalNamespace"] = {"enabled": True}
            # HNS buckets require uniform bucket-level access.
            kwargs["iamConfiguration"] = {"uniformBucketLevelAccess": {"enabled": True}}
            # When uniformBucketLevelAccess is enabled, ACLs cannot be used.
            # We must explicitly set them to None to prevent the parent
            # method from using default ACLs.
            kwargs["acl"] = None
            kwargs["default_acl"] = None

        bucket, key, _ = self.split_path(path)
        # If the key is empty, the path refers to a bucket, not an object.
        # Defer to the parent method to handle bucket creation.
        if not key:
            return await super()._mkdir(path, create_parents=create_parents, **kwargs)

        is_hns = False
        # If creating an HNS bucket, check for its existence first.
        if create_parents and enable_hierarchical_namespace:
            if not await self._exists(bucket):
                await super()._mkdir(bucket, create_parents=True, **kwargs)
                is_hns = True  # Skip HNS check since we just created it.

        if not is_hns:
            # If the bucket was not created above, we need to check its type.
            is_hns = await self._is_bucket_hns_enabled(bucket)

        if not is_hns:
            return await super()._mkdir(path, create_parents=create_parents, **kwargs)

        logger.info(f"Using HNS-aware mkdir for '{path}'.")
        parent = f"projects/_/buckets/{bucket}"
        folder_id = key.rstrip("/")
        request = storage_control_v2.CreateFolderRequest(
            parent=parent,
            folder_id=folder_id,
            recursive=create_parents,
        )
        try:
            logger.debug(f"create_folder request: {request}")
            await self._get_control_plane_client()
            await self._storage_control_client.create_folder(request=request)
            # Instead of invalidating the parent cache, update it to add the new entry.
            parent_path = self._parent(path)
            if parent_path in self.dircache:
                new_entry = {
                    "Key": key.rstrip("/"),
                    "Size": 0,
                    "name": path,
                    "size": 0,
                    "type": "directory",
                    "storageClass": "DIRECTORY",
                }
                self.dircache[parent_path].append(new_entry)
        except api_exceptions.Conflict as e:
            logger.debug(f"Directory already exists: {path}: {e}")
        except api_exceptions.FailedPrecondition as e:
            # This error can occur if create_parents=False and the parent dir doesn't exist.
            # Translate it to FileNotFoundError for fsspec compatibility.
            raise FileNotFoundError(
                f"mkdir for '{path}' failed due to a precondition error: {e}"
            ) from e

    mkdir = asyn.sync_wrapper(_mkdir)

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

                request = storage_control_v2.GetFolderRequest(name=folder_resource_name)

                # Verify existence using get_folder API
                response = await self._storage_control_client.get_folder(
                    request=request
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
            logger.info(f"Using HNS-aware rmdir for '{path}'.")
            folder_name = f"projects/_/buckets/{bucket}/folders/{key.rstrip('/')}"
            request = storage_control_v2.DeleteFolderRequest(name=folder_name)

            logger.debug(f"delete_folder request: {request}")
            await self._storage_control_client.delete_folder(request=request)

            # Remove the directory from the cache and from its parent's listing.
            self.dircache.pop(path, None)
            parent = self._parent(path)
            if parent in self.dircache:
                # Remove the deleted directory entry from the parent's listing.
                self.dircache[parent] = [
                    e for e in self.dircache[parent] if e.get("name") != path
                ]
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
        self, path, recursive=False, maxdepth=None, detail=False
    ):
        if maxdepth is not None and maxdepth < 1:
            raise ValueError("maxdepth must be at least 1")

        if isinstance(path, str):
            out = await self._expand_path_with_details(
                [path], recursive, maxdepth, detail=detail
            )
        else:
            out = {} if detail else set()
            path = [self._strip_protocol(p) for p in path]

            for p in path:
                if has_magic(p):
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
        deletion in HNS-enabled buckets. For non-HNS buckets, it falls back
        to the parent implementation.

        For HNS buckets, it first expands the path to get a list of all files
        and directories, then categorizes them. It deletes files in batches
        and then deletes directories individually.

        Args:
            path (str or list): The path(s) to delete.
            recursive (bool): If True, deletes directories and their contents.
            maxdepth (int, optional): The maximum depth to traverse for deletion.
            batchsize (int): The number of files to delete in a single batch request.
        """
        if isinstance(path, list):
            # For HNS check, we can check for bucket type from the first path.
            bucket, _, _ = self.split_path(path[0]) if path else (None, None, None)
        else:
            bucket, _, _ = self.split_path(path)

        is_hns = await self._is_bucket_hns_enabled(bucket)

        if not is_hns:
            # Fall back to the parent's async rm implementation for non-HNS buckets.
            return await super()._rm(
                path, recursive=recursive, maxdepth=maxdepth, batchsize=batchsize
            )

        paths = await self._expand_path_with_details(
            path, recursive=recursive, maxdepth=maxdepth, detail=True
        )

        # Separate files and directories based on their type.
        # Directories must be deleted from the deepest first.
        files = list({p["name"] for p in paths if p["type"] == "file"})
        dirs = sorted(
            list({p["name"] for p in paths if p["type"] == "directory"}),
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
        files_task = self.loop.create_task(
            super()._find(
                path,
                withdirs=False,  # Fetch files only
                detail=True,  # Get full details for merging and populating cache
                prefix=prefix,
                versions=versions,
                maxdepth=maxdepth,
                update_cache=False,  # Defer caching until merging files and folders
                **kwargs,
            )
        )

        # 2. Fetch all folders recursively. This is necessary to find all folders,
        # especially empty ones.
        folders_task = self.loop.create_task(
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
        folders = []
        base_path = self.split_path(path)[1].rstrip("/")
        full_prefix = f"{base_path}/{prefix}".strip("/") if base_path else prefix

        folder_id = full_prefix
        if folder_id and not folder_id.endswith("/"):
            folder_id += "/"
        parent = f"projects/_/buckets/{bucket}"
        request = storage_control_v2.ListFoldersRequest(parent=parent, prefix=folder_id)
        logger.debug(f"list_folders request: {request}")

        async for folder in await self._storage_control_client.list_folders(
            request=request
        ):
            folders.append(self._create_folder_entry(bucket, folder))

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
            finalize_on_close = kwargs.get("finalize_on_close", False)
            await writer.close(finalize_on_close=finalize_on_close)

        self.invalidate_cache(self._parent(rpath))

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
            for i in range(0, len(data), chunksize):
                await writer.append(data[i : i + chunksize])
        finally:
            finalize_on_close = kwargs.get("finalize_on_close", False)
            await writer.close(finalize_on_close=finalize_on_close)

        self.invalidate_cache(self._parent(path))

    async def _get_file(self, rpath, lpath, callback=None, **kwargs):
        bucket, key, generation = self.split_path(rpath)
        if not await self._is_zonal_bucket(bucket):
            return await super()._get_file(
                rpath,
                lpath,
                callback=callback,
                **kwargs,
            )

        if os.path.isdir(lpath):
            return
        callback = callback or NoOpCallback()

        mrd = None
        try:
            await self._get_grpc_client()
            mrd = await AsyncMultiRangeDownloader.create_mrd(
                self.grpc_client, bucket, key, generation
            )

            size = mrd.persisted_size
            if size is None:
                logger.warning(
                    "AsyncMultiRangeDownloader (MRD) has no 'persisted_size'. "
                    "Falling back to _info() to get the file size. "
                    "This may result in incorrect behavior for unfinalized objects."
                )
                size = (await self._info(rpath))["size"]
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
        except Exception as e:
            # Clean up the corrupted file before raising error
            if os.path.exists(lpath):
                os.remove(lpath)
            raise e
        finally:
            if mrd:
                await mrd.close()

    async def _do_list_objects(
        self,
        path,
        max_results=None,
        delimiter="/",
        prefix="",
        versions=False,
        **kwargs,
    ):
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


async def upload_chunk(fs, location, data, offset, size, content_type):
    """
    Uploads a chunk of data using AsyncAppendableObjectWriter for zonal buckets.
    Finalizes the upload when the total uploaded data size reaches the specified size.
    Delegates to core upload_chunk implementaion for Non-Zonal buckets.
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
        await location.close(finalize_on_close=False)
        raise

    if (location.offset or 0) >= size:
        logger.debug("Uploaded data is equal or greater than size. Finalizing upload.")
        await location.close(finalize_on_close=True)


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
    Delegates to core initiate_upload implementaion for Non-Zonal buckets.

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
    Delegates to core simple_upload implementaion for Non-Zonal buckets.

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
        finalize_on_close = kwargs.get("finalize_on_close", False)
        await writer.close(finalize_on_close=finalize_on_close)
