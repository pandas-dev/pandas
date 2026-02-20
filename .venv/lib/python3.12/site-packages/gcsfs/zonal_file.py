import logging

from fsspec import asyn
from google.cloud.storage.asyncio.async_multi_range_downloader import (
    AsyncMultiRangeDownloader,
)

from gcsfs import zb_hns_utils
from gcsfs.core import DEFAULT_BLOCK_SIZE, GCSFile

from .caching import (  # noqa: F401 Unused import to register GCS-Specific caches, Please do not remove it.
    ReadAheadChunked,
)

logger = logging.getLogger("gcsfs.zonal_file")


class ZonalFile(GCSFile):
    """
    ZonalFile is subclass of GCSFile and handles data operations from
    Zonal buckets only using a high-performance gRPC path.
    """

    def __init__(
        self,
        gcsfs,
        path,
        mode="rb",
        block_size=DEFAULT_BLOCK_SIZE,
        autocommit=True,
        cache_type="readahead_chunked",
        cache_options=None,
        acl=None,
        consistency="md5",
        metadata=None,
        content_type=None,
        timeout=None,
        fixed_key_metadata=None,
        generation=None,
        kms_key_name=None,
        finalize_on_close=False,
        **kwargs,
    ):
        """
        Initializes the ZonalFile object.

        For Zonal buckets, `finalize_on_close` is set to `False` by default to optimize
        for write throughput and keep the file appendable. This means that when exiting
        a `with` block or closing, the file will not be automatically finalized. To
        ensure the write is finalized, `.commit()` must be called explicitly or
        `finalize_on_close` must be set to `True` when opening the file.
        """
        bucket, key, generation = gcsfs._split_path(path)
        if not key:
            raise OSError("Attempt to open a bucket")
        self.mrd = None
        self.finalize_on_close = finalize_on_close
        self.finalized = False
        self.mode = mode
        self.gcsfs = gcsfs
        object_size = None
        if "r" in self.mode:
            self.mrd = asyn.sync(
                self.gcsfs.loop, self._init_mrd, bucket, key, generation
            )
            object_size = self.mrd.persisted_size
            if object_size is None:
                logger.warning(
                    "AsyncMultiRangeDownloader (MRD) exists but has no 'persisted_size'. "
                    "This may result in incorrect behavior for unfinalized objects."
                )
        elif "w" or "a" in self.mode:
            self.aaow = asyn.sync(
                self.gcsfs.loop,
                self._init_aaow,
                bucket,
                key,
                generation,
            )
        else:
            raise NotImplementedError(
                "Only read, write and append operations are currently supported for Zonal buckets."
            )

        super().__init__(
            gcsfs,
            path,
            mode,
            block_size,
            autocommit,
            cache_type,
            cache_options,
            acl,
            consistency,
            metadata,
            content_type,
            timeout,
            fixed_key_metadata,
            generation,
            kms_key_name,
            # Zonal buckets support append; this prevents GCSFile from forcing 'w' mode
            supports_append="a" in mode,
            # pass persisted_size here so that Cache is initialized with correct object size
            size=object_size,
            **kwargs,
        )

    async def _init_mrd(self, bucket_name, object_name, generation=None):
        """
        Initializes the AsyncMultiRangeDownloader.
        """
        await self.gcsfs._get_grpc_client()
        return await AsyncMultiRangeDownloader.create_mrd(
            self.gcsfs.grpc_client, bucket_name, object_name, generation
        )

    async def _init_aaow(self, bucket_name, object_name, generation=None):
        """
        Initializes the AsyncAppendableObjectWriter.
        """
        # generation is needed while creating aaow to append to existing objects
        if "a" in self.mode and generation is None:
            try:
                # self.path might not be set yet, so reconstruct full path
                info = await self.gcsfs._info(f"{bucket_name}/{object_name}")
                generation = info.get("generation")
            except FileNotFoundError:
                # if file doesn't exist, we don't need generation
                pass
        await self.gcsfs._get_grpc_client()
        return await zb_hns_utils.init_aaow(
            self.gcsfs.grpc_client, bucket_name, object_name, generation
        )

    def _fetch_range(self, start=None, end=None, chunk_lengths=None):
        """
        Overrides the default _fetch_range to implement the gRPC read path.

        See super() class for documentation.
        """
        if end is not None and chunk_lengths is not None:
            raise ValueError(
                "The end and chunk_lengths arguments are mutually exclusive and cannot be used together."
            )

        try:
            if chunk_lengths is not None:
                return asyn.sync(
                    self.fs.loop,
                    self.gcsfs._fetch_range_split,
                    self.path,
                    start=start,
                    chunk_lengths=chunk_lengths,
                    size=self.size,
                    mrd=self.mrd,
                )
            return self.gcsfs.cat_file(self.path, start=start, end=end, mrd=self.mrd)
        except RuntimeError as e:
            if "not satisfiable" in str(e):
                return b"" if chunk_lengths is None else [b""]
            raise

    def write(self, data):
        """
        Writes data using AsyncAppendableObjectWriter.

        For more details, see the documentation for AsyncAppendableObjectWriter:
        https://github.com/googleapis/python-storage/blob/9e6fefdc24a12a9189f7119bc9119e84a061842f/google/cloud/storage/_experimental/asyncio/async_appendable_object_writer.py#L38
        """
        if self.closed:
            raise ValueError("I/O operation on closed file.")
        if not self.writable():
            raise ValueError("File not in write mode")
        if self.forced:
            raise ValueError("This file has been force-flushed, can only close")

        asyn.sync(self.gcsfs.loop, self.aaow.append, data)

    def flush(self, force=False):
        """
        Flushes the AsyncAppendableObjectWriter, sending all buffered data
        to the server.
        """
        if self.closed:
            raise ValueError("Flush on closed file.")
        if force and self.forced:
            raise ValueError("Force flush cannot be called more than once.")
        if self.finalized:
            logger.warning("File is already finalized. Ignoring flush call.")
            return
        if force:
            self.forced = True

        if self.readable():
            # no-op to flush on read-mode
            return

        asyn.sync(self.gcsfs.loop, self.aaow.flush)

    def commit(self):
        """
        Commits the write by finalizing the AsyncAppendableObjectWriter.
        """
        if not self.writable():  # No-op
            logger.warning("File not in write mode.")
            return
        if self.finalized:  # No-op
            logger.warning("This file has already been finalized.")
            return
        asyn.sync(self.gcsfs.loop, self.aaow.finalize)
        self.finalized = True
        # File is already finalized, avoid finalizing again on close
        self.finalize_on_close = False

    def discard(self):
        """Discard is not applicable for Zonal Buckets. Log a warning instead."""
        logger.warning(
            "Discard is not applicable for Zonal Buckets. \
            Data is uploaded via streaming and cannot be cancelled."
        )

    def _initiate_upload(self):
        """Initiates the upload for Zonal buckets using gRPC."""
        from gcsfs.extended_gcsfs import initiate_upload

        self.location = asyn.sync(
            self.gcsfs.loop,
            initiate_upload,
            self.gcsfs,
            self.bucket,
            self.key,
            self.content_type,
            self.metadata,
            self.fixed_key_metadata,
            mode="create" if "x" in self.mode else "overwrite",
            kms_key_name=self.kms_key_name,
            timeout=self.timeout,
        )

    def _simple_upload(self):
        """Performs a simple upload for Zonal buckets using gRPC."""
        from gcsfs.extended_gcsfs import simple_upload

        self.buffer.seek(0)
        data = self.buffer.read()
        asyn.sync(
            self.gcsfs.loop,
            simple_upload,
            self.gcsfs,
            self.bucket,
            self.key,
            data,
            self.metadata,
            self.consistency,
            self.content_type,
            self.fixed_key_metadata,
            mode="create" if "x" in self.mode else "overwrite",
            kms_key_name=self.kms_key_name,
            timeout=self.timeout,
            finalize_on_close=self.finalize_on_close,
        )

    def _upload_chunk(self, final=False):
        raise NotImplementedError(
            "_upload_chunk is not implemented yet for ZonalFile. Please use write() instead."
        )

    def close(self):
        """
        Closes the ZonalFile and the underlying AsyncMultiRangeDownloader and AsyncAppendableObjectWriter.
        If in write mode, finalizes the write if finalize_on_close is True.
        """
        if self.closed:
            return
        # super is closed before aaow since flush may need aaow
        super().close()
        if hasattr(self, "mrd") and self.mrd:
            asyn.sync(self.gcsfs.loop, self.mrd.close)
        # Only close aaow if the stream is open
        if hasattr(self, "aaow") and self.aaow and self.aaow._is_stream_open:
            asyn.sync(
                self.gcsfs.loop,
                self.aaow.close,
                finalize_on_close=self.finalize_on_close,
            )
