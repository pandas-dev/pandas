import errno
import logging
from contextlib import contextmanager, AsyncExitStack
from botocore.exceptions import ClientError


logger = logging.getLogger("s3fs")


@contextmanager
def ignoring(*exceptions):
    try:
        yield
    except exceptions:
        pass


class S3BucketRegionCache:
    # See https://github.com/aio-libs/aiobotocore/issues/866
    # for details.

    def __init__(self, session, **client_kwargs):
        self._session = session
        self._stack = AsyncExitStack()
        self._client = None
        self._client_kwargs = client_kwargs
        self._buckets = {}
        self._regions = {}

    async def get_bucket_client(self, bucket_name=None):
        if bucket_name in self._buckets:
            return self._buckets[bucket_name]

        general_client = await self.get_client()
        if bucket_name is None:
            return general_client

        try:
            response = await general_client.head_bucket(Bucket=bucket_name)
        except ClientError as e:
            logger.debug("RC: HEAD_BUCKET call for %r has failed", bucket_name)
            response = e.response

        region = (
            response["ResponseMetadata"]
            .get("HTTPHeaders", {})
            .get("x-amz-bucket-region")
        )

        if not region:
            logger.debug(
                "RC: No region in HEAD_BUCKET call response for %r, returning the general client",
                bucket_name,
            )
            return general_client

        if region not in self._regions:
            logger.debug(
                "RC: Creating a new regional client for %r on the region %r",
                bucket_name,
                region,
            )
            self._regions[region] = await self._stack.enter_async_context(
                self._session.create_client(
                    "s3", region_name=region, **self._client_kwargs
                )
            )

        client = self._buckets[bucket_name] = self._regions[region]
        return client

    async def get_client(self):
        if not self._client:
            self._client = await self._stack.enter_async_context(
                self._session.create_client("s3", **self._client_kwargs)
            )
        return self._client

    async def clear(self):
        logger.debug("RC: discarding all clients")
        self._buckets.clear()
        self._regions.clear()
        self._client = None
        await self._stack.aclose()

    async def __aenter__(self):
        return self

    async def __aexit__(self, *exc_args):
        await self.clear()


class FileExpired(IOError):
    """
    Is raised, when the file content has been changed from a different process after
    opening the file. Reading the file would lead to invalid or inconsistent output.
    This can also be triggered by outdated file-information inside the directory cache.
    In this case ``S3FileSystem.invalidate_cache`` can be used to force an update of
    the file-information when opening the file.
    """

    def __init__(self, filename: str, e_tag: str):
        super().__init__(
            errno.EBUSY,
            "The remote file corresponding to filename %s and Etag %s no longer exists."
            % (filename, e_tag),
        )


def title_case(string):
    """
    TitleCases a given string.

    Parameters
    ----------
    string : underscore separated string
    """
    return "".join(x.capitalize() for x in string.split("_"))


class ParamKwargsHelper(object):
    """
    Utility class to help extract the subset of keys that an s3 method is
    actually using

    Parameters
    ----------
    s3 : boto S3FileSystem
    """

    _kwarg_cache = {}

    def __init__(self, s3):
        self.s3 = s3

    def _get_valid_keys(self, model_name):
        if model_name not in self._kwarg_cache:
            model = self.s3.meta.service_model.operation_model(model_name)
            valid_keys = (
                set(model.input_shape.members.keys())
                if model.input_shape is not None
                else set()
            )
            self._kwarg_cache[model_name] = valid_keys
        return self._kwarg_cache[model_name]

    def filter_dict(self, method_name, d):
        model_name = title_case(method_name)
        valid_keys = self._get_valid_keys(model_name)
        if isinstance(d, SSEParams):
            d = d.to_kwargs()
        return {k: v for k, v in d.items() if k in valid_keys}


class SSEParams(object):
    def __init__(
        self,
        server_side_encryption=None,
        sse_customer_algorithm=None,
        sse_customer_key=None,
        sse_kms_key_id=None,
    ):
        self.ServerSideEncryption = server_side_encryption
        self.SSECustomerAlgorithm = sse_customer_algorithm
        self.SSECustomerKey = sse_customer_key
        self.SSEKMSKeyId = sse_kms_key_id

    def to_kwargs(self):
        return {k: v for k, v in self.__dict__.items() if v is not None}


def _get_brange(size, block):
    """
    Chunk up a file into zero-based byte ranges

    Parameters
    ----------
    size : file size
    block : block size
    """
    for offset in range(0, size, block):
        yield offset, min(offset + block - 1, size - 1)
