import io

from botocore.httpchecksum import (
    _CHECKSUM_CLS,
    AwsChunkedWrapper,
    FlexibleChecksumError,
    _apply_request_header_checksum,
    _register_checksum_algorithm_feature_id,
    base64,
    conditionally_calculate_md5,
    determine_content_length,
    logger,
)

from aiobotocore._helpers import resolve_awaitable
from aiobotocore.response import HttpxStreamingBody, StreamingBody

try:
    import httpx
except ImportError:
    httpx = None


class AioAwsChunkedWrapper(AwsChunkedWrapper):
    async def read(self, size=None):
        # Normalize "read all" size values to None
        if size is not None and size <= 0:
            size = None

        # If the underlying body is done and we have nothing left then
        # end the stream
        if self._complete and not self._remaining:
            return b""

        # While we're not done and want more bytes
        want_more_bytes = size is None or size > len(self._remaining)
        while not self._complete and want_more_bytes:
            self._remaining += await self._make_chunk()
            want_more_bytes = size is None or size > len(self._remaining)

        # If size was None, we want to return everything
        if size is None:
            size = len(self._remaining)

        # Return a chunk up to the size asked for
        to_return = self._remaining[:size]
        self._remaining = self._remaining[size:]
        return to_return

    async def _make_chunk(self):
        # NOTE: Chunk size is not deterministic as read could return less. This
        # means we cannot know the content length of the encoded aws-chunked
        # stream ahead of time without ensuring a consistent chunk size

        raw_chunk = await resolve_awaitable(self._raw.read(self._chunk_size))
        hex_len = hex(len(raw_chunk))[2:].encode("ascii")
        self._complete = not raw_chunk

        if self._checksum:
            self._checksum.update(raw_chunk)

        if self._checksum and self._complete:
            name = self._checksum_name.encode("ascii")
            checksum = self._checksum.b64digest().encode("ascii")
            return b"0\r\n%s:%s\r\n\r\n" % (name, checksum)

        return b"%s\r\n%s\r\n" % (hex_len, raw_chunk)

    def __aiter__(self):
        return self

    async def __anext__(self):
        while not self._complete:
            return await self._make_chunk()
        raise StopAsyncIteration()


# unfortunately we can't inherit from botocore's StreamingChecksumBody due to
# subclassing
class StreamingChecksumBody(StreamingBody):
    def __init__(self, raw_stream, content_length, checksum, expected):
        super().__init__(raw_stream, content_length)
        self._checksum = checksum
        self._expected = expected

    async def read(self, amt=None):
        chunk = await super().read(amt=amt)
        self._checksum.update(chunk)
        if amt is None or (not chunk and amt > 0):
            self._validate_checksum()
        return chunk

    async def readinto(self, b: bytearray):
        chunk = await self.__wrapped__.content.read(len(b))
        amount_read = len(chunk)
        b[:amount_read] = chunk

        if amount_read == len(b):
            view = b
        else:
            view = memoryview(b)[:amount_read]

        self._checksum.update(view)
        if amount_read == 0 and len(b) > 0:
            self._validate_checksum()
        return amount_read

    def _validate_checksum(self):
        if self._checksum.digest() != base64.b64decode(self._expected):
            error_msg = (
                f"Expected checksum {self._expected} did not match calculated "
                f"checksum: {self._checksum.b64digest()}"
            )
            raise FlexibleChecksumError(error_msg=error_msg)


# TODO: fix inheritance? read & _validate_checksum are the exact same as above
# only diff is super class and how to call __init__
class HttpxStreamingChecksumBody(HttpxStreamingBody):
    def __init__(self, raw_stream, content_length, checksum, expected):
        # HttpxStreamingbody doesn't use content_length
        super().__init__(raw_stream)
        self._checksum = checksum
        self._expected = expected

    # TODO: this class is largely (or possibly entirely) untested. The tests need to be
    # more thoroughly rewritten wherever they directly create Streamingbody,
    # StreamingChecksumBody, etc.

    async def read(self, amt=None):
        chunk = await super().read(amt=amt)
        self._checksum.update(chunk)
        if amt is None or (not chunk and amt > 0):
            self._validate_checksum()
        return chunk

    async def readinto(self, b: bytearray):
        chunk = await self.__wrapped__.content.read(len(b))
        amount_read = len(chunk)
        b[:amount_read] = chunk

        if amount_read == len(b):
            view = b
        else:
            view = memoryview(b)[:amount_read]

        self._checksum.update(view)
        if amount_read == 0 and len(b) > 0:
            self._validate_checksum()
        return amount_read

    def _validate_checksum(self):
        if self._checksum.digest() != base64.b64decode(self._expected):
            error_msg = (
                f"Expected checksum {self._expected} did not match calculated "
                f"checksum: {self._checksum.b64digest()}"
            )
            raise FlexibleChecksumError(error_msg=error_msg)


async def handle_checksum_body(
    http_response, response, context, operation_model
):
    headers = response["headers"]
    checksum_context = context.get("checksum", {})
    algorithms = checksum_context.get("response_algorithms")

    if not algorithms:
        return

    for algorithm in algorithms:
        header_name = f"x-amz-checksum-{algorithm}"
        # If the header is not found, check the next algorithm
        if header_name not in headers:
            continue

        # If a - is in the checksum this is not valid Base64. S3 returns
        # checksums that include a -# suffix to indicate a checksum derived
        # from the hash of all part checksums. We cannot wrap this response
        if "-" in headers[header_name]:
            continue

        if operation_model.has_streaming_output:
            response["body"] = _handle_streaming_response(
                http_response, response, algorithm
            )
        else:
            response["body"] = await _handle_bytes_response(
                http_response, response, algorithm
            )

        # Expose metadata that the checksum check actually occurred
        checksum_context = response["context"].get("checksum", {})
        checksum_context["response_algorithm"] = algorithm
        response["context"]["checksum"] = checksum_context
        return

    logger.debug(
        'Skipping checksum validation. Response did not contain one of the following algorithms: %s.',
        algorithms,
    )


def _handle_streaming_response(http_response, response, algorithm):
    checksum_cls = _CHECKSUM_CLS.get(algorithm)
    header_name = f"x-amz-checksum-{algorithm}"
    if httpx is not None and isinstance(http_response.raw, httpx.Response):
        streaming_cls = HttpxStreamingChecksumBody
    else:
        streaming_cls = StreamingChecksumBody
    return streaming_cls(
        http_response.raw,
        response["headers"].get("content-length"),
        checksum_cls(),
        response["headers"][header_name],
    )


async def _handle_bytes_response(http_response, response, algorithm):
    body = await http_response.content
    header_name = f"x-amz-checksum-{algorithm}"
    checksum_cls = _CHECKSUM_CLS.get(algorithm)
    checksum = checksum_cls()
    checksum.update(body)
    expected = response["headers"][header_name]
    if checksum.digest() != base64.b64decode(expected):
        error_msg = (
            f"Expected checksum {expected} did not match calculated "
            f"checksum: {checksum.b64digest()}"
        )
        raise FlexibleChecksumError(error_msg=error_msg)
    return body


def apply_request_checksum(request):
    checksum_context = request.get("context", {}).get("checksum", {})
    algorithm = checksum_context.get("request_algorithm")

    if not algorithm:
        return

    if algorithm == "conditional-md5":
        # Special case to handle the http checksum required trait
        conditionally_calculate_md5(request)
    elif algorithm["in"] == "header":
        _apply_request_header_checksum(request)
    elif algorithm["in"] == "trailer":
        _apply_request_trailer_checksum(request)
    else:
        raise FlexibleChecksumError(
            error_msg="Unknown checksum variant: {}".format(algorithm["in"])
        )
    if "request_algorithm_header" in checksum_context:
        request_algorithm_header = checksum_context["request_algorithm_header"]
        request["headers"][request_algorithm_header["name"]] = (
            request_algorithm_header["value"]
        )


def _apply_request_trailer_checksum(request):
    checksum_context = request.get("context", {}).get("checksum", {})
    algorithm = checksum_context.get("request_algorithm")
    location_name = algorithm["name"]
    checksum_cls = _CHECKSUM_CLS.get(algorithm["algorithm"])

    headers = request["headers"]
    body = request["body"]

    if location_name in headers:
        # If the header is already set by the customer, skip calculation
        return

    # Cannot set this as aiohttp complains
    headers["Transfer-Encoding"] = "chunked"
    if "Content-Encoding" in headers:
        # We need to preserve the existing content encoding and add
        # aws-chunked as a new content encoding.
        headers["Content-Encoding"] += ",aws-chunked"
    else:
        headers["Content-Encoding"] = "aws-chunked"
    headers["X-Amz-Trailer"] = location_name
    _register_checksum_algorithm_feature_id(algorithm)

    content_length = determine_content_length(body)
    if content_length is not None:
        # Send the decoded content length if we can determine it. Some
        # services such as S3 may require the decoded content length
        headers["X-Amz-Decoded-Content-Length"] = str(content_length)

        if "Content-Length" in headers:
            del headers["Content-Length"]
            logger.debug(
                "Removing the Content-Length header since 'chunked' is specified for Transfer-Encoding."
            )

    if isinstance(body, (bytes, bytearray)):
        body = io.BytesIO(body)

    request["body"] = AioAwsChunkedWrapper(
        body,
        checksum_cls=checksum_cls,
        checksum_name=location_name,
    )
