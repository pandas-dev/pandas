from io import BytesIO

from google.cloud.storage.asyncio.async_appendable_object_writer import (
    AsyncAppendableObjectWriter,
)


async def download_range(offset, length, mrd):
    """
    Downloads a byte range from the file asynchronously.
    """
    # If length = 0, mrd returns till end of file, so handle that case here
    if length == 0:
        return b""
    buffer = BytesIO()
    await mrd.download_ranges([(offset, length, buffer)])
    return buffer.getvalue()


async def init_aaow(grpc_client, bucket_name, object_name, generation=None):
    """
    Creates and opens the AsyncAppendableObjectWriter.
    """

    writer = AsyncAppendableObjectWriter(
        client=grpc_client,
        bucket_name=bucket_name,
        object_name=object_name,
        generation=generation,
    )
    await writer.open()
    return writer
