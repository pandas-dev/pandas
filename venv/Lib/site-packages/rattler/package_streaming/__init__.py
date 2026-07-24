from os import PathLike
from typing import Optional, Tuple

from rattler.networking.client import Client
from rattler.rattler import download_bytes as py_download_bytes
from rattler.rattler import download_to_path as py_download_to_path
from rattler.rattler import download_to_writer as py_download_to_writer
from rattler.rattler import download_and_extract as py_download_and_extract
from rattler.rattler import extract as py_extract
from rattler.rattler import extract_tar_bz2 as py_extract_tar_bz2
from rattler.rattler import fetch_raw_package_file_from_url as py_fetch_raw_package_file_from_url


def extract(path: PathLike[str], dest: PathLike[str]) -> Tuple[bytes, bytes]:
    """Extract a file to a destination."""
    return py_extract(path, dest)


def extract_tar_bz2(path: PathLike[str], dest: PathLike[str]) -> Tuple[bytes, bytes]:
    """Extract a tar.bz2 file to a destination."""
    return py_extract_tar_bz2(path, dest)


async def download_to_path(client: Client, url: str, dest: PathLike[str]) -> None:
    """
    Stream a package archive from a URL to a destination path.

    This method does not buffer the whole response in Python memory. Response
    bytes are fetched incrementally and written directly to `dest`.
    """
    await py_download_to_path(client._client, url, dest)


async def download_bytes(client: Client, url: str) -> bytes:
    """
    Download a package archive from a URL into memory.

    This is a convenience API. The full response body is buffered before the
    `bytes` object is returned, so peak memory use scales with the artifact
    size.
    """
    return await py_download_bytes(client._client, url)


async def download_to_writer(client: Client, url: str, writer: object) -> None:
    """
    Stream a package archive from a URL into a Python writer.

    The response body is fetched incrementally. For each chunk, `writer.write`
    is called with a `bytes` object. The writer must provide a synchronous
    `write(bytes)` method, for example `io.BytesIO()` or an open binary file.
    """
    await py_download_to_writer(client._client, url, writer)


async def download_and_extract(
    client: Client, url: str, dest: PathLike[str], expected_sha: Optional[bytes] = None
) -> Tuple[bytes, bytes]:
    """Download a file from a URL and extract it to a destination."""
    return await py_download_and_extract(client._client, url, dest, expected_sha)


async def fetch_raw_package_file_from_url(client: Client, url: str, path: str) -> bytes:
    """
    Fetch raw bytes for a file inside a remote `.conda` package using sparse
    range requests.
    """
    return await py_fetch_raw_package_file_from_url(client._client, url, path)
