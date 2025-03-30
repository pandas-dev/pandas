import base64
from base64 import b64encode
from hashlib import md5
from typing import Optional

from .retry import ChecksumError

try:
    import crcmod
except ImportError:
    crcmod = None


class ConsistencyChecker:
    def __init__(self):
        pass

    def update(self, data: bytes):
        pass

    def validate_json_response(self, gcs_object):
        pass

    def validate_headers(self, headers):
        pass

    def validate_http_response(self, r):
        pass


class MD5Checker(ConsistencyChecker):
    def __init__(self):
        self.md = md5()

    def update(self, data):
        self.md.update(data)

    def validate_json_response(self, gcs_object):
        mdback = gcs_object["md5Hash"]
        if b64encode(self.md.digest()) != mdback.encode():
            raise ChecksumError("MD5 checksum failed")

    def validate_headers(self, headers):
        if headers is not None and "X-Goog-Hash" in headers:

            dig = [
                bit.split("=")[1]
                for bit in headers["X-Goog-Hash"].split(",")
                if bit and bit.strip().startswith("md5=")
            ]
            if dig:
                if b64encode(self.md.digest()).decode().rstrip("=") != dig[0]:
                    raise ChecksumError("Checksum failure")
            else:
                raise NotImplementedError(
                    "No md5 checksum available to do consistency check. GCS does "
                    "not provide md5 sums for composite objects."
                )

    def validate_http_response(self, r):
        return self.validate_headers(r.headers)


class SizeChecker(ConsistencyChecker):
    def __init__(self):
        self.size = 0

    def update(self, data: bytes):
        self.size += len(data)

    def validate_json_response(self, gcs_object):
        assert int(gcs_object["size"]) == self.size, "Size mismatch"

    def validate_http_response(self, r):
        assert r.content_length == self.size


class Crc32cChecker(ConsistencyChecker):
    def __init__(self):
        self.crc32c = crcmod.Crc(0x11EDC6F41, initCrc=0, xorOut=0xFFFFFFFF)

    def update(self, data: bytes):
        self.crc32c.update(data)

    def validate_json_response(self, gcs_object):
        # docs for gcs_object: https://cloud.google.com/storage/docs/json_api/v1/objects
        digest = self.crc32c.digest()
        digest_b64 = base64.b64encode(digest).decode()
        expected = gcs_object["crc32c"]

        if digest_b64 != expected:
            raise ChecksumError(f'Expected "{expected}". Got "{digest_b64}"')

    def validate_headers(self, headers):
        if headers is not None:
            hasher = headers.get("X-Goog-Hash", "")
            crc = [h.split("=", 1)[1] for h in hasher.split(",") if "crc32c" in h]
            if not crc:
                raise NotImplementedError("No crc32c checksum was provided by google!")
            if crc[0] != b64encode(self.crc32c.digest()).decode():
                raise ChecksumError()

    def validate_http_response(self, r):
        return self.validate_headers(r.headers)


def get_consistency_checker(consistency: Optional[str]) -> ConsistencyChecker:
    if consistency == "size":
        return SizeChecker()
    elif consistency == "md5":
        return MD5Checker()
    elif consistency == "crc32c":
        if crcmod is None:
            raise ImportError(
                "The python package `crcmod` is required for `consistency='crc32c'`. "
                "This can be installed with `pip install gcsfs[crc]`"
            )
        else:
            return Crc32cChecker()
    else:
        return ConsistencyChecker()
