from typing import Final

from authlib.jose.rfc7516 import JWEZipAlgorithm

GZIP_HEAD: Final[bytes]
MAX_SIZE: Final = 256000

class DeflateZipAlgorithm(JWEZipAlgorithm):
    name: str
    description: str
    def compress(self, s: bytes) -> bytes: ...
    def decompress(self, s: bytes) -> bytes: ...

def register_jwe_rfc7518() -> None: ...
