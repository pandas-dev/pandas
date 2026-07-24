from .ec_key import ECKey as ECKey
from .jwe_algs import AESAlgorithm as AESAlgorithm, ECDHESAlgorithm as ECDHESAlgorithm, u32be_len_input as u32be_len_input
from .jwe_encs import CBCHS2EncAlgorithm as CBCHS2EncAlgorithm
from .oct_key import OctKey as OctKey
from .rsa_key import RSAKey as RSAKey

__all__ = [
    "register_jws_rfc7518",
    "register_jwe_rfc7518",
    "OctKey",
    "RSAKey",
    "ECKey",
    "u32be_len_input",
    "AESAlgorithm",
    "ECDHESAlgorithm",
    "CBCHS2EncAlgorithm",
]

def register_jws_rfc7518(cls) -> None: ...
def register_jwe_rfc7518(cls) -> None: ...
