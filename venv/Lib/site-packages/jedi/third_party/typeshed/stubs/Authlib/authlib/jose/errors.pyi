from _typeshed import Incomplete

from authlib.common.errors import AuthlibBaseError

class JoseError(AuthlibBaseError): ...

class DecodeError(JoseError):
    error: str

class MissingAlgorithmError(JoseError):
    error: str

class UnsupportedAlgorithmError(JoseError):
    error: str

class BadSignatureError(JoseError):
    error: str
    result: Incomplete
    def __init__(self, result) -> None: ...

class InvalidHeaderParameterNameError(JoseError):
    error: str
    def __init__(self, name) -> None: ...

class InvalidCritHeaderParameterNameError(JoseError):
    error: str
    def __init__(self, name: str) -> None: ...

class InvalidEncryptionAlgorithmForECDH1PUWithKeyWrappingError(JoseError):
    error: str
    def __init__(self) -> None: ...

class InvalidAlgorithmForMultipleRecipientsMode(JoseError):
    error: str
    def __init__(self, alg) -> None: ...

class KeyMismatchError(JoseError):
    error: str
    description: str

class MissingEncryptionAlgorithmError(JoseError):
    error: str
    description: str

class UnsupportedEncryptionAlgorithmError(JoseError):
    error: str
    description: str

class UnsupportedCompressionAlgorithmError(JoseError):
    error: str
    description: str

class InvalidUseError(JoseError):
    error: str
    description: str

class InvalidClaimError(JoseError):
    error: str
    claim_name: Incomplete
    def __init__(self, claim) -> None: ...

class MissingClaimError(JoseError):
    error: str
    def __init__(self, claim) -> None: ...

class InsecureClaimError(JoseError):
    error: str
    def __init__(self, claim) -> None: ...

class ExpiredTokenError(JoseError):
    error: str
    description: str

class InvalidTokenError(JoseError):
    error: str
    description: str
