class FFIError(Exception):
    __module__: str

class CDefError(Exception):
    __module__: str

class VerificationError(Exception):
    __module__: str

class VerificationMissing(Exception):
    __module__: str

class PkgConfigError(Exception):
    __module__: str
