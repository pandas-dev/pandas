from typing import TypedDict, type_check_only

@type_check_only
class _VersionDict(TypedDict):
    version: str

@type_check_only
class _OptionalVersionDict(TypedDict):
    version: str | None

@type_check_only
class _PlatformDict(TypedDict):
    system: str
    release: str

@type_check_only
class _ImplementationDict(_VersionDict):
    name: str

@type_check_only
class _PyOpenSSLDict(_OptionalVersionDict):
    openssl_version: str

@type_check_only
class _InfoDict(TypedDict):
    platform: _PlatformDict
    implementation: _ImplementationDict
    system_ssl: _VersionDict
    using_pyopenssl: bool
    using_charset_normalizer: bool
    pyOpenSSL: _PyOpenSSLDict
    urllib3: _VersionDict
    chardet: _OptionalVersionDict
    charset_normalizer: _OptionalVersionDict
    cryptography: _VersionDict
    idna: _VersionDict
    requests: _VersionDict

def info() -> _InfoDict: ...
def main() -> None: ...
