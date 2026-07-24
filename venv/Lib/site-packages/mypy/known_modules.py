"""Known Python module names for fuzzy matching import suggestions.

This module provides a curated list of popular Python package import names
for suggesting corrections when a user mistypes an import statement.

Sources:
- Python standard library (typeshed/stdlib/VERSIONS)
- Top 100 PyPI packages by downloads (https://github.com/hugovk/top-pypi-packages)

Note: These are import names, not PyPI package names.
"""

from __future__ import annotations

from typing import Final

from mypy.modulefinder import StdlibVersions

POPULAR_THIRD_PARTY_MODULES: Final[frozenset[str]] = frozenset(
    {
        # Cloud
        "boto3",
        "botocore",
        "aiobotocore",
        "s3transfer",
        "s3fs",
        "awscli",
        # HTTP / Networking
        "urllib3",
        "requests",
        "certifi",
        "idna",
        "charset_normalizer",
        "httpx",
        "httpcore",
        "aiohttp",
        "yarl",
        "multidict",
        "requests_oauthlib",
        "oauthlib",
        "h11",
        # Typing / Extensions
        "typing_extensions",
        "annotated_types",
        "typing_inspection",
        # Core Utilities
        "setuptools",
        "packaging",
        "pip",
        "wheel",
        "virtualenv",
        "platformdirs",
        "filelock",
        "zipp",
        "importlib_metadata",
        # Data Science / Numerical
        "numpy",
        "pandas",
        "scipy",
        "pyarrow",
        # Serialization / Config
        "yaml",
        "pydantic",
        "pydantic_core",
        "attrs",
        "tomli",
        "jsonschema",
        "jsonschema_specifications",
        "jmespath",
        # Cryptography / Security
        "cryptography",
        "cffi",
        "pycparser",
        "rsa",
        "pyjwt",
        "pyasn1",
        "pyasn1_modules",
        # Date / Time
        "dateutil",
        "pytz",
        "tzdata",
        # Google / gRPC
        "google",
        "grpc",
        "grpc_status",
        "grpc_tools",
        "protobuf",
        "googleapis_common_protos",
        # Testing
        "pytest",
        "pluggy",
        "iniconfig",
        # CLI / Terminal
        "click",
        "colorama",
        "rich",
        "tqdm",
        # Web Frameworks
        "starlette",
        # Templates / Markup
        "jinja2",
        "markupsafe",
        "pygments",
        "markdown_it",
        "mdurl",
        # Async
        "anyio",
        "greenlet",
        "aiosignal",
        "aiohappyeyeballs",
        "frozenlist",
        # Database
        "sqlalchemy",
        # Parsing / XML
        "pyparsing",
        "et_xmlfile",
        # OpenTelemetry
        "opentelemetry",
        # Other Popular Modules
        "six",
        "fsspec",
        "wrapt",
        "propcache",
        "rpds",
        "pathspec",
        "PIL",
        "psutil",
        "referencing",
        "trove_classifiers",
        "openpyxl",
        "dotenv",
        "yandexcloud",
        "cachetools",
    }
)


_known_modules_cache: frozenset[str] | None = None


def reset_known_modules_cache() -> None:
    global _known_modules_cache
    _known_modules_cache = None


def get_stdlib_modules(
    stdlib_versions: StdlibVersions, python_version: tuple[int, int] | None = None
) -> frozenset[str]:
    modules: set[str] = set()
    for module, (min_ver, max_ver) in stdlib_versions.items():
        if python_version is not None:
            if python_version < min_ver:
                continue
            if max_ver is not None and python_version > max_ver:
                continue
        top_level = module.split(".")[0]
        # Skip private and very short modules to avoid false positives and noise
        if top_level.startswith("_") or len(top_level) <= 2:
            continue
        modules.add(top_level)
    return frozenset(modules)


def get_known_modules(
    stdlib_versions: StdlibVersions | None = None, python_version: tuple[int, int] | None = None
) -> frozenset[str]:
    global _known_modules_cache
    if _known_modules_cache is not None:
        return _known_modules_cache
    modules: set[str] = set(POPULAR_THIRD_PARTY_MODULES)
    if stdlib_versions is not None:
        modules = modules.union(get_stdlib_modules(stdlib_versions, python_version))
    _known_modules_cache = frozenset(modules)
    return _known_modules_cache
