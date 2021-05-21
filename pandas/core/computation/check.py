from pandas.compat._optional import (
    get_version,
    import_optional_dependency,
)

ne = import_optional_dependency("numexpr", errors="warn")
NUMEXPR_INSTALLED = ne is not None
NUMEXPR_VERSION = get_version(ne) if ne else None

__all__ = ["NUMEXPR_INSTALLED", "NUMEXPR_VERSION"]
