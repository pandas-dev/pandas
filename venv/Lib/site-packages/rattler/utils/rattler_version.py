try:
    from rattler.rattler import get_rattler_version as _get_rattler_version

    rattler_version_string = _get_rattler_version()
except ImportError:
    # this is only useful for documentation
    import warnings

    warnings.warn("rattler binary missing!", stacklevel=2)
    rattler_version_string = ""


def get_rattler_version() -> str:
    """
    Return the version of the Python Rattler package as a string.

    If the Rattler binary is missing, returns an empty string.
    """
    return rattler_version_string
