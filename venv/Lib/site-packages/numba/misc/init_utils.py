"""Collection of miscellaneous initialization utilities."""

from collections import namedtuple

version_info = namedtuple('version_info',
                          ('major minor patch short full '
                           'string tuple git_revision'))


def generate_version_info(version):
    """Process a version string into a structured version_info object.

    Parameters
    ----------
    version: str
        a string describing the current version

    Returns
    -------
    version_info: tuple
        structured version information

    See also
    --------
    Look at the definition of 'version_info' in this module for details.

    """
    parts = version.split('.')

    def try_int(x):
        try:
            return int(x)
        except ValueError:
            return None
    major = try_int(parts[0]) if len(parts) >= 1 else None
    minor = try_int(parts[1]) if len(parts) >= 2 else None
    patch = try_int(parts[2]) if len(parts) >= 3 else None
    short = (major, minor)
    full = (major, minor, patch)
    string = version
    tup = tuple(parts)
    git_revision = tup[3] if len(tup) >= 4 else None
    return version_info(major, minor, patch, short, full, string, tup,
                        git_revision)
