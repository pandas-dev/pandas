"""The Sphinx documentation toolchain."""

# Keep this file executable as-is in Python 3!
# (Otherwise getting the version out of it when packaging is impossible.)

from __future__ import annotations

from sphinx.util._pathlib import _StrPath

TYPE_CHECKING = False
if TYPE_CHECKING:
    from typing import Final

__version__: Final = '9.1.0'
__display_version__: Final = __version__  # used for command line version

#: Version info for better programmatic use.
#:
#: A tuple of five elements; for Sphinx version 1.2.1 beta 3 this would be
#: ``(1, 2, 1, 'beta', 3)``. The fourth element can be one of: ``alpha``,
#: ``beta``, ``rc``, ``final``. ``final`` always has 0 as the last element.
#:
#: .. versionadded:: 1.2
#:    Before version 1.2, check the string ``sphinx.__version__``.
version_info: Final = (9, 1, 0, 'final', 0)

package_dir: Final = _StrPath(__file__).resolve().parent
del _StrPath

_in_development = False
if _in_development:
    # Only import subprocess if needed
    import subprocess

    try:
        if ret := subprocess.run(
            ('git', 'rev-parse', '--short', 'HEAD'),
            capture_output=True,
            check=False,
            cwd=package_dir,
            encoding='utf-8',
            errors='surrogateescape',
        ).stdout:
            __display_version__ += f'+/{ret.strip()}'  # type: ignore[misc]
        del ret
    finally:
        del subprocess
del _in_development
