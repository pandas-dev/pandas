from . import version
from .cli import VersionCLI

__version_numeric__ = version.version
__version__ = str(version)


if __name__ == "__main__":
    VersionCLI(version).parse_args()
