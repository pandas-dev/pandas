from .versioner import Versioner

version = Versioner("../../VERSION", __file__)
__version_numeric__ = version.version
__version__ = str(version)


if __name__ == "__main__":
    from .cli import VersionCLI
    VersionCLI(version).parse_args()
