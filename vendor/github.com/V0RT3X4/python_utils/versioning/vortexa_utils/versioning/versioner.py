"""Class to track the version of a package."""
import os
from dataclasses import dataclass, field, InitVar
from collections import namedtuple
from .utils import cached_property


@dataclass
class Versioner(object):
    VERSION_FILE: str = "VERSION"
    MODULE_FILE: InitVar[str] = None
    VERSION_SEP: str = "."
    VERSION_SPEC: str = "Major.minor.patch"
    __version__: namedtuple = field(default=None, init=False)
    __version_file__: namedtuple = field(default=None, init=False)

    def __post_init__(self, MODULE_FILE):
        parts = []
        if MODULE_FILE is not None:
            dir = os.path.dirname(os.path.abspath(MODULE_FILE))
            parts.append(dir)
        parts.append(self.VERSION_FILE)
        path = os.path.join(*parts)
        self.__version_file__ = os.path.abspath(path)

    @cached_property
    def SemanticVersion(self):
        version_type = namedtuple(
            "SemanticVersion",
            self.VERSION_SPEC.lower().split(self.VERSION_SEP)
        )
        return version_type

    def init_version(self):
        fields = self.SemanticVersion._fields
        version = (
            1 if i == len(fields) - 1 else 0
            for i, field in enumerate(fields)
        )
        self.version = self.SemanticVersion(*version)
        self.write()
        return self.version

    def new_version(self, spec_flags):
        bumped = False
        for spec, ver in zip(spec_flags, self.version):
            if bumped:
                yield 0
            elif spec:
                bumped = True
                yield ver + 1
            else:
                yield ver

    def update_version(self, spec_flags):
        version = self.SemanticVersion(*self.new_version(spec_flags))
        self.version = version
        self.write()
        return version

    def read(self):
        try:
            with open(self.__version_file__, "r") as file:
                version_string = file.readline().strip()
        except FileNotFoundError:
            version = self.init_version()
        else:
            if version_string == "":
                version = self.init_version()
            else:
                version = self.parse_verion(version_string)
        self.version = version
        return version

    def write(self):
        with open(self.__version_file__, "w") as file:
            file.write(str(self))

    @property
    def version(self):
        if self.__version__ is None:
            self.read()
        return self.__version__

    @version.setter
    def version(self, version):
        if isinstance(version, str):
            version = self.parse_verion(version)
        if isinstance(version, self.SemanticVersion):
            self.__version__ = version
        else:
            raise TypeError("Version is not str or self.SemanticVersion")

    def parse_verion(self, version: str):
        parts = (int(v) for v in version.split(self.VERSION_SEP))
        return self.SemanticVersion(*parts)

    def __str__(self):
        return self.VERSION_SEP.join(str(v) for v in self.version)
