from __future__ import annotations
import os
from typing import List, Optional, Tuple, TYPE_CHECKING
from rattler.lock.channel import LockChannel
from rattler.lock.environment import Environment
from rattler.repo_data.record import RepoDataRecord

from rattler.rattler import PyLockFile

if TYPE_CHECKING:
    from rattler.lock.platform import LockPlatform


class LockFile:
    """
    Represents a lock-file for both Conda packages and Pypi packages.
    Lock-files can store information for multiple platforms and for multiple environments.
    """

    _lock_file: PyLockFile

    def __init__(self, platforms: List[LockPlatform]) -> None:
        """
        Create a new rattler-lock file with the given platforms.

        Packages can be added using the `add_conda_package` and `add_pypi_package` methods.
        Channels can be set using the `set_channels` method.

        Args:
            platforms: The list of platforms this lock file supports.
        """
        self._lock_file = PyLockFile([p._inner for p in platforms])

    @staticmethod
    def from_path(path: os.PathLike[str]) -> LockFile:
        """
        Parses a rattler-lock file from a file.

        Examples
        --------
        ```python
        >>> lock_file = LockFile.from_path("./pixi.lock")
        >>> lock_file
        LockFile()
        >>>
        ```
        """
        return LockFile._from_py_lock_file(PyLockFile.from_path(path))

    def to_path(self, path: os.PathLike[str]) -> None:
        """
        Writes the rattler-lock to a file.

        Examples
        --------
        ```python
        >>> import tempfile
        >>> lock_file = LockFile.from_path("./pixi.lock")
        >>> with tempfile.NamedTemporaryFile() as fp:
        ...     lock_file.to_path(fp.name)
        >>>
        ```
        """
        return self._lock_file.to_path(path)

    def environments(self) -> List[Tuple[str, Environment]]:
        """
        Returns an iterator over all environments defined in the lock-file.

        Examples
        --------
        ```python
        >>> lock_file = LockFile.from_path("./pixi.lock")
        >>> sorted(lock_file.environments())
        [('default', Environment()), ('docs', Environment()), ('repl', Environment()), ('test', Environment())]
        >>>
        ```
        """
        return [(name, Environment._from_py_environment(e)) for (name, e) in self._lock_file.environments()]

    def environment(self, name: str) -> Optional[Environment]:
        """
        Returns the environment with the given name.

        Examples
        --------
        ```python
        >>> lock_file = LockFile.from_path("./pixi.lock")
        >>> lock_file.environment("default")
        Environment()
        >>> lock_file.environment("doesnt-exist")
        >>>
        ```
        """
        if env := self._lock_file.environment(name):
            return Environment._from_py_environment(env)
        return None

    def default_environment(self) -> Optional[Environment]:
        """
        Returns the environment with the default name as defined by [`DEFAULT_ENVIRONMENT_NAME`].

        Examples
        --------
        ```python
        >>> lock_file = LockFile.from_path("./pixi.lock")
        >>> lock_file.default_environment()
        Environment()
        >>>
        ```
        """
        return Environment._from_py_environment(self._lock_file.default_environment())

    def platforms(self) -> List[LockPlatform]:
        """
        Returns all platforms defined in the lock-file.

        Examples
        --------
        ```python
        >>> from rattler import LockFile, LockPlatform
        >>> lock_file = LockFile([LockPlatform("linux-64"), LockPlatform("osx-arm64")])
        >>> lock_file.platforms()
        [LockPlatform(name='linux-64'), LockPlatform(name='osx-arm64')]
        >>>
        ```
        """
        from rattler.lock.platform import LockPlatform

        return [LockPlatform._from_py_lock_platform(p) for p in self._lock_file.platforms()]

    def set_channels(self, environment: str, channels: List[LockChannel]) -> None:
        """
        Sets the channels for the given environment.

        Args:
            environment: The name of the environment.
            channels: The list of channels to set.

        Examples
        --------
        ```python
        >>> from rattler import LockFile, LockPlatform, LockChannel
        >>> lock_file = LockFile([LockPlatform("linux-64")])
        >>> lock_file.set_channels("default", [LockChannel("https://conda.anaconda.org/conda-forge/")])
        >>>
        ```
        """
        self._lock_file.set_channels(environment, [c._channel for c in channels])

    def add_conda_package(self, environment: str, platform: LockPlatform, record: RepoDataRecord) -> None:
        """
        Adds a conda package to the lock file for the given environment and platform.

        The platform must be one of the platforms specified when creating the lock file.

        Args:
            environment: The name of the environment.
            platform: The platform to add the package for.
            record: The repo data record of the package to add.

        Raises:
            Exception: If the platform is not one of the platforms in this lock file.

        Examples
        --------
        ```python
        >>> from rattler import LockFile, LockPlatform
        >>> lock_file = LockFile([LockPlatform("linux-64")])
        >>> # lock_file.add_conda_package("default", LockPlatform("linux-64"), some_record)
        >>>
        ```
        """
        self._lock_file.add_conda_package(environment, platform._inner, record._record)

    def add_pypi_package(
        self, environment: str, platform: LockPlatform, name: str, version: str, location: str
    ) -> None:
        """
        Adds a pypi package to the lock file for the given environment and platform.

        The platform must be one of the platforms specified when creating the lock file.

        Args:
            environment: The name of the environment.
            platform: The platform to add the package for.
            name: The name of the package.
            version: The version of the package.
            location: The URL or path where the package can be found.

        Raises:
            Exception: If the platform is not one of the platforms in this lock file.

        Examples
        --------
        ```python
        >>> from rattler import LockFile, LockPlatform
        >>> lock_file = LockFile([LockPlatform("linux-64")])
        >>> # lock_file.add_pypi_package("default", LockPlatform("linux-64"), "requests", "2.28.0", "https://...")
        >>>
        ```
        """
        self._lock_file.add_pypi_package(environment, platform._inner, name, version, location)

    @classmethod
    def _from_py_lock_file(cls, py_lock_file: PyLockFile) -> LockFile:
        """
        Construct Rattler LockFile from FFI PyLockFile object.
        """
        lock_file = cls.__new__(cls)
        lock_file._lock_file = py_lock_file
        return lock_file

    def __repr__(self) -> str:
        """
        Returns a representation of the LockFile.
        """
        return "LockFile()"
