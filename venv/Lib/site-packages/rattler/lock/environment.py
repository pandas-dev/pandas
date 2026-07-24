from __future__ import annotations
from typing import Dict, List, Optional, Union

from rattler.channel import Channel
from rattler.lock.channel import LockChannel
from rattler.lock.package import LockedPackage, PypiLockedPackage
from rattler.platform.platform import Platform

from rattler.rattler import PyEnvironment
from rattler.repo_data.record import RepoDataRecord

from rattler.lock.platform import LockPlatform


class Environment:
    """
    Information about a specific environment in the lock-file.
    """

    _env: PyEnvironment

    def __init__(
        self, name: str, requirements: Dict[Platform, List[RepoDataRecord]], channels: List[Union[Channel, LockChannel]]
    ) -> None:
        """
        Create a new environment.
        """
        # Convert LockChannel to Channel for compatibility
        py_channels = []
        for channel in channels:
            if isinstance(channel, LockChannel):
                py_channels.append(Channel(str(channel))._channel)
            else:
                py_channels.append(channel._channel)

        self._env = PyEnvironment(
            name=name,
            # TODO: move this logic to rust
            records={
                platform._inner: [record._record for record in records] for (platform, records) in requirements.items()
            },
            channels=py_channels,
        )

    def platforms(self) -> List[LockPlatform]:
        """
        Returns all the platforms for which we have a locked-down environment.

        Examples
        --------
        ```python
        >>> from rattler import LockFile
        >>> lock_file = LockFile.from_path("../test-data/test.lock")
        >>> env = lock_file.default_environment()
        >>> env.platforms()
        [...]
        >>>
        ```
        """
        return [LockPlatform._from_py_lock_platform(p) for p in self._env.platforms()]

    def channels(self) -> List[LockChannel]:
        """
        Returns the channels that are used by this environment.
        Note that the order of the channels is significant.
        The first channel is the highest priority channel.

        Examples
        --------
        ```python
        >>> from rattler import LockFile
        >>> lock_file = LockFile.from_path("../test-data/test.lock")
        >>> env = lock_file.default_environment()
        >>> env.channels()
        [LockChannel(url="https://conda.anaconda.org/conda-forge/")]
        >>>
        ```
        """
        return [LockChannel._from_py_lock_channel(c) for c in self._env.channels()]

    def packages(self, platform: LockPlatform) -> Optional[List[LockedPackage]]:
        """
        Returns all the packages for a specific platform in this environment.

        Examples
        --------
        ```python
        >>> from rattler import LockFile
        >>> lock_file = LockFile.from_path("../test-data/test.lock")
        >>> env = lock_file.default_environment()
        >>> platform = env.platforms()[0]
        >>> env.packages(platform)[0]
        CondaLockedBinaryPackage(name='tzdata',location='https://conda.anaconda.org/conda-forge/noarch/tzdata-2024a-h0c530f3_0.conda')
        >>>
        ```
        """
        if packages := self._env.packages(platform._inner):
            return [LockedPackage._from_py_locked_package(p) for p in packages]
        return None

    def packages_by_platform(self) -> Dict[LockPlatform, List[LockedPackage]]:
        """
        Returns a list of all packages and platforms defined for this environment.

        Examples
        --------
        ```python
        >>> from rattler import LockFile
        >>> lock_file = LockFile.from_path("../test-data/test.lock")
        >>> env = lock_file.default_environment()
        >>> pkgs = env.packages_by_platform()
        >>> list(pkgs.keys())
        [LockPlatform(...)]
        >>>
        ```
        """
        return {
            LockPlatform._from_py_lock_platform(platform): [LockedPackage._from_py_locked_package(p) for p in packages]
            for (platform, packages) in self._env.packages_by_platform()
        }

    def pypi_packages(
        self,
    ) -> Dict[str, List[PypiLockedPackage]]:
        """
        Returns all pypi packages for all platforms.

        The keys are platform names (e.g., "linux-64", "osx-arm64").

        Examples
        --------
        ```python
        >>> from rattler import LockFile
        >>> lock_file = LockFile.from_path("../test-data/test.lock")
        >>> env = lock_file.default_environment()
        >>> pypi_packages = env.pypi_packages()
        >>> pypi_packages["osx-arm64"][0]
        PypiLockedPackage(name='charset-normalizer',location='https://files.pythonhosted.org/packages/3a/52/9f9d17c3b54dc238de384c4cb5a2ef0e27985b42a0e5cc8e8a31d918d48d/charset_normalizer-3.3.2-cp312-cp312-macosx_11_0_arm64.whl#sha256=55086ee1064215781fff39a1af09518bc9255b50d6333f2e4c74ca09fac6a8f6')
        >>>
        ```
        """
        return {
            platform_name: [PypiLockedPackage._from_py_locked_package(pypi) for pypi in pypi_tup]
            for (platform_name, pypi_tup) in self._env.pypi_packages().items()
        }

    def conda_repodata_records(self) -> Dict[str, List[RepoDataRecord]]:
        """
        Returns all conda packages for all platforms.

        The keys are platform names (e.g., "linux-64", "osx-arm64").

        Examples
        --------
        ```python
        >>> from rattler import LockFile
        >>> lock_file = LockFile.from_path("../test-data/test.lock")
        >>> env = lock_file.default_environment()
        >>> env.conda_repodata_records()
        {'osx-arm64': [RepoDataRecord(...), ...]}
        >>>
        ```
        """
        return {
            platform_name: [RepoDataRecord._from_py_record(r) for r in records]
            for (platform_name, records) in self._env.conda_repodata_records().items()
        }

    def conda_repodata_records_for_platform(self, platform: LockPlatform) -> Optional[List[RepoDataRecord]]:
        """
        Takes all the conda packages, converts them to [`RepoDataRecord`] and returns them or
        returns an error if the conversion failed. Returns `None` if the specified platform is not
        defined for this environment.

        Examples
        --------
        ```python
        >>> from rattler import LockFile
        >>> lock_file = LockFile.from_path("../test-data/test.lock")
        >>> env = lock_file.default_environment()
        >>> platform = env.platforms()[0]
        >>> rdr = env.conda_repodata_records_for_platform(platform)
        >>> rdr
        [...]
        >>> rdr[0]
        RepoDataRecord(...)
        >>>
        ```
        """
        if records := self._env.conda_repodata_records_for_platform(platform._inner):
            return [RepoDataRecord._from_py_record(r) for r in records]
        return None

    def pypi_packages_for_platform(self, platform: LockPlatform) -> Optional[List[PypiLockedPackage]]:
        """
        Returns all the pypi packages and their associated environment data for the specified
        platform. Returns `None` if the platform is not defined for this environment.

        Examples
        --------
        ```python
        >>> from rattler import LockFile
        >>> lock_file = LockFile.from_path("../test-data/test.lock")
        >>> env = lock_file.default_environment()
        >>> platform = env.platforms()[0]
        >>> osx_pypi_pkgs = env.pypi_packages_for_platform(platform)
        >>> osx_pypi_pkgs
        [...]
        >>> osx_pypi_pkgs[0]
        PypiLockedPackage(name='charset-normalizer',location='https://files.pythonhosted.org/packages/3a/52/9f9d17c3b54dc238de384c4cb5a2ef0e27985b42a0e5cc8e8a31d918d48d/charset_normalizer-3.3.2-cp312-cp312-macosx_11_0_arm64.whl#sha256=55086ee1064215781fff39a1af09518bc9255b50d6333f2e4c74ca09fac6a8f6')
        >>>
        ```
        """
        if data := self._env.pypi_packages_for_platform(platform._inner):
            return [PypiLockedPackage._from_py_locked_package(pkg) for pkg in data]
        return None

    @classmethod
    def _from_py_environment(cls, py_environment: PyEnvironment) -> Environment:
        """
        Construct Rattler Environment from FFI PyEnvironment object.
        """
        env = cls.__new__(cls)
        env._env = py_environment
        return env

    def __repr__(self) -> str:
        """
        Returns a representation of the Environment.
        """
        return "Environment()"
