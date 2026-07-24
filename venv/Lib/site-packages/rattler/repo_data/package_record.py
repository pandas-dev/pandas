from __future__ import annotations
import os
from typing import Dict, List, Optional, TYPE_CHECKING
import datetime

from rattler import VersionWithSource
from rattler.match_spec.match_spec import MatchSpec
from rattler.package.no_arch_type import NoArchType, NoArchLiteral
from rattler.package.package_name import PackageName
from rattler.platform.platform import Platform
from rattler.rattler import PyRecord, ParsePlatformError

if TYPE_CHECKING:
    import networkx as nx


class PackageRecord:
    """
    A single record in the Conda repodata. A single
    record refers to a single binary distribution
    of a package on a Conda channel.
    """

    _record: PyRecord

    def matches(self, spec: MatchSpec) -> bool:
        """
        Match a [`PackageRecord`] against a [`MatchSpec`].

        Examples
        --------
        ```python
        >>> from rattler import MatchSpec
        >>> record = PackageRecord.from_index_json(
        ...     "../test-data/conda-meta/pysocks-1.7.1-pyh0701188_6.json"
        ... )
        >>> spec = MatchSpec("pysocks")
        >>> record.matches(spec)
        True
        >>> spec = MatchSpec("pysocks>=1.7")
        >>> record.matches(spec)
        True
        >>> spec = MatchSpec("pysocks<1.7")
        >>> record.matches(spec)
        False
        >>>
        ```
        """
        return spec.matches(self)

    @staticmethod
    def from_index_json(
        path: os.PathLike[str],
        size: Optional[int] = None,
        sha256: Optional[str] = None,
        md5: Optional[str] = None,
    ) -> PackageRecord:
        """
        Builds a PackageRecord from an `index.json`.
        These can be found in `info` directory inside an
        extracted package archive.

        Examples
        --------
        ```python
        >>> record = PackageRecord.from_index_json(
        ...     "../test-data/conda-meta/pysocks-1.7.1-pyh0701188_6.json"
        ... )
        >>> assert isinstance(record, PackageRecord)
        >>> record.name
        PackageName("pysocks")
        >>> record.version
        VersionWithSource(version="1.7.1", source="1.7.1")
        >>> record.build
        'pyh0701188_6'
        >>>
        ```
        """
        return PackageRecord._from_py_record(PyRecord.from_index_json(path, size, sha256, md5))

    @staticmethod
    def sort_topologically(records: List[PackageRecord]) -> List[PackageRecord]:
        """
        Sorts the records topologically.
        This function is deterministic, meaning that it will return the same result
        regardless of the order of records and of the depends vector inside the records.
        Note that this function only works for packages with unique names.

        Examples
        --------
        ```python
        >>> from os import listdir
        >>> from os.path import isfile, join
        >>> from rattler import PrefixRecord
        >>> records = [
        ...     PrefixRecord.from_path(join("../test-data/conda-meta/", f))
        ...     for f in listdir("../test-data/conda-meta")
        ...     if isfile(join("../test-data/conda-meta", f))
        ... ]
        >>> sorted = PackageRecord.sort_topologically(records)
        >>> sorted[0].name
        PackageName("python_abi")
        >>> # Verify it's deterministic by sorting again
        >>> sorted2 = PackageRecord.sort_topologically(records)
        >>> [str(r) for r in sorted] == [str(r) for r in sorted2]
        True
        >>>
        ```
        """
        return [PackageRecord._from_py_record(p) for p in PyRecord.sort_topologically(records)]

    @staticmethod
    def to_graph(records: List[PackageRecord]) -> nx.DiGraph:  # type: ignore[type-arg]
        """
        Converts a list of PackageRecords to a DAG (`networkx.DiGraph`).
        The nodes in the graph are the PackageRecords and the edges are the dependencies.

        Note: Virtual packages (starting with `__`) are skipped.
        """
        try:
            import networkx as nx
        except ImportError:
            raise ImportError("networkx is not installed") from None

        names_to_records = {record.name: record for record in records}

        graph = nx.DiGraph()  # type: ignore[var-annotated]
        for record in records:
            graph.add_node(record)
            for dep in record.depends:
                name = dep.split(" ")[0]
                if name.startswith("__"):
                    # this is a virtual package, so we just skip it
                    continue
                graph.add_edge(record, names_to_records[PackageName(name)])

        return graph

    @staticmethod
    def validate(records: List[PackageRecord]) -> None:
        """
        Validate that the given package records are valid w.r.t. 'depends' and 'constrains'.

        This function will return nothing if all records form a valid environment, i.e., all dependencies
        of each package are satisfied by the other packages in the list.
        If there is a dependency that is not satisfied, this function will raise an exception.

        Examples
        --------
        ```python
        >>> from os import listdir
        >>> from os.path import isfile, join
        >>> from rattler import PrefixRecord
        >>> from rattler.exceptions import ValidatePackageRecordsError
        >>> records = [
        ...     PrefixRecord.from_path(join("../test-data/conda-meta/", f))
        ...     for f in sorted(listdir("../test-data/conda-meta"))
        ...     if isfile(join("../test-data/conda-meta", f))
        ... ]
        >>> try:
        ...     PackageRecord.validate(records)
        ... except ValidatePackageRecordsError as e:
        ...     print(e)
        package 'libsqlite=3.40.0=hcfcfb64_0' has dependency 'ucrt >=10.0.20348.0', which is not in the environment
        >>>
        ```
        """
        return PyRecord.validate(records)

    @classmethod
    def _from_py_record(cls, py_record: PyRecord) -> PackageRecord:
        """
        Construct Rattler PackageRecord from FFI PyRecord object.
        """

        # quick sanity check
        assert py_record.is_package_record
        record = cls.__new__(cls)
        record._record = py_record
        return record

    def __init__(
        self,
        name: str | PackageName,
        version: str | VersionWithSource,
        build: str,
        build_number: int,
        subdir: str | Platform,
        arch: Optional[str] = None,
        platform: Optional[str] = None,
        noarch: Optional[NoArchType | NoArchLiteral] = None,
        depends: Optional[List[str]] = None,
        constrains: Optional[List[str]] = None,
        sha256: Optional[bytes] = None,
        md5: Optional[bytes] = None,
        size: Optional[int] = None,
        features: Optional[List[str]] = None,
        legacy_bz2_md5: Optional[bytes] = None,
        legacy_bz2_size: Optional[int] = None,
        license: Optional[str] = None,
        license_family: Optional[str] = None,
        python_site_packages_path: Optional[str] = None,
        extra_depends: Optional[Dict[str, List[str]]] = None,
    ) -> None:
        if isinstance(subdir, str):
            try:
                subdir = Platform(subdir)
            except ParsePlatformError:
                # if the string is not a valid platform, we just keep it as a string
                pass

        if isinstance(subdir, Platform):
            if arch is None:
                subdir_arch = subdir.arch
                arch = str(subdir_arch) if subdir_arch is not None else arch
            platform = subdir.only_platform if platform is None else platform

        # convert str to PackageName
        if isinstance(name, str):
            name = PackageName(name)

        # convert str to VersionWithSource
        if isinstance(version, str):
            version = VersionWithSource(version)

        if not isinstance(noarch, NoArchType):
            noarch = NoArchType(noarch)

        self._record = PyRecord.create(
            name._name,
            (version._version, version._source),
            build,
            build_number,
            str(subdir),
            str(arch) if arch is not None else None,
            platform,
            noarch._noarch,
            python_site_packages_path,
        )

        if constrains is not None:
            self._record.constrains = constrains
        if depends is not None:
            self._record.depends = depends
        if sha256 is not None:
            self._record.sha256 = sha256
        if md5 is not None:
            self._record.md5 = md5
        if size is not None:
            self._record.size = size
        if features is not None:
            self._record.features = features
        if legacy_bz2_md5 is not None:
            self._record.legacy_bz2_md5 = legacy_bz2_md5
        if legacy_bz2_size is not None:
            self._record.legacy_bz2_size = legacy_bz2_size
        if license is not None:
            self._record.license = license
        if license_family is not None:
            self._record.license_family = license_family
        if extra_depends is not None:
            self._record.extra_depends = extra_depends

    @property
    def arch(self) -> Optional[str]:
        """
        Optionally the architecture the package supports.

        Examples
        --------
        ```python
        >>> from rattler import PrefixRecord
        >>> record = PrefixRecord.from_path(
        ...     "../test-data/conda-meta/libsqlite-3.40.0-hcfcfb64_0.json"
        ... )
        >>> record.arch
        'x86_64'
        >>> record.arch = "arm64"
        >>> record.arch
        'arm64'
        >>> record.arch = None
        >>> record.arch is None
        True
        >>>
        ```
        """
        return self._record.arch

    @arch.setter
    def arch(self, value: Optional[str]) -> None:
        self._record.arch = value

    @property
    def build(self) -> str:
        """
        The build string of the package.

        Examples
        --------
        ```python
        >>> from rattler import PrefixRecord
        >>> record = PrefixRecord.from_path(
        ...     "../test-data/conda-meta/libsqlite-3.40.0-hcfcfb64_0.json"
        ... )
        >>> record.build
        'hcfcfb64_0'
        >>> record.build = "new_build_1"
        >>> record.build
        'new_build_1'
        >>>
        ```
        """
        return self._record.build

    @build.setter
    def build(self, value: str) -> None:
        self._record.build = value

    @property
    def build_number(self) -> int:
        """
        The build number of the package.

        Examples
        --------
        ```python
        >>> from rattler import PrefixRecord
        >>> record = PrefixRecord.from_path(
        ...     "../test-data/conda-meta/libsqlite-3.40.0-hcfcfb64_0.json"
        ... )
        >>> record.build_number
        0
        >>> record.build_number = 42
        >>> record.build_number
        42
        >>>
        ```
        """
        return self._record.build_number

    @build_number.setter
    def build_number(self, value: int) -> None:
        self._record.build_number = value

    @property
    def constrains(self) -> List[str]:
        """
        Additional constraints on packages.
        Constrains are different from depends in that packages
        specified in depends must be installed next to this package,
        whereas packages specified in constrains are not required to
        be installed, but if they are installed they must follow
        these constraints.

        Examples
        --------
        ```python
        >>> from rattler import PrefixRecord
        >>> record = PrefixRecord.from_path(
        ...     "../test-data/conda-meta/libsqlite-3.40.0-hcfcfb64_0.json"
        ... )
        >>> record.constrains
        []
        >>> record.constrains = ["python >=3.6"]
        >>> record.constrains
        ['python >=3.6']
        >>>
        ```
        """
        return self._record.constrains

    @constrains.setter
    def constrains(self, value: List[str]) -> None:
        self._record.constrains = value

    @property
    def depends(self) -> List[str]:
        """
        Specification of packages this package depends on.

        Examples
        --------
        ```python
        >>> from rattler import PrefixRecord
        >>> record = PrefixRecord.from_path(
        ...     "../test-data/conda-meta/libsqlite-3.40.0-hcfcfb64_0.json"
        ... )
        >>> record.depends
        ['ucrt >=10.0.20348.0', 'vc >=14.2,<15', 'vs2015_runtime >=14.29.30139']
        >>> record.depends = ["python >=3.6"]
        >>> record.depends
        ['python >=3.6']
        >>>
        ```
        """
        return self._record.depends

    @depends.setter
    def depends(self, value: List[str]) -> None:
        self._record.depends = value

    @property
    def extra_depends(self) -> Dict[str, List[str]]:
        """
        Conditional or optional dependencies. Maps a condition name (e.g. an
        extra/feature name) to a list of dependency specifications that are
        required when that condition is active.

        Examples
        --------
        ```python
        >>> record = PackageRecord(
        ...     name="requests",
        ...     version="2.28.0",
        ...     build="py3-none-any",
        ...     build_number=0,
        ...     subdir="noarch",
        ... )
        >>> record.extra_depends
        {}
        >>> record.extra_depends = {"security": ["cryptography >=3.0"]}
        >>> record.extra_depends
        {'security': ['cryptography >=3.0']}
        >>>
        ```
        """
        return self._record.extra_depends

    @extra_depends.setter
    def extra_depends(self, value: Dict[str, List[str]]) -> None:
        self._record.extra_depends = value

    @property
    def features(self) -> Optional[str]:
        """
        Features are a deprecated way to specify different feature
        sets for the conda solver. This is not supported anymore and
        should not be used. Instead, mutex packages should be used
        to specify mutually exclusive features.

        Examples
        --------
        ```python
        >>> record = PackageRecord.from_index_json(
        ...     "../test-data/conda-meta/pysocks-1.7.1-pyh0701188_6.json"
        ... )
        >>> record.features is None
        True
        >>> record.features = "feature1"
        >>> record.features
        'feature1'
        >>> record.features = None
        >>> record.features is None
        True
        >>>
        ```
        """
        return self._record.features

    @features.setter
    def features(self, value: Optional[str]) -> None:
        self._record.features = value

    @property
    def legacy_bz2_md5(self) -> Optional[bytes]:
        """
        A deprecated md5 hash.

        Examples
        --------
        ```python
        >>> record = PackageRecord.from_index_json(
        ...     "../test-data/conda-meta/pysocks-1.7.1-pyh0701188_6.json"
        ... )
        >>> record.legacy_bz2_md5 is None
        True
        >>> record.legacy_bz2_md5 = bytes.fromhex("2ddbbaf3a82b46ac7214681262e3d746")
        >>> record.legacy_bz2_md5.hex()
        '2ddbbaf3a82b46ac7214681262e3d746'
        >>> record.legacy_bz2_md5 = None
        >>> record.legacy_bz2_md5 is None
        True
        >>>
        ```
        """
        return self._record.legacy_bz2_md5

    @legacy_bz2_md5.setter
    def legacy_bz2_md5(self, value: Optional[bytes]) -> None:
        self._record.legacy_bz2_md5 = value

    @property
    def legacy_bz2_size(self) -> Optional[int]:
        """
        A deprecated package archive size.

        Examples
        --------
        ```python
        >>> record = PackageRecord.from_index_json(
        ...     "../test-data/conda-meta/pysocks-1.7.1-pyh0701188_6.json"
        ... )
        >>> record.legacy_bz2_size is None
        True
        >>> record.legacy_bz2_size = 42
        >>> record.legacy_bz2_size
        42
        >>> record.legacy_bz2_size = None
        >>> record.legacy_bz2_size is None
        True
        >>>
        ```
        """
        return self._record.legacy_bz2_size

    @legacy_bz2_size.setter
    def legacy_bz2_size(self, value: Optional[int]) -> None:
        self._record.legacy_bz2_size = value

    @property
    def license(self) -> Optional[str]:
        """
        The specific license of the package.

        Examples
        --------
        ```python
        >>> from rattler import PrefixRecord
        >>> record = PrefixRecord.from_path(
        ...     "../test-data/conda-meta/libsqlite-3.40.0-hcfcfb64_0.json"
        ... )
        >>> record.license
        'Unlicense'
        >>> record.license = "MIT"
        >>> record.license
        'MIT'
        >>> record.license = None
        >>> record.license is None
        True
        >>>
        ```
        """
        return self._record.license

    @license.setter
    def license(self, value: Optional[str]) -> None:
        self._record.license = value

    @property
    def license_family(self) -> Optional[str]:
        """
        The license family.

        Examples
        --------
        ```python
        >>> from rattler import PrefixRecord
        >>> record = PrefixRecord.from_path(
        ...     "../test-data/conda-meta/pip-23.0-pyhd8ed1ab_0.json"
        ... )
        >>> record.license_family
        'MIT'
        >>> record.license_family = "BSD"
        >>> record.license_family
        'BSD'
        >>> record.license_family = None
        >>> record.license_family is None
        True
        >>>
        ```

        """
        return self._record.license_family

    @license_family.setter
    def license_family(self, value: Optional[str]) -> None:
        self._record.license_family = value

    @property
    def md5(self) -> Optional[bytes]:
        """
        Optionally a MD5 hash of the package archive.

        Examples
        --------
        ```python
        >>> from rattler import PackageRecord
        >>> record = PackageRecord(
        ...     name="foo", version="1.0", build="bar", build_number=0, subdir="linux-64"
        ... )
        >>> record.md5 = bytes.fromhex("5e5a97795de72f8cc3baf3d9ea6327a2")
        >>> record.md5.hex()
        '5e5a97795de72f8cc3baf3d9ea6327a2'
        >>> record.md5 = None
        >>> record.md5 is None
        True
        >>>
        ```
        """
        return self._record.md5

    @md5.setter
    def md5(self, value: Optional[bytes]) -> None:
        self._record.md5 = value

    @property
    def name(self) -> PackageName:
        """
        The name of the package.

        Examples
        --------
        ```python
        >>> from rattler import PrefixRecord, PackageName
        >>> record = PrefixRecord.from_path(
        ...     "../test-data/conda-meta/libsqlite-3.40.0-hcfcfb64_0.json"
        ... )
        >>> record.name
        PackageName("libsqlite")
        >>> record.name = PackageName("newname")
        >>> record.name
        PackageName("newname")
        >>>
        ```
        """
        return PackageName._from_py_package_name(self._record.name)

    @name.setter
    def name(self, value: PackageName) -> None:
        self._record.name = value._name

    @property
    def noarch(self) -> NoArchType:
        """
        The noarch type of the package.

        Examples
        --------
        ```python
        >>> from rattler import PrefixRecord
        >>> record = PrefixRecord.from_path(
        ...     "../test-data/conda-meta/libsqlite-3.40.0-hcfcfb64_0.json"
        ... )
        >>> record.noarch
        NoArchType(None)
        >>> record = PrefixRecord.from_path(
        ...     "../test-data/conda-meta/pip-23.0-pyhd8ed1ab_0.json"
        ... )
        >>> record.noarch
        NoArchType("python")
        >>> record.noarch = NoArchType("generic")
        >>> record.noarch
        NoArchType("generic")
        >>> record.noarch = NoArchType(None)
        >>> record.noarch.none
        True
        >>>
        ```
        """
        return NoArchType._from_py_no_arch_type(self._record.noarch)

    @noarch.setter
    def noarch(self, value: NoArchType) -> None:
        self._record.noarch = value._noarch

    @property
    def platform(self) -> Optional[str]:
        """
        Optionally the platform the package supports.

        Examples
        --------
        ```python
        >>> from rattler import PrefixRecord
        >>> record = PrefixRecord.from_path(
        ...     "../test-data/conda-meta/libsqlite-3.40.0-hcfcfb64_0.json"
        ... )
        >>> record.platform
        'win32'
        >>> record.platform = "linux"
        >>> record.platform
        'linux'
        >>> record.platform = None
        >>> record.platform is None
        True
        >>>
        ```
        """
        return self._record.platform

    @platform.setter
    def platform(self, value: Optional[str]) -> None:
        self._record.platform = value

    @property
    def sha256(self) -> Optional[bytes]:
        """
        Optionally a SHA256 hash of the package archive.

        Examples
        --------
        ```python
        >>> from rattler import PackageRecord
        >>> record = PackageRecord(
        ...     name="foo", version="1.0", build="bar", build_number=0, subdir="linux-64"
        ... )
        >>> record.sha256 = bytes.fromhex("edd7dd24fc070fad8ca690a920d94b6312a376faa96b47c657f9ef5fe5a97dd1")
        >>> record.sha256.hex()
        'edd7dd24fc070fad8ca690a920d94b6312a376faa96b47c657f9ef5fe5a97dd1'
        >>> record.sha256 = None
        >>> record.sha256 is None
        True
        >>>
        ```
        """
        return self._record.sha256

    @sha256.setter
    def sha256(self, value: Optional[bytes]) -> None:
        self._record.sha256 = value

    @property
    def size(self) -> Optional[int]:
        """
        Optionally the size of the package archive in bytes.

        Examples
        --------
        ```python
        >>> from rattler import PrefixRecord
        >>> record = PrefixRecord.from_path(
        ...     "../test-data/conda-meta/libsqlite-3.40.0-hcfcfb64_0.json"
        ... )
        >>> record.size
        669941
        >>> record.size = 42
        >>> record.size
        42
        >>> record.size = None
        >>> record.size is None
        True
        >>>
        ```
        """
        return self._record.size

    @size.setter
    def size(self, value: Optional[int]) -> None:
        self._record.size = value

    @property
    def subdir(self) -> str:
        """
        The subdirectory where the package can be found.

        Examples
        --------
        ```python
        >>> from rattler import PrefixRecord
        >>> record = PrefixRecord.from_path(
        ...     "../test-data/conda-meta/libsqlite-3.40.0-hcfcfb64_0.json"
        ... )
        >>> record.subdir
        'win-64'
        >>> record.subdir = "linux-64"
        >>> record.subdir
        'linux-64'
        >>>
        ```
        """
        return self._record.subdir

    @subdir.setter
    def subdir(self, value: str) -> None:
        self._record.subdir = value

    @property
    def timestamp(self) -> Optional[datetime.datetime]:
        """
        The date this entry was created.

        Examples
        --------
        ```python
        >>> from rattler import PrefixRecord
        >>> import datetime
        >>> record = PrefixRecord.from_path(
        ...     "../test-data/conda-meta/libsqlite-3.40.0-hcfcfb64_0.json"
        ... )
        >>> record.timestamp
        datetime.datetime(2022, 11, 17, 15, 7, 19, 781000, tzinfo=datetime.timezone.utc)
        >>> new_time = datetime.datetime(2023, 1, 1, tzinfo=datetime.timezone.utc)
        >>> record.timestamp = new_time
        >>> record.timestamp
        datetime.datetime(2023, 1, 1, 0, 0, tzinfo=datetime.timezone.utc)
        >>> record.timestamp = None
        >>> record.timestamp is None
        True
        >>>
        ```
        """
        if self._record.timestamp:
            return datetime.datetime.fromtimestamp(self._record.timestamp / 1000.0, tz=datetime.timezone.utc)

        return self._record.timestamp

    @timestamp.setter
    def timestamp(self, value: Optional[datetime.datetime]) -> None:
        if value is not None:
            self._record.timestamp = int(value.timestamp() * 1000)
        else:
            self._record.timestamp = None

    @property
    def track_features(self) -> List[str]:
        """
        Track features are nowadays only used to downweight
        packages (ie. give them less priority).
        To that effect, the number of track features is
        counted (number of commas) and the package is downweighted
        by the number of track_features.

        Examples
        --------
        ```python
        >>> from rattler import PrefixRecord
        >>> record = PrefixRecord.from_path(
        ...     "../test-data/conda-meta/libsqlite-3.40.0-hcfcfb64_0.json"
        ... )
        >>> record.track_features
        []
        >>> record.track_features = ["feature1", "feature2"]
        >>> record.track_features
        ['feature1', 'feature2']
        >>>
        ```
        """
        return self._record.track_features

    @track_features.setter
    def track_features(self, value: List[str]) -> None:
        self._record.track_features = value

    @property
    def version(self) -> VersionWithSource:
        """
        The version of the package.

        Examples
        --------
        ```python
        >>> from rattler import PrefixRecord, VersionWithSource
        >>> record = PrefixRecord.from_path(
        ...     "../test-data/conda-meta/libsqlite-3.40.0-hcfcfb64_0.json"
        ... )
        >>> record.version
        VersionWithSource(version="3.40.0", source="3.40.0")
        >>> record.version = VersionWithSource("1.0.0")
        >>> record.version
        VersionWithSource(version="1.0.0", source="1.0.0")
        >>>
        ```
        """
        return VersionWithSource._from_py_version(*self._record.version)

    @version.setter
    def version(self, value: VersionWithSource) -> None:
        self._record.version = (value._version, value._source)

    @property
    def python_site_packages_path(self) -> Optional[str]:
        """
        Optionally a path within the environment of the site-packages directory. This field is only
        present for python interpreter packages.
        This field was introduced with <https://github.com/conda/ceps/blob/main/cep-17.md>.

        Examples
        --------
        ```python
        >>> from rattler import PrefixRecord
        >>> record = PrefixRecord.from_path(
        ...     "../test-data/conda-meta/python-3.11.9-h932a869_0_cpython.json"
        ... )
        >>> record.python_site_packages_path
        'lib/python3.11/site-packages'
        >>>
        ```
        """
        return self._record.python_site_packages_path

    @python_site_packages_path.setter
    def python_site_packages_path(self, value: Optional[str]) -> None:
        """
        Sets the optional path within the environment of the site-packages directory.
        Examples
        --------
        ```python
        >>> from rattler import PrefixRecord
        >>> record = PrefixRecord.from_path(
        ...     "../test-data/conda-meta/libsqlite-3.40.0-hcfcfb64_0.json"
        ... )
        >>> record.python_site_packages_path
        None
        >>> record.python_site_packages_path = "lib/something"
        >>> record.python_site_packages_path
        'lib/something'
        >>>
        ```
        """
        self._record.python_site_packages_path = value

    def __eq__(self, other: object) -> bool:
        """
        Returns True if the two records are equal.

        Examples
        --------
        ```python
        >>> a = PackageRecord.from_index_json(
        ...     "../test-data/conda-meta/pysocks-1.7.1-pyh0701188_6.json"
        ... )
        >>> b = PackageRecord.from_index_json(
        ...     "../test-data/conda-meta/pysocks-1.7.1-pyh0701188_6.json"
        ... )
        >>> a == b
        True
        >>>
        ```
        """
        if not isinstance(other, PackageRecord):
            return NotImplemented
        return self._record == other._record

    def __ne__(self, other: object) -> bool:
        """
        Returns True if the two records are not equal.

        Examples
        --------
        ```python
        >>> a = PackageRecord.from_index_json(
        ...     "../test-data/conda-meta/pysocks-1.7.1-pyh0701188_6.json"
        ... )
        >>> b = PackageRecord.from_index_json(
        ...     "../test-data/conda-meta/pysocks-1.7.1-pyh0701188_6.json"
        ... )
        >>> a != b
        False
        >>>
        ```
        """
        if not isinstance(other, PackageRecord):
            return NotImplemented
        return self._record != other._record

    def __lt__(self, other: object) -> bool:
        """
        Returns True if this record is less than the other.

        Ordering is defined by package name, track features, version,
        build number, and timestamp.

        Examples
        --------
        ```python
        >>> a = PackageRecord("foo", "1.0", "build_0", 0, "noarch")
        >>> b = PackageRecord("foo", "2.0", "build_0", 0, "noarch")
        >>> a < b
        True
        >>>
        ```
        """
        if not isinstance(other, PackageRecord):
            return NotImplemented
        return self._record < other._record

    def __le__(self, other: object) -> bool:
        """
        Returns True if this record is less than or equal to the other.

        Ordering is defined by package name, track features, version,
        build number, and timestamp.

        Examples
        --------
        ```python
        >>> a = PackageRecord("foo", "1.0", "build_0", 0, "noarch")
        >>> b = PackageRecord("foo", "2.0", "build_0", 0, "noarch")
        >>> a <= b
        True
        >>>
        ```
        """
        if not isinstance(other, PackageRecord):
            return NotImplemented
        return self._record <= other._record

    def __gt__(self, other: object) -> bool:
        """
        Returns True if this record is greater than the other.

        Ordering is defined by package name, track features, version,
        build number, and timestamp.

        Examples
        --------
        ```python
        >>> a = PackageRecord("foo", "2.0", "build_0", 0, "noarch")
        >>> b = PackageRecord("foo", "1.0", "build_0", 0, "noarch")
        >>> a > b
        True
        >>>
        ```
        """
        if not isinstance(other, PackageRecord):
            return NotImplemented
        return self._record > other._record

    def __ge__(self, other: object) -> bool:
        """
        Returns True if this record is greater than or equal to the other.

        Ordering is defined by package name, track features, version,
        build number, and timestamp.

        Examples
        --------
        ```python
        >>> a = PackageRecord("foo", "2.0", "build_0", 0, "noarch")
        >>> b = PackageRecord("foo", "1.0", "build_0", 0, "noarch")
        >>> a >= b
        True
        >>>
        ```
        """
        if not isinstance(other, PackageRecord):
            return NotImplemented
        return self._record >= other._record

    def __hash__(self) -> int:
        """
        Returns the hash of the record.

        Examples
        --------
        ```python
        >>> a = PackageRecord.from_index_json(
        ...     "../test-data/conda-meta/pysocks-1.7.1-pyh0701188_6.json"
        ... )
        >>> isinstance(hash(a), int)
        True
        >>>
        ```
        """
        return hash(self._record)

    def __str__(self) -> str:
        """
        Returns the string representation of the PackageRecord.

        Examples
        --------
        ```python
        >>> record = PackageRecord.from_index_json(
        ...     "../test-data/conda-meta/pysocks-1.7.1-pyh0701188_6.json"
        ... )
        >>> str(record)
        'pysocks=1.7.1=pyh0701188_6'
        >>>
        ```
        """
        return self._record.as_str()

    def __repr__(self) -> str:
        """
        Returns a representation of the PackageRecord.

        Examples
        --------
        ```python
        >>> record = PackageRecord.from_index_json(
        ...     "../test-data/conda-meta/pysocks-1.7.1-pyh0701188_6.json"
        ... )
        >>> record
        PackageRecord("pysocks=1.7.1=pyh0701188_6")
        >>>
        ```
        """
        return f'PackageRecord("{self.__str__()}")'

    def to_json(self) -> str:
        """
        Convert the PackageRecord to a JSON string.

        Examples
        --------
        ```python
        >>> import json
        >>> record = PackageRecord.from_index_json(
        ...     "../test-data/conda-meta/pysocks-1.7.1-pyh0701188_6.json"
        ... )
        >>> json_data = record.to_json()
        >>> isinstance(json_data, str)
        True
        >>> as_dict = json.loads(json_data)
        >>> as_dict["name"]
        'pysocks'
        >>> as_dict["version"]
        '1.7.1'
        >>>
        ```
        """
        return self._record.to_json()
