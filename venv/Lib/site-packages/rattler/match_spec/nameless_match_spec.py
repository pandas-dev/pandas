from __future__ import annotations
from typing import TYPE_CHECKING, Optional
from rattler.channel.channel import Channel

from rattler.rattler import PyNamelessMatchSpec

if TYPE_CHECKING:
    from rattler.match_spec import MatchSpec
    from rattler.repo_data import PackageRecord


class NamelessMatchSpec:
    """
    Similar to a `MatchSpec` but does not include the package name.
    This is useful in places where the package name is already known
    (e.g. `foo = "3.4.1 *cuda"`).
    """

    def __init__(
        self,
        spec: str,
        strict: bool = False,
        extras: bool = True,
        conditionals: bool = True,
        flags: bool = True,
    ) -> None:
        """
        Create a new version spec.

        When `strict` is `True`, some ambiguous version specs are rejected.

        When `extras` is `True`, extras syntax (`[extras=[foo,bar]]`) is
        allowed. When `conditionals` is `True`, conditionals syntax
        (`>=1.0[when="python >=3.6"]`) is allowed. When `flags` is `True`,
        flags syntax (`[flags=[cuda]]`) is allowed.

        ```python
        >>> NamelessMatchSpec(">=24.0")
        NamelessMatchSpec(">=24.0")
        >>> NamelessMatchSpec("24")
        NamelessMatchSpec("==24")
        >>>
        ```
        """
        if isinstance(spec, str):
            self._nameless_match_spec = PyNamelessMatchSpec(
                spec,
                strict,
                extras,
                conditionals,
                flags,
            )
        else:
            raise TypeError(
                "NamelessMatchSpec constructor received unsupported type"
                f" {type(spec).__name__!r} for the 'spec' parameter"
            )

    @property
    def version(self) -> Optional[str]:
        """
        The version spec of the package (e.g. `1.2.3`, `>=1.2.3`, `1.2.*`)
        """
        return self._nameless_match_spec.version

    @property
    def build(self) -> Optional[str]:
        """
        The build string of the package (e.g. `py37_0`, `py37h6de7cb9_0`, `py*`)
        """
        return self._nameless_match_spec.build

    @property
    def build_number(self) -> Optional[str]:
        """
        The build number of the package.
        """
        return self._nameless_match_spec.build_number

    @property
    def file_name(self) -> Optional[str]:
        """
        Match the specific filename of the package.
        """
        return self._nameless_match_spec.file_name

    @property
    def channel(self) -> Optional[Channel]:
        """
        The channel of the package.
        """
        channel = self._nameless_match_spec.channel
        return channel and Channel._from_py_channel(channel)

    @property
    def subdir(self) -> Optional[str]:
        """
        The subdir of the channel.
        """
        return self._nameless_match_spec.subdir

    @property
    def namespace(self) -> Optional[str]:
        """
        The namespace of the package.
        """
        return self._nameless_match_spec.namespace

    @property
    def md5(self) -> Optional[bytes]:
        """
        The md5 hash of the package.
        """
        return self._nameless_match_spec.md5

    @property
    def sha256(self) -> Optional[bytes]:
        """
        The sha256 hash of the package.
        """
        return self._nameless_match_spec.sha256

    @property
    def extras(self) -> Optional[list[str]]:
        """The extras (optional dependencies) of the package."""
        return self._nameless_match_spec.extras

    @property
    def condition(self) -> Optional[str]:
        """The condition under which this match spec applies."""
        return self._nameless_match_spec.condition

    def matches(self, package_record: PackageRecord) -> bool:
        """
        Match a MatchSpec against a PackageRecord
        """
        return self._nameless_match_spec.matches(package_record._record)

    @classmethod
    def _from_py_nameless_match_spec(cls, py_nameless_match_spec: PyNamelessMatchSpec) -> NamelessMatchSpec:
        """
        Construct py-rattler NamelessMatchSpec from PyNamelessMatchSpec FFI object.
        """
        nameless_match_spec = cls.__new__(cls)
        nameless_match_spec._nameless_match_spec = py_nameless_match_spec

        return nameless_match_spec

    @classmethod
    def from_match_spec(cls, spec: MatchSpec) -> NamelessMatchSpec:
        """
        Constructs a NamelessMatchSpec from a MatchSpec.

        Examples
        --------
        ```python
        >>> from rattler import MatchSpec
        >>> NamelessMatchSpec.from_match_spec(MatchSpec("foo ==3.4"))
        NamelessMatchSpec("==3.4")
        >>>
        ```
        """
        return cls._from_py_nameless_match_spec(PyNamelessMatchSpec.from_match_spec(spec._match_spec))

    def __str__(self) -> str:
        """
        Returns a string representation of the NamelessMatchSpec.

        Examples
        --------
        ```python
        >>> str(NamelessMatchSpec("3.4"))
        '==3.4'
        >>>
        ```
        """
        return self._nameless_match_spec.as_str()

    def __repr__(self) -> str:
        """
        Returns a representation of the NamelessMatchSpec.

        Examples
        --------
        ```python
        >>> NamelessMatchSpec("3.4")
        NamelessMatchSpec("==3.4")
        >>>
        ```
        """
        return f'NamelessMatchSpec("{self._nameless_match_spec.as_str()}")'
