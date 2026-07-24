from __future__ import annotations


class Timeout(TimeoutError):  # ruff:ignore[error-suffix-on-exception-name]  # public exception name; renaming breaks the API
    """Raised when the lock could not be acquired in *timeout* seconds."""

    def __init__(self, lock_file: str) -> None:
        super().__init__()
        self._lock_file = lock_file

    def __reduce__(self) -> tuple[type[Timeout], tuple[str]]:
        # __init__ needs lock_file, so pickle must restore it as a constructor arg
        return self.__class__, (self._lock_file,)

    def __str__(self) -> str:  # pragma: needs hard-link
        return f"The file lock '{self._lock_file}' could not be acquired."

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}({self.lock_file!r})"

    @property
    def lock_file(self) -> str:
        """The path of the file lock."""
        return self._lock_file


class SoftFileLockLifetimeWarning(DeprecationWarning):
    """The configured soft-lock lifetime permits overlapping live holders after expiry."""


class LeaseSettingsMismatch(ValueError):  # ruff:ignore[error-suffix-on-exception-name]  # public exception name; renaming breaks the API
    """A lease contender disagrees with the published claim about how long the lease lasts."""


class SoftFileLockProtocolError(OSError):
    """Raised when strict soft-lock state cannot be interpreted without risking overlap."""

    def __init__(self, lock_file: str, claim_name: str | None, reason: str) -> None:
        self._lock_file = lock_file
        self._claim_name = claim_name
        self._reason = reason
        super().__init__(self.__str__())

    def __reduce__(
        self,
    ) -> tuple[type[SoftFileLockProtocolError], tuple[str, str | None, str]]:  # pragma: needs hard-link
        return self.__class__, (self._lock_file, self._claim_name, self._reason)

    def __str__(self) -> str:
        location = self._lock_file if self._claim_name is None else f"{self._lock_file}: claim {self._claim_name!r}"
        return f"Invalid strict soft-lock state at {location}: {self._reason}"

    @property
    def lock_file(self) -> str:  # pragma: needs hard-link
        """The requested lock path."""
        return self._lock_file

    @property
    def claim_name(self) -> str | None:  # pragma: needs hard-link
        """The claim that caused the error, if scanning identified one."""
        return self._claim_name

    @property
    def reason(self) -> str:  # pragma: needs hard-link
        """The protocol validation failure."""
        return self._reason


__all__ = [
    "LeaseSettingsMismatch",
    "SoftFileLockLifetimeWarning",
    "SoftFileLockProtocolError",
    "Timeout",
]
