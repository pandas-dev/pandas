from __future__ import annotations
from enum import Enum

from typing import Iterable, Optional
from pathlib import Path
import os
from rattler.platform.platform import Platform, PlatformLiteral

from rattler.rattler import (
    PyActivationVariables,
    PyActivator,
    PyShellEnum,
    PyActivationResult,
)


class PathModificationBehavior(Enum):
    """
    The behavior to use when modifying the PATH environment variable.
    Prepend will add the new path to the beginning of the PATH variable.
    Append will add the new path to the end of the PATH variable.
    Replace will replace the entire PATH variable with the new path.
    """

    Prepend = "prepend"
    Append = "append"
    Replace = "replace"


class ActivationVariables:
    """An object that holds the state of the current environment."""

    def __init__(
        self,
        current_prefix: Optional[os.PathLike[str]] = None,
        current_path: Optional[Iterable[str] | Iterable[os.PathLike[str]]] | None = None,
        path_modification_behavior: PathModificationBehavior = PathModificationBehavior.Prepend,
    ) -> None:
        """
        Construct a new ActivationVariables object.

        current_prefix: The current activated conda prefix (usually
            `os.environ["CONDA_PREFIX"]`). This prefix is going to be deactivated.
        current_path: The current PATH environment variable (usually
            `os.environ["PATH"].split(os.pathsep)`).
        path_modification_behavior: The behavior to use when modifying the PATH
            environment variable. One of "Prepend", "Append", or "Replace".
            Defaults to "Prepend".
        """
        self._activation_variables = PyActivationVariables(
            current_prefix,
            current_path or os.environ.get("PATH", "").split(os.pathsep),
            path_modification_behavior.value,
        )

    def __str__(self) -> str:
        return self._activation_variables.as_str()


class ActivationResult:
    """An object that holds the result of activating a conda environment."""

    _py_activation_result: PyActivationResult

    @classmethod
    def _from_py_activation_result(cls, py_activation_result: PyActivationResult) -> ActivationResult:
        """Construct Rattler version from FFI PyActivationResult object."""
        activation_result = cls.__new__(cls)
        activation_result._py_activation_result = py_activation_result
        return activation_result

    @property
    def path(self) -> Path:
        """The new PATH environment variable."""
        return self._py_activation_result.path

    @property
    def script(self) -> str:
        """The script to run to activate the environment."""
        return self._py_activation_result.script


class Shell:
    """An enum of supported shells."""

    bash = PyShellEnum.Bash
    zsh = PyShellEnum.Zsh
    fish = PyShellEnum.Fish
    xonsh = PyShellEnum.Xonsh
    powershell = PyShellEnum.PowerShell
    cmd_exe = PyShellEnum.CmdExe


def activate(
    prefix: Path,
    activation_variables: ActivationVariables,
    shell: Optional[Shell] = None,
    platform: Optional[Platform | PlatformLiteral] = None,
) -> ActivationResult:
    """
    Return an ActivationResult object that contains the new PATH environment variable
    and the script to run to activate the environment.

    Arguments:
        prefix: The path to the conda prefix to activate.
        activation_variables: The current activation variables.
        shell: The shell to generate the activation script for.
        platform: The platform to generate the activation script for.
                  If None, the current platform is used.

    Returns:
        An ActivationResult object containing the new PATH environment variable
        and the script to run to activate the environment.

    Examples
    --------
    ```python
    >>> from rattler.shell import Shell, activate, ActivationVariables
    >>> from rattler.platform import Platform
    >>> from pathlib import Path
    >>> import sys
    >>> p = Path("/path/to/conda/prefix")
    >>> actvars = ActivationVariables()
    >>> a = activate(p, actvars, Shell.xonsh)
    >>> print(a)
    <rattler.shell.shell.ActivationResult object at ...>
    >>>
    ```
    """
    platform = Platform(platform) if isinstance(platform, str) else platform or Platform.current()
    shell = shell or Shell.bash
    return ActivationResult._from_py_activation_result(
        PyActivator.activate(prefix, activation_variables._activation_variables, platform._inner, shell)
    )
