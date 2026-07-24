from __future__ import annotations

import os
import shutil
from argparse import SUPPRESS, ArgumentDefaultsHelpFormatter, ArgumentParser, Namespace
from collections import OrderedDict
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from argparse import Action
    from collections.abc import Mapping, Sequence

from virtualenv.config.convert import get_type
from virtualenv.config.env_var import get_env_var
from virtualenv.config.ini import IniConfig


class VirtualEnvOptions(Namespace):
    def __init__(self, **kwargs: Any) -> None:  # ruff:ignore[any-type]
        super().__init__(**kwargs)
        self._src: str | None = None
        self._sources: dict[str, str] = {}

    def set_src(self, key: str, value: Any, src: str) -> None:  # ruff:ignore[any-type]
        """Set an option value and record where it came from.

        :param key: the option name
        :param value: the option value
        :param src: the source of the value (e.g. ``"cli"``, ``"env var"``, ``"default"``)

        """
        setattr(self, key, value)
        if src.startswith("env var"):
            src = "env var"
        self._sources[key] = src

    def __setattr__(self, key: str, value: Any) -> None:  # ruff:ignore[any-type]
        if (src := getattr(self, "_src", None)) is not None:
            self._sources[key] = src
        super().__setattr__(key, value)

    def get_source(self, key: str) -> str | None:
        """Return the source that provided a given option value.

        :param key: the option name

        :returns: the source string (e.g. ``"cli"``, ``"env var"``, ``"default"``), or ``None`` if not tracked

        """
        return self._sources.get(key)

    @property
    def verbosity(self) -> int | None:
        """The verbosity level, computed as ``verbose - quiet``, clamped to zero.

        :returns: the verbosity level, or ``None`` if neither ``--verbose`` nor ``--quiet`` has been parsed yet

        """
        if not hasattr(self, "verbose") and not hasattr(self, "quiet"):
            return None
        return max(self.verbose - self.quiet, 0)

    def __repr__(self) -> str:
        return f"{type(self).__name__}({', '.join(f'{k}={v}' for k, v in vars(self).items() if not k.startswith('_'))})"


class VirtualEnvConfigParser(ArgumentParser):
    """Custom option parser which updates its defaults by checking the configuration files and environmental vars."""

    def __init__(
        self,
        options: VirtualEnvOptions | None = None,
        env: Mapping[str, str] | None = None,
        *args: Any,  # ruff:ignore[any-type]
        **kwargs: Any,  # ruff:ignore[any-type]
    ) -> None:
        env = os.environ if env is None else env
        self.file_config = IniConfig(env)
        self.epilog_list = []
        self.env = env
        kwargs["epilog"] = self.file_config.epilog
        kwargs["add_help"] = False
        kwargs["formatter_class"] = HelpFormatter
        kwargs["prog"] = "virtualenv"
        super().__init__(*args, **kwargs)
        self._fixed = set()
        if options is not None and not isinstance(options, VirtualEnvOptions):
            msg = "options must be of type VirtualEnvOptions"
            raise TypeError(msg)
        self.options = VirtualEnvOptions() if options is None else options
        self._interpreter = None
        self._app_data = None

    def _fix_defaults(self) -> None:
        for action in self._actions:
            action_id = id(action)
            if action_id not in self._fixed:
                self._fix_default(action)
                self._fixed.add(action_id)

    def _fix_default(self, action: Action) -> None:
        if hasattr(action, "default") and hasattr(action, "dest") and action.default != SUPPRESS:
            as_type = get_type(action)
            names = OrderedDict((i.lstrip("-").replace("-", "_"), None) for i in action.option_strings)
            outcome = None
            for name in names:
                outcome = get_env_var(name, as_type, self.env)
                if outcome is not None:
                    break
            if outcome is None and self.file_config:
                for name in names:
                    outcome = self.file_config.get(name, as_type)
                    if outcome is not None:
                        break
            if outcome is not None:
                action.default, default_source = outcome
                vars(action)["default_source"] = default_source
            else:
                outcome = action.default, "default"
            self.options.set_src(action.dest, *outcome)

    def enable_help(self) -> None:
        self._fix_defaults()
        self.add_argument("-h", "--help", action="help", default=SUPPRESS, help="show this help message and exit")

    def parse_known_args(  # ty: ignore[invalid-method-override]
        self, args: Sequence[str] | None = None, namespace: VirtualEnvOptions | None = None
    ) -> tuple[VirtualEnvOptions, list[str]]:
        if namespace is None:
            namespace = self.options
        elif namespace is not self.options:
            msg = "can only pass in parser.options"
            raise ValueError(msg)
        self._fix_defaults()
        self.options._src = "cli"  # ruff:ignore[private-member-access]
        try:
            namespace.env = self.env
            return super().parse_known_args(args, namespace=namespace)
        finally:
            self.options._src = None  # ruff:ignore[private-member-access]


class HelpFormatter(ArgumentDefaultsHelpFormatter):
    def __init__(self, prog: str, **kwargs: Any) -> None:  # ruff:ignore[any-type]
        super().__init__(prog, max_help_position=32, width=shutil.get_terminal_size().columns, **kwargs)

    def _get_help_string(self, action: Action) -> str | None:
        text = super()._get_help_string(action)
        if text is not None and hasattr(action, "default_source"):
            default = " (default: %(default)s)"
            if text.endswith(default):
                text = f"{text[: -len(default)]} (default: %(default)s -> from %(default_source)s)"
        return text


__all__ = [
    "HelpFormatter",
    "VirtualEnvConfigParser",
    "VirtualEnvOptions",
]
