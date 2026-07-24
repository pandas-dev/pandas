from __future__ import annotations

from collections import OrderedDict, defaultdict
from typing import TYPE_CHECKING, NamedTuple

from virtualenv.create.describe import Describe
from virtualenv.create.via_global_ref.builtin.builtin_way import VirtualenvBuiltin

from .base import ComponentBuilder

if TYPE_CHECKING:
    from collections.abc import Sequence

    from python_discovery import PythonInfo

    from virtualenv.config.cli.parser import VirtualEnvConfigParser, VirtualEnvOptions
    from virtualenv.create.creator import Creator, CreatorMeta


class CreatorInfo(NamedTuple):
    key_to_class: dict[str, type[Creator]]
    key_to_meta: dict[str, CreatorMeta]
    describe: type[Describe] | None
    builtin_key: str


class CreatorSelector(ComponentBuilder):
    def __init__(self, interpreter: PythonInfo, parser: VirtualEnvConfigParser) -> None:
        creators, self.key_to_meta, self.describe, self.builtin_key = self.for_interpreter(interpreter)
        super().__init__(interpreter, parser, "creator", creators)  # ty: ignore[invalid-argument-type]

    @classmethod
    def for_interpreter(cls, interpreter: PythonInfo) -> CreatorInfo:
        key_to_class, key_to_meta, builtin_key, describe = OrderedDict(), {}, None, None
        errors = defaultdict(list)
        for key, creator_class in cls.options("virtualenv.create").items():
            if key == "builtin":
                msg = "builtin creator is a reserved name"
                raise RuntimeError(msg)
            meta = creator_class.can_create(interpreter)  # ty: ignore[unresolved-attribute]
            if meta:
                if meta.error:
                    errors[meta.error].append(creator_class)
                else:
                    if "builtin" not in key_to_class and issubclass(creator_class, VirtualenvBuiltin):
                        builtin_key = key
                        key_to_class["builtin"] = creator_class
                        key_to_meta["builtin"] = meta
                    key_to_class[key] = creator_class
                    key_to_meta[key] = meta
            if describe is None and issubclass(creator_class, Describe) and creator_class.can_describe(interpreter):
                describe = creator_class
        if not key_to_meta:
            if errors:
                rows = [f"{k} for creators {', '.join(i.__name__ for i in v)}" for k, v in errors.items()]
                raise RuntimeError("\n".join(rows))
            msg = f"No virtualenv implementation for {interpreter}"
            raise RuntimeError(msg)
        return CreatorInfo(
            key_to_class=key_to_class,
            key_to_meta=key_to_meta,
            describe=describe,
            builtin_key=builtin_key or "",
        )

    def add_selector_arg_parse(self, name: str, choices: Sequence[str]) -> None:
        # prefer the built-in venv if present, otherwise fallback to first defined type
        choices = sorted(choices, key=lambda a: 0 if a == "builtin" else 1)
        default_value = self._get_default(choices)
        self.parser.add_argument(
            f"--{name}",
            choices=choices,
            default=default_value,
            required=False,
            help=f"create environment via{'' if self.builtin_key is None else f' (builtin = {self.builtin_key})'}",
        )

    @staticmethod
    def _get_default(choices: list[str]) -> str:
        return next(iter(choices))

    def populate_selected_argparse(self, selected: str, app_data: object) -> None:
        self.parser.description = f"options for {self.name} {selected}"
        assert self._impl_class is not None  # ruff:ignore[assert]  # Set by handle_selected_arg_parse
        self._impl_class.add_parser_arguments(self.parser, self.interpreter, self.key_to_meta[selected], app_data)  # ty: ignore[unresolved-attribute]

    def create(self, options: VirtualEnvOptions) -> Creator:
        options.meta = self.key_to_meta[getattr(options, self.name)]
        assert self._impl_class is not None  # ruff:ignore[assert]  # Set by handle_selected_arg_parse
        if not issubclass(self._impl_class, Describe):
            options.describe = self.describe(options, self.interpreter)  # ty: ignore[call-non-callable, invalid-argument-type]
        return super().create(options)  # ty: ignore[invalid-return-type]


__all__ = [
    "CreatorInfo",
    "CreatorSelector",
]
