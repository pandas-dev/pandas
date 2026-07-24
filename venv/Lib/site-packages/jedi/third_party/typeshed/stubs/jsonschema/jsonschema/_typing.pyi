from collections.abc import Callable, Iterable
from typing import Any, Protocol
from typing_extensions import TypeAlias

from jsonschema.protocols import Validator
from referencing.jsonschema import Schema

class SchemaKeywordValidator(Protocol):
    def __call__(self, validator: Validator, value: Any, instance: Any, schema: Schema) -> None: ...

id_of: TypeAlias = Callable[[Schema], str | None]  # noqa: Y042

ApplicableValidators: TypeAlias = Callable[[Schema], Iterable[tuple[str, Any]]]
