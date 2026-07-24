# Licensed to the Apache Software Foundation (ASF) under one
# or more contributor license agreements.  See the NOTICE file
# distributed with this work for additional information
# regarding copyright ownership.  The ASF licenses this file
# to you under the Apache License, Version 2.0 (the
# "License"); you may not use this file except in compliance
# with the License.  You may obtain a copy of the License at
#
#   http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing,
# software distributed under the License is distributed on an
# "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
# KIND, either express or implied.  See the License for the
# specific language governing permissions and limitations
# under the License.

from __future__ import annotations

from abc import ABC, abstractmethod
from collections.abc import Callable, Iterable, Sequence
from functools import cached_property
from typing import Any, TypeAlias
from typing import Literal as TypingLiteral

from pydantic import ConfigDict, Field, SerializeAsAny, model_validator
from pydantic_core.core_schema import ValidatorFunctionWrapHandler

from pyiceberg.expressions.literals import AboveMax, BelowMin, Literal, literal
from pyiceberg.schema import Accessor, Schema
from pyiceberg.typedef import IcebergBaseModel, IcebergRootModel, L, LiteralValue, StructProtocol
from pyiceberg.types import DoubleType, FloatType, NestedField
from pyiceberg.utils.singleton import Singleton


def _to_unbound_term(term: str | UnboundTerm) -> UnboundTerm:
    return Reference(term) if isinstance(term, str) else term


def _to_literal_set(values: Iterable[L] | Iterable[Literal[L]]) -> set[Literal[L]]:
    return {_to_literal(v) for v in values}


def _to_literal(value: L | Literal[L]) -> Literal[L]:
    if isinstance(value, Literal):
        return value
    else:
        return literal(value)


class BooleanExpression(IcebergBaseModel, ABC):
    """An expression that evaluates to a boolean."""

    @abstractmethod
    def __invert__(self) -> BooleanExpression:
        """Transform the Expression into its negated version."""

    def __and__(self, other: BooleanExpression) -> BooleanExpression:
        """Perform and operation on another expression."""
        if not isinstance(other, BooleanExpression):
            raise ValueError(f"Expected BooleanExpression, got: {other}")

        return And(self, other)

    def __or__(self, other: BooleanExpression) -> BooleanExpression:
        """Perform or operation on another expression."""
        if not isinstance(other, BooleanExpression):
            raise ValueError(f"Expected BooleanExpression, got: {other}")

        return Or(self, other)

    @model_validator(mode="wrap")
    @classmethod
    def handle_primitive_type(cls, v: Any, handler: ValidatorFunctionWrapHandler) -> BooleanExpression:
        """Apply custom deserialization logic before validation."""
        # Already a BooleanExpression? return as-is so we keep the concrete subclass.
        if isinstance(v, BooleanExpression):
            return v

        # Handle different input formats
        if isinstance(v, bool):
            return AlwaysTrue() if v is True else AlwaysFalse()

        if isinstance(v, dict) and (field_type := v.get("type")):
            # Unary
            if field_type == "is-null":
                return IsNull(**v)
            elif field_type == "not-null":
                return NotNull(**v)
            elif field_type == "is-nan":
                return IsNaN(**v)
            elif field_type == "not-nan":
                return NotNaN(**v)

            # Literal
            elif field_type == "lt":
                return LessThan(**v)
            elif field_type == "lt-eq":
                return LessThanOrEqual(**v)
            elif field_type == "gt":
                return GreaterThan(**v)
            elif field_type == "gt-eq":
                return GreaterThanOrEqual(**v)
            elif field_type == "eq":
                return EqualTo(**v)
            elif field_type == "not-eq":
                return NotEqualTo(**v)
            elif field_type == "starts-with":
                return StartsWith(**v)
            elif field_type == "not-starts-with":
                return NotStartsWith(**v)

            # Set
            elif field_type == "in":
                return In(**v)
            elif field_type == "not-in":
                return NotIn(**v)

            # Other
            elif field_type == "and":
                return And(**v)
            elif field_type == "or":
                return Or(**v)
            elif field_type == "not":
                return Not(**v)

        return handler(v)


SerializableBooleanExpression: TypeAlias = SerializeAsAny["BooleanExpression"]


def _build_balanced_tree(
    operator_: Callable[[BooleanExpression, BooleanExpression], BooleanExpression], items: Sequence[BooleanExpression]
) -> BooleanExpression:
    """
    Recursively constructs a balanced binary tree of BooleanExpressions using the provided binary operator.

    This function is a safer and more scalable alternative to:
        reduce(operator_, items)

    Using `reduce` creates a deeply nested, unbalanced tree (e.g., operator_(a, operator_(b, operator_(c, ...)))),
    which grows linearly with the number of items. This can lead to RecursionError exceptions in Python
    when the number of expressions is large (e.g., >1000).

    In contrast, this function builds a balanced binary tree with logarithmic depth (O(log n)),
    helping avoid recursion issues and ensuring that expression trees remain stable, predictable,
    and safe to traverse — especially in tools like PyIceberg that operate on large logical trees.

    Parameters:
        operator_ (Callable): A binary operator function (e.g., pyiceberg.expressions.Or, And) that takes two
            BooleanExpressions and returns a combined BooleanExpression.
        items (Sequence[BooleanExpression]): A sequence of BooleanExpression objects to combine.

    Returns:
        BooleanExpression: The balanced combination of all input BooleanExpressions.

    Raises:
        ValueError: If the input sequence is empty.
    """
    if not items:
        raise ValueError("No expressions to combine")
    if len(items) == 1:
        return items[0]
    mid = len(items) // 2

    left = _build_balanced_tree(operator_, items[:mid])
    right = _build_balanced_tree(operator_, items[mid:])
    return operator_(left, right)


class Term:
    """A simple expression that evaluates to a value."""


class Bound:
    """Represents a bound value expression."""


class Unbound(ABC):
    """Represents an unbound value expression."""

    @abstractmethod
    def bind(self, schema: Schema, case_sensitive: bool = True) -> Bound | BooleanExpression: ...

    @property
    @abstractmethod
    def as_bound(self) -> type[Bound]: ...


class BoundTerm(Term, Bound, ABC):
    """Represents a bound term."""

    @abstractmethod
    def ref(self) -> BoundReference:
        """Return the bound reference."""

    @abstractmethod
    def eval(self, struct: StructProtocol) -> Any:  # pylint: disable=W0613
        """Return the value at the referenced field's position in an object that abides by the StructProtocol."""


class BoundReference(BoundTerm):
    """A reference bound to a field in a schema.

    Args:
        field (NestedField): A referenced field in an Iceberg schema.
        accessor (Accessor): An Accessor object to access the value at the field's position.
    """

    field: NestedField
    accessor: Accessor

    def __init__(self, field: NestedField, accessor: Accessor):
        self.field = field
        self.accessor = accessor

    def eval(self, struct: StructProtocol) -> Any:
        """Return the value at the referenced field's position in an object that abides by the StructProtocol.

        Args:
            struct (StructProtocol): A row object that abides by the StructProtocol and returns values given a position.
        Returns:
            Any: The value at the referenced field's position in `struct`.
        """
        return self.accessor.get(struct)

    def __eq__(self, other: Any) -> bool:
        """Return the equality of two instances of the BoundReference class."""
        return self.field == other.field if isinstance(other, BoundReference) else False

    def __repr__(self) -> str:
        """Return the string representation of the BoundReference class."""
        return f"BoundReference(field={repr(self.field)}, accessor={repr(self.accessor)})"

    def ref(self) -> BoundReference:
        return self

    def __hash__(self) -> int:
        """Return hash value of the BoundReference class."""
        return hash(str(self))


class UnboundTerm(Term, Unbound, ABC):
    """Represents an unbound term."""

    @abstractmethod
    def bind(self, schema: Schema, case_sensitive: bool = True) -> BoundTerm: ...


class Reference(UnboundTerm, IcebergRootModel[str]):
    """A reference not yet bound to a field in a schema.

    Args:
        name (str): The name of the field.

    Note:
        An unbound reference is sometimes referred to as a "named" reference.
    """

    root: str = Field()

    def __init__(self, name: str) -> None:
        super().__init__(name)

    def __repr__(self) -> str:
        """Return the string representation of the Reference class."""
        return f"Reference(name={repr(self.root)})"

    def __str__(self) -> str:
        """Return the string representation of the Reference class."""
        return f"Reference(name={repr(self.root)})"

    def bind(self, schema: Schema, case_sensitive: bool = True) -> BoundReference:
        """Bind the reference to an Iceberg schema.

        Args:
            schema (Schema): An Iceberg schema.
            case_sensitive (bool): Whether to consider case when binding the reference to the field.

        Raises:
            ValueError: If an empty name is provided.

        Returns:
            BoundReference: A reference bound to the specific field in the Iceberg schema.
        """
        field = schema.find_field(name_or_id=self.name, case_sensitive=case_sensitive)
        accessor = schema.accessor_for_field(field.field_id)
        return self.as_bound(field=field, accessor=accessor)

    @property
    def name(self) -> str:
        return self.root

    @property
    def as_bound(self) -> type[BoundReference]:
        return BoundReference


class And(BooleanExpression):
    """AND operation expression - logical conjunction."""

    model_config = ConfigDict(arbitrary_types_allowed=True)

    type: TypingLiteral["and"] = Field(default="and", alias="type")
    left: SerializableBooleanExpression = Field()
    right: SerializableBooleanExpression = Field()

    def __init__(self, left: BooleanExpression, right: BooleanExpression, *rest: BooleanExpression, **_: Any) -> None:
        if isinstance(self, And) and not hasattr(self, "left") and not hasattr(self, "right"):
            super().__init__(left=left, right=right)

    def __new__(cls, left: BooleanExpression, right: BooleanExpression, *rest: BooleanExpression, **_: Any) -> BooleanExpression:
        if rest:
            return _build_balanced_tree(And, (left, right, *rest))
        if left is AlwaysFalse() or right is AlwaysFalse():
            return AlwaysFalse()
        elif left is AlwaysTrue():
            return right
        elif right is AlwaysTrue():
            return left
        else:
            return super().__new__(cls)

    def __eq__(self, other: Any) -> bool:
        """Return the equality of two instances of the And class."""
        return self.left == other.left and self.right == other.right if isinstance(other, And) else False

    def __str__(self) -> str:
        """Return the string representation of the And class."""
        return f"And(left={str(self.left)}, right={str(self.right)})"

    def __repr__(self) -> str:
        """Return the string representation of the And class."""
        return f"And(left={repr(self.left)}, right={repr(self.right)})"

    def __invert__(self) -> BooleanExpression:
        """Transform the Expression into its negated version."""
        # De Morgan's law: not (A and B) = (not A) or (not B)
        return Or(~self.left, ~self.right)

    def __getnewargs__(self) -> tuple[BooleanExpression, BooleanExpression]:
        """Pickle the And class."""
        return (self.left, self.right)


class Or(BooleanExpression):
    """OR operation expression - logical disjunction."""

    model_config = ConfigDict(arbitrary_types_allowed=True)

    type: TypingLiteral["or"] = Field(default="or", alias="type")
    left: SerializableBooleanExpression = Field()
    right: SerializableBooleanExpression = Field()

    def __init__(self, left: BooleanExpression, right: BooleanExpression, *rest: BooleanExpression, **_: Any) -> None:
        if isinstance(self, Or) and not hasattr(self, "left") and not hasattr(self, "right"):
            super().__init__(left=left, right=right)

    def __new__(cls, left: BooleanExpression, right: BooleanExpression, *rest: BooleanExpression, **_: Any) -> BooleanExpression:
        if rest:
            return _build_balanced_tree(Or, (left, right, *rest))
        if left is AlwaysTrue() or right is AlwaysTrue():
            return AlwaysTrue()
        elif left is AlwaysFalse():
            return right
        elif right is AlwaysFalse():
            return left
        else:
            return super().__new__(cls)

    def __str__(self) -> str:
        """Return the string representation of the Or class."""
        return f"{str(self.__class__.__name__)}(left={repr(self.left)}, right={repr(self.right)})"

    def __eq__(self, other: Any) -> bool:
        """Return the equality of two instances of the Or class."""
        return self.left == other.left and self.right == other.right if isinstance(other, Or) else False

    def __repr__(self) -> str:
        """Return the string representation of the Or class."""
        return f"Or(left={repr(self.left)}, right={repr(self.right)})"

    def __invert__(self) -> BooleanExpression:
        """Transform the Expression into its negated version."""
        # De Morgan's law: not (A or B) = (not A) and (not B)
        return And(~self.left, ~self.right)

    def __getnewargs__(self) -> tuple[BooleanExpression, BooleanExpression]:
        """Pickle the Or class."""
        return (self.left, self.right)


class Not(BooleanExpression):
    """NOT operation expression - logical negation."""

    model_config = ConfigDict(arbitrary_types_allowed=True)

    type: TypingLiteral["not"] = Field(default="not")
    child: SerializableBooleanExpression = Field()

    def __init__(self, child: BooleanExpression, **_: Any) -> None:
        super().__init__(child=child)

    def __new__(cls, child: BooleanExpression, **_: Any) -> BooleanExpression:
        if child is AlwaysTrue():
            return AlwaysFalse()
        elif child is AlwaysFalse():
            return AlwaysTrue()
        elif isinstance(child, Not):
            return child.child
        else:
            return super().__new__(cls)

    def __str__(self) -> str:
        """Return the string representation of the Not class."""
        return f"Not(child={self.child})"

    def __repr__(self) -> str:
        """Return the string representation of the Not class."""
        return f"Not(child={repr(self.child)})"

    def __eq__(self, other: Any) -> bool:
        """Return the equality of two instances of the Not class."""
        return self.child == other.child if isinstance(other, Not) else False

    def __invert__(self) -> BooleanExpression:
        """Transform the Expression into its negated version."""
        return self.child

    def __getnewargs__(self) -> tuple[BooleanExpression]:
        """Pickle the Not class."""
        return (self.child,)


class AlwaysTrue(BooleanExpression, Singleton, IcebergRootModel[bool]):
    """TRUE expression."""

    root: bool = True

    def __invert__(self) -> AlwaysFalse:
        """Transform the Expression into its negated version."""
        return AlwaysFalse()

    def __str__(self) -> str:
        """Return the string representation of the AlwaysTrue class."""
        return "AlwaysTrue()"

    def __repr__(self) -> str:
        """Return the string representation of the AlwaysTrue class."""
        return "AlwaysTrue()"


class AlwaysFalse(BooleanExpression, Singleton, IcebergRootModel[bool]):
    """FALSE expression."""

    root: bool = False

    def __invert__(self) -> AlwaysTrue:
        """Transform the Expression into its negated version."""
        return AlwaysTrue()

    def __str__(self) -> str:
        """Return the string representation of the AlwaysFalse class."""
        return "AlwaysFalse()"

    def __repr__(self) -> str:
        """Return the string representation of the AlwaysFalse class."""
        return "AlwaysFalse()"


class BoundPredicate(Bound, BooleanExpression, ABC):
    model_config = ConfigDict(arbitrary_types_allowed=True)

    term: BoundTerm

    def __init__(self, term: BoundTerm, **kwargs: Any) -> None:
        super().__init__(term=term, **kwargs)

    def __eq__(self, other: Any) -> bool:
        """Return the equality of two instances of the BoundPredicate class."""
        if isinstance(other, self.__class__):
            return self.term == other.term
        return False

    def __str__(self) -> str:
        """Return the string representation of the BoundPredicate class."""
        return f"{self.__class__.__name__}(term={str(self.term)})"

    @property
    @abstractmethod
    def as_unbound(self) -> type[UnboundPredicate]: ...


class UnboundPredicate(Unbound, BooleanExpression, ABC):
    model_config = ConfigDict(arbitrary_types_allowed=True)

    term: UnboundTerm

    def __init__(self, term: str | UnboundTerm, **kwargs: Any) -> None:
        super().__init__(term=_to_unbound_term(term), **kwargs)

    def __eq__(self, other: Any) -> bool:
        """Return the equality of two instances of the UnboundPredicate class."""
        return self.term == other.term if isinstance(other, self.__class__) else False

    @abstractmethod
    def bind(self, schema: Schema, case_sensitive: bool = True) -> BooleanExpression: ...

    @property
    @abstractmethod
    def as_bound(self) -> type[BoundPredicate]: ...


class UnaryPredicate(UnboundPredicate, ABC):
    type: TypingLiteral["is-null", "not-null", "is-nan", "not-nan"] = Field()

    model_config = {"arbitrary_types_allowed": True}

    def __init__(self, term: str | UnboundTerm, **_: Any) -> None:
        unbound = _to_unbound_term(term)
        super().__init__(term=unbound)

    def __str__(self) -> str:
        """Return the string representation of the UnaryPredicate class."""
        # Sort to make it deterministic
        return f"{str(self.__class__.__name__)}(term={str(self.term)})"

    def bind(self, schema: Schema, case_sensitive: bool = True) -> BoundUnaryPredicate:
        bound_term = self.term.bind(schema, case_sensitive)
        bound_type = self.as_bound
        return bound_type(bound_term)  # type: ignore[misc]

    def __repr__(self) -> str:
        """Return the string representation of the UnaryPredicate class."""
        return f"{str(self.__class__.__name__)}(term={repr(self.term)})"

    @property
    @abstractmethod
    def as_bound(self) -> type[BoundUnaryPredicate]: ...  # type: ignore


class BoundUnaryPredicate(BoundPredicate, ABC):
    def __repr__(self) -> str:
        """Return the string representation of the BoundUnaryPredicate class."""
        return f"{str(self.__class__.__name__)}(term={repr(self.term)})"

    @property
    @abstractmethod
    def as_unbound(self) -> type[UnaryPredicate]: ...

    def __getnewargs__(self) -> tuple[BoundTerm]:
        """Pickle the BoundUnaryPredicate class."""
        return (self.term,)


class BoundIsNull(BoundUnaryPredicate):
    def __new__(cls, term: BoundTerm) -> BooleanExpression:  # pylint: disable=W0221
        if term.ref().field.required:
            return AlwaysFalse()
        return super().__new__(cls)

    def __invert__(self) -> BoundNotNull:
        """Transform the Expression into its negated version."""
        return BoundNotNull(self.term)

    @property
    def as_unbound(self) -> type[IsNull]:
        return IsNull


class BoundNotNull(BoundUnaryPredicate):
    def __new__(cls, term: BoundTerm) -> BooleanExpression:  # pylint: disable=W0221
        if term.ref().field.required:
            return AlwaysTrue()
        return super().__new__(cls)

    def __invert__(self) -> BoundIsNull:
        """Transform the Expression into its negated version."""
        return BoundIsNull(self.term)

    @property
    def as_unbound(self) -> type[NotNull]:
        return NotNull


class IsNull(UnaryPredicate):
    type: TypingLiteral["is-null"] = Field(default="is-null")

    def __invert__(self) -> NotNull:
        """Transform the Expression into its negated version."""
        return NotNull(self.term)

    @property
    def as_bound(self) -> type[BoundIsNull]:  # type: ignore
        return BoundIsNull


class NotNull(UnaryPredicate):
    type: TypingLiteral["not-null"] = Field(default="not-null")

    def __invert__(self) -> IsNull:
        """Transform the Expression into its negated version."""
        return IsNull(self.term)

    @property
    def as_bound(self) -> type[BoundNotNull]:  # type: ignore
        return BoundNotNull


class BoundIsNaN(BoundUnaryPredicate):
    def __new__(cls, term: BoundTerm) -> BooleanExpression:  # pylint: disable=W0221
        bound_type = term.ref().field.field_type
        if isinstance(bound_type, (FloatType, DoubleType)):
            return super().__new__(cls)
        return AlwaysFalse()

    def __invert__(self) -> BoundNotNaN:
        """Transform the Expression into its negated version."""
        return BoundNotNaN(self.term)

    @property
    def as_unbound(self) -> type[IsNaN]:
        return IsNaN


class BoundNotNaN(BoundUnaryPredicate):
    def __new__(cls, term: BoundTerm) -> BooleanExpression:  # pylint: disable=W0221
        bound_type = term.ref().field.field_type
        if isinstance(bound_type, (FloatType, DoubleType)):
            return super().__new__(cls)
        return AlwaysTrue()

    def __invert__(self) -> BoundIsNaN:
        """Transform the Expression into its negated version."""
        return BoundIsNaN(self.term)

    @property
    def as_unbound(self) -> type[NotNaN]:
        return NotNaN


class IsNaN(UnaryPredicate):
    type: TypingLiteral["is-nan"] = Field(default="is-nan")

    def __invert__(self) -> NotNaN:
        """Transform the Expression into its negated version."""
        return NotNaN(self.term)

    @property
    def as_bound(self) -> type[BoundIsNaN]:  # type: ignore
        return BoundIsNaN


class NotNaN(UnaryPredicate):
    type: TypingLiteral["not-nan"] = Field(default="not-nan")

    def __invert__(self) -> IsNaN:
        """Transform the Expression into its negated version."""
        return IsNaN(self.term)

    @property
    def as_bound(self) -> type[BoundNotNaN]:  # type: ignore
        return BoundNotNaN


class SetPredicate(UnboundPredicate, ABC):
    model_config = ConfigDict(arbitrary_types_allowed=True)

    type: TypingLiteral["in", "not-in"] = Field(default="in")
    literals: set[LiteralValue] = Field(alias="values")

    def __init__(
        self, term: str | UnboundTerm, literals: Iterable[Any] | Iterable[LiteralValue] | None = None, **kwargs: Any
    ) -> None:
        if literals is None and "values" in kwargs:
            literals = kwargs["values"]

        if literals is None:
            literal_set: set[LiteralValue] = set()
        else:
            literal_set = _to_literal_set(literals)
        super().__init__(term=_to_unbound_term(term), values=literal_set)

    def bind(self, schema: Schema, case_sensitive: bool = True) -> BoundSetPredicate:
        bound_term = self.term.bind(schema, case_sensitive)
        literal_set = self.literals
        return self.as_bound(bound_term, {lit.to(bound_term.ref().field.field_type) for lit in literal_set})  # type: ignore

    def __str__(self) -> str:
        """Return the string representation of the SetPredicate class."""
        # Sort to make it deterministic
        literals_str = ", ".join(sorted([str(literal) for literal in self.literals]))
        return f"{str(self.__class__.__name__)}({str(self.term)}, {{{literals_str}}})"

    def __repr__(self) -> str:
        """Return the string representation of the SetPredicate class."""
        # Sort to make it deterministic
        literals_repr = ", ".join(sorted([repr(literal) for literal in self.literals]))
        return f"{str(self.__class__.__name__)}({repr(self.term)}, {{{literals_repr}}})"

    def __eq__(self, other: Any) -> bool:
        """Return the equality of two instances of the SetPredicate class."""
        return self.term == other.term and self.literals == other.literals if isinstance(other, self.__class__) else False

    def __getnewargs__(self) -> tuple[UnboundTerm, set[Any]]:
        """Pickle the SetPredicate class."""
        return (self.term, self.literals)

    @property
    @abstractmethod
    def as_bound(self) -> type[BoundSetPredicate]:  # type: ignore
        return BoundSetPredicate


class BoundSetPredicate(BoundPredicate, ABC):
    literals: set[LiteralValue]

    def __init__(self, term: BoundTerm, literals: set[LiteralValue]) -> None:
        literal_set = _to_literal_set(literals)
        super().__init__(term=term, literals=literal_set)

    @cached_property
    def value_set(self) -> set[Any]:
        return {lit.value for lit in self.literals}

    def __str__(self) -> str:
        """Return the string representation of the BoundSetPredicate class."""
        # Sort to make it deterministic
        literals_str = ", ".join(sorted([str(literal) for literal in self.literals]))
        return f"{str(self.__class__.__name__)}({str(self.term)}, {{{literals_str}}})"

    def __repr__(self) -> str:
        """Return the string representation of the BoundSetPredicate class."""
        # Sort to make it deterministic
        literals_repr = ", ".join(sorted([repr(literal) for literal in self.literals]))
        return f"{str(self.__class__.__name__)}({repr(self.term)}, {{{literals_repr}}})"

    def __eq__(self, other: Any) -> bool:
        """Return the equality of two instances of the BoundSetPredicate class."""
        return self.term == other.term and self.literals == other.literals if isinstance(other, self.__class__) else False

    def __getnewargs__(self) -> tuple[BoundTerm, set[LiteralValue]]:
        """Pickle the BoundSetPredicate class."""
        return (self.term, self.literals)

    @property
    @abstractmethod
    def as_unbound(self) -> type[SetPredicate]: ...


class BoundIn(BoundSetPredicate):
    def __new__(cls, term: BoundTerm, literals: set[LiteralValue]) -> BooleanExpression:  # pylint: disable=W0221
        count = len(literals)
        if count == 0:
            return AlwaysFalse()
        elif count == 1:
            return BoundEqualTo(term, next(iter(literals)))
        else:
            return super().__new__(cls)

    def __invert__(self) -> BoundNotIn:
        """Transform the Expression into its negated version."""
        return BoundNotIn(self.term, self.literals)

    def __eq__(self, other: Any) -> bool:
        """Return the equality of two instances of the BoundIn class."""
        return self.term == other.term and self.literals == other.literals if isinstance(other, self.__class__) else False

    @property
    def as_unbound(self) -> type[In]:
        return In


class BoundNotIn(BoundSetPredicate):
    def __new__(  # pylint: disable=W0221
        cls,
        term: BoundTerm,
        literals: set[LiteralValue],
    ) -> BooleanExpression:
        count = len(literals)
        if count == 0:
            return AlwaysTrue()
        elif count == 1:
            return BoundNotEqualTo(term, next(iter(literals)))
        else:
            return super().__new__(cls)

    def __invert__(self) -> BoundIn:
        """Transform the Expression into its negated version."""
        return BoundIn(self.term, self.literals)

    @property
    def as_unbound(self) -> type[NotIn]:
        return NotIn


class In(SetPredicate):
    type: TypingLiteral["in"] = Field(default="in", alias="type")

    def __new__(  # pylint: disable=W0221
        cls, term: str | UnboundTerm, literals: Iterable[Any] | Iterable[LiteralValue] | None = None, **kwargs: Any
    ) -> In:
        if literals is None and "values" in kwargs:
            literals = kwargs["values"]

        if literals is None:
            literals_set: set[LiteralValue] = set()
        else:
            literals_set = _to_literal_set(literals)
        count = len(literals_set)
        if count == 0:
            return AlwaysFalse()
        elif count == 1:
            return EqualTo(term, next(iter(literals_set)))
        else:
            return super().__new__(cls)

    def __invert__(self) -> NotIn:
        """Transform the Expression into its negated version."""
        return NotIn(self.term, self.literals)

    @property
    def as_bound(self) -> type[BoundIn]:  # type: ignore
        return BoundIn


class NotIn(SetPredicate, ABC):
    type: TypingLiteral["not-in"] = Field(default="not-in", alias="type")

    def __new__(  # pylint: disable=W0221
        cls, term: str | UnboundTerm, literals: Iterable[Any] | Iterable[LiteralValue] | None = None, **kwargs: Any
    ) -> NotIn:
        if literals is None and "values" in kwargs:
            literals = kwargs["values"]

        if literals is None:
            literals_set: set[LiteralValue] = set()
        else:
            literals_set = _to_literal_set(literals)
        count = len(literals_set)
        if count == 0:
            return AlwaysTrue()
        elif count == 1:
            return NotEqualTo(term, next(iter(literals_set)))
        else:
            return super().__new__(cls)

    def __invert__(self) -> In:
        """Transform the Expression into its negated version."""
        return In(self.term, self.literals)

    @property
    def as_bound(self) -> type[BoundNotIn]:  # type: ignore
        return BoundNotIn


class LiteralPredicate(UnboundPredicate, ABC):
    type: TypingLiteral["lt", "lt-eq", "gt", "gt-eq", "eq", "not-eq", "starts-with", "not-starts-with"] = Field(alias="type")
    term: UnboundTerm
    value: LiteralValue = Field()
    model_config = ConfigDict(populate_by_name=True, frozen=True, arbitrary_types_allowed=True)

    def __init__(self, term: str | UnboundTerm, literal: Any | None = None, **kwargs: Any) -> None:
        if literal is None and "value" in kwargs:
            literal = kwargs["value"]

        super().__init__(term=_to_unbound_term(term), value=_to_literal(literal))

    @property
    def literal(self) -> LiteralValue:
        return self.value

    def bind(self, schema: Schema, case_sensitive: bool = True) -> BoundLiteralPredicate:
        bound_term = self.term.bind(schema, case_sensitive)
        lit = self.literal.to(bound_term.ref().field.field_type)

        if isinstance(lit, AboveMax):
            if isinstance(self, (LessThan, LessThanOrEqual, NotEqualTo)):
                return AlwaysTrue()
            elif isinstance(self, (GreaterThan, GreaterThanOrEqual, EqualTo)):
                return AlwaysFalse()
        elif isinstance(lit, BelowMin):
            if isinstance(self, (GreaterThan, GreaterThanOrEqual, NotEqualTo)):
                return AlwaysTrue()
            elif isinstance(self, (LessThan, LessThanOrEqual, EqualTo)):
                return AlwaysFalse()

        return self.as_bound(bound_term, lit)  # type: ignore

    def __eq__(self, other: Any) -> bool:
        """Return the equality of two instances of the LiteralPredicate class."""
        if isinstance(other, self.__class__):
            return self.term == other.term and self.literal == other.literal
        return False

    def __str__(self) -> str:
        """Return the string representation of the LiteralPredicate class."""
        return f"{str(self.__class__.__name__)}(term={repr(self.term)}, literal={repr(self.literal)})"

    def __repr__(self) -> str:
        """Return the string representation of the LiteralPredicate class."""
        return f"{str(self.__class__.__name__)}(term={repr(self.term)}, literal={repr(self.literal)})"

    @property
    @abstractmethod
    def as_bound(self) -> type[BoundLiteralPredicate]: ...  # type: ignore


class BoundLiteralPredicate(BoundPredicate, ABC):
    literal: LiteralValue

    def __init__(self, term: BoundTerm, literal: LiteralValue):  # pylint: disable=W0621
        super().__init__(term=term, literal=literal)

    def __eq__(self, other: Any) -> bool:
        """Return the equality of two instances of the BoundLiteralPredicate class."""
        if isinstance(other, self.__class__):
            return self.term == other.term and self.literal == other.literal
        return False

    def __str__(self) -> str:
        """Return the string representation of the BoundLiteralPredicate class."""
        return f"{self.__class__.__name__}(term={str(self.term)}, literal={repr(self.literal)})"

    def __repr__(self) -> str:
        """Return the string representation of the BoundLiteralPredicate class."""
        return f"{str(self.__class__.__name__)}(term={repr(self.term)}, literal={repr(self.literal)})"

    @property
    @abstractmethod
    def as_unbound(self) -> type[LiteralPredicate]: ...


class BoundEqualTo(BoundLiteralPredicate):
    def __invert__(self) -> BoundNotEqualTo:
        """Transform the Expression into its negated version."""
        return BoundNotEqualTo(self.term, self.literal)

    @property
    def as_unbound(self) -> type[EqualTo]:
        return EqualTo


class BoundNotEqualTo(BoundLiteralPredicate):
    def __invert__(self) -> BoundEqualTo:
        """Transform the Expression into its negated version."""
        return BoundEqualTo(self.term, self.literal)

    @property
    def as_unbound(self) -> type[NotEqualTo]:
        return NotEqualTo


class BoundGreaterThanOrEqual(BoundLiteralPredicate):
    def __invert__(self) -> BoundLessThan:
        """Transform the Expression into its negated version."""
        return BoundLessThan(self.term, self.literal)

    @property
    def as_unbound(self) -> type[GreaterThanOrEqual]:
        return GreaterThanOrEqual


class BoundGreaterThan(BoundLiteralPredicate):
    def __invert__(self) -> BoundLessThanOrEqual:
        """Transform the Expression into its negated version."""
        return BoundLessThanOrEqual(self.term, self.literal)

    @property
    def as_unbound(self) -> type[GreaterThan]:
        return GreaterThan


class BoundLessThan(BoundLiteralPredicate):
    def __invert__(self) -> BoundGreaterThanOrEqual:
        """Transform the Expression into its negated version."""
        return BoundGreaterThanOrEqual(self.term, self.literal)

    @property
    def as_unbound(self) -> type[LessThan]:
        return LessThan


class BoundLessThanOrEqual(BoundLiteralPredicate):
    def __invert__(self) -> BoundGreaterThan:
        """Transform the Expression into its negated version."""
        return BoundGreaterThan(self.term, self.literal)

    @property
    def as_unbound(self) -> type[LessThanOrEqual]:
        return LessThanOrEqual


class BoundStartsWith(BoundLiteralPredicate):
    def __invert__(self) -> BoundNotStartsWith:
        """Transform the Expression into its negated version."""
        return BoundNotStartsWith(self.term, self.literal)

    @property
    def as_unbound(self) -> type[StartsWith]:
        return StartsWith


class BoundNotStartsWith(BoundLiteralPredicate):
    def __invert__(self) -> BoundStartsWith:
        """Transform the Expression into its negated version."""
        return BoundStartsWith(self.term, self.literal)

    @property
    def as_unbound(self) -> type[NotStartsWith]:
        return NotStartsWith


class EqualTo(LiteralPredicate):
    type: TypingLiteral["eq"] = Field(default="eq", alias="type")

    def __invert__(self) -> NotEqualTo:
        """Transform the Expression into its negated version."""
        return NotEqualTo(self.term, self.literal)

    @property
    def as_bound(self) -> type[BoundEqualTo]:  # type: ignore
        return BoundEqualTo


class NotEqualTo(LiteralPredicate):
    type: TypingLiteral["not-eq"] = Field(default="not-eq", alias="type")

    def __invert__(self) -> EqualTo:
        """Transform the Expression into its negated version."""
        return EqualTo(self.term, self.literal)

    @property
    def as_bound(self) -> type[BoundNotEqualTo]:  # type: ignore
        return BoundNotEqualTo


class LessThan(LiteralPredicate):
    type: TypingLiteral["lt"] = Field(default="lt", alias="type")

    def __invert__(self) -> GreaterThanOrEqual:
        """Transform the Expression into its negated version."""
        return GreaterThanOrEqual(self.term, self.literal)

    @property
    def as_bound(self) -> type[BoundLessThan]:  # type: ignore
        return BoundLessThan


class GreaterThanOrEqual(LiteralPredicate):
    type: TypingLiteral["gt-eq"] = Field(default="gt-eq", alias="type")

    def __invert__(self) -> LessThan:
        """Transform the Expression into its negated version."""
        return LessThan(self.term, self.literal)

    @property
    def as_bound(self) -> type[BoundGreaterThanOrEqual]:  # type: ignore
        return BoundGreaterThanOrEqual


class GreaterThan(LiteralPredicate):
    type: TypingLiteral["gt"] = Field(default="gt", alias="type")

    def __invert__(self) -> LessThanOrEqual:
        """Transform the Expression into its negated version."""
        return LessThanOrEqual(self.term, self.literal)

    @property
    def as_bound(self) -> type[BoundGreaterThan]:  # type: ignore
        return BoundGreaterThan


class LessThanOrEqual(LiteralPredicate):
    type: TypingLiteral["lt-eq"] = Field(default="lt-eq", alias="type")

    def __invert__(self) -> GreaterThan:
        """Transform the Expression into its negated version."""
        return GreaterThan(self.term, self.literal)

    @property
    def as_bound(self) -> type[BoundLessThanOrEqual]:  # type: ignore
        return BoundLessThanOrEqual


class StartsWith(LiteralPredicate):
    type: TypingLiteral["starts-with"] = Field(default="starts-with", alias="type")

    def __invert__(self) -> NotStartsWith:
        """Transform the Expression into its negated version."""
        return NotStartsWith(self.term, self.literal)

    @property
    def as_bound(self) -> type[BoundStartsWith]:  # type: ignore
        return BoundStartsWith


class NotStartsWith(LiteralPredicate):
    type: TypingLiteral["not-starts-with"] = Field(default="not-starts-with", alias="type")

    def __invert__(self) -> StartsWith:
        """Transform the Expression into its negated version."""
        return StartsWith(self.term, self.literal)

    @property
    def as_bound(self) -> type[BoundNotStartsWith]:  # type: ignore
        return BoundNotStartsWith
