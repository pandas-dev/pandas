# This file is dual licensed under the terms of the Apache License, Version
# 2.0, and the BSD License. See the LICENSE file in the root of this repository
# for complete details.

from __future__ import annotations

import builtins
import dataclasses
import enum
import sys
import types
import typing

if sys.version_info < (3, 11):
    import typing_extensions

    LiteralString = typing_extensions.LiteralString
else:
    LiteralString = typing.LiteralString

from cryptography.hazmat.bindings._rust import declarative_asn1
from cryptography.hazmat.bindings._rust import x509 as rust_x509

if sys.version_info < (3, 10):
    NoneType = type(None)
else:
    NoneType = types.NoneType  # type: ignore[nonetype-type]

T = typing.TypeVar("T", covariant=True)
U = typing.TypeVar("U")
Tag = typing.TypeVar("Tag", bound=LiteralString)


@dataclasses.dataclass(frozen=True)
class Variant(typing.Generic[U, Tag]):
    """
    A tagged variant for CHOICE fields with the same underlying type.

    Use this when you have multiple CHOICE alternatives with the same type
    and need to distinguish between them:

        foo: (
            Annotated[Variant[int, typing.Literal["IntA"]], Implicit(0)]
            | Annotated[Variant[int, typing.Literal["IntB"]], Implicit(1)]
        )

    Usage:
        example = Example(foo=Variant(5, "IntA"))
        decoded.foo.value  # The int value
        decoded.foo.tag    # "IntA" or "IntB"
    """

    value: U
    tag: str


decode_der = declarative_asn1.decode_der
encode_der = declarative_asn1.encode_der


_X509_TYPES = (
    rust_x509.Certificate,
    rust_x509.CertificateSigningRequest,
    rust_x509.CertificateRevocationList,
)


def _check_x509_field_annotations(
    field_type: typing.Any,
    annotation: declarative_asn1.Annotation,
    field_name: str,
) -> None:
    if field_type in _X509_TYPES and isinstance(annotation.encoding, Implicit):
        raise TypeError(
            f"field '{field_name}' has an IMPLICIT annotation, but "
            "IMPLICIT annotations are not supported for X.509 types."
        )


def _is_union(field_type: type) -> bool:
    # NOTE: types.UnionType for `T | U`, typing.Union for `Union[T, U]`.
    # TODO: Drop the `hasattr()` once the minimum supported Python version
    # is >= 3.10.
    union_types = (
        (types.UnionType, typing.Union)
        if hasattr(types, "UnionType")
        else (typing.Union,)
    )
    return typing.get_origin(field_type) in union_types


def _resolve_type_aliases(field_type: typing.Any) -> typing.Any:
    # Recursively resolve PEP 695 (`type X = ...`) type aliases (Python
    # 3.12+) to their underlying value, so that the rest of the
    # normalization logic never encounters an alias. Aliases can refer
    # to other aliases and can appear at any level of nesting (e.g.
    # inside `Annotated[...]`, unions, or `list[...]`).
    if sys.version_info < (3, 12):
        return field_type

    while isinstance(field_type, typing.TypeAliasType):
        field_type = field_type.__value__

    args = typing.get_args(field_type)
    resolved_args = tuple(_resolve_type_aliases(arg) for arg in args)
    if resolved_args == args:
        # No aliases anywhere inside: return the type unchanged.
        return field_type

    if _is_union(field_type):
        # `X | Y` unions can't be rebuilt through their origin like
        # other generics below: `typing.get_origin` returns
        # `types.UnionType` for them, which is only subscriptable on
        # Python 3.14+. Rebuilding through `typing.Union` also
        # flattens any nested union introduced by an alias of a union
        # (e.g. `Time | int` where `type Time = UTCTime |
        # GeneralizedTime`), just like `typing.Union` would have done
        # if the alias had been written inline.
        return typing.Union[resolved_args]

    # An alias appeared inside a generic (e.g. `Annotated[Time, ...]`,
    # `list[MyInt]`, or `SetOf[MyInt]`): re-parameterize the generic
    # with the resolved arguments. Subscripting with a tuple is
    # equivalent to subscripting with multiple arguments.
    return typing.get_origin(field_type)[resolved_args]


def _extract_annotation(
    metadata: tuple, field_name: str
) -> declarative_asn1.Annotation:
    default = None
    encoding = None
    size = None
    for raw_annotation in metadata:
        if isinstance(raw_annotation, Default):
            if default is not None:
                raise TypeError(
                    f"multiple DEFAULT annotations found in field "
                    f"'{field_name}'"
                )
            default = raw_annotation.value
        elif isinstance(raw_annotation, declarative_asn1.Encoding):
            if encoding is not None:
                raise TypeError(
                    f"multiple IMPLICIT/EXPLICIT annotations found in field "
                    f"'{field_name}'"
                )
            encoding = raw_annotation
        elif isinstance(raw_annotation, declarative_asn1.Size):
            if size is not None:
                raise TypeError(
                    f"multiple SIZE annotations found in field '{field_name}'"
                )
            size = raw_annotation
        else:
            raise TypeError(f"unsupported annotation: {raw_annotation}")

    return declarative_asn1.Annotation(
        default=default, encoding=encoding, size=size
    )


def _normalize_field_type(
    field_type: typing.Any, field_name: str
) -> declarative_asn1.AnnotatedType:
    field_type = _resolve_type_aliases(field_type)

    # Strip the `Annotated[...]` off, and populate the annotation
    # from it if it exists.
    if typing.get_origin(field_type) is typing.Annotated:
        annotation = _extract_annotation(field_type.__metadata__, field_name)
        field_type, *_ = typing.get_args(field_type)
    else:
        annotation = declarative_asn1.Annotation()

    if annotation.size is not None and (
        typing.get_origin(field_type) not in (builtins.list, SetOf)
        and field_type
        not in (
            builtins.bytes,
            builtins.str,
            BitString,
            IA5String,
            PrintableString,
        )
    ):
        raise TypeError(
            f"field '{field_name}' has a SIZE annotation, but SIZE "
            "annotations are only supported for fields of types: "
            "[SEQUENCE OF, SET OF, BIT STRING, OCTET STRING, UTF8String, "
            "PrintableString, IA5String]"
        )

    if field_type is TLV:
        if isinstance(annotation.encoding, Implicit):
            raise TypeError(
                f"field '{field_name}' has an IMPLICIT annotation, but "
                "IMPLICIT annotations are not supported for TLV types."
            )
        elif annotation.default is not None:
            raise TypeError(
                f"field '{field_name}' has a DEFAULT annotation, but "
                "DEFAULT annotations are not supported for TLV types."
            )

    _check_x509_field_annotations(field_type, annotation, field_name)

    if hasattr(field_type, "__asn1_root__"):
        root_type = field_type.__asn1_root__
        if not isinstance(
            root_type,
            (
                declarative_asn1.Type.Sequence,
                declarative_asn1.Type.Set,
                declarative_asn1.Type.ValueSet,
            ),
        ):
            raise TypeError(f"unsupported root type: {root_type}")
        return declarative_asn1.AnnotatedType(
            typing.cast(declarative_asn1.Type, root_type), annotation
        )
    elif _is_union(field_type):
        union_args = typing.get_args(field_type)
        if len(union_args) == 2 and NoneType in union_args:
            # A Union between a type and None is an OPTIONAL
            optional_type = (
                union_args[0] if union_args[1] is type(None) else union_args[1]
            )
            if optional_type is TLV:
                raise TypeError(
                    "optional TLV types (`TLV | None`) are not "
                    "currently supported"
                )
            # For optional types, the annotation is associated with the
            # union, so we check it against the inner type here.
            _check_x509_field_annotations(
                optional_type, annotation, field_name
            )
            annotated_type = _normalize_field_type(optional_type, field_name)

            if not annotated_type.annotation.is_empty():
                raise TypeError(
                    "optional (`X | None`) types cannot have `X` "
                    "annotated: annotations must apply to the union "
                    "(i.e: `Annotated[X | None, annotation]`)"
                )

            if annotation.default is not None:
                raise TypeError(
                    "optional (`X | None`) types should not have a DEFAULT "
                    "annotation"
                )

            rust_field_type = declarative_asn1.Type.Option(annotated_type)

        else:
            # Otherwise, the Union is a CHOICE
            if isinstance(annotation.encoding, Implicit):
                # CHOICEs cannot be IMPLICIT. See X.680 section 31.2.9.
                raise TypeError(
                    "CHOICE (`X | Y | ...`) types should not have an IMPLICIT "
                    "annotation"
                )
            variants = [
                _type_to_variant(arg, field_name)
                for arg in union_args
                if arg is not type(None)
            ]

            # Union types should either be all Variants
            # (`Variant[..] | Variant[..] | etc`) or all non Variants
            are_union_types_tagged = variants[0].tag_name is not None
            if any(
                (v.tag_name is not None) != are_union_types_tagged
                for v in variants
            ):
                raise TypeError(
                    "When using `asn1.Variant` in a union, all the other "
                    "types in the union must also be `asn1.Variant`"
                )

            if are_union_types_tagged:
                tags = {v.tag_name for v in variants}
                if len(variants) != len(tags):
                    raise TypeError(
                        "When using `asn1.Variant` in a union, the tags used "
                        "must be unique"
                    )

            rust_choice_type = declarative_asn1.Type.Choice(variants)
            # If None is part of the union types, this is an OPTIONAL CHOICE
            rust_field_type = (
                declarative_asn1.Type.Option(
                    declarative_asn1.AnnotatedType(
                        rust_choice_type, declarative_asn1.Annotation()
                    )
                )
                if NoneType in union_args
                else rust_choice_type
            )

    elif typing.get_origin(field_type) is builtins.list:
        inner_type = _normalize_field_type(
            typing.get_args(field_type)[0], field_name
        )
        rust_field_type = declarative_asn1.Type.SequenceOf(inner_type)
    elif typing.get_origin(field_type) is SetOf:
        inner_type = _normalize_field_type(
            typing.get_args(field_type)[0], field_name
        )
        rust_field_type = declarative_asn1.Type.SetOf(inner_type)
    else:
        rust_field_type = declarative_asn1.non_root_python_to_rust(field_type)

    return declarative_asn1.AnnotatedType(rust_field_type, annotation)


# Convert a type to a Variant. Used with types inside Union
# annotations (T1, T2, etc in `Union[T1, T2, ...]`).
def _type_to_variant(
    t: typing.Any, field_name: str
) -> declarative_asn1.Variant:
    is_annotated = typing.get_origin(t) is typing.Annotated
    inner_type = typing.get_args(t)[0] if is_annotated else t

    # Check if this is a Variant[T, Tag] type
    if typing.get_origin(inner_type) is Variant:
        value_type, tag_literal = typing.get_args(inner_type)
        if typing.get_origin(tag_literal) is not typing.Literal:
            raise TypeError(
                "When using `asn1.Variant` in a type annotation, the second "
                "type parameter must be a `typing.Literal` type. E.g: "
                '`Variant[int, typing.Literal["MyInt"]]`.'
            )
        tag_name = typing.get_args(tag_literal)[0]

        if hasattr(value_type, "__asn1_root__"):
            rust_type = value_type.__asn1_root__
        else:
            rust_type = declarative_asn1.non_root_python_to_rust(value_type)

        if is_annotated:
            ann_type = declarative_asn1.AnnotatedType(
                rust_type,
                _extract_annotation(t.__metadata__, field_name),
            )
        else:
            ann_type = declarative_asn1.AnnotatedType(
                rust_type,
                declarative_asn1.Annotation(),
            )

        return declarative_asn1.Variant(Variant, ann_type, tag_name)
    else:
        # Plain type (not a tagged Variant)
        return declarative_asn1.Variant(
            inner_type,
            _normalize_field_type(t, field_name),
            None,
        )


def _annotate_fields(
    raw_fields: dict[str, type],
) -> dict[str, declarative_asn1.AnnotatedType]:
    fields = {}
    for field_name, field_type in raw_fields.items():
        # Recursively normalize the field type into something that the
        # Rust code can understand.
        annotated_field_type = _normalize_field_type(field_type, field_name)
        fields[field_name] = annotated_field_type

    return fields


def _register_asn1_sequence(cls: type[U]) -> None:
    raw_fields = typing.get_type_hints(cls, include_extras=True)
    root = declarative_asn1.Type.Sequence(cls, _annotate_fields(raw_fields))

    setattr(cls, "__asn1_root__", root)


def _register_asn1_set(cls: type[U]) -> None:
    raw_fields = typing.get_type_hints(cls, include_extras=True)
    root = declarative_asn1.Type.Set(cls, _annotate_fields(raw_fields))

    setattr(cls, "__asn1_root__", root)


# Due to https://github.com/python/mypy/issues/19731, we can't define an alias
# for `dataclass_transform` that conditionally points to `typing` or
# `typing_extensions` depending on the Python version. We work around it by
# making the whole decorated class conditional on the Python version.
if sys.version_info < (3, 11):

    @typing_extensions.dataclass_transform(kw_only_default=True)
    def sequence(cls: type[U]) -> type[U]:
        # We use `dataclasses.dataclass` to add an __init__ method
        # to the class with keyword-only parameters.
        if sys.version_info >= (3, 10):
            dataclass_cls = dataclasses.dataclass(
                repr=False,
                eq=False,
                # `match_args` was added in Python 3.10 and defaults
                # to True
                match_args=False,
                # `kw_only` was added in Python 3.10 and defaults to
                # False
                kw_only=True,
            )(cls)
        else:
            dataclass_cls = dataclasses.dataclass(
                repr=False,
                eq=False,
            )(cls)
        _register_asn1_sequence(dataclass_cls)
        return dataclass_cls

    @typing_extensions.dataclass_transform(kw_only_default=True)
    def set(cls: type[U]) -> type[U]:
        # We use `dataclasses.dataclass` to add an __init__ method
        # to the class with keyword-only parameters.
        if sys.version_info >= (3, 10):
            dataclass_cls = dataclasses.dataclass(
                repr=False,
                eq=False,
                # `match_args` was added in Python 3.10 and defaults
                # to True
                match_args=False,
                # `kw_only` was added in Python 3.10 and defaults to
                # False
                kw_only=True,
            )(cls)
        else:
            dataclass_cls = dataclasses.dataclass(
                repr=False,
                eq=False,
            )(cls)
        _register_asn1_set(dataclass_cls)
        return dataclass_cls

else:

    @typing.dataclass_transform(kw_only_default=True)
    def sequence(cls: type[U]) -> type[U]:
        # Only add an __init__ method, with keyword-only
        # parameters.
        dataclass_cls = dataclasses.dataclass(
            repr=False,
            eq=False,
            match_args=False,
            kw_only=True,
        )(cls)
        _register_asn1_sequence(dataclass_cls)
        return dataclass_cls

    @typing.dataclass_transform(kw_only_default=True)
    def set(cls: type[U]) -> type[U]:
        # Only add an __init__ method, with keyword-only
        # parameters.
        dataclass_cls = dataclasses.dataclass(
            repr=False,
            eq=False,
            match_args=False,
            kw_only=True,
        )(cls)
        _register_asn1_set(dataclass_cls)
        return dataclass_cls


def value_set(
    value_type: type,
) -> typing.Callable[[type[U]], type[U]]:
    """
    A class decorator that registers an `enum.Enum` subclass as an
    ASN.1 value set of the given underlying type. All the member
    values must be instances of `value_type`. Members are encoded as
    their value; decoding fails if the decoded value does not match
    any member.
    """
    rust_type = declarative_asn1.non_root_python_to_rust(value_type)

    def decorator(cls: type[U]) -> type[U]:
        if not issubclass(cls, enum.Enum):
            raise TypeError(
                "value sets can only be defined from enum.Enum subclasses"
            )
        members = list(cls)
        if not members:
            raise TypeError(
                f"value set '{cls.__name__}' must have at least one member"
            )
        for member in members:
            if not isinstance(member.value, value_type):
                raise TypeError(
                    f"member '{member.name}' of value set '{cls.__name__}' "
                    f"must have a value of type "
                    f"'{value_type.__name__}', got: "
                    f"'{type(member.value).__name__}'"
                )
        inner = declarative_asn1.AnnotatedType(
            rust_type, declarative_asn1.Annotation()
        )
        # Map from member value to member, used for O(1) lookups when
        # decoding. This requires the member values to be hashable.
        value_map = {member.value: member for member in members}
        root = declarative_asn1.Type.ValueSet(cls, inner, value_map)

        setattr(cls, "__asn1_root__", root)
        return cls

    return decorator


# TODO: replace with `Default[U]` once the min Python version is >= 3.12
@dataclasses.dataclass(frozen=True)
class Default(typing.Generic[U]):
    value: U


SetOf = declarative_asn1.SetOf

Explicit = declarative_asn1.Encoding.Explicit
Implicit = declarative_asn1.Encoding.Implicit
Size = declarative_asn1.Size

PrintableString = declarative_asn1.PrintableString
IA5String = declarative_asn1.IA5String
UTCTime = declarative_asn1.UTCTime
GeneralizedTime = declarative_asn1.GeneralizedTime
BitString = declarative_asn1.BitString
TLV = declarative_asn1.Tlv
Null = declarative_asn1.Null
