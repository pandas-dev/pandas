#######################################################################
# Copyright (c) 2019-present, Blosc Development Team <blosc@blosc.org>
# All rights reserved.
#
# This source code is licensed under a BSD-style license (found in the
# LICENSE file in the root directory of this source tree)
#######################################################################

"""Row-level validation via an internally-generated Pydantic model.

All Pydantic-specific logic is isolated here.  CTable and the rest of the
schema layer never import from Pydantic directly.
"""

from __future__ import annotations

import math
from dataclasses import MISSING
from typing import Any

import numpy as np
from pydantic import BaseModel, Field, ValidationError, create_model

from blosc2.list_array import _coerce_struct_item, coerce_list_cell
from blosc2.schema import ListSpec, NDArraySpec, StructSpec
from blosc2.schema_compiler import CompiledSchema  # noqa: TC001


def build_validator_model(schema: CompiledSchema) -> type[BaseModel]:
    """Return (and cache) a Pydantic model class for *schema*.

    Built once per schema; subsequent calls return the cached class.
    The model enforces all constraints declared in each column's
    :class:`~blosc2.schema.SchemaSpec` (``ge``, ``le``, ``gt``, ``lt``,
    ``max_length``, ``min_length``, ``pattern``).

    Nullable columns (those with a ``null_value``) are typed as
    ``Optional[T]`` with ``default=None`` so that null sentinels can be
    passed as ``None`` and bypass constraint validation entirely — no
    placeholder guessing required.
    """
    if schema.validator_model is not None:
        return schema.validator_model

    field_definitions: dict[str, Any] = {}
    for col in schema.columns:
        pydantic_kwargs = col.spec.to_pydantic_kwargs()
        is_nullable = getattr(col.spec, "null_value", None) is not None or bool(
            getattr(col.spec, "nullable", False)
        )
        if isinstance(col.spec, ListSpec):
            item_type = col.spec.item_spec.python_type
            py_type = list[item_type] | None if is_nullable else list[item_type]
        else:
            py_type = col.py_type | None if is_nullable else col.py_type

        if col.default is MISSING:
            default = None if is_nullable else MISSING
            if default is MISSING:
                field_definitions[col.name] = (py_type, Field(**pydantic_kwargs))
            else:
                field_definitions[col.name] = (py_type, Field(default=default, **pydantic_kwargs))
        else:
            field_definitions[col.name] = (py_type, Field(default=col.default, **pydantic_kwargs))

    cls_name = schema.row_cls.__name__ if schema.row_cls is not None else "Unknown"
    model_cls = create_model(f"_Validator_{cls_name}", **field_definitions)
    schema.validator_model = model_cls
    return model_cls


def _is_null_value(val, null_value) -> bool:
    """Return True if *val* equals the null sentinel, handling NaN correctly."""
    import math

    if null_value is None:
        return False
    try:
        if isinstance(null_value, (float, np.floating)) and math.isnan(null_value):
            return isinstance(val, (float, np.floating)) and math.isnan(val)
    except TypeError:
        pass
    return val == null_value


def _mask_nulls(schema: CompiledSchema, row: dict[str, Any]) -> tuple[dict[str, Any], dict[str, Any]]:
    """Replace null sentinel values with ``None`` so Pydantic skips constraint checks.

    Nullable columns are declared as ``Optional[T]`` in the validator model,
    so passing ``None`` is always valid regardless of ``ge``/``le``/``pattern``
    constraints.  The original sentinel is stashed in *nulled* and restored
    after validation.

    Returns (masked_row, nulled) where nulled maps column name → sentinel value.
    """
    masked = dict(row)
    nulled: dict[str, Any] = {}
    for col in schema.columns:
        nv = getattr(col.spec, "null_value", None)
        if nv is None:
            continue
        val = row.get(col.name)
        if isinstance(col.spec, NDArraySpec):
            try:
                arr = np.asarray(val, dtype=col.spec.dtype)
                is_null = arr.shape == col.spec.item_shape and bool(
                    np.isnan(arr).all()
                    if isinstance(nv, (float, np.floating)) and math.isnan(nv)
                    else (arr == nv).all()
                )
            except Exception:
                is_null = val is None
        else:
            is_null = _is_null_value(val, nv)
        if is_null:
            nulled[col.name] = val
            masked[col.name] = None
    return masked, nulled


def validate_row(schema: CompiledSchema, row: dict[str, Any]) -> dict[str, Any]:
    """Validate a single row dict and return the coerced values.

    Parameters
    ----------
    schema:
        Compiled schema for the table.
    row:
        ``{column_name: value}`` mapping for one row.

    Returns
    -------
    dict
        Validated (and Pydantic-coerced) values ready for storage.

    Raises
    ------
    ValueError
        If any constraint is violated.  The message includes the column
        name and the violated constraint.
    """
    model_cls = build_validator_model(schema)
    normalized = dict(row)
    for col in schema.columns:
        if isinstance(col.spec, ListSpec) and col.name in normalized:
            normalized[col.name] = coerce_list_cell(col.spec, normalized[col.name])
        elif (
            isinstance(col.spec, StructSpec) and col.name in normalized and normalized[col.name] is not None
        ):
            normalized[col.name] = _coerce_struct_item(col.spec, normalized[col.name])
    masked_row, nulled = _mask_nulls(schema, normalized)
    try:
        instance = model_cls(**masked_row)
    except ValidationError as exc:
        # Re-raise as a plain ValueError so callers don't need to import Pydantic.
        raise ValueError(str(exc)) from exc
    result = instance.model_dump()
    result.update(nulled)
    return result


def validate_rows_rowwise(schema: CompiledSchema, rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Validate a list of row dicts.  Returns a list of validated dicts.

    Parameters
    ----------
    schema:
        Compiled schema for the table.
    rows:
        List of ``{column_name: value}`` mappings.

    Raises
    ------
    ValueError
        On the first row that violates a constraint, with the row index
        and the Pydantic error details.
    """
    model_cls = build_validator_model(schema)
    result = []
    for i, row in enumerate(rows):
        normalized = dict(row)
        for col in schema.columns:
            if isinstance(col.spec, ListSpec) and col.name in normalized:
                normalized[col.name] = coerce_list_cell(col.spec, normalized[col.name])
            elif (
                isinstance(col.spec, StructSpec)
                and col.name in normalized
                and normalized[col.name] is not None
            ):
                normalized[col.name] = _coerce_struct_item(col.spec, normalized[col.name])
        masked_row, nulled = _mask_nulls(schema, normalized)
        try:
            instance = model_cls(**masked_row)
        except ValidationError as exc:
            raise ValueError(f"Row {i}: {exc}") from exc
        validated = instance.model_dump()
        validated.update(nulled)
        result.append(validated)
    return result
