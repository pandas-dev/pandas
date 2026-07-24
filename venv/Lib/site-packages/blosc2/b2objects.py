#######################################################################
# Copyright (c) 2019-present, Blosc Development Team <blosc@blosc.org>
# All rights reserved.
#
# SPDX-License-Identifier: BSD-3-Clause
#######################################################################

from __future__ import annotations

import inspect
import pathlib
import textwrap
from dataclasses import asdict
from typing import Any

import numpy as np

import blosc2
from blosc2.dsl_kernel import DSLKernel, kernel_from_source

_B2OBJECT_META_KEY = "b2o"
_B2OBJECT_VERSION = 1
_B2OBJECT_DSL_VERSION = 1
_B2OBJECT_USER_VLMETA_KEY = "_b2o_user_vlmeta"


def make_b2object_carrier(
    kind: str,
    shape,
    dtype,
    *,
    chunks=None,
    blocks=None,
    **kwargs,
):
    meta = dict(kwargs.pop("meta", {}))
    meta[_B2OBJECT_META_KEY] = {"kind": kind, "version": _B2OBJECT_VERSION}
    kwargs["meta"] = meta
    return blosc2.empty(shape=shape, dtype=dtype, chunks=chunks, blocks=blocks, **kwargs)


def write_b2object_payload(array, payload: dict[str, Any]) -> None:
    array.schunk.vlmeta[_B2OBJECT_META_KEY] = payload


def write_b2object_user_vlmeta(array, user_vlmeta: dict[str, Any]) -> None:
    array.schunk.vlmeta[_B2OBJECT_USER_VLMETA_KEY] = user_vlmeta


def read_b2object_user_vlmeta(obj) -> dict[str, Any]:
    schunk = getattr(obj, "schunk", obj)
    if _B2OBJECT_USER_VLMETA_KEY not in schunk.vlmeta:
        return {}
    return schunk.vlmeta[_B2OBJECT_USER_VLMETA_KEY]


def encode_operand_reference(obj):
    return blosc2.Ref.from_object(obj).to_dict()


def decode_operand_reference(payload, *, base_path=None):
    if (
        payload.get("kind") in {"urlpath", "dictstore_key"}
        and base_path is not None
        and not pathlib.Path(payload["urlpath"]).is_absolute()
    ):
        payload = dict(payload)
        payload["urlpath"] = (base_path / payload["urlpath"]).as_posix()
    ref = blosc2.Ref.from_dict(payload)
    return ref.open()


def encode_b2object_payload(obj) -> dict[str, Any] | None:
    if isinstance(obj, blosc2.C2Array):
        return blosc2.Ref.c2array_ref(obj.path, obj.urlbase).to_dict()
    if isinstance(obj, blosc2.LazyExpr):
        expression = obj.expression_tosave if hasattr(obj, "expression_tosave") else obj.expression
        operands = obj.operands_tosave if hasattr(obj, "operands_tosave") else obj.operands
        return {
            "kind": "lazyexpr",
            "version": _B2OBJECT_VERSION,
            "expression": expression,
            "operands": {key: encode_operand_reference(value) for key, value in operands.items()},
        }
    if isinstance(obj, blosc2.LazyUDF):
        if not isinstance(obj.func, DSLKernel):
            raise TypeError("Structured Blosc2 msgpack payload only supports LazyUDF backed by DSLKernel")
        udf_func = obj.func.func
        udf_name = getattr(udf_func, "__name__", obj.func.__name__)
        try:
            udf_source = textwrap.dedent(inspect.getsource(udf_func)).lstrip()
        except Exception:
            udf_source = obj.func.dsl_source
        if udf_source is None:
            raise ValueError("Structured LazyUDF msgpack payload requires recoverable DSL kernel source")
        kwargs = {}
        for key, value in obj.kwargs.items():
            if key in {"dtype", "shape"}:
                continue
            if isinstance(value, blosc2.CParams | blosc2.DParams):
                kwargs[key] = asdict(value)
            else:
                kwargs[key] = value
        return {
            "kind": "lazyudf",
            "version": _B2OBJECT_VERSION,
            "function_kind": "dsl",
            "dsl_version": _B2OBJECT_DSL_VERSION,
            "name": udf_name,
            "udf_source": udf_source,
            "dtype": np.dtype(obj.dtype).str,
            "shape": list(obj.shape),
            "operands": {f"o{i}": encode_operand_reference(value) for i, value in enumerate(obj.inputs)},
            "kwargs": kwargs,
        }
    return None


def decode_b2object_payload(payload: dict[str, Any], *, carrier_path=None):
    kind = payload.get("kind")
    version = payload.get("version")
    if version != _B2OBJECT_VERSION:
        raise ValueError(f"Unsupported persisted Blosc2 object version: {version!r}")
    if kind == "c2array":
        ref = blosc2.Ref.from_dict(payload)
        return ref.open()
    if kind == "lazyexpr":
        return decode_structured_lazyexpr(payload, carrier_path=carrier_path)
    if kind == "lazyudf":
        return decode_structured_lazyudf(payload, carrier_path=carrier_path)
    raise ValueError(f"Unsupported persisted Blosc2 object kind: {kind!r}")


def decode_structured_lazyexpr(payload, *, carrier_path=None):
    expression = payload.get("expression")
    if not isinstance(expression, str):
        raise TypeError("Structured LazyExpr payload requires a string 'expression'")
    operands_payload = payload.get("operands")
    if not isinstance(operands_payload, dict):
        raise TypeError("Structured LazyExpr payload requires a mapping 'operands'")
    operands, missing_ops = decode_operand_mapping(operands_payload, base_path=carrier_path)
    if missing_ops:
        exc = blosc2.exceptions.MissingOperands(expression, missing_ops)
        exc.expr = expression
        exc.missing_ops = missing_ops
        raise exc
    return blosc2.lazyexpr(expression, operands=operands)


def decode_operand_mapping(operands_payload, *, base_path=None):
    operands = {}
    missing_ops = {}
    for key, value in operands_payload.items():
        try:
            operands[key] = decode_operand_reference(value, base_path=base_path)
        except FileNotFoundError:
            ref = blosc2.Ref.from_dict(value)
            if ref.kind in {"urlpath", "dictstore_key"}:
                missing_ops[key] = pathlib.Path(ref.urlpath)
            else:
                raise
    return operands, missing_ops


def decode_structured_lazyudf(payload, *, carrier_path=None):
    function_kind = payload.get("function_kind")
    if function_kind != "dsl":
        raise ValueError(f"Unsupported structured LazyUDF function kind: {function_kind!r}")
    dsl_version = payload.get("dsl_version")
    if dsl_version != _B2OBJECT_DSL_VERSION:
        raise ValueError(f"Unsupported structured LazyUDF DSL version: {dsl_version!r}")
    udf_source = payload.get("udf_source")
    if not isinstance(udf_source, str):
        raise TypeError("Structured LazyUDF payload requires a string 'udf_source'")
    name = payload.get("name")
    if not isinstance(name, str):
        raise TypeError("Structured LazyUDF payload requires a string 'name'")
    dtype = payload.get("dtype")
    if not isinstance(dtype, str):
        raise TypeError("Structured LazyUDF payload requires a string 'dtype'")
    shape_payload = payload.get("shape")
    if not isinstance(shape_payload, list):
        raise TypeError("Structured LazyUDF payload requires a list 'shape'")
    operands_payload = payload.get("operands")
    if not isinstance(operands_payload, dict):
        raise TypeError("Structured LazyUDF payload requires a mapping 'operands'")
    kwargs = payload.get("kwargs", {})
    if not isinstance(kwargs, dict):
        raise TypeError("Structured LazyUDF payload requires a mapping 'kwargs'")

    func = kernel_from_source(udf_source, name)
    ordered_operands_payload = {f"o{n}": operands_payload[f"o{n}"] for n in range(len(operands_payload))}
    operands, missing_ops = decode_operand_mapping(ordered_operands_payload, base_path=carrier_path)
    if missing_ops:
        exc = blosc2.exceptions.MissingOperands(name, missing_ops)
        exc.expr = name
        exc.missing_ops = missing_ops
        raise exc
    return blosc2.lazyudf(
        func, tuple(operands.values()), dtype=np.dtype(dtype), shape=tuple(shape_payload), **kwargs
    )


def read_b2object_marker(obj) -> dict[str, Any] | None:
    schunk = getattr(obj, "schunk", obj)
    if _B2OBJECT_META_KEY not in schunk.meta:
        return None
    return schunk.meta[_B2OBJECT_META_KEY]


def read_b2object_payload(obj) -> dict[str, Any]:
    schunk = getattr(obj, "schunk", obj)
    return schunk.vlmeta[_B2OBJECT_META_KEY]


def open_b2object(obj):
    marker = read_b2object_marker(obj)
    if marker is None:
        return None

    payload = read_b2object_payload(obj)
    if marker.get("version") != _B2OBJECT_VERSION:
        raise ValueError(f"Unsupported persisted Blosc2 object version: {marker.get('version')!r}")
    if marker.get("kind") != payload.get("kind"):
        raise ValueError("Persisted Blosc2 object marker/payload kind mismatch")
    carrier_path = None
    schunk = getattr(obj, "schunk", obj)
    if getattr(schunk, "urlpath", None) is not None:
        carrier_path = pathlib.Path(schunk.urlpath).parent
    opened = decode_b2object_payload(payload, carrier_path=carrier_path)
    if isinstance(opened, blosc2.LazyExpr | blosc2.LazyUDF):
        opened.array = obj
        opened.schunk = schunk
        opened._set_user_vlmeta(read_b2object_user_vlmeta(obj), sync=False)
    return opened
