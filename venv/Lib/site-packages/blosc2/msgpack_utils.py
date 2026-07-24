#######################################################################
# Copyright (c) 2019-present, Blosc Development Team <blosc@blosc.org>
# All rights reserved.
#
# SPDX-License-Identifier: BSD-3-Clause
#######################################################################

from __future__ import annotations

from msgpack import ExtType, packb, unpackb

from blosc2 import blosc2_ext
from blosc2.b2objects import decode_b2object_payload, encode_b2object_payload
from blosc2.ref import Ref

# Msgpack extension type codes are application-defined.  Reserve code 42 in
# python-blosc2 for values serialized as Blosc2 CFrames via ``to_cframe()`` and
# reconstructed with ``blosc2.from_cframe()``.  Keep this stable for backward
# compatibility with persisted msgpack payloads produced by this package.
_BLOSC2_EXT_CODE = 42
# Reserve code 43 for structured Blosc2 reference objects that are not naturally
# serialized as CFrames.  The payload is a msgpack-encoded mapping with a
# stable ``kind`` and ``version`` envelope.
_BLOSC2_STRUCTURED_EXT_CODE = 43
_BLOSC2_STRUCTURED_VERSION = 1


def _encode_structured_reference(obj):
    import blosc2

    if isinstance(obj, blosc2.Ref):
        payload = {"kind": "ref", "version": _BLOSC2_STRUCTURED_VERSION, "ref": obj.to_dict()}
        return ExtType(_BLOSC2_STRUCTURED_EXT_CODE, packb(payload, use_bin_type=True))
    payload = encode_b2object_payload(obj)
    if payload is not None:
        return ExtType(_BLOSC2_STRUCTURED_EXT_CODE, packb(payload, use_bin_type=True))
    return None


def _decode_structured_reference(data):
    payload = unpackb(data)
    if not isinstance(payload, dict):
        raise TypeError("Structured Blosc2 msgpack payload must decode to a mapping")

    version = payload.get("version")
    if version != _BLOSC2_STRUCTURED_VERSION:
        raise ValueError(f"Unsupported structured Blosc2 msgpack payload version: {version!r}")

    kind = payload.get("kind")
    if kind == "ref":
        ref_payload = payload.get("ref")
        return Ref.from_dict(ref_payload)
    if kind in {"c2array", "lazyexpr", "lazyudf"}:
        return decode_b2object_payload(payload)
    raise ValueError(f"Unsupported structured Blosc2 msgpack payload kind: {kind!r}")


def _encode_msgpack_ext(obj):
    import blosc2

    if isinstance(
        obj, blosc2.NDArray | blosc2.SChunk | blosc2.ObjectArray | blosc2.BatchArray | blosc2.EmbedStore
    ):
        return ExtType(_BLOSC2_EXT_CODE, obj.to_cframe())
    structured = _encode_structured_reference(obj)
    if structured is not None:
        return structured
    return blosc2_ext.encode_tuple(obj)


def msgpack_packb(value):
    return packb(value, default=_encode_msgpack_ext, strict_types=True, use_bin_type=True)


def decode_tuple_list_hook(obj):
    if obj and isinstance(obj[0], str) and obj[0] == "__tuple__":
        return tuple(obj[1:])
    return obj


def _decode_msgpack_ext(code, data):
    import blosc2

    if code == _BLOSC2_EXT_CODE:
        return blosc2.from_cframe(data, copy=True)
    if code == _BLOSC2_STRUCTURED_EXT_CODE:
        return _decode_structured_reference(data)
    return ExtType(code, data)


def msgpack_unpackb(payload):
    return unpackb(payload, list_hook=decode_tuple_list_hook, ext_hook=_decode_msgpack_ext)
