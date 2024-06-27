"""Utilities to manipulate JSON objects."""

# NOTE: this is a copy of ipykernel/jsonutils.py (+blackified)

# Copyright (c) IPython Development Team.
# Distributed under the terms of the Modified BSD License.
from __future__ import annotations

import math
import numbers
import re
import types
from binascii import b2a_base64
from datetime import datetime
from typing import Any

# -----------------------------------------------------------------------------
# Globals and constants
# -----------------------------------------------------------------------------

# timestamp formats
ISO8601 = "%Y-%m-%dT%H:%M:%S.%f"
ISO8601_PAT = re.compile(
    r"^(\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2})(\.\d{1,6})?Z?([\+\-]\d{2}:?\d{2})?$"
)

# holy crap, strptime is not threadsafe.
# Calling it once at import seems to help.
datetime.strptime("1", "%d")

# -----------------------------------------------------------------------------
# Classes and functions
# -----------------------------------------------------------------------------


# constants for identifying png/jpeg data
PNG = b"\x89PNG\r\n\x1a\n"
# front of PNG base64-encoded
PNG64 = b"iVBORw0KG"
JPEG = b"\xff\xd8"
# front of JPEG base64-encoded
JPEG64 = b"/9"
# constants for identifying gif data
GIF_64 = b"R0lGODdh"
GIF89_64 = b"R0lGODlh"
# front of PDF base64-encoded
PDF64 = b"JVBER"


def encode_images(format_dict: dict[str, str]) -> dict[str, str]:
    """b64-encodes images in a displaypub format dict

    Perhaps this should be handled in json_clean itself?

    Parameters
    ----------

    format_dict : dict
        A dictionary of display data keyed by mime-type

    Returns
    -------

    format_dict : dict
        A copy of the same dictionary,
        but binary image data ('image/png', 'image/jpeg' or 'application/pdf')
        is base64-encoded.

    """
    return format_dict


def json_clean(obj: Any) -> Any:
    """Clean an object to ensure it's safe to encode in JSON.

    Atomic, immutable objects are returned unmodified.  Sets and tuples are
    converted to lists, lists are copied and dicts are also copied.

    Note: dicts whose keys could cause collisions upon encoding (such as a dict
    with both the number 1 and the string '1' as keys) will cause a ValueError
    to be raised.

    Parameters
    ----------
    obj : any python object

    Returns
    -------
    out : object

      A version of the input which will not cause an encoding error when
      encoded as JSON.  Note that this function does not *encode* its inputs,
      it simply sanitizes it so that there will be no encoding errors later.

    """
    # types that are 'atomic' and ok in json as-is.
    atomic_ok = (str, type(None))

    # containers that we need to convert into lists
    container_to_list = (tuple, set, types.GeneratorType)

    # Since bools are a subtype of Integrals, which are a subtype of Reals,
    # we have to check them in that order.

    if isinstance(obj, bool):
        return obj

    if isinstance(obj, numbers.Integral):
        # cast int to int, in case subclasses override __str__ (e.g. boost enum, #4598)
        return int(obj)

    if isinstance(obj, numbers.Real):
        # cast out-of-range floats to their reprs
        if math.isnan(obj) or math.isinf(obj):
            return repr(obj)
        return float(obj)

    if isinstance(obj, atomic_ok):
        return obj

    if isinstance(obj, bytes):
        return b2a_base64(obj).decode("ascii")

    if isinstance(obj, container_to_list) or (
        hasattr(obj, "__iter__") and hasattr(obj, "__next__")
    ):
        obj = list(obj)

    if isinstance(obj, list):
        return [json_clean(x) for x in obj]

    if isinstance(obj, dict):
        # First, validate that the dict won't lose data in conversion due to
        # key collisions after stringification.  This can happen with keys like
        # True and 'true' or 1 and '1', which collide in JSON.
        nkeys = len(obj)
        nkeys_collapsed = len(set(map(str, obj)))
        if nkeys != nkeys_collapsed:
            raise ValueError(
                "dict cannot be safely converted to JSON: "
                "key collision would lead to dropped values"
            )
        # If all OK, proceed by making the new dict that will be json-safe
        out = {}
        for k, v in iter(obj.items()):
            out[str(k)] = json_clean(v)
        return out
    if isinstance(obj, datetime):
        return obj.strftime(ISO8601)

    # we don't understand it, it's probably an unserializable object
    raise ValueError("Can't clean for JSON: %r" % obj)
