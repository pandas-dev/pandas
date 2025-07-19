# Copyright (c) 2005-2011, NumPy Developers.
# All rights reserved.

# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:

#     * Redistributions of source code must retain the above copyright
#        notice, this list of conditions and the following disclaimer.

#     * Redistributions in binary form must reproduce the above
#        copyright notice, this list of conditions and the following
#        disclaimer in the documentation and/or other materials provided
#        with the distribution.

#     * Neither the name of the NumPy Developers nor the names of any
#        contributors may be used to endorse or promote products derived
#        from this software without specific prior written permission.

# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
from __future__ import annotations

from typing import Any

try:
    # requires numpy>=2.0
    from numpy import isdtype  # type: ignore[attr-defined,unused-ignore]

    HAS_STRING_DTYPE = True
except ImportError:
    import numpy as np
    from numpy.typing import DTypeLike

    kind_mapping = {
        "bool": np.bool_,
        "signed integer": np.signedinteger,
        "unsigned integer": np.unsignedinteger,
        "integral": np.integer,
        "real floating": np.floating,
        "complex floating": np.complexfloating,
        "numeric": np.number,
    }

    def isdtype(
        dtype: np.dtype[Any] | type[Any], kind: DTypeLike | tuple[DTypeLike, ...]
    ) -> bool:
        kinds = kind if isinstance(kind, tuple) else (kind,)
        str_kinds = {k for k in kinds if isinstance(k, str)}
        type_kinds = {k.type for k in kinds if isinstance(k, np.dtype)}

        if unknown_kind_types := set(kinds) - str_kinds - type_kinds:
            raise TypeError(
                f"kind must be str, np.dtype or a tuple of these, got {unknown_kind_types}"
            )
        if unknown_kinds := {k for k in str_kinds if k not in kind_mapping}:
            raise ValueError(
                f"unknown kind: {unknown_kinds}, must be a np.dtype or one of {list(kind_mapping)}"
            )

        # verified the dtypes already, no need to check again
        translated_kinds = {kind_mapping[k] for k in str_kinds} | type_kinds
        if isinstance(dtype, np.generic):
            return isinstance(dtype, translated_kinds)
        else:
            return any(np.issubdtype(dtype, k) for k in translated_kinds)

    HAS_STRING_DTYPE = False
