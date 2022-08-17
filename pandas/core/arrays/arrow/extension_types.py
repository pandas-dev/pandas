from __future__ import annotations

import json
import warnings

import pyarrow

from pandas._typing import IntervalInclusiveType
from pandas.util._decorators import deprecate_kwarg
from pandas.util._exceptions import find_stack_level

from pandas.core.arrays.interval import VALID_CLOSED


class ArrowPeriodType(pyarrow.ExtensionType):
    def __init__(self, freq) -> None:
        # attributes need to be set first before calling
        # super init (as that calls serialize)
        self._freq = freq
        pyarrow.ExtensionType.__init__(self, pyarrow.int64(), "pandas.period")

    @property
    def freq(self):
        return self._freq

    def __arrow_ext_serialize__(self) -> bytes:
        metadata = {"freq": self.freq}
        return json.dumps(metadata).encode()

    @classmethod
    def __arrow_ext_deserialize__(cls, storage_type, serialized) -> ArrowPeriodType:
        metadata = json.loads(serialized.decode())
        return ArrowPeriodType(metadata["freq"])

    def __eq__(self, other):
        if isinstance(other, pyarrow.BaseExtensionType):
            return type(self) == type(other) and self.freq == other.freq
        else:
            return NotImplemented

    def __hash__(self) -> int:
        return hash((str(self), self.freq))

    def to_pandas_dtype(self):
        import pandas as pd

        return pd.PeriodDtype(freq=self.freq)


# register the type with a dummy instance
_period_type = ArrowPeriodType("D")
pyarrow.register_extension_type(_period_type)


class ArrowIntervalType(pyarrow.ExtensionType):
    @deprecate_kwarg(old_arg_name="closed", new_arg_name="inclusive")
    def __init__(self, subtype, inclusive: IntervalInclusiveType) -> None:
        # attributes need to be set first before calling
        # super init (as that calls serialize)
        assert inclusive in VALID_CLOSED
        self._inclusive: IntervalInclusiveType = inclusive
        if not isinstance(subtype, pyarrow.DataType):
            subtype = pyarrow.type_for_alias(str(subtype))
        self._subtype = subtype

        storage_type = pyarrow.struct([("left", subtype), ("right", subtype)])
        pyarrow.ExtensionType.__init__(self, storage_type, "pandas.interval")

    @property
    def subtype(self):
        return self._subtype

    @property
    def inclusive(self) -> IntervalInclusiveType:
        return self._inclusive

    @property
    def closed(self) -> IntervalInclusiveType:
        warnings.warn(
            "Attribute `closed` is deprecated in favor of `inclusive`.",
            FutureWarning,
            stacklevel=find_stack_level(),
        )
        return self._inclusive

    def __arrow_ext_serialize__(self) -> bytes:
        metadata = {"subtype": str(self.subtype), "inclusive": self.inclusive}
        return json.dumps(metadata).encode()

    @classmethod
    def __arrow_ext_deserialize__(cls, storage_type, serialized) -> ArrowIntervalType:
        metadata = json.loads(serialized.decode())
        subtype = pyarrow.type_for_alias(metadata["subtype"])
        inclusive = metadata["inclusive"]
        return ArrowIntervalType(subtype, inclusive)

    def __eq__(self, other):
        if isinstance(other, pyarrow.BaseExtensionType):
            return (
                type(self) == type(other)
                and self.subtype == other.subtype
                and self.inclusive == other.inclusive
            )
        else:
            return NotImplemented

    def __hash__(self) -> int:
        return hash((str(self), str(self.subtype), self.inclusive))

    def to_pandas_dtype(self):
        import pandas as pd

        return pd.IntervalDtype(self.subtype.to_pandas_dtype(), self.inclusive)


# register the type with a dummy instance
_interval_type = ArrowIntervalType(pyarrow.int64(), "left")
pyarrow.register_extension_type(_interval_type)
