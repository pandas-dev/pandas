"""Define signature overrides for Cython classes used in the documentation."""

from datetime import timedelta
import inspect

import pandas as pd

_OFFSET_PARAMETERS = {
    "n": inspect.Parameter(
        "n",
        inspect.Parameter.POSITIONAL_OR_KEYWORD,
        default=1,
    ),
    "normalize": inspect.Parameter(
        "normalize",
        inspect.Parameter.POSITIONAL_OR_KEYWORD,
        default=False,
    ),
    "weekmask": inspect.Parameter(
        "weekmask",
        inspect.Parameter.POSITIONAL_OR_KEYWORD,
        default="Mon Tue Wed Thu Fri",
    ),
    "holidays": inspect.Parameter(
        "holidays",
        inspect.Parameter.POSITIONAL_OR_KEYWORD,
        default=None,
    ),
    "calendar": inspect.Parameter(
        "calendar",
        inspect.Parameter.POSITIONAL_OR_KEYWORD,
        default=None,
    ),
    "offset": inspect.Parameter(
        "offset",
        inspect.Parameter.POSITIONAL_OR_KEYWORD,
        default=timedelta(0),
    ),
}

_OFFSET_SIGNATURES = {
    pd.tseries.offsets.BusinessDay: (
        "n",
        "normalize",
        "offset",
    ),
    pd.tseries.offsets.CustomBusinessDay: (
        "n",
        "normalize",
        "weekmask",
        "holidays",
        "calendar",
        "offset",
    ),
}


def apply_signature_overrides() -> None:
    for offset_class, parameter_names in _OFFSET_SIGNATURES.items():
        signature = inspect.Signature(
            parameters=[_OFFSET_PARAMETERS[name] for name in parameter_names]
        )

        assert offset_class() == offset_class(
            **{name: param.default for name, param in signature.parameters.items()}
        )

        offset_class.__signature__ = signature
