from typing import Literal

import optype.numpy as onp

from ._resampling import PermutationMethod, PermutationTestResult

def bws_test(
    x: onp.ToComplex1D,
    y: onp.ToComplex1D,
    *,
    alternative: Literal["two-sided", "less", "greater"] = "two-sided",
    method: PermutationMethod | None = None,
) -> PermutationTestResult: ...
