"""
Dispatch in dask.array.

Also see backends.py
"""

from __future__ import annotations

from dask.utils import Dispatch

concatenate_lookup = Dispatch("concatenate")
take_lookup = Dispatch("take")
tensordot_lookup = Dispatch("tensordot")
einsum_lookup = Dispatch("einsum")
empty_lookup = Dispatch("empty")
divide_lookup = Dispatch("divide")
percentile_lookup = Dispatch("percentile")
numel_lookup = Dispatch("numel")
nannumel_lookup = Dispatch("nannumel")
to_numpy_dispatch = Dispatch("to_numpy_dispatch")
to_cupy_dispatch = Dispatch("to_cupy_dispatch")
