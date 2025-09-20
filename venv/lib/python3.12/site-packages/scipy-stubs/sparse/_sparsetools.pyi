# NOTE: The scipy source is generated in `scipy/sparse/_generate_sparsetools.py`

# Type codes used:
#
# - i:  integer scalar
# - I:  integer array
# - T:  data array
# - V:  std::vector<integer>*
# - W:  std::vector<data>*
# - *:  indicates that the next argument is an output argument
# - v:  void
# - l:  64-bit integer scalar
#
# See `scipy/sparse/sparsetools/sparsetools.cxx` for more details.

from collections.abc import Sequence
from typing import Never, TypeAlias

import numpy as np
import optype.numpy as onp
import optype.numpy.compat as npc

_i: TypeAlias = int | np.int32 | np.int64  # noqa: PYI042
_l: TypeAlias = int | np.int64  # noqa: PYI042
_I: TypeAlias = onp.ArrayND[np.int32 | np.int64]
_T: TypeAlias = onp.ArrayND[npc.number | np.bool_]
_V: TypeAlias = Sequence[_i]
_W: TypeAlias = Sequence[complex | npc.number | np.bool_]

###

# bsr.h
def bsr_diagonal(a0: _i, a1: _i, a2: _i, a3: _i, a4: _i, a5: _I, a6: _I, a7: _T, a8: _T, /) -> None: ...  # v iiiiiIIT*T
def bsr_tocsr(a0: _i, a1: _i, a2: _i, a3: _i, a4: _I, a5: _I, a6: _T, a7: _I, a8: _I, a9: _T, /) -> None: ...  # v iiiiIIT*I*I*T
def bsr_scale_columns(a0: _i, a1: _i, a2: _i, a3: _i, a4: _I, a5: _I, a6: _T, a7: _T, /) -> None: ...  # v iiiiII*TT
def bsr_scale_rows(a0: _i, a1: _i, a2: _i, a3: _i, a4: _I, a5: _I, a6: _T, a7: _T, /) -> None: ...  # v iiiiII*TT
def bsr_sort_indices(a0: _i, a1: _i, a2: _i, a3: _i, a4: _I, a5: _I, a6: _T, /) -> None: ...  # v iiii*I*I*T
def bsr_transpose(
    a0: _i, a1: _i, a2: _i, a3: _i, a4: _I, a5: _I, a6: _T, a7: _I, a8: _I, a9: _T, /
) -> None: ...  # v iiiiIIT*I*I*T
def bsr_matmat(
    a0: _i, a1: _i, a2: _i, a3: _i, a4: _i, a5: _i, a6: _I, a7: _I, a8: _T, a9: _I, a10: _I, a11: _T, a12: _I, a13: _I, a14: _T, /
) -> None: ...  # v iiiiiiIITIIT*I*I*T
def bsr_matvec(a0: _i, a1: _i, a2: _i, a3: _i, a4: _I, a5: _I, a6: _T, a7: _T, a8: _T, /) -> None: ...  # v iiiiIITT*T
def bsr_matvecs(a0: _i, a1: _i, a2: _i, a3: _i, a4: _i, a5: _I, a6: _I, a7: _T, a8: _T, a9: _T, /) -> None: ...  # v iiiiiIITT*T
def bsr_elmul_bsr(
    a0: _i, a1: _i, a2: _i, a3: _i, a4: _I, a5: _I, a6: _T, a7: _I, a8: _I, a9: _T, a10: _I, a11: _I, a12: _T, /
) -> None: ...  # v iiiiIITIIT*I*I*T
def bsr_eldiv_bsr(
    a0: _i, a1: _i, a2: _i, a3: _i, a4: _I, a5: _I, a6: _T, a7: _I, a8: _I, a9: _T, a10: _I, a11: _I, a12: _T, /
) -> None: ...  # v iiiiIITIIT*I*I*T
def bsr_plus_bsr(
    a0: _i, a1: _i, a2: _i, a3: _i, a4: _I, a5: _I, a6: _T, a7: _I, a8: _I, a9: _T, a10: _I, a11: _I, a12: _T, /
) -> None: ...  # v iiiiIITIIT*I*I*T
def bsr_minus_bsr(
    a0: _i, a1: _i, a2: _i, a3: _i, a4: _I, a5: _I, a6: _T, a7: _I, a8: _I, a9: _T, a10: _I, a11: _I, a12: _T, /
) -> None: ...  # v iiiiIITIIT*I*I*T
def bsr_maximum_bsr(
    a0: _i, a1: _i, a2: _i, a3: _i, a4: _I, a5: _I, a6: _T, a7: _I, a8: _I, a9: _T, a10: _I, a11: _I, a12: _T, /
) -> None: ...  # v iiiiIITIIT*I*I*T
def bsr_minimum_bsr(
    a0: _i, a1: _i, a2: _i, a3: _i, a4: _I, a5: _I, a6: _T, a7: _I, a8: _I, a9: _T, a10: _I, a11: _I, a12: _T, /
) -> None: ...  # v iiiiIITIIT*I*I*T
def bsr_ne_bsr(
    a0: _i, a1: _i, a2: _i, a3: _i, a4: _I, a5: _I, a6: _T, a7: _I, a8: _I, a9: _T, a10: _I, a11: _I, a12: _T, /
) -> None: ...  # v iiiiIITIIT*I*I*T
def bsr_lt_bsr(
    a0: _i, a1: _i, a2: _i, a3: _i, a4: _I, a5: _I, a6: _T, a7: _I, a8: _I, a9: _T, a10: _I, a11: _I, a12: _T, /
) -> None: ...  # v iiiiIITIIT*I*I*T
def bsr_gt_bsr(
    a0: _i, a1: _i, a2: _i, a3: _i, a4: _I, a5: _I, a6: _T, a7: _I, a8: _I, a9: _T, a10: _I, a11: _I, a12: _T, /
) -> None: ...  # v iiiiIITIIT*I*I*T
def bsr_le_bsr(
    a0: _i, a1: _i, a2: _i, a3: _i, a4: _I, a5: _I, a6: _T, a7: _I, a8: _I, a9: _T, a10: _I, a11: _I, a12: _T, /
) -> None: ...  # v iiiiIITIIT*I*I*T
def bsr_ge_bsr(
    a0: _i, a1: _i, a2: _i, a3: _i, a4: _I, a5: _I, a6: _T, a7: _I, a8: _I, a9: _T, a10: _I, a11: _I, a12: _T, /
) -> None: ...  # v iiiiIITIIT*I*I*T

# csc.h
def csc_matvec(a0: _i, a1: _i, a2: _I, a3: _I, a4: _T, a5: _T, a6: _T, /) -> None: ...  # v iiIITT*T
def csc_matvecs(a0: _i, a1: _i, a2: _i, a3: _I, a4: _I, a5: _T, a6: _T, a7: _T, /) -> None: ...  # v iiiIITT*T

# csr.h
def csr_matmat_maxnnz(a0: _i, a1: _i, a2: _I, a3: _I, a4: _I, a5: _I, /) -> np.int64: ...  # l iiIIII
def csr_matmat(
    a0: _i, a1: _i, a2: _I, a3: _I, a4: _T, a5: _I, a6: _I, a7: _T, a8: _I, a9: _I, a10: _T, /
) -> None: ...  # v iiIITIIT*I*I*T
def csr_diagonal(a0: _i, a1: _i, a2: _i, a3: _I, a4: _I, a5: _T, a6: _T, /) -> None: ...  # v iiiIIT*T
def csr_tocsc(a0: _i, a1: _i, a2: _I, a3: _I, a4: _T, a5: _I, a6: _I, a7: _T, /) -> None: ...  # v iiIIT*I*I*T
def csr_tobsr(a0: _i, a1: _i, a2: _i, a3: _i, a4: _I, a5: _I, a6: _T, a7: _I, a8: _I, a9: _T, /) -> None: ...  # v iiiiIIT*I*I*T
def csr_todense(a0: _i, a1: _i, a2: _I, a3: _I, a4: _I, a5: _T, /) -> None: ...  # v iiIIT*T
def csr_matvec(a0: _i, a1: _i, a2: _I, a3: _I, a4: _T, a5: _T, a6: _T, /) -> None: ...  # v iiIITT*T
def csr_matvecs(a0: _i, a1: _i, a2: _i, a3: _I, a4: _I, a5: _T, a6: _T, a7: _T, /) -> None: ...  # v iiiIITT*T
def csr_elmul_csr(
    a0: _i, a1: _i, a2: _I, a3: _I, a4: _T, a5: _I, a6: _I, a7: _T, a8: _I, a9: _I, a10: _T, /
) -> None: ...  # v iiIITIIT*I*I*T
def csr_eldiv_csr(
    a0: _i, a1: _i, a2: _I, a3: _I, a4: _T, a5: _I, a6: _I, a7: _T, a8: _I, a9: _I, a10: _T, /
) -> None: ...  # v iiIITIIT*I*I*T
def csr_plus_csr(
    a0: _i, a1: _i, a2: _I, a3: _I, a4: _T, a5: _I, a6: _I, a7: _T, a8: _I, a9: _I, a10: _T, /
) -> None: ...  # v iiIITIIT*I*I*T
def csr_minus_csr(
    a0: _i, a1: _i, a2: _I, a3: _I, a4: _T, a5: _I, a6: _I, a7: _T, a8: _I, a9: _I, a10: _T, /
) -> None: ...  # v iiIITIIT*I*I*T
def csr_maximum_csr(
    a0: _i, a1: _i, a2: _I, a3: _I, a4: _T, a5: _I, a6: _I, a7: _T, a8: _I, a9: _I, a10: _T, /
) -> None: ...  # v iiIITIIT*I*I*T
def csr_minimum_csr(
    a0: _i, a1: _i, a2: _I, a3: _I, a4: _T, a5: _I, a6: _I, a7: _T, a8: _I, a9: _I, a10: _T, /
) -> None: ...  # v iiIITIIT*I*I*T
def csr_ne_csr(
    a0: _i, a1: _i, a2: _I, a3: _I, a4: _T, a5: _I, a6: _I, a7: _T, a8: _I, a9: _I, a10: _T, /
) -> None: ...  # v iiIITIIT*I*I*T
def csr_lt_csr(
    a0: _i, a1: _i, a2: _I, a3: _I, a4: _T, a5: _I, a6: _I, a7: _T, a8: _I, a9: _I, a10: _T, /
) -> None: ...  # v iiIITIIT*I*I*T
def csr_gt_csr(
    a0: _i, a1: _i, a2: _I, a3: _I, a4: _T, a5: _I, a6: _I, a7: _T, a8: _I, a9: _I, a10: _T, /
) -> None: ...  # v iiIITIIT*I*I*T
def csr_le_csr(
    a0: _i, a1: _i, a2: _I, a3: _I, a4: _T, a5: _I, a6: _I, a7: _T, a8: _I, a9: _I, a10: _T, /
) -> None: ...  # v iiIITIIT*I*I*T
def csr_ge_csr(
    a0: _i, a1: _i, a2: _I, a3: _I, a4: _T, a5: _I, a6: _I, a7: _T, a8: _I, a9: _I, a10: _T, /
) -> None: ...  # v iiIITIIT*I*I*T
def csr_scale_rows(a0: _i, a1: _i, a2: _I, a3: _I, a4: _T, a5: _T, /) -> None: ...  # v iiII*TT
def csr_scale_columns(a0: _i, a1: _i, a2: _I, a3: _I, a4: _T, a5: _T, /) -> None: ...  # v iiII*TT
def csr_sort_indices(a0: _i, a1: _I, a2: _I, a3: _T, /) -> None: ...  # v iI*I*T
def csr_eliminate_zeros(a0: _i, a1: _i, a2: _I, a3: _I, a4: _T, /) -> None: ...  # v ii*I*I*T
def csr_sum_duplicates(a0: _i, a1: _i, a2: _I, a3: _I, a4: _T, /) -> None: ...  # v ii*I*I*T
def get_csr_submatrix(
    a0: _i, a1: _i, a2: _I, a3: _I, a4: _T, a5: _i, a6: _i, a7: _i, a8: _i, a9: _V, a10: _V, a11: _W, /
) -> None: ...  # v iiIITiiii*V*V*W
def csr_row_index(a0: _i, a1: _I, a2: _I, a3: _I, a4: _T, a5: _I, a6: _T, /) -> None: ...  # v iIIIT*I*T
def csr_row_slice(a0: _i, a1: _i, a2: _i, a3: _I, a4: _I, a5: _T, a6: _I, a7: _T, /) -> None: ...  # v iiiIIT*I*T
def csr_column_index1(a0: _i, a1: _i, a2: _i, a3: _I, a4: _I, a5: _T, a6: _I, a7: _T, /) -> None: ...  # v iiiIIT*I*T
def csr_column_index2(a0: _I, a1: _I, a2: _i, a3: _I, a4: _T, a5: _I, a6: _T, /) -> None: ...  # v IIiIT*I*T
def csr_sample_values(a0: _i, a1: _i, a2: _I, a3: _I, a4: _T, a5: _i, a6: _I, a7: _I, a8: _T, /) -> None: ...  # v iiIITiII*T
def csr_count_blocks(a0: _i, a1: _i, a2: _i, a3: _i, a4: _I, a5: _I, /) -> int: ...  # i iiiiII
def csr_sample_offsets(a0: _i, a1: _i, a2: _I, a3: _I, a4: _i, a5: _I, a6: _I, a7: _I, /) -> int: ...  # i iiIIiII*I
def csr_hstack(a0: _i, a1: _i, a2: _I, a3: _I, a4: _I, a5: _T, a6: _I, a7: _I, a8: _T, /) -> None: ...  # v iiIIIT*I*I*T
def csr_has_sorted_indices(a0: _i, a1: _I, a2: _I, /) -> int: ...  # i iII
def csr_has_canonical_format(a0: _i, a1: _I, a2: _I, /) -> int: ...  # i iII
def expandptr(a0: _i, a1: _I, a2: _I, /) -> None: ...  # v iI*I
def test_throw_error() -> Never: ...  # i  # raises `MemoryError` if called

# coo.h
def coo_tocsr(a0: _i, a1: _i, a2: _i, a3: _I, a4: _I, a5: _T, a6: _I, a7: _I, a8: _T, /) -> None: ...  # v iiiIIT*I*I*T
def coo_todense(a0: _i, a1: _i, a2: _l, a3: _I, a4: _I, a5: _T, a6: _T, a7: _i, /) -> None: ...  # v iilIIT*Ti
def coo_todense_nd(a0: _I, a1: _l, a2: _l, a3: _I, a4: _T, a5: _T, a6: _i, /) -> None: ...  # v IllIT*Ti
def coo_matvec(a0: _l, a1: _I, a2: _I, a3: _T, a4: _T, a5: _T, /) -> None: ...  # v lIITT*T
def coo_matvec_nd(a0: _l, a1: _l, a2: _I, a3: _I, a4: _T, a5: _T, a6: _T, /) -> None: ...  # v llIITT*T
def coo_matmat_dense(a0: _l, a1: _l, a2: _I, a3: _I, a4: _T, a5: _T, a6: _T, /) -> None: ...  # v llIITT*T
def coo_matmat_dense_nd(a0: _l, a1: _l, a2: _l, a3: _I, a4: _I, a5: _I, a6: _T, a7: _T, a8: _T, /) -> None: ...  # v lllIIITT*T

# dia.h
def dia_matmat(
    a0: _i, a1: _i, a2: _i, a3: _i, a4: _I, a5: _T, a6: _i, a7: _i, a8: _i, a9: _I, a10: _T, a11: _V, a12: _W, /
) -> None: ...  # v iiiiITiiiIT*V*W
def dia_matvec(a0: _i, a1: _i, a2: _i, a3: _i, a4: _T, a5: _T, a6: _T, a7: _T, /) -> None: ...  # v iiiiITT*T
def dia_matvecs(a0: _i, a1: _i, a2: _i, a3: _i, a4: _I, a5: _T, a6: _i, a7: _T, a8: _T, /) -> None: ...  # v iiiiITiT*T

# csgraph.h
def cs_graph_components(a0: _i, a1: _I, a2: _I, a3: _I, /) -> int: ...  # i iII*I
