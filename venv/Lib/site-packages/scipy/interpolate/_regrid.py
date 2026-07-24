"""
Regrid (2-D smoothing B-splines via separable 1-D FITPACK kernels)
==================================================================

1) Overview
-----------
This module fits a bivariate tensor-product B-spline surface to gridded data
`Z[i, j]` sampled on strictly increasing coordinates `x[i]`, `y[j]`. It mirrors
the adaptive spirit of FITPACK's REGRID/fpgrre (knot growth + smoothing
parameter search to meet a target residual `s`) while keeping the control flow
explicit in Python and returning a 2-D `NdBSpline`.

**Key conventions**
- **Penalty scaling is 1/p** - *same as FITPACK*. The augmented system stacks
  `(D / p)` under the data matrix. Smaller `p` -> stronger smoothing; larger `p`
  -> weaker smoothing (approaching interpolation).
- `p == -1` is used **as a sentinel for `p = ∞`** (interpolatory limit): when
  seen, the solver omits penalty rows entirely.

2) Mathematics (and how this differs from FITPACK's REGRID)
-----------------------------------------------------------
Let `A_x`, `A_y` be banded 1-D design matrices and `D_x`, `D_y` the banded
1-D roughness (difference) matrices returned by FITPACK/Dierckx's 1-D APIs
(`data_matrix`, `disc`). The 2-D smoothing objective is

    minimize  ||A c - z||^2 + (1/p) ||D c||^2,            (1)

implemented by the **augmented** least squares system

    [ A ] c ~ [ z ]
    [D/p]     [ 0 ].

**Same as FITPACK:** Equation (1) uses the *1/p* convention.
**Different in this module:** We realize the 2-D problem by composing **1-D**
banded operators in two separable passes (x then y), instead of calling a
monolithic 2-D routine. This makes the mathematics transparent while producing
the same normal equations structure that REGRID targets internally.

3) Residual energy: definition and use
--------------------------------------
After solving on the current knots (initially at the interpolatory limit,
`p = ∞`), we evaluate the surface and form

    R = Z - Zhat,                         fp = sum(R[i, j]^2).

We then **project** residual energy to knot spans along each axis:

- Row energy  `row_energy[i] = sum(R[i, j]^2, j)`  -> accumulated into `fpintx`.
- Column energy `col_energy[j] = sum(R[i, j]^2, i)` -> accumulated into `fpinty`.

These per-span energies guide **adaptive knot insertion**: the algorithm picks
high-energy spans (skipping zero-length ones) and inserts data-aligned knots
(e.g., median sample within the span), with simple batch sizing and headroom
limits. This is conceptually the same idea as REGRID's `fpint` arrays; here it
is implemented explicitly and vectorized in Python for clarity.

4) Solver used here vs. FITPACK's REGRID
----------------------------------------
**This module (separable 1-D composition):**
- Uses **1-D FITPACK/Dierckx kernels** (`data_matrix`, `disc`, `qr_reduce`,
  `fpback`) to solve 2-D via two passes:
  1) augment/QR/backsolve along **x**, producing an intermediate;
  2) augment/QR/backsolve along **y**, producing the coefficient grid `C`.
- The augmented rows are stacked as `[A; D/p]` (1/p scaling, **same as FITPACK**).
- If `fp > s` after knot growth, performs a **scalar search in `p`** using a
  ratio-of-roots routine; `p == -1` is treated as `p = ∞` for the interpolatory
  reference without penalty rows.

**FITPACK REGRID (monolithic 2-D routine):**
- Implements the same mathematical objective and **1/p** penalty scaling,
  but inside a specialized, fully 2-D Fortran routine with in-place Givens/QR
  and a rational update for `p`.
- Handles residual partitioning, knot insertion, and `p` updates internally.

**Practical difference:** We build the 2-D solve from **1-D building blocks**,
which keeps each stage observable/testable and uses the same low-level kernels
as FITPACK, but without relying on the single monolithic REGRID entry point.

5) Execution flow (who calls whom)
----------------------------------

**Top-level**
1. `_regrid`
   - Validates shapes/monotonicity and normalizes `bbox`.
   - Dispatches to `_regrid_fitpack(...)`.

**Core driver**
2. `_regrid_fitpack`
   - Clips inputs to `bbox` via `_apply_bbox_grid` -> `(x_fit, y_fit, Z_fit)`.
   - If `s == 0` (interpolation path):
     - Build initial not-a-knot vectors: `tx = _not_a_knot(x_fit, kx)`,
       `ty = _not_a_knot(y_fit, ky)`.
     - Build design matrices: `(Ax, Ay, Q) = build_design_matrices(...)`.
     - Solve once at `p = -1` (inf): `C0, fp = _solve_2d_fitpack(...)`.
     - Return `return_NdBSpline(...)`.
   - Else (`s > 0`, smoothing with adaptive knots):
     - Initialize clamped no-interior-knot vectors with
       `_initialise_knots` (for x and y).
     - **Knot-growth loop** (up to `len(x_fit)+len(y_fit)` iterations):
       1) Build design matrices: `(Ax, Ay, Q) = build_design_matrices(...)`.
       2) Interpolatory reference solve (`p = -1`): `C0, fp = _solve_2d_fitpack(...)`.
       3) If both `tx`, `ty` are at minimal size, record `fp0 = fp`
          (used to bracket the p-search).
       4) If `fp < s`: **stop growth** and proceed to p-search / finalize.
       5) Compute residuals for knot insertion:
          - Convert packed to CSR once per axis for evaluation:
            `_Ax = Ax.tocsr(...)`, `_Ay = Ay.tocsr(...)`
          - Evaluate `Z0 = _Ax @ C0 @ _Ay.T`, residual `R = Z_fit - Z0`.
       6) Insert knots on alternating axes using energy heuristics:
          - If last axis was `y`, grow `tx` with
            `_add_knots(..., residuals=np.sum(R**2, axis=1), ...)`.
          - Else, grow `ty` with `_add_knots(..., residuals=np.sum(R**2, axis=0), ...)`.
          - Update `fpold`, `nplus{ x|y }`, `last_axis`.
     - If growth ended with both axes still minimal: return `return_NdBSpline(...)`.
     - **Finite-p smoothing (if needed):**
       - Build 1-D penalty operators: `(Drx, _, ncx) = disc(tx, kx)`,
         `(Dry, _, ncy) = disc(ty, ky)`, wrap as `PackedMatrix`.
       - Rebuild design matrices on the final `(tx, ty)`.
       - Run `_p_search_hit_s(...)` to find `p*` such that `fp(p*) ~ s`:
         - Internally constructs `F(...)` which evaluates `fp(p)`
           by calling `_solve_2d_fitpack(...)`.
         - Uses `root_rati` on `g(p) = fp(p) - s` with the interpolatory
           reference at `p = inf` and `fp0` bracket.
       - Return `return_NdBSpline(...)`.

**Separable solver (used by both interpolation and p-search)**
3. `_solve_2d_fitpack`
   - Forms augmented banded systems via `_stack_augmented_fitpack`:
     - X-side: `Ax_aug = [Ax; Dx/p]` if `p != -1`, else just `Ax`.
     - Y-side: `Ay_aug = [Ay; Dy/p]` if `p != -1`, else just `Ay`.
     - Pads RHS `Q` with zeros when penalties are stacked.
   - **Stage X**: `_dierckx.qr_reduce(Ax_aug, ...)` then
     `_dierckx.fpback(...)` -> intermediate `T`.
   - Transpose `T` to feed Y solves column-wise; pad if needed.
   - **Stage Y**: `_dierckx.qr_reduce(Ay_aug, ...)` then
     `_dierckx.fpback(...)` -> coefficients `C`.
   - Build CSR design matrices once: `_Ax = Ax.tocsr(...)`, `_Ay = Ay.tocsr(...)`.
   - Evaluate `zhat = _Ax @ C.T @ _Ay.T`, compute `fp = ||Z_fit - zhat||^2`.
   - Return `C.T` (in `(nx_coef, ny_coef)` layout) and `fp`.

**Helpers**
- `_apply_bbox_grid(...)` - slices to `(x_fit, y_fit, Z_fit)`.
- `build_design_matrices(...)` - wraps `_dierckx.data_matrix` and
  returns `PackedMatrix` wrappers + `Q`.
- `_initialise_knots(...)`, `_add_knots(...)` - FITPACK-style knot bookkeeping/growth.
- `disc(...)` - 1-D roughness operators (packed band form).
- `F` + `root_rati` - maps `p -> fp(p)` and finds `p*` with `fp(p*) ~ s`.
- `return_NdBSpline(...)` - packs `(tx, ty, C)` into an `NdBSpline`.
"""
import numpy as np
from scipy.interpolate._ndbspline import NdBSpline
from scipy.interpolate._fitpack_repro import (
    root_rati, disc, add_knot, _not_a_knot)
from . import _dierckx      # type: ignore[attr-defined]
from scipy.sparse import csr_array
from scipy._lib._util import _validate_int

def _ndbspline_call_like_bivariate(ndbs, x, y, dx=0, dy=0, grid=True):
    """
    Evaluate a 2D `NdBSpline` like a classical bivariate API.

    Parameters
    ----------
    ndbs : NdBSpline
        A 2D spline object (``len(ndbs.t) == 2``).
    x, y : array_like
        Sample locations. If ``grid=True``, these must be 1-D strictly
        increasing vectors. If ``grid=False``, they can be broadcastable
        arrays of the same shape.
    dx, dy : int, optional
        Derivative orders along `x` and `y` respectively, by default 0.
    grid : bool, optional
        If True, evaluate on the cartesian product of `x` and `y`;
        otherwise treat `(x, y)` as paired coordinates, by default True.

    Returns
    -------
    ndarray or (ndarray, dict)
        Evaluated values with shape:
        - ``(len(x), len(y), ...)`` if ``grid=True``.
        - ``x.shape + ...`` if ``grid=False``.

    Raises
    ------
    ValueError
        If `ndbs` is not 2D, derivatives are negative, or monotonicity checks fail.

    Notes
    -----
    This is a thin convenience wrapper around ``NdBSpline.__call__`` with input
    validation and optional profiling.
    """
    if len(ndbs.t) != 2:
        raise ValueError("ndbs must be a 2D NdBSpline (len(t) == 2).")

    dx = _validate_int(dx, 'dx')
    dy = _validate_int(dy, 'dy')
    if dx < 0 or dy < 0:
        raise ValueError("order of derivative must be positive or zero")

    trailing = ndbs.c.shape[2:]
    x = np.asarray(x)
    y = np.asarray(y)

    if grid:
        if x.size == 0 or y.size == 0:
            vals = np.zeros((x.size, y.size) + trailing, dtype=ndbs.c.dtype)
            return vals

        if (x.size >= 2) and (not np.all(np.diff(x) >= 0.0)):
            raise ValueError("x must be strictly increasing when `grid` is True")
        if (y.size >= 2) and (not np.all(np.diff(y) >= 0.0)):
            raise ValueError("y must be strictly increasing when `grid` is True")

        X, Y = np.meshgrid(x, y, indexing="ij")
        xi = np.stack((X, Y), axis=-1)  # (len(x), len(y), 2)

        vals = ndbs(xi, nu=(dx, dy), extrapolate=ndbs.extrapolate)

        return vals
    else:
        if x.shape != y.shape:
            x, y = np.broadcast_arrays(x, y)

        if x.size == 0:
            return np.zeros(x.shape + trailing, dtype=ndbs.c.dtype)
        xi = np.stack((x.ravel(), y.ravel()), axis=-1)
        vals = ndbs(xi, nu=(dx, dy), extrapolate=ndbs.extrapolate)
        return vals.reshape(x.shape + trailing)


def return_NdBSpline(fp, tck, degrees):
    """
    Build a 2D ``NdBSpline`` from knot vectors and a coefficient grid.

    Parameters
    ----------
    fp : float
        Residual sum of squares of the produced fit (kept for upstream use).
    tck : tuple
        Tuple ``(tx, ty, C)`` where ``tx``, ``ty`` are knot vectors and ``C``
        is a coefficient array with shape ``(nx - kx - 1, ny - ky - 1)`` or
        a compatible shape that can be reshaped to that.
    degrees : tuple of int
        Degrees ``(kx, ky)`` along x and y.

    Returns
    -------
    NdBSpline
        The constructed 2D spline.

    Notes
    -----
    Only repacks the coefficient grid; ``fp`` is not used internally here.
    """
    nx, ny = len(tck[0]), len(tck[1])
    kx, ky = degrees
    c = tck[2].reshape(nx - kx - 1, ny - ky - 1)
    return NdBSpline((tck[0], tck[1]), c, degrees)


class PackedMatrix:
    """A simplified CSR format for when non-zeros in each row are consecutive.

    Assuming that each row of an `(m, nc)` matrix 1) only has `nz` non-zeros, and
    2) these non-zeros are consecutive, we only store an `(m, nz)` matrix of
    non-zeros and a 1D array of row offsets. This way, a row `i` of the original
    matrix A is ``A[i, offset[i]: offset[i] + nz]``.

    """
    def __init__(self, a, offset, nc):
        self.a = a
        self.offset = offset
        self.nc = nc

    @property
    def shape(self):
        return self.a.shape[0], self.nc

    def todense(self):
        out = np.zeros(self.shape)
        nelem = self.a.shape[1]
        for i in range(out.shape[0]):
            nel = min(self.nc - self.offset[i], nelem)
            out[i, self.offset[i]:self.offset[i] + nel] = self.a[i, :nel]
        return out

    def tocsr(self, k, m, len_t):
        # Inlined from https://github.com/scipy/scipy/blob/v1.16.2/scipy/interpolate/_bsplines.py#L455-L466
        data = self.a.ravel()

        # Convert from per-row offsets to the CSR indices/indptr format
        indices = np.repeat(self.offset, k+1).reshape(-1, k+1)
        indices = indices + np.arange(k+1, dtype=self.offset.dtype)
        indices = indices.ravel()

        indptr = np.arange(0, (m + 1) * (k + 1), k + 1,
                           dtype=self.offset.dtype)

        return csr_array(
            (data, indices, indptr),
            shape=(m, len_t - k - 1)
    )


def _stack_augmented_fitpack(A, D, nc, k, p):
    """
    Builds augmented banded matrix.

    Parameters
    ----------
    A : PackedMatrix
        Banded data/design matrix for one axis (from `_dierckx.data_matrix`).
    D : PackedMatrix
        Banded roughness (difference) penalty matrix for the same axis
        (from `disc`).
    nc : int
        Number of top (data) rows from `A` to include.
    k : int
        Spline degree (used only for sizing in the current implementation).
    p : float
        Smoothing parameter. The effective penalization term is scaled as **1/p**:
        larger `p` means *less* smoothing (approaching interpolation).
        If `p == -1`, it signals *p -> inf*, i.e. a pure interpolatory system
        with **no** penalty rows appended.

    Returns
    -------
    AA : ndarray
        Augmented banded matrix with `A` stacked over `(D / p)` when `p != -1`.
    offset : ndarray
        Concatenated band offsets for the augmented matrix.
    nc : int
        Returned unchanged for downstream convenience.
    """
    if p == -1:
        return A.a.copy(), A.offset.copy(), nc

    nz = k + 1
    AA = np.zeros((nc + D.shape[0], k + 2), dtype=float)
    AA[:nc, :nz] = A.a[:nc, :]
    AA[nc:, :] = D.a / p
    offset = np.concatenate((A.offset, D.offset))
    return AA, offset, nc

def _solve_2d_fitpack(Ax, Ay, Q, p,
                      kx, tx, x_x,
                      ky, ty, x_y, z,
                      Dx=None, Dy=None):
    """
    Solve the 2-D tensor-product spline system using separable banded QR.

    ================================================================
    Mathematical model (step by step, plain text)
    ================================================================

    Shapes:
        Z      : (mx, my)  -> original data
        Ax, Ay : design matrices for x and y
        Dx, Dy : roughness penalty matrices for x and y
        C      : (nx, ny)  -> spline coefficients to solve for

    Surface approximation:
        Zhat = Ax * C * Ay^T

    Objective (smoothing formulation):
        minimize ||Ax*C*Ay^T - Z||^2 + (1/p)*(||Dx*C||^2 + ||C*Dy^T||^2)

    In practice (FITPACK-style separable approach), we solve this in two stages:

    --------------------------------------------------------
    Stage 1 (x-direction solve for all y-columns together):
    --------------------------------------------------------

        For each column of Z:
            minimize ||Ax*T - Z||^2 + (1/p)*||Dx*T||^2

        This is equivalent to the augmented least-squares system:
            [Ax]       [Z]
            [Dx/p] * T = [0]

        i.e.  minimize || [Ax; Dx/p]*T - [Z; 0] ||^2

        The solution T is obtained by QR reduction and back-substitution.

    --------------------------------------------------------
    Stage 2 (y-direction solve using transposed result):
    --------------------------------------------------------
        Now treat T^T as the new RHS for the y-direction:
            minimize ||Ay*C^T - T^T||^2 + (1/p)*||Dy*C^T||^2

        Equivalent to augmented system:
            [Ay]       [T^T]
            [Dy/p] * C^T = [0]

        i.e.  minimize || [Ay; Dy/p]*C^T - [T^T; 0] ||^2

        Solving this gives C^T (then transposed back to C).

    --------------------------------------------------------
    Interpolation limit:
    --------------------------------------------------------
        If p == -1, penalties are omitted (Dx, Dy are not stacked).
        The solver behaves as a near-interpolating system.

    --------------------------------------------------------
    Residual computation:
    --------------------------------------------------------
        Zhat = Ax * C * Ay^T
        R    = Z - Zhat
        fp   = sum(R^2)

    Parameters
    ----------
    Ax, Ay : PackedMatrix
        Banded data matrices for the x and y axes.
    Q : ndarray, shape (mx, my)
        RHS data grid (copied from `Z`).
    p : float
        Smoothing parameter. The penalty term is scaled as **1/p**.
        Setting `p == -1` signals *p -> inf* (interpolation, omit penalty).
    kx, ky : int
        Spline degrees along x and y.
    tx, ty : ndarray
        Knot vectors along x and y.
    x_x, x_y : ndarray
        Sample coordinates.
    z : ndarray
        Original data grid for residual evaluation.
    Dx, Dy : ndarray
        Banded roughness penalty matrices for x and y.
        Optional, Only needed when ``p != -1``.

    Returns
    -------
    C : ndarray
        2-D B-spline coefficient grid.
    fp : float
        Residual sum of squares between fitted surface and `z`.
    R : ndarray, shape (mx, my)
        Residual matrix ``z - zhat``, where ``zhat = Ax @ C @ Ay.T``.

    Notes
    -----
    This performs two separable QR solves (x then y), each augmented by
    `(D / p)` when `p != -1`.  Setting `p = -1` skips all penalty rows,
    yielding an interpolatory surface.  The resulting coefficients and residual
    follow the same conventions as FITPACK's `fpgrre`.
    """
    # Dummy unit weights for FITPACK fpback APIs.
    w_x = np.ones_like(x_x)
    w_y = np.ones_like(x_y)

    # https://github.com/scipy/scipy/blob/v1.16.2/scipy/interpolate/fitpack/fpgrre.f#L97-L105
    # Build the augmented banded matrix for x:
    #   - If p != -1, stack (Dx / p) under Ax for FITPACK-style smoothing.
    #   - If p == -1, _stack_augmented_fitpack omits the penalty part entirely.
    # Returns:
    #   Ax_aug      : augmented banded matrix (data [+ penalty]).
    #   offset_aug_x: band offsets compatible with Ax_aug.
    #   nc_augx     : number of top data rows within Ax_aug (== ncx).
    Ax_aug, offset_aug_x, nc_augx = _stack_augmented_fitpack(
        Ax, Dx, Ax.shape[0], kx, p)
    nc_x = Ax.nc

    # Same for y: build Ay_aug with (Dy / p) stacked if p != -1.
    Ay_aug, offset_aug_y, nc_augy = _stack_augmented_fitpack(
        Ay, Dy, Ay.shape[0], ky, p)
    nc_y = Ay.nc

    # If we stacked penalty rows on the x side, the RHS must be padded with zeros
    # to match the augmented row count for the QR reduction call.
    if p != -1: # https://github.com/scipy/scipy/blob/v1.16.2/scipy/interpolate/fitpack/fpgrre.f#L97
        # Dx.shape[0] is the number of penalty rows; add that many zero rows
        # so Ax_aug and Q have compatible leading dimensions for in-place QR.
        # https://github.com/scipy/scipy/blob/v1.16.2/scipy/interpolate/fitpack/fpgrre.f#L110-L118
        Q = np.vstack([Q, np.zeros((Dx.shape[0], Q.shape[1]), dtype=float)])

    # https://github.com/scipy/scipy/blob/v1.16.2/scipy/interpolate/fitpack/fpgrre.f#L106-L175
    # Perform in-place banded QR reduction of the x-augmented system:
    # This orthogonalizes/eliminates along x for all RHS columns in Q simultaneously.
    # After this, fpback can do x-direction back-substitution to
    # get c^T (partial coeffs).
    _dierckx.qr_reduce(Ax_aug, offset_aug_x, nc_augx, Q)

    # https://github.com/scipy/scipy/blob/v1.16.2/scipy/interpolate/fitpack/fpgrre.f#L246-L253
    # Back-substitute along x to solve the reduced system:
    #   cT has shape (ncoef_x, num_y_data) in this calling pattern, i.e. per y-column.
    # The API uses:
    #   Ax_aug, nc_augx: reduced upper structure
    #   x_x, tx, kx, w_x: x-sample grid, knot vector, degree, and (unit) weights
    #   Q: RHS (current)
    # fpback returns shape (nc_x, num_y_data)
    T, _, _ = _dierckx.fpback(
        Ax_aug, nc_x, x_x,
        Q, tx, kx, w_x,
        Q, False
    )

    # We now want to treat the *y*-direction solve with these as the new RHS.
    # Transpose so each column corresponds to a y-solve RHS consistently.
    Q = np.ascontiguousarray(T.T)

    # If we stacked penalty rows on the y side, pad RHS with zeros to match Ay_aug.
    if p != -1: # https://github.com/scipy/scipy/blob/v1.16.2/scipy/interpolate/fitpack/fpgrre.f#L97
        # https://github.com/scipy/scipy/blob/v1.16.2/scipy/interpolate/fitpack/fpgrre.f#L110-L118
        Q = np.vstack([Q, np.zeros((Dy.shape[0], Q.shape[1]), dtype=float)])

    # https://github.com/scipy/scipy/blob/v1.16.2/scipy/interpolate/fitpack/fpgrre.f#L176-L245
    # Perform in-place banded QR reduction along y for all columns of Q.
    _dierckx.qr_reduce(Ay_aug, offset_aug_y, nc_augy, Q)

    # https://github.com/scipy/scipy/blob/v1.16.2/scipy/interpolate/fitpack/fpgrre.f#L254-L269
    # Final back-substitution along y:
    # Returns:
    #   C  : coefficient matrix with shape (nc_y, nc_x)
    #   fp : FITPACK's internal residual metric from the y-solve
    #        (we recompute below anyway)
    C, _, fp = _dierckx.fpback(
        Ay_aug, nc_y,
        x_y, Q, ty, ky, w_y,    # y-grid, y-knots, degree, weights
        Q,                      # RHS -> solution becomes coefficients along y
        False
    )

    # Build explicit design matrices to evaluate the fitted surface:
    # _Ax: [mx × nx_coef], _Ay: [my × ny_coef]
    # Note: We call PackedMatrix.tocsr here because matrix multiplication
    # with the packed banded format (returned by _dierckx.data_matrix)
    # is not implemented. PackedMatrix.tocsr returns the design matrix,
    # in CSR format, that supports standard @ operations for residual
    # evaluation and diagnostics.
    _Ax = Ax.tocsr(kx, x_x.shape[0], len(tx))
    _Ay = Ay.tocsr(ky, x_y.shape[0], len(ty))

    # Evaluate the fitted surface: zhat = Ax * C^T * Ay^T
    # Note: C currently aligns so that C.T matches x-first multiplication order.
    zhat = _Ax @ C.T @ _Ay.T

    # Compute residual matrix R and fp in one pass; return R so callers that
    # need per-span energy for knot placement can reuse it directly instead of
    # recomputing zhat a second time.
    R = z - zhat
    fp = np.sum(np.square(R))

    # Return coefficients in the conventional (nx_coef, ny_coef) orientation,
    # fp, and the residual matrix R.
    return C.T, fp, R

class F:
    """
    Callable wrapper for computing `fp(p)` for a fixed spline configuration.

    Parameters
    ----------
    Ax, Ay : PackedMatrix
        Banded data matrices.
    Dx, Dy : PackedMatrix
        Banded penalty matrices.
    kx, ky : int
        Degrees along x and y.
    tx, ty : ndarray
        Knot vectors along x and y.
    x_x, x_y : ndarray
        Sample coordinates.
    w_x, w_y : ndarray
        Weights (usually ones).
    z : ndarray
        Data grid for computing the residual.

    Attributes
    ----------
    C : ndarray
        Coefficient matrix from the most recent solve.
    fp : float
        Residual value from the most recent solve.

    Notes
    -----
    The penalty is applied as **1/p**, so smaller `p` values yield heavier
    smoothing. Setting `p == -1` corresponds to *p = inf*, i.e. interpolation.
    Intended for use by `_p_search_hit_s` to iteratively evaluate `fp(p)`.
    """

    def __init__(self, Ax, Dx, Ay, Dy, Q,
                 kx, tx, x_x, ky, ty,
                 x_y, z):
        self.Ax = Ax
        self.Dx = Dx
        self.Ay = Ay
        self.Dy = Dy
        self.Q = Q
        self.kx = kx
        self.tx = tx
        self.x_x = x_x
        self.ky = ky
        self.ty = ty
        self.x_y = x_y
        self.z = z

    def __call__(self, p):
        # _stack_augmented_fitpack always allocates fresh arrays (np.zeros for
        # p != -1, .copy() for p == -1), so qr_reduce never writes back into
        # the original PackedMatrix.a buffers.  The four copies below are
        # therefore unnecessary and are removed to reduce memory pressure.
        C, fp, _ = _solve_2d_fitpack(
            self.Ax, self.Ay, self.Q.copy(),
            p, self.kx, self.tx, self.x_x,
            self.ky, self.ty, self.x_y, self.z,
            Dx=self.Dx, Dy=self.Dy)
        self.C = C
        self.fp = fp
        return fp

# https://github.com/scipy/scipy/blob/v1.16.2/scipy/interpolate/fitpack/fpregr.f#L301-L367
def _p_search_hit_s(
    Ax, Dx, Ay, Dy, Q, kx,
    tx, x_x, ky, ty, x_y, z, s, fp0, *,
    p_init=1.0, tol_rel=1e-3, maxit=40):
    """
    Search for a smoothing parameter `p` such that `fp(p) ~ s`.

    Parameters
    ----------
    Ax, Ay : PackedMatrix
        Banded data matrices.
    Dx, Dy : PackedMatrix
        Banded penalty matrices.
    Q : ndarray
        RHS data grid (copy of `Z`).
    kx, ky : int
        Spline degrees.
    tx, ty : ndarray
        Knot vectors.
    x_x, x_y : ndarray
        Sample coordinates.
    w_x, w_y : ndarray
        Sample weights.
    z : ndarray
        Original data grid for residuals.
    s : float
        Target smoothing residual (`fp` target).
    fp0 : float or None
        Residual at `p = inf` (interpolatory limit,
                               represented by `p == -1`).
    p_init : float, optional
        Starting guess for the finite `p` search, default 1.0.
    tol_rel : float, optional
        Relative tolerance for matching `fp(p)` to `s`.
    maxit : int, optional
        Maximum iterations for the root search.

    Returns
    -------
    p_star : float
        Smoothing parameter for which `fp(p_star)` ~ `s`.
    C_star : ndarray
        Coefficient grid corresponding to `p_star`.
    fp_star : float
        Residual at `p_star`.

    Notes
    -----
    The solver treats `p == -1` as *p = inf* (interpolatory, no penalty).
    For finite `p`, the penalty scales as **1/p** - smaller `p` increases
    smoothing. A ratio-of-roots search (`root_rati`) iteratively adjusts `p`
    until the residual `fp(p)` matches the target `s` within tolerance.
    """

    fp_at = F(Ax, Dx, Ay, Dy, Q, kx,
              tx, x_x, ky, ty, x_y, z)

    def g(p):
        return fp_at(p) - s

    fpms = g(-1)

    bracket = ((0.0, fp0 - s), (np.inf, fpms))
    ftol = max(s * tol_rel, 1e-12)

    r = root_rati(g, p_init, bracket, ftol, maxit=maxit)
    p_star = r.root
    fp_star = fp_at(p_star)
    C_star = fp_at.C

    return p_star, C_star, fp_star


def _apply_bbox_grid(x, y, Z, bbox):
    """
    Restrict (x, y, Z) to a rectangular bounding box.

    Parameters
    ----------
    x, y : ndarray
        Monotonic sample coordinates.
    Z : ndarray, shape (len(x), len(y))
        Data grid.
    bbox : sequence of 4 scalars or None
        ``(xb, xe, yb, ye)``; any element may be None to skip clipping.

    Returns
    -------
    x_fit, y_fit, Z_fit : ndarray
        Sliced arrays restricted to bbox.
    ix, iy : slice or ndarray
        Indexers mapping from full arrays to the restricted ones.

    Raises
    ------
    ValueError
        If bbox is invalid or excludes all samples along an axis.
    """
    if all([bboxi is None for bboxi in bbox]):
        return x, y, Z, slice(None), slice(None)

    xb, xe, yb, ye = bbox
    if not (xb < xe and yb < ye):
        raise ValueError("bbox must satisfy xb < xe and yb < ye")

    ix = np.where((x >= xb) & (x <= xe))[0]
    iy = np.where((y >= yb) & (y <= ye))[0]
    if ix.size == 0 or iy.size == 0:
        raise ValueError("bbox excludes all samples in x or y.")

    return x[ix], y[iy], Z[np.ix_(ix, iy)], np.s_[ix], np.s_[iy]


def _build_design_matrices(x, y, z, tx, ty, kx, ky):

    w_x = np.ones_like(x)
    w_y = np.ones_like(y)

    Ax, offset_x, nc_x = _dierckx.data_matrix(x, tx, kx, w_x)
    Ay, offset_y, nc_y = _dierckx.data_matrix(y, ty, ky, w_y)
    Q = z.copy()

    return (PackedMatrix(Ax, offset_x, nc_x),
            PackedMatrix(Ay, offset_y, nc_y),
            Q)


def _initialise_knots(m, xb, xe, k, nest=None):
    """
    Initialize a non-periodic knot vector.

    Parameters
    ----------
    m : int
        Number of data points (equivalent to len(x) if x were provided).
    xb, xe : float
        Domain endpoints used to seed the initial knot vector with no internal knots.
    k : int
        Spline degree.
    nest : int, optional
        Storage cap for knots. If None, defaults to max(m + k + 1, 2*k + 3).
        Must satisfy nest >= 2*(k + 1); otherwise a ValueError is raised.

    Returns
    -------
    t : 1-D ndarray
        Initial knot vector with no internal knots: [xb]*(k+1) + [xe]*(k+1).
    nest : int
        The finalized storage cap for knots.
    nmin : int
        Lower bound on knot count.
    nmax : int
        Upper bound on knot count.

    What this does
    --------------
    - Computes defaults and bounds used by FITPACK-style knot growth:
        * nest: storage cap for knots (defaults to max(m + k + 1, 2*k + 3))
        * nmin: minimal knot count (2*(k+1))
        * nmax: maximal knot count (m + k + 1)
    - Returns an initial knot vector with no internal knots:
        t = [xb]*(k+1) + [xe]*(k+1)
    """
    if nest is None:
        nest = max(m + k + 1, 2*k + 3)
    else:
        if nest < 2*(k + 1):
            raise ValueError(f"`nest` too small: {nest = } < 2*(k+1) = {2*(k+1)}.")

    nmin = 2*(k + 1)    # the number of knots for an LSQ polynomial approximation
    nmax = m + k + 1  # the number of knots for the spline interpolation

    # start from no internal knots
    t = np.asarray([xb]*(k+1) + [xe]*(k+1))

    return t, nest, nmin, nmax


TOL = 0.001

def _add_knots(x, k, s, t, nmin, nmax,
               nest, fp, fpold,
               residuals, nplus):
    """
    Knot-growth helper for knot-finding loop (non-periodic).

    Parameters
    ----------
    x : 1-D ndarray
        Strictly increasing sample coordinates.
    k : int
        Spline degree.
    s : float
        Target smoothing.
    t : 1-D ndarray
        Current knot vector to be grown.
    nmin, nmax : int
        Lower/upper bounds on knot count (from initialisation).
    nest : int
        Storage cap for total knots.
    fp, fpold : float
        Current and previous residual sums of squares. Used to update nplus.
    residuals : 1-D ndarray
        Most recent residual signal used by `add_knot` to decide placement.
    nplus : int
        Previous iteration's proposed number of knots; used to update the next nplus.

    Returns
    -------
    t_new, nplus : tuple
        Updated knot vector and the nplus chosen for this step.
        If n >= nmax, t_new is a not-a-knot layout. If n >= nest, t_new is the
        current vector respecting the storage cap.

    What this function does
    -----------------------
    - Assumes the caller has already decided to GROW (i.e., checks
      like |fp - s| < acc or fp < s has FAILED).
    - Updates nplus (how many knots to add next) using the FITPACK heuristic
      based on the previous improvement (delta = fpold - fp).
    - Inserts up to nplus new internal knots using `add_knot(x, t, k, residuals)`.
    - Stops early if storage or interpolation caps are reached:
        * if n >= nmax: switch to interpolation layout (not-a-knot) and return
        * if n >= nest: return current t respecting storage cap

    How it compares with _fitpack_repro.py::_generate_knots_impl
    -------------------------------------------------------------
    Similarities:
      1) Same growth logic for nplus:
         - Use delta = fpold - fp with ratio fpms/delta
         - Apply min/max caps (doubling and halving behavior)
      3) Same storage guard:
         - If n reaches nest, stop and return current t
      4) Same end behavior at the "interpolating" cap:
         - When n >= nmax, switch to not-a-knot layout and return

    Differences:
      1) API style:
         - _generate_knots_impl is a generator that yields trial knot vectors and
           recomputes residuals/fp internally on each iteration.
         - `_add_knots` is a stateful helper that only grows knots; it expects the
           caller to handle residual computation and fp/fpold updates between calls.
      2) Periodicity:
         - _generate_knots_impl supports periodic=True.
         - `_add_knots` is non-periodic only; it uses not-a-knot when n >= nmax.
      3) Residual computation:
         - _generate_knots_impl calls an internal residual routine each iteration.
         - `_add_knots` does not compute residuals; the caller must supply:
             residuals (used by add_knot), fp, fpold.
      4) Return values:
         - _generate_knots_impl yields multiple t's and eventually returns None.
         - `_add_knots` returns:
             * (t_new, nplus) after inserting knots,
             * (not_a_knot_t, nplus) if n >= nmax,
             * (t, nplus) if n >= nest (storage cap).
    """

    acc = s * TOL
    n = t.size
    fpms = fp - s

    # ### c  increase the number of knots. ###

    # c  determine the number of knots nplus we are going to add.
    # https://github.com/scipy/scipy/blob/v1.16.2/scipy/interpolate/fitpack/fpregr.f#L248-L261
    if n == nmin:
        nplus = 1
    else:
        delta = fpold - fp
        npl1 = int(nplus * fpms / delta) if delta > acc else nplus*2
        nplus = min(nplus*2, max(npl1, nplus//2, 1))

    # actually add knots
    # https://github.com/scipy/scipy/blob/v1.16.2/scipy/interpolate/fitpack/fpregr.f#L271-L281
    # https://github.com/scipy/scipy/blob/v1.16.2/scipy/interpolate/fitpack/fpregr.f#L288-L295
    for j in range(nplus):
        t = add_knot(x, t, k, residuals)

        # check if we have enough knots already

        n = t.shape[0]
        # c  if n = nmax, sinf(x) is an interpolating spline.
        # c  if n=nmax we locate the knots as for interpolation.
        # https://github.com/scipy/scipy/blob/v1.16.2/scipy/interpolate/fitpack/fpregr.f#L276-L279
        if n >= nmax:
            # https://github.com/scipy/scipy/blob/v1.16.2/scipy/interpolate/fitpack/fpregr.f#L93-L109
            # https://github.com/scipy/scipy/blob/v1.16.2/scipy/interpolate/fitpack/fpregr.f#L114-L131
            return _not_a_knot(x, k), nplus

        # c  if n=nest we cannot increase the number of knots because of
        # c  the storage capacity limitation.
        # https://github.com/scipy/scipy/blob/v1.16.2/scipy/interpolate/fitpack/fpregr.f#L280
        # https://github.com/scipy/scipy/blob/v1.16.2/scipy/interpolate/fitpack/fpregr.f#L294
        if n >= nest:
            return t, nplus

    return t, nplus


def _regrid_fitpack(
    x, y, Z, *, kx=3, ky=3, s=0.0,
    maxit=50, nestx=None, nesty=None,
    bbox=None):
    """
    Core adaptive bivariate spline fitter using the 1/p-penalty convention.

    Parameters
    ----------
    x, y : array_like
        Strictly increasing coordinate vectors.
    Z : array_like, shape (len(x), len(y))
        Data grid.
    kx, ky : int, optional
        Spline degrees along x and y, default 3 (cubic).
    s : float, optional
        Target residual (`fp` target). `s = 0` requests an interpolatory
        surface; `s > 0` triggers smoothing with penalty weight **1/p**.
    maxit : int, optional
        Maximum iterations for the `p`-search when smoothing, default 50.
    nestx, nesty : int or None
        Max coefficient counts per axis (nesting limits).
    bbox : sequence of 4 scalars
        Optional domain limits `(xb, xe, yb, ye)`. Use `None` entries to skip.

    Returns
    -------
    NdBSpline
        Fitted 2-D spline surface.

    Notes
    -----
    The internal smoothing parameter `p` follows the **inverse**-penalty
    rule: penalty term is 1/p.  Hence, larger `p` -> weaker smoothing
    (approaching interpolation), while smaller `p` -> stronger smoothing.
    A sentinel value `p == -1` is interpreted as *p = inf*, corresponding to
    an exact (interpolatory) fit.

    The iterative process adaptively grows knot vectors based on residual
    energy and optionally performs a 1-D search over `p` to satisfy `fp ~ s`.
    """
    x_fit, y_fit, Z_fit, _, _ = _apply_bbox_grid(x, y, Z, bbox)

    if x_fit.size < (kx + 1) or y_fit.size < (ky + 1):
        raise ValueError(
            f"Not enough samples inside bbox for degrees (kx={kx}, ky={ky}). "
            f"Need at least k+1 per axis: ({kx+1}, {ky+1}). "
            f"Got ({x_fit.size}, {y_fit.size})."
        )

    xb = float(x_fit[0] if bbox[0] is None else bbox[0])
    xe = float(x_fit[-1] if bbox[1] is None else bbox[1])
    yb = float(y_fit[0] if bbox[2] is None else bbox[2])
    ye = float(y_fit[-1] if bbox[3] is None else bbox[3])

    p = -1

    if s == 0.0:
        if nestx is not None or nesty is not None:
            raise ValueError("s == 0 is interpolation only")
        # For special-case k=1 (e.g., Lyche and Morken, Eq.(2.16)),
        # _not_a_knot produces desired knot vector
        tx = _not_a_knot(x_fit, kx)
        ty = _not_a_knot(y_fit, ky)
        (Ax, Ay, Q) = _build_design_matrices(
             x_fit, y_fit, Z, tx, ty, kx, ky)
        C0, fp, _ = _solve_2d_fitpack(Ax, Ay, Q, p,
                                     kx, tx, x_fit, ky, ty,
                                     y_fit, Z_fit)
        return return_NdBSpline(fp, (tx, ty, C0), (kx, ky))

    tx, nestx, nminx, nmaxx = _initialise_knots(x_fit.size, xb, xe, kx, nest=nestx)
    ty, nesty, nminy, nmaxy = _initialise_knots(y_fit.size, yb, ye, ky, nest=nesty)

    fpold = None
    last_axis = "y"
    mpm = len(x) + len(y)
    fp0 = None
    nplusx = None
    nplusy = None

    # https://github.com/scipy/scipy/blob/v1.16.2/scipy/interpolate/fitpack/fpregr.f#L51-L300
    for _ in range(mpm):

        (Ax, Ay, Q) = _build_design_matrices(
             x_fit, y_fit, Z, tx, ty, kx, ky)
        # _solve_2d_fitpack now returns R = Z_fit - zhat alongside C0 and fp,
        # so we can reuse it directly for knot placement instead of
        # recomputing zhat a second time via tocsr + dense matmul.
        C0, fp, R = _solve_2d_fitpack(Ax, Ay, Q, p,
                                      kx, tx, x_fit,
                                      ky, ty, y_fit,
                                      Z_fit)

        # https://github.com/scipy/scipy/blob/v1.16.2/scipy/interpolate/fitpack/fpregr.f#L190
        # https://github.com/scipy/scipy/blob/v1.16.2/scipy/interpolate/fitpack/fpregr.f#L224
        if len(tx) == nminx and len(ty) == nminy:
            fp0 = fp

        if fp < s:
            break

        # https://github.com/scipy/scipy/blob/v1.16.2/scipy/interpolate/fitpack/fpregr.f#L265-L295
        if last_axis == "y":
            tx, nplusx = _add_knots(
                x_fit, kx, s, tx, nmin=nminx, nmax=nmaxx,
                nest=nestx, fp=fp, fpold=fpold,
                residuals=np.sum(R**2, axis=1),
                nplus=nplusx)
            last_axis = "x"
        else:
            ty, nplusy = _add_knots(
                y_fit, ky, s, ty, nmin=nminy, nmax=nmaxy,
                nest=nesty, fp=fp, fpold=fpold,
                residuals=np.sum(R**2, axis=0),
                nplus=nplusy)
            last_axis = "y"

        # When both knot vectors have reached their maximum size no further
        # growth is possible; break to avoid redundant loop iterations.
        if len(tx) >= nmaxx and len(ty) >= nmaxy:
            break

        fpold = fp

    if len(tx) == nminx and len(ty) == nminy:
        return return_NdBSpline(fp, (tx, ty, C0), (kx, ky))

    p = 1
    Drx, offset_dx, nc_dx = disc(tx, kx)
    Dry, offset_dy, nc_dy = disc(ty, ky)
    Drx = PackedMatrix(Drx, offset_dx, nc_dx)
    Dry = PackedMatrix(Dry, offset_dy, nc_dy)
    (Ax, Ay, Q) = _build_design_matrices(
        x_fit, y_fit, Z, tx, ty, kx, ky)
    _, C_sm, fp_sm = _p_search_hit_s(Ax, Drx, Ay, Dry, Q,
                                     kx, tx, x_fit, ky,
                                     ty, y_fit, Z_fit, s,
                                     fp0, maxit=maxit, p_init=p)
    return return_NdBSpline(fp_sm, (tx, ty, C_sm), (kx, ky))


def _regrid(x, y, z, *, bbox=None, kx=3, ky=3, s=0.0, maxit=50):
    """
    Interface for 2-D smoothing B-spline fitting (1/p penalty form).

    Parameters
    ----------
    x, y : array_like
        Strictly increasing 1-D coordinate vectors.
    z : array_like, shape (len(x), len(y))
        Data grid.
    bbox : sequence of 4 scalars
        Optional bounding box ``(xb, xe, yb, ye)``; use ``None`` entries to disable.
    kx, ky : int, optional
        Spline degrees along `x` and `y`, default cubic (3).
    s : float, optional
        Target smoothing residual (`fp` target). Must satisfy ``s >= 0``.
        The underlying formulation uses a **1/p** penalty, meaning:
        - small `p` -> heavy smoothing,
        - large `p` -> light smoothing (approaching interpolation).
        Setting `p == -1` internally denotes *p = inf*, i.e. a pure interpolant.
    maxit : int, optional
        Maximum iterations for `p`-search if invoked.

    Returns
    -------
    NdBSpline
        Fitted bivariate spline surface.
    """
    if bbox is None:
        bbox = [None]*4

    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    z = np.asarray(z, dtype=float)
    bbox = np.ravel(bbox)
    s = float(s)

    if not np.all(np.diff(x) > 0.0):
        raise ValueError("x must be strictly increasing")
    if not np.all(np.diff(y) > 0.0):
        raise ValueError("y must be strictly increasing")
    if x.size != z.shape[0]:
        raise ValueError("x dimension of z must have same number of elements as x")
    if y.size != z.shape[1]:
        raise ValueError("y dimension of z must have same number of elements as y")
    if (s < 0.0):
        raise ValueError("s should be s >= 0.0")
    if not bbox.shape == (4,):
        raise ValueError(f"bbox shape should be (4,), found: {bbox.shape}")

    return _regrid_fitpack(
        x, y, z, kx=kx, ky=ky, s=s, maxit=maxit,
        nestx=None, nesty=None, bbox=bbox)
