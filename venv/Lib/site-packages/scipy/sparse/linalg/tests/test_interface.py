"""Test functions for the sparse.linalg._interface module
"""

import contextlib
from functools import partial
from itertools import product
import operator
from types import ModuleType
from typing import NamedTuple

import pytest
from pytest import raises as assert_raises, warns
from numpy.testing import assert_, assert_equal

import numpy as np
from scipy._lib._array_api import (
    SCIPY_ARRAY_API, SCIPY_DEVICE, is_torch, xp_assert_close, is_lazy_array,
    xp_assert_equal, xp_ravel, is_numpy, make_xp_pytest_marks, is_jax_array,
)
from scipy._external import array_api_extra as xpx
import scipy.sparse as sparse

import scipy.sparse.linalg._interface as interface
from scipy.sparse.linalg import LinearOperator, aslinearoperator
from scipy.sparse._sputils import matrix
from scipy._lib._gcutils import assert_deallocated

pytestmark = make_xp_pytest_marks(
    (LinearOperator, "__init__"),
    (LinearOperator, "matvec"),
    (LinearOperator, "matmat"),
    (LinearOperator, "dot"),
    (LinearOperator, "rdot"),
    (LinearOperator, "rmatvec"),
    (LinearOperator, "rmatmat"),
    aslinearoperator
)

def generate_broadcastable_shapes(nshapes, *, ndim=2, min=0, max=10, rng=None):
    rng = np.random.default_rng(rng)
    min = np.broadcast_to(min, ndim)  # so min and max can be scalars or array-like 
    max = np.broadcast_to(max, ndim)
    batch_shape = tuple(rng.integers(min_, max_+1) for min_, max_ in zip(min, max))
    shapes = np.repeat([batch_shape], nshapes, axis=0)

    # make some elements of some shapes 1 (while preserving overall batch shape)
    for column in shapes.T:
        column[rng.integers(1, nshapes):] = 1
    # permute elements between shapes (while preserving overall batch shape)
    shapes = list(rng.permuted(shapes, axis=0))
    # potentially trim preceeding 1s from a shape
    for i in range(len(shapes)):
        shape = shapes[i]
        j = np.where(shape != 1)[0][0] if np.any(shape != 1) else ndim
        if rng.random() < 0.25:
            shapes[i] = shape[rng.integers(j+1):]
            break

    assert np.broadcast_shapes(*shapes) == batch_shape
    return [tuple(int(el) for el in shape) for shape in shapes]


class TestLinearOperator:
    def setup_method(self):
        self.A = np.array([[1,2,3],
                           [4,5,6]])
        self.B = np.array([[1,2],
                           [3,4],
                           [5,6]])
        self.C = np.array([[1,2],
                           [3,4]])

    def test_matvec(self, xp):
        def get_matvecs(A):
            return [{
                'shape': A.shape,
                'dtype': xp.float64,
                'matvec': lambda x: A @ x,
                'rmatvec': lambda x: xp.conj(A.T) @ x,
                'rmatmat': lambda x: xp.conj(A.T) @ x,
                'matmat': lambda x: A @ x
            }]
        
        _asarray = partial(xp.asarray, dtype=xp.complex128)

        for matvecs in get_matvecs(_asarray(self.A)):
            A = interface.LinearOperator(**matvecs, xp=xp)

            assert A.args == ()

            xp_assert_equal(A.matvec(_asarray([1,2,3])), _asarray([14,32]))
            with pytest.warns(FutureWarning, match="column vectors"):
                xp_assert_equal(
                    A.matvec(_asarray([[1],[2],[3]])),
                    _asarray([[14],[32]]),
                )
            xp_assert_equal(A @ _asarray([1,2,3]), _asarray([14,32]))
            xp_assert_equal(A @ _asarray([[1],[2],[3]]), _asarray([[14],[32]]))
            xp_assert_equal(A.dot(_asarray([1,2,3])), _asarray([14,32]))
            xp_assert_equal(A.dot(_asarray([[1],[2],[3]])), _asarray([[14],[32]]))

            if not SCIPY_ARRAY_API:
                xp_assert_equal(
                    A.matvec(matrix([[1],[2],[3]])),
                    _asarray([[14],[32]]),
                )
                xp_assert_equal(A @ matrix([[1],[2],[3]]), _asarray([[14],[32]]))
                xp_assert_equal(A.dot(matrix([[1],[2],[3]])), _asarray([[14],[32]]))

            xp_assert_equal((2*A) @ _asarray([1,1,1]), _asarray([12,30]))
            xp_assert_equal((2 * A).rmatvec(_asarray([1, 1])), _asarray([10, 14, 18]))
            xp_assert_equal((2*A).H.matvec(_asarray([1,1])), _asarray([10, 14, 18]))
            xp_assert_equal(
                (2*A).adjoint().matvec(_asarray([1,1])),
                _asarray([10, 14, 18]),
            )
            xp_assert_equal((2*A) @ _asarray([[1],[1],[1]]), _asarray([[12],[30]]))
            xp_assert_equal(
                (2 * A).matmat(_asarray([[1], [1], [1]])),
                _asarray([[12], [30]]),
            )
            xp_assert_equal((A*2) @ _asarray([1,1,1]), _asarray([12,30]))
            xp_assert_equal((A*2) @ _asarray([[1],[1],[1]]), _asarray([[12],[30]]))
            xp_assert_equal((2j*A) @ _asarray([1,1,1]), _asarray([12j,30j]))
            xp_assert_equal((A+A) @ _asarray([1,1,1]), _asarray([12, 30]))
            xp_assert_equal((A + A).rmatvec(_asarray([1, 1])), _asarray([10, 14, 18]))
            xp_assert_equal((A+A).H.matvec(_asarray([1,1])), _asarray([10, 14, 18]))
            xp_assert_equal(
                (A + A).adjoint().matvec(_asarray([1,1])),
                _asarray([10, 14, 18]),
            )
            xp_assert_equal((A+A) @ _asarray([[1],[1],[1]]), _asarray([[12], [30]]))
            xp_assert_equal(
                (A+A).matmat(_asarray([[1],[1],[1]])),
                _asarray([[12], [30]])
            )
            xp_assert_equal((-A) @ _asarray([1,1,1]), _asarray([-6,-15]))
            xp_assert_equal((-A) @ _asarray([[1],[1],[1]]), _asarray([[-6],[-15]]))
            xp_assert_equal((A-A) @ _asarray([1,1,1]), _asarray([0,0]))
            xp_assert_equal((A - A) @ _asarray([[1], [1], [1]]), _asarray([[0], [0]]))

            A_ = xp.asarray(self.A, dtype=xp.complex128)
            X = xp.asarray([[1, 2], [3, 4]], dtype=xp.complex128)
            xp_assert_equal((2 * A).rmatmat(X), (2 * A_).T @ X)
            xp_assert_equal((A * 2).rmatmat(X), (A_ * 2).T @ X)
            xp_assert_equal((2j * A).rmatmat(X), xp.conj((2j * A_).T) @ X)
            xp_assert_equal((A * 2j).rmatmat(X), xp.conj((A_ * 2j).T) @ X)
            xp_assert_equal((A + A).rmatmat(X), (A_ + A_).T @ X)
            xp_assert_equal((A + 2j * A).rmatmat(X), xp.conj((A_ + 2j * A_).T) @ X)
            xp_assert_equal((-A).rmatmat(X), (-A_).T @ X)
            xp_assert_equal((A - A).rmatmat(X), (A_ - A_).T @ X)
            xp_assert_equal((2j * A).rmatmat(2j * X), xp.conj((2j * A_).T) @ (2j * X))

            z = A+A
            assert len(z.args) == 2 and z.args[0] is A and z.args[1] is A
            z = 2*A
            assert len(z.args) == 2
            if not is_lazy_array(xp.empty(0)):
                # lazy_xp_function breaks object identity
                assert z.args[0] is A and z.args[1] == 2

            array_object = type(xp.empty(0))
            assert isinstance(A.matvec(_asarray([1, 2, 3])), array_object)
            ctx = (
                pytest.warns(FutureWarning, match="column vectors")
                if not is_jax_array(xp.empty(0))
                else contextlib.nullcontext()
            )
            with ctx:
                assert isinstance(A.matvec(_asarray([[1],[2],[3]])), array_object)
            assert isinstance(A @ _asarray([1,2,3]), array_object)
            assert isinstance(A @ _asarray([[1],[2],[3]]), array_object)
            assert isinstance(A.dot(_asarray([1,2,3])), array_object)
            assert isinstance(A.dot(_asarray([[1],[2],[3]])), array_object)

            if not SCIPY_ARRAY_API:
                assert isinstance(A.matvec(matrix([[1],[2],[3]])), array_object)
                assert isinstance(A @ matrix([[1],[2],[3]]), array_object)
                assert isinstance(A.dot(matrix([[1],[2],[3]])), array_object)

            assert isinstance(2*A, interface._ScaledLinearOperator)
            assert isinstance(2j*A, interface._ScaledLinearOperator)
            assert isinstance(A+A, interface._SumLinearOperator)
            assert isinstance(-A, interface._ScaledLinearOperator)
            assert isinstance(A-A, interface._SumLinearOperator)
            assert isinstance(A/2, interface._ScaledLinearOperator)
            assert isinstance(A/2j, interface._ScaledLinearOperator)

            if not is_lazy_array(xp.empty(0)):
                # lazy_xp_function breaks object identity
                assert ((A * 3) / 3).args[0] is A  # check for simplification

            # Test that prefactor is of _ScaledLinearOperator is not mutated
            # when the operator is multiplied by a number
            result = A @ _asarray([1, 2, 3])
            B = A * 3
            C = A / 5
            xp_assert_equal(A @ _asarray([1, 2, 3]), result)

            assert (2j * A).dtype == xp.complex128

            # Test division by non-scalar
            msg = "Can only divide a linear operator by a scalar."
            with assert_raises(ValueError, match=msg):
                _ = A / xp.asarray([1, 2])

            assert_raises(ValueError, A.matvec, xp.asarray([1,2]))
            assert_raises(ValueError, A.matvec, xp.asarray([1,2,3,4]))
            assert_raises(ValueError, A.matvec, xp.asarray([[1],[2]]))
            assert_raises(ValueError, A.matvec, xp.asarray([[1],[2],[3],[4]]))

            assert_raises(ValueError, lambda: A@A)
            assert_raises(ValueError, lambda: A**2)

        for matvecsA, matvecsB in product(get_matvecs(_asarray(self.A)),
                                          get_matvecs(_asarray(self.B))):
            A = interface.LinearOperator(**matvecsA, xp=xp)
            B = interface.LinearOperator(**matvecsB, xp=xp)
            AtimesB = _asarray(self.A @ self.B)
            X = _asarray([[1, 2], [3, 4]])

            xp_assert_equal((A @ B).rmatmat(X), AtimesB.T @ X)
            xp_assert_equal((2j * A @ B).rmatmat(X), xp.conj((2j * AtimesB).T) @ X)

            xp_assert_equal((A @ B) @ _asarray([1,1]), _asarray([50,113]))
            xp_assert_equal((A @ B) @ _asarray([[1],[1]]), _asarray([[50],[113]]))
            xp_assert_equal((A @ B).matmat(_asarray([[1],[1]])), _asarray([[50],[113]]))

            xp_assert_equal((A @ B).rmatvec(_asarray([1, 1])), _asarray([71, 92]))
            xp_assert_equal((A @ B).H.matvec(_asarray([1, 1])), _asarray([71, 92]))
            xp_assert_equal(
                (A @ B).adjoint().matvec(_asarray([1, 1])),
                _asarray([71, 92]),
            )

            assert isinstance(A@B, interface._ProductLinearOperator)

            assert_raises(ValueError, lambda: A+B)
            assert_raises(ValueError, lambda: A**2)

            z = A@B
            assert len(z.args) == 2
            if not is_lazy_array(xp.empty(0)):
                # lazy_xp_function breaks object identity
                assert z.args[0] is A and z.args[1] is B

        for matvecsC in get_matvecs(_asarray(self.C)):
            C = interface.LinearOperator(**matvecsC, xp=xp)
            X = _asarray([[1, 2], [3, 4]])
            C_ = _asarray(self.C)

            xp_assert_equal(C.rmatmat(X), (C_).T @ X)
            xp_assert_equal((C**2).rmatmat(X), (C_ @ C_).T @ X)

            xp_assert_equal((C**2)@_asarray([1,1]), _asarray([17,37]))
            xp_assert_equal((C**2).rmatvec(_asarray([1, 1])), _asarray([22, 32]))
            xp_assert_equal((C**2).H.matvec(_asarray([1, 1])), _asarray([22, 32]))
            xp_assert_equal(
                (C**2).adjoint().matvec(_asarray([1, 1])),
                _asarray([22, 32]),
            )
            xp_assert_equal((C**2).matmat(_asarray([[1],[1]])), _asarray([[17],[37]]))

            assert isinstance(C**2, interface._PowerLinearOperator)

    @pytest.mark.skip_xp_backends(
        "array_api_strict",
        reason="https://github.com/data-apis/array-api-strict/issues/188",
    )
    def test_matmul(self, xp):
        _asarray = partial(xp.asarray, dtype=xp.complex128)
        A_ = _asarray(self.A)
        D = {'shape': A_.shape,
             'dtype': xp.complex128,
             'matvec': lambda x: xp.reshape(A_ @ x, (A_.shape[0],)),
             'rmatvec': lambda x: xp.reshape(xp.conj(A_.T) @ x, (A_.shape[1],)),
             'rmatmat': lambda x: xp.conj(A_.T) @ x,
             'matmat': lambda x: A_ @ x}
        A = interface.LinearOperator(**D, xp=xp)
        B = _asarray([[1 + 1j, 2, 3],
                      [4, 5, 6],
                      [7, 8, 9]])
        b = B[0, ...]

        xp_assert_equal(operator.matmul(A, b), A * b)
        xp_assert_equal(
            operator.matmul(A, xp.reshape(b, (-1, 1))),
            A * xp.reshape(b, (-1, 1)),
        )
        xp_assert_equal(operator.matmul(A, B), A @ B)
        xp_assert_equal(operator.matmul(b, A.H), b * A.H)
        xp_assert_equal(operator.matmul(b, A.adjoint()), b * A.adjoint())
        xp_assert_equal(
            operator.matmul(xp.reshape(b, (1, -1)), A.H),
            xp.reshape(b, (1, -1)) * A.H,
        )
        xp_assert_equal(
            operator.matmul(b.reshape(1, -1), A.adjoint()),
            b.reshape(1, -1) * A.adjoint(),
        )
        xp_assert_equal(operator.matmul(B, A.H), B @ A.H)
        xp_assert_equal(operator.matmul(B, A.adjoint()), B @ A.adjoint())
        assert_raises(ValueError, operator.matmul, A, 2)
        assert_raises(ValueError, operator.matmul, 2, A)

    def test_dimension_mismatch(self, xp):
        msg = "(D|d)imension mismatch"
        A = interface.aslinearoperator(xp.ones((3, 4)))
        # matvec
        A.matvec(xp.ones(4))
        with pytest.raises(ValueError, match=msg):
            A.matvec(xp.ones(3))
        # rmatvec
        A.rmatvec(xp.ones(3))
        with pytest.raises(ValueError, match=msg):
            A.rmatvec(xp.ones(4))
        # matmat
        A.matmat(xp.ones((4, 3)))
        with pytest.raises(ValueError, match=msg):
            A.matmat(xp.ones((3, 3)))
        # rmatmat
        A.rmatmat(xp.ones((3, 3)))
        with pytest.raises(ValueError, match=msg):
            A.rmatmat(xp.ones((4, 4)))
        # dot
        A.dot(xp.ones((4, 3)))
        with pytest.raises(ValueError, match=msg):
            A.dot(xp.ones((3, 3)))
        # rdot
        A.rdot(xp.ones((3, 3)))
        with pytest.raises(ValueError, match=msg):
            A.rdot(xp.ones((4, 4)))
        

@pytest.mark.skip_xp_backends("dask.array", reason="https://github.com/dask/dask/issues/11711")
class TestDotTests:
    """
    This class aims to help ensure correctness of the LinearOperator
    interface, by verifying correctness properties based on equivalent
    computations using 'forward' and 'adjoint' modes.
    """
    class OperatorArgs(NamedTuple):
        """
        shape: (core) shape of the operator
        op_dtype: dtype of the operator
        data_dtype: real dtype corresponding to op_dtype for data generation
        complex: the operator has a complex dtype
        """
        shape: tuple[int, ...]
        op_dtype: str
        data_dtype: str
        complex: bool

    real_square_args: OperatorArgs = OperatorArgs(
        (12, 12), "float64", "float64", False
    )
    integer_square_args: OperatorArgs = OperatorArgs(
        (9, 9), "int32", "float32", False
    )
    complex_square_args: OperatorArgs = OperatorArgs(
        (13, 13), "complex64", "float32", True
    )
    real_overdetermined_args: OperatorArgs = OperatorArgs(
        (17, 11), "float64", "float64", False
    )
    complex_overdetermined_args: OperatorArgs = OperatorArgs(
        (17, 11), "complex128", "float64", True
    )
    real_underdetermined_args: OperatorArgs = OperatorArgs(
        (5, 9), "float64", "float64", False
    )

    square_args_list: list[OperatorArgs] = [
        real_square_args, integer_square_args, complex_square_args
    ]
    all_args_list: list[OperatorArgs]  = square_args_list + [
        real_overdetermined_args, complex_overdetermined_args, real_underdetermined_args
    ]

    def check_matvec(
        self, xp: ModuleType, op: interface.LinearOperator, data_dtype: str,
        complex_data: bool = False,
    ):
        """
        This check verifies the equivalence of the forward and adjoint computation,
        using `matvec` and `rmatvec` respectively, on randomised data.

        Data is generated with the real dtype `data_dtype` and operated on by the
        linear operator `op`.

        If `complex_data` is set to `True`, complex data is instead generated
        by combining randomised real and imaginary components, each generated
        with `data_dtype`.

        If `check_operators` is set to `True`, equivalence is checked between
        `matvec` and `*` and `@`,
        and between `rmatvec` and the composition of `.H` and `*`.

        If `check_dot` is set to `True`, equivalence is checked between
        `matvec` and `.dot`,
        and between `rmatvec` and the composition of `.H` and `.dot`.
        """
        rng = np.random.default_rng(42)

        dtype = np.dtype(data_dtype)
        *batch_shape, M, N = op.shape

        u_broadcast, v_broadcast = generate_broadcastable_shapes(
            2, ndim=3, min=0, max=5, rng=63
        )

        for u_shape, v_shape in (
            # test broadcasting for unbatched RHS
            ((N,), (M,)),
            # as well as with the batch shape of `op` + 3-D broadcast dims
            ((*u_broadcast, *batch_shape, N), (*v_broadcast, *batch_shape, M))
        ):
            u = rng.standard_normal(u_shape, dtype=dtype)
            v = rng.standard_normal(v_shape, dtype=dtype)
            if complex_data:
                u = u + (1j * rng.standard_normal(u_shape, dtype=dtype))
                v = v + (1j * rng.standard_normal(v_shape, dtype=dtype))
            u = xp.asarray(u)
            v = xp.asarray(v)
    
            op_u = op.matvec(u)
            opH_v = op.rmatvec(v)
    
            op_u_H_v = xp.vecdot(op_u, v, axis=-1)
            uH_opH_v = xp.vecdot(u, opH_v, axis=-1)
    
            rtol = 1e-12 if np.finfo(data_dtype).eps < 1e-8 else 1e-5
            atol = 2e-15 if np.finfo(data_dtype).eps < 1e-8 else 1e-5
            xp_assert_close(op_u_H_v, uH_opH_v, rtol=rtol, atol=atol)

    def check_matmat(
        self, xp, op: interface.LinearOperator, data_dtype: str,
        complex_data: bool = False, check_operators: bool = False,
        check_dot: bool = False,
    ):
        """
        This check verifies the equivalence of the forward and adjoint computation,
        using `matmat` and `rmatmat` respectively, on randomised data.

        Data is generated with the real dtype `data_dtype` and operated on by the
        linear operator `op`.

        If `complex_data` is set to `True`, complex data is instead generated
        by combining randomised real and imaginary components, each generated
        with `data_dtype`.

        If `check_operators` is set to `True`, equivalence is checked between
        `matmat` and `*` and `@`,
        and between `rmatvec` and the composition of `.H` and `@`.

        If `check_dot` is set to `True`, equivalence is checked between
        `matmat` and `.dot`,
        and between `rmatmat` and the composition of `.H` and `.dot`.
        """
        rng = np.random.default_rng(42)
        k = rng.integers(2, 100)

        dtype = np.dtype(data_dtype)
        *batch_shape, M, N = op.shape

        # Test `U` and `V` with the same batch shape as `op`
        U = rng.standard_normal(size=(*batch_shape, N, k), dtype=dtype)
        V = rng.standard_normal(size=(*batch_shape, M, k), dtype=dtype)
        if complex_data:
            U = U + (1j * rng.standard_normal(size=(*batch_shape, N, k), dtype=dtype))
            V = V + (1j * rng.standard_normal(size=(*batch_shape, M, k), dtype=dtype))
        U = xp.asarray(U)
        V = xp.asarray(V)

        op_U = op.matmat(U)
        opH_V = op.rmatmat(V)

        if check_operators:
            xp_assert_close(op_U, op * U)
            xp_assert_close(op_U, op @ U)
            xp_assert_close(opH_V, op.H * V)
            xp_assert_close(opH_V, op.H @ V)

        if check_dot:
            xp_assert_close(op_U, op.dot(U))
            xp_assert_close(opH_V, op.H.dot(V))

        op_U_H = xp.conj(op_U).mT
        UH = xp.conj(U).mT

        op_U_H_V = xp.matmul(op_U_H, V)
        UH_opH_V = xp.matmul(UH, opH_V)

        rtol = 1e-12 if np.finfo(data_dtype).eps < 1e-8 else 1e-5
        atol = 2e-15 if np.finfo(data_dtype).eps < 1e-8 else 1e-5
        xp_assert_close(op_U_H_V, UH_opH_V, rtol=rtol, atol=atol)

    @pytest.mark.parametrize("batch_shape", [(), (3,), (3, 4, 5,), (0,)])
    @pytest.mark.parametrize("args", square_args_list)
    def test_identity_square(
        self, args: OperatorArgs, batch_shape: tuple[int, ...], xp
    ):
        """
        Simple identity operator on square matrices.
        Tests batches of RHS via `args.batch_shape`.
        """
        def identity(x):
            shape = (*xpx.broadcast_shapes(batch_shape, x.shape[:-1]), x.shape[-1])
            return xp.broadcast_to(x, shape)

        shape = batch_shape + args.shape
        dtype = getattr(xp, args.op_dtype)
        op = interface.LinearOperator(  # type:ignore[call-arg]
            shape=shape, dtype=dtype,
            matvec=identity, rmatvec=identity, xp=xp,
        )

        self.check_matvec(xp, op, data_dtype=args.data_dtype, complex_data=args.complex)
        self.check_matmat(xp, op, data_dtype=args.data_dtype, complex_data=args.complex)

    @pytest.mark.parametrize("batch_shape", [(), (3,), (3, 4, 5,), (0,)])
    @pytest.mark.parametrize("args", all_args_list)
    def test_identity_nonsquare(
        self, args: OperatorArgs, batch_shape: tuple[int, ...], xp
    ):
        """
        Identity operator with zero-padding on non-square matrices.
        Tests batches of RHS via `args.batch_shape`.
        """
        M, N = args.shape

        def pad_core_dim(x, target_size):
            pad_widths = [(0, 0)] * x.ndim # no padding for batch dims
            pad_widths[-1] = (0, target_size - x.shape[-1]) # padding for core dim
            return xpx.pad(x, pad_widths, mode="constant")

        def mv(x):
            match np.sign(x.shape[-1] - M):
                case 0:  # square
                    pass
                case 1:  # crop x to size
                    x = x[..., :M]
                case -1:  # pad with zeros
                    x = pad_core_dim(x, target_size=M)

            shape = (*xpx.broadcast_shapes(batch_shape, x.shape[:-1]), x.shape[-1])
            return xp.broadcast_to(x, shape)

        def rmv(x):            
            match np.sign(N - x.shape[-1]):
                case 0:  # square
                    pass
                case 1:  # pad with zeros
                    x = pad_core_dim(x, target_size=N)
                case -1:  # crop x to size
                    x = x[..., :N]

            shape = (*xpx.broadcast_shapes(batch_shape, x.shape[:-1]), x.shape[-1])
            return xp.broadcast_to(x, shape)

        shape = batch_shape + args.shape
        dtype = getattr(xp, args.op_dtype)
        op = interface.LinearOperator(  # type:ignore[call-arg]
            shape=shape, dtype=dtype, matvec=mv, rmatvec=rmv, xp=xp
        )
        
        self.check_matvec(xp, op, data_dtype=args.data_dtype, complex_data=args.complex)
        self.check_matmat(xp, op, data_dtype=args.data_dtype, complex_data=args.complex)

    @pytest.mark.parametrize("batch_shape", [(), (3,), (3, 4, 5,), (0,)])
    @pytest.mark.parametrize("args", square_args_list)
    def test_scaling_square(self, args: OperatorArgs, batch_shape: tuple[int, ...], xp):
        """
        Simple (complex) scaling operator on square matrices.
        Tests batches of RHS via `args.batch_shape`.
        """
        def scale(x):
            scale = (3 + 2j) / abs(3 + 2j)
            shape = (*xpx.broadcast_shapes(batch_shape, x.shape[:-1]), x.shape[-1])
            return xp.broadcast_to(scale * x, shape)

        def r_scale(x):
            scale = (3 - 2j) / abs(3 - 2j)
            shape = (*xpx.broadcast_shapes(batch_shape, x.shape[:-1]), x.shape[-1])
            return xp.broadcast_to(scale * x, shape)

        shape = batch_shape + args.shape
        dtype = getattr(xp, args.op_dtype)
        op = interface.LinearOperator(  # type:ignore[call-arg]
            shape=shape, dtype=dtype, matvec=scale, rmatvec=r_scale, xp=xp
        )
        self.check_matvec(xp, op, data_dtype=args.data_dtype, complex_data=args.complex)
        self.check_matmat(xp, op, data_dtype=args.data_dtype, complex_data=args.complex)

    @pytest.mark.parametrize("batch_shape", [(), (3,), (3, 4, 5,), (0,)])
    def test_subclass_matmat(self, batch_shape: tuple[int, ...], xp):
        """
        Simple rotation operator defined by `matmat` and `adjoint`,
        subclassing `LinearOperator`.
        Tests batches of RHS via `batch_shape`.
        """
        class RotOp(interface.LinearOperator):

            def __init__(self, dtype, shape, theta):
                self._theta = theta
                super().__init__(dtype, shape, xp=xp)

            def _matmat(self, X):
                theta = self._theta
                R = xp.asarray([
                    [np.cos(theta), -np.sin(theta)],
                    [np.sin(theta),  np.cos(theta)]
                ])
                B = R @ X
                shape = xpx.broadcast_shapes(batch_shape, X.shape[:-2]) + X.shape[-2:]
                return xp.broadcast_to(B, shape)

            def _adjoint(self):
                negative_theta = -self._theta
                return RotOp(self.dtype, self.shape, negative_theta)

        theta = xp.pi / 2
        data_dtype = "float64"
        op = RotOp(shape=(*batch_shape, 2, 2), dtype=xp.float64, theta=theta)

        self.check_matvec(xp, op, data_dtype=data_dtype, complex_data=False)
        self.check_matmat(xp, op, data_dtype=data_dtype, complex_data=False)

    @pytest.mark.parametrize("batch_shape", [(), (3,), (3, 4, 5,), (0,)])
    @pytest.mark.skip_xp_backends("dask.array", reason="https://github.com/dask/dask/issues/11711")
    def test_aslinearop_dense(
        self, batch_shape: tuple[int, ...], xp
    ):
        """
        Test operators coming from `aslinearoperator`,
        *including batched LHS*.
        """
        rng = np.random.default_rng(42)
        matrix = xp.asarray(
            rng.standard_normal((*batch_shape, 4, 4)),
            dtype=xp.complex128,
        )
        op = interface.aslinearoperator(matrix)
        data_dtype = "float64"
        self.check_matvec(xp, op, data_dtype=data_dtype, complex_data=True)
        self.check_matmat(xp, op, data_dtype=data_dtype, complex_data=True)

    @pytest.mark.parametrize("batch_shape", [(), (3,), (3, 4, 5,), (0,)])
    @pytest.mark.skip_xp_backends(np_only=True)
    def test_aslinearop_sparse(self, batch_shape: tuple[int, ...], xp):
        """
        Test operators coming from `aslinearoperator` with sparse arrays,
        *including batched LHS*
        """
        matrix = sparse.random_array((*batch_shape, 4, 4), dtype=xp.complex128)
        op = interface.aslinearoperator(matrix)
        data_dtype = "float64"
        self.check_matvec(xp, op, data_dtype=data_dtype, complex_data=True)
        self.check_matmat(xp, op, data_dtype=data_dtype, complex_data=True)


class TestAsLinearOperator:
    def make_cases(self, original, dtype, xp=None):
        cases = []

        # Test default implementations of _adjoint and _rmatvec, which
        # refer to each other.
        _xp = xp if xp is not None else np
        def mv(x):
            y = original @ x
            if len(x.shape) == 2:
                y = _xp.reshape(y, (-1, 1))
            return y

        def rmv(x):
            return _xp.conj(original.T) @ x

        class BaseMatlike(interface.LinearOperator):
            args = ()

            def __init__(self, dtype):
                super().__init__(
                    dtype=dtype,
                    shape=original.shape,
                    xp=xp,
                )

            def _matvec(self, x):
                return mv(x)

        class HasRmatvec(BaseMatlike):
            args = ()

            def _rmatvec(self,x):
                return rmv(x)

        class HasAdjoint(BaseMatlike):
            args = ()

            def _adjoint(self):
                shape = self.shape[1], self.shape[0]
                matvec = partial(rmv)
                rmatvec = partial(mv)
                return interface.LinearOperator(matvec=matvec,
                                                rmatvec=rmatvec,
                                                dtype=dtype,
                                                shape=shape,
                                                xp=xp,)

        class HasRmatmat(HasRmatvec):
            def _matmat(self, x):
                return original @ x

            def _rmatmat(self, x):
                return _xp.conj(original.T) @ x

        if xp:
            cases.append((original, original))

            cases.append((HasRmatvec(dtype), original))
            cases.append((HasAdjoint(dtype), original))
            cases.append((HasRmatmat(dtype), original))
            
            cases.append((interface.aslinearoperator(original.T).T, original))
            cases.append((
                interface.aslinearoperator(original.T).H,
                xp.conj(original),
            ))
            cases.append((
                interface.aslinearoperator(original.T).adjoint(),
                xp.conj(original)),
            )

            return cases
        
        cases.append((matrix(original, dtype=dtype), original))
        cases.append((np.array(original, dtype=dtype), original))
        cases.append((sparse.csr_array(original, dtype=dtype), original))

        cases.append((HasRmatvec(dtype), original))
        cases.append((HasAdjoint(dtype), original))
        cases.append((HasRmatmat(dtype), original))
        return cases
    
    def setup_method(self):
        self.cases = []
        make_cases = self.make_cases
        original = np.array([[1,2,3], [4,5,6]])
        self.cases += make_cases(original, np.int32)
        self.cases += make_cases(original, np.float32)
        self.cases += make_cases(original, np.float64)
        self.cases += [(interface.aslinearoperator(M).T, A.T)
                       for M, A in make_cases(original.T, np.float64)]
        self.cases += [(interface.aslinearoperator(M).H, A.T.conj())
                       for M, A in make_cases(original.T, np.float64)]
        self.cases += [(interface.aslinearoperator(M).adjoint(), A.T.conj())
                       for M, A in make_cases(original.T, np.float64)]

        original = np.array([[1, 2j, 3j], [4j, 5j, 6]])
        self.cases += make_cases(original, np.complex128)
        self.cases += [(interface.aslinearoperator(M).T, A.T)
                       for M, A in make_cases(original.T, np.complex128)]
        self.cases += [(interface.aslinearoperator(M).H, A.T.conj())
                       for M, A in make_cases(original.T, np.complex128)]
        self.cases += [(interface.aslinearoperator(M).adjoint(), A.T.conj())
                       for M, A in make_cases(original.T, np.complex128)]

    @pytest.mark.skip_xp_backends(np_only=True)
    def test_basic(self, xp):

        for M, A_array in self.cases:
            A = interface.aslinearoperator(M)

            xs = [np.array([1, 2, 3]),
                  np.array([[1], [2], [3]])]
            ys = [np.array([1, 2]), np.array([[1], [2]])]

            if A.dtype == np.complex128:
                xs += [np.array([1, 2j, 3j]),
                       np.array([[1], [2j], [3j]])]
                ys += [np.array([1, 2j]), np.array([[1], [2j]])]

            x2 = np.array([[1, 4], [2, 5], [3, 6]])

            for x in xs:
                ctx = (
                    pytest.warns(FutureWarning, match="column vectors")
                    if x.shape[-1] == 1
                    else contextlib.nullcontext()
                )
                with ctx:
                    assert_equal(A.matvec(x), A_array.dot(x))

                assert_equal(A @ x, A_array.dot(x))

            assert_equal(A.matmat(x2), A_array.dot(x2))
            assert_equal(A @ x2, A_array.dot(x2))

            for y in ys:
                ctx = (
                    pytest.warns(FutureWarning, match="column vectors")
                    if y.shape[-1] == 1
                    else contextlib.nullcontext()
                )
                with ctx:
                    assert_equal(A.rmatvec(y), A_array.T.conj().dot(y))
                    assert_equal(A.T.matvec(y), A_array.T.dot(y))
                    assert_equal(A.H.matvec(y), A_array.T.conj().dot(y))
                    assert_equal(A.adjoint().matvec(y), A_array.T.conj().dot(y))

            for y in ys:
                if y.ndim < 2:
                    continue
                assert_equal(A.rmatmat(y), A_array.T.conj().dot(y))
                assert_equal(A.T.matmat(y), A_array.T.dot(y))
                assert_equal(A.H.matmat(y), A_array.T.conj().dot(y))
                assert_equal(A.adjoint().matmat(y), A_array.T.conj().dot(y))

            if hasattr(M,'dtype'):
                assert_equal(A.dtype, M.dtype)

            assert_(hasattr(A, 'args'))

    @pytest.mark.skip_xp_backends(np_only=True)
    def test_dot(self, xp):

        for M, A_array in self.cases:
            A = interface.aslinearoperator(M)

            x0 = np.array([1, 2, 3])
            x1 = np.array([[1], [2], [3]])
            x2 = np.array([[1, 4], [2, 5], [3, 6]])

            assert_equal(A.dot(x0), A_array.dot(x0))
            assert_equal(A.dot(x1), A_array.dot(x1))
            assert_equal(A.dot(x2), A_array.dot(x2))

    @pytest.mark.parametrize("dtype", ["int64", "float64", "complex128"])
    def test_xp(self, dtype, xp):
        if dtype == "int64" and is_torch(xp) and SCIPY_DEVICE != "cpu":
            pytest.skip("\"addmm_cuda\" not implemented for 'Long'")
        dtype = getattr(xp, dtype)
        original = xp.asarray([[1, 2, 3], [4, 5, 6]], dtype=dtype)
        for M, A_array in self.make_cases(original, dtype, xp=xp):
            A = interface.aslinearoperator(M)

            xs = [xp.asarray([1, 2, 3]),
                  xp.asarray([[1], [2], [3]])]
            ys = [xp.asarray([1, 2]), xp.asarray([[1], [2]])]

            if A.dtype == xp.complex128:
                xs += [xp.asarray([1, 2j, 3j]),
                       xp.asarray([[1], [2j], [3j]])]
                ys += [xp.asarray([1, 2j]), xp.asarray([[1], [2j]])]

            x2 = xp.asarray([[1, 4], [2, 5], [3, 6]], dtype=dtype)

            for x in xs:
                x = xp.astype(x, dtype)
                xp_assert_equal(A @ x, A_array @ x)

            xp_assert_equal(A.matmat(x2), A_array @ x2)
            xp_assert_equal(A @ x2, A_array @ x2)

            for y in ys:
                if y.ndim < 2:
                    continue
                y = xp.astype(y, dtype)
                xp_assert_equal(A.rmatmat(y), xp.conj(A_array.T) @ y)
                xp_assert_equal(A.T.matmat(y), xp.conj(A_array.T) @ y)
                xp_assert_equal(A.H.matmat(y), xp.conj(A_array.T) @ y)
                xp_assert_equal(A.adjoint().matmat(y), xp.conj(A_array.T) @ y)

            if hasattr(M, 'dtype'):
                assert_equal(A.dtype, M.dtype)

            assert hasattr(A, 'args')
            

def test_repr(xp):
    A = interface.LinearOperator(shape=(1, 1), matvec=lambda x: xp.asarray([1]), xp=xp)
    repr_A = repr(A)
    assert 'unspecified dtype' not in repr_A, repr_A


def test_identity(xp):
    ident = interface.IdentityOperator((3, 3), xp=xp)
    xp_assert_equal(ident @ xp.asarray([1, 2, 3]), xp.asarray([1, 2, 3]))
    xp_assert_equal(xp_ravel(ident.dot(xp.reshape(xp.arange(9), (3, 3)))), xp.arange(9))

    assert_raises(ValueError, ident.matvec, xp.asarray([1, 2, 3, 4]))


def test_attributes(xp):
    A = interface.aslinearoperator(xp.reshape(xp.arange(16, dtype=xp.float64), (4, 4)))

    def always_four_ones(x):
        x = xp.asarray(x)
        assert x.shape == (3,) or x.shape == (3, 1)
        return xp.ones(4)

    B = interface.LinearOperator(shape=(4, 3), matvec=always_four_ones, xp=xp)

    ops = [A, B, A * B, A @ B, A.H, A.adjoint(), A + A, B + B, A**4]
    for op in ops:
        assert hasattr(op, "dtype")
        assert hasattr(op, "shape")
        assert hasattr(op, "_matvec")


def matvec_for_pickle(x):
    """ Needed for test_pickle as local functions are not pickleable """
    return x


@pytest.mark.skip_xp_backends(
    "array_api_strict",
    reason="pickle-ability is not guaranteed by the standard"
)
def test_pickle(xp):
    import pickle

    protocol_min = 0 if is_numpy(xp) else 2
    for protocol in range(protocol_min, pickle.HIGHEST_PROTOCOL + 1):
        A = interface.LinearOperator((3, 3), matvec_for_pickle, xp=xp)
        s = pickle.dumps(A, protocol=protocol)
        B = pickle.loads(s)

        for k in A.__dict__:
            assert getattr(A, k) == getattr(B, k)


def test_inheritance(xp):
    class Empty(interface.LinearOperator):
        pass

    with warns(RuntimeWarning, match="should implement at least"):
        assert_raises(TypeError, Empty)

    class Identity(interface.LinearOperator):
        def __init__(self, n):
            super().__init__(dtype=None, shape=(n, n), xp=xp)

        def _matvec(self, x):
            return x

    id3 = Identity(3)
    xp_assert_equal(id3.matvec(xp.asarray([1, 2, 3])), xp.asarray([1, 2, 3]))
    assert_raises(NotImplementedError, id3.rmatvec, xp.asarray([4, 5, 6]))

    class MatmatOnly(interface.LinearOperator):
        def __init__(self, A):
            super().__init__(A.dtype, A.shape, xp=xp)
            self.A = A

        def _matmat(self, x):
            return self.A @ x

    mm = MatmatOnly(xp.asarray(np.random.randn(5, 3)))
    assert mm.matvec(xp.asarray(np.random.randn(3))).shape == (5,)


def test_dtypes_of_operator_sum(xp):
    # gh-6078

    mat_complex = xp.asarray(np.random.rand(2,2) + 1j * np.random.rand(2,2))
    mat_real = xp.asarray(np.random.rand(2,2))

    complex_operator = interface.aslinearoperator(mat_complex)
    real_operator = interface.aslinearoperator(mat_real)

    sum_complex = complex_operator + complex_operator
    sum_real = real_operator + real_operator

    assert sum_real.dtype == xp.float64
    assert sum_complex.dtype == xp.complex128


def test_no_double_init(xp):
    call_count = [0]

    def matvec(v):
        call_count[0] += 1
        return v

    # It should call matvec exactly once (in order to determine the
    # operator dtype)
    interface.LinearOperator((2, 2), matvec=matvec, xp=xp)
    assert_equal(call_count[0], 1)


INT_DTYPES = (np.int8, np.int16, np.int32, np.int64)
REAL_DTYPES = (np.float32, np.float64, np.longdouble)
COMPLEX_DTYPES = (np.complex64, np.complex128, np.clongdouble)
INEXACTDTYPES = REAL_DTYPES + COMPLEX_DTYPES
ALLDTYPES = INT_DTYPES + INEXACTDTYPES


@pytest.mark.parametrize("test_dtype", ALLDTYPES)
def test_determine_lo_dtype_from_matvec(test_dtype, xp):
    if "longdouble" in test_dtype.__name__ and not is_numpy(xp):
        pytest.skip("longdoubles are only tested for `np`")
    # gh-19209
    scalar = xp.asarray(np.array(1, dtype=test_dtype))
    def mv(v):
        return xp.stack([scalar * v[0], v[1]])

    lo = interface.LinearOperator((2, 2), matvec=mv, xp=xp)
    # expected dtype depends on if mixed exact-inexact promotion is defined
    # since dtype determination follows the following procedure:
    # - take the dtype from calling `matvec` on an `int8`
    # - unless that fails (e.g. via overflow or no mixed exact-inexact promotion),
    #   in which case use the default integral dtype
    expected = scalar.dtype
    if xp.isdtype(expected, ("real floating", "complex floating")):
        try:
            xp.asarray(2) + xp.asarray(2.0)
        except TypeError:
            expected = xpx.default_dtype(xp, kind="integral")
    assert lo.dtype == expected


def test_determine_lo_dtype_for_int(xp):
    # gh-19209
    # test Python int larger than int8 max cast to some int
    def mv(v):
        return xp.asarray([128 * v[0], v[1]])

    lo = interface.LinearOperator((2, 2), matvec=mv, xp=xp)
    assert xp.isdtype(lo.dtype, "integral")


def test_adjoint_conjugate(xp):
    X = xp.asarray([[1j]], dtype=xp.complex128)
    A = interface.aslinearoperator(X)

    B = 1j * A
    Y = 1j * X

    v = xp.asarray([1], dtype=xp.complex128)

    xp_assert_equal(B.dot(v), xp.vecdot(Y, v))
    xp_assert_equal(B.H.dot(v), xp.vecdot(xp.conj(Y.T), v))
    xp_assert_equal(B.adjoint().dot(v), xp.vecdot(xp.conj(Y.T), v))


def test_ndim(xp):
    X = xp.asarray([[1]])
    A = interface.aslinearoperator(X)
    assert A.ndim == 2


def test_transpose_noconjugate(xp):
    X = xp.asarray([[1j]], dtype=xp.complex128)
    A = interface.aslinearoperator(X)

    B = 1j * A
    Y = 1j * X

    v = xp.asarray([1], dtype=xp.complex128)

    xp_assert_equal(B.dot(v), xp.vecdot(Y, v))
    xp_assert_equal(B.T.dot(v), xp.vecdot(Y.T, v))


@pytest.mark.skip_xp_backends("dask.array", reason="https://github.com/dask/dask/issues/11711")
@pytest.mark.skip_xp_backends(
    "array_api_strict",
    reason="https://github.com/data-apis/array-api-strict/issues/188"
)
def test_transpose_multiplication(xp):
    _asarray = partial(xp.asarray, dtype=xp.complex128)
    class MyMatrix(interface.LinearOperator):
        def __init__(self, A):
            super().__init__(A.dtype, A.shape, xp=xp)
            self.A = A
        def _matmat(self, other): return self.A @ other
        def _rmatmat(self, other): return self.A.mT @ other

    A = MyMatrix(_asarray([[1, 2], [3, 4]]))
    X = _asarray([1, 2])
    X_T = X
    B = _asarray([[10, 20], [30, 40]])
    X2 = xp.reshape(X, (-1, 1))
    Y = _asarray([[1, 2], [3, 4]])

    xp_assert_equal(A @ B, Y @ B)
    xp_assert_equal(B.T @ A, B.T @ Y)
    xp_assert_equal(A.T @ B, Y.mT @ B)
    xp_assert_equal(A @ X, Y @ X)
    xp_assert_equal(X_T @ A, X_T @ Y)
    xp_assert_equal(A.T @ X, Y.mT @ X)
    xp_assert_equal(A @ X2, Y @ X2)
    xp_assert_equal(X2.mT @ A, X2.mT @ Y)
    xp_assert_equal(A.T @ X2, Y.mT @ X2)


@pytest.mark.skip_xp_backends(np_only=True)
def test_sparse_matmat_exception():
    B = sparse.eye_array(2)
    # well defined matmat via `aslinearoperator`
    A = interface.aslinearoperator(sparse.eye_array(2))
    assert isinstance(A @ B, sparse.sparray)
    assert isinstance(B @ A, sparse.sparray)
    xp_assert_equal((A @ B).toarray(), np.eye(2))
    xp_assert_equal((B @ A).toarray(), np.eye(2))
    # ill-defined matmat via default fallback to matvec
    A = interface.LinearOperator((2, 2), matvec=lambda x: x)
    msg = "Try wrapping the matrix with `aslinearoperator` first."
    with assert_raises(TypeError, match=msg):
        A @ B
    with assert_raises(TypeError, match=msg):
        B @ A
    # after using `aslinearoperator`
    B = interface.aslinearoperator(B)
    assert isinstance(A @ B, interface.LinearOperator)
    assert isinstance(B @ A, interface.LinearOperator)
    xp_assert_equal((A @ B).matvec(np.ones(2)), np.ones(2))


@pytest.mark.skip_xp_backends(
    "jax.numpy", reason="remaining reference (due to lazy_xp_function?)"
)
def test_MatrixLinearOperator_refcycle(xp):
    # gh-10634
    # Test that MatrixLinearOperator can be automatically garbage collected
    A = xp.eye(2)
    with assert_deallocated(interface.MatrixLinearOperator, A, xp) as op:
        op.adjoint()
        del op


@pytest.mark.parametrize("left", [False, True])
@pytest.mark.parametrize("operator_definition", ["subclass_matvec", "subclass_matmat",
                                                 "__init__matvec", "__init__matmat",
                                                 "aslinearoperator"])
@pytest.mark.parametrize("batch_A", [(), (5,), (0,)])
@pytest.mark.parametrize("batch_x", [(), (6, 1), (0, 1)])
@pytest.mark.parametrize("dtype", [np.float64, np.complex128])
@pytest.mark.skip_xp_backends("dask.array", reason="https://github.com/dask/dask/issues/11711")
@pytest.mark.skip_xp_backends(
    "array_api_strict",
    reason="https://github.com/data-apis/array-api-strict/issues/188"
)
def test_batch(left, operator_definition, batch_A, batch_x, dtype, xp):
    # TODO ideas:
    # - test lower-precision types
    # - test `transpose`, `adjoint`, `__mul__`, etc.
    # - test composite LinearOperators
    rng = np.random.default_rng(41981392342349823)

    m, n, k = 4, 3, 2
    A_ = rng.random((*batch_A, m, n))
    x_row = rng.random((*batch_x, n if left else m))
    x_col = rng.random((*batch_x, n if left else m, 1))
    x_mat = rng.random((*batch_x, n if left else m, k))

    if dtype == np.complex128:
        A_ = A_ + 1j * rng.random(A_.shape)
        x_row = x_row + 1j * rng.random(x_row.shape)
        x_col = x_col + 1j * rng.random(x_col.shape)
        x_mat = x_mat + 1j * rng.random(x_mat.shape)
    
    A_, x_row, x_col, x_mat = (xp.asarray(x) for x in (A_, x_row, x_col, x_mat))

    def matvec(A, x):
        if not batch_x:
            assert x.ndim == 1
            return A @ x
        else:
            # It might make it easier on the author of LinearOperators to
            # "ravel" all batch dimensions into one before calling their `matvec`
            # implementation. That way, they only have to think about vectorizing
            # w.r.t. one batch dimension rather than an arbitrary number of batch
            # dimensions. In this case, `x.ndim` would be exactly 2.
            assert x.ndim >= 2
            return (A @ (x[..., xp.newaxis]))[..., 0]

    def matmat(A, X):
        if not batch_x:
            assert X.ndim == 2
        else:
            # Similar to above, we could ravel all batch dimensions before calling
            # the user's `matmat`. Then `x.ndim` would be exactly 3.
            assert X.ndim >= 3
        return A @ X

    if operator_definition == "aslinearoperator":
        A = interface.aslinearoperator(A_)
    elif operator_definition == "__init__matvec":
        A = interface.LinearOperator(
            shape=A_.shape, dtype=A_.dtype,
            matvec=lambda x: matvec(A_, x),
            rmatvec=lambda x: matvec(xp.conj(A_.mT), x),
            xp=xp,
        )
    elif operator_definition == "__init__matmat":
        # A = interface.LinearOperator(shape=A_.shape, dtype=A_.dtype,
        #                              matmat=lambda X: matmat(A_, X))
        pytest.skip("should work but doesn't - see gh-24510")
    elif operator_definition == "subclass_matvec":
        class MyLinearOperator(interface.LinearOperator):
            def _matvec(self, x):
                return matvec(A_, x)
            def _rmatvec(self, x):
                return matvec(xp.conj(A_.mT), x)
        A = MyLinearOperator(shape=A_.shape, dtype=A_.dtype, xp=xp)
    elif operator_definition == "subclass_matmat":
        class MyLinearOperator(interface.LinearOperator):
            def _matmat(self, X):
                return matmat(A_, X)
            def _rmatmat(self, X):
                return matmat(xp.conj(A_.mT), X)
        A = MyLinearOperator(shape=A_.shape, dtype=A_.dtype, xp=xp)

    # Test matvec
    # a. with row vector (or batch of row vectors)
    xp_assert_close(
        A.matvec(x_row) if left else A.rmatvec(x_row),
        matvec(A_, x_row) if left else matvec(xp.conj(A_.mT), x_row)
    )
    # b. with column vector (or batch of column vectors)
    if batch_x:
        message = "Dimension mismatch:..."
        with pytest.raises(ValueError, match=message):
            A.matvec(x_col) if left else A.rmatvec(x_col)
    else:
        message = ("Calling `matvec` on 'column vectors'..." if left
                   else "Calling `rmatvec` on 'column vectors'...")
        with pytest.warns(FutureWarning, match=message):
            xp_assert_close(
                A.matvec(x_col) if left else A.rmatvec(x_col),
                A_ @ x_col if left else xp.conj(A_.mT) @ x_col,
            )
    # c. with matrix (or batch of matrices)
    with pytest.raises(ValueError, match="Dimension mismatch:..."):
        A.matvec(x_mat) if left else A.rmatvec(x_mat)

    # Test matmat
    # a. with row vector (or batch of row vectors)
    message = "Dimension mismatch:..." if batch_x else "Expected at least 2-d..."
    with pytest.raises(ValueError, match=message):
        A.matmat(x_row) if left else A.rmatmat(x_row)
    # b. with column vector (or batch of column vectors)
    xp_assert_close(
        A.matmat(x_col) if left else A.rmatmat(x_col),
        A_ @ x_col if left else xp.conj(A_.mT) @ x_col,
    )
    # c. with matrix (or batch of matrices)
    xp_assert_close(
        A.matmat(x_mat) if left else A.rmatmat(x_mat),
        A_ @ x_mat if left else xp.conj(A_.mT) @ x_mat,
    )

    # test __matmul__ (via `@`), `__call__`, and `dot`
    if left:
        # a. with row vector (or batch of row vectors)
        if batch_x:
            with pytest.raises(ValueError, match="Dimension mismatch:..."):
                A @ x_row
        else:
            xp_assert_close(A @ x_row, A_ @ x_row)
            xp_assert_close(A(x_row), A_ @ x_row)
            xp_assert_close(A.dot(x_row), A_ @ x_row)
        # b. with column vector (or batch of column vectors)
        xp_assert_close(A @ x_col, A_ @ x_col)
        xp_assert_close(A(x_col), A_ @ x_col)
        xp_assert_close(A.dot(x_col), A_ @ x_col)
        # c. with matrix (or batch of matrices)
        xp_assert_close(A @ x_mat, A_ @ x_mat)
        xp_assert_close(A(x_mat), A_ @ x_mat)
        xp_assert_close(A.dot(x_mat), A_ @ x_mat)
    else:
        # a. with row vector (or batch of row vectors)
        broadcastable = True
        try:
            # `batch_x[-1]` becomes part of the core shape for matrix
            # multiplication with 1-D row vector
            # success thus depends on whether `batch_A` can broadcast
            # with `batch_x[:-1]`
            np.broadcast_shapes(batch_A, batch_x[:-1])
        except ValueError:
            broadcastable = False
        if not broadcastable:
            # Maybe should be "Dimension mismatch..."?
            msg = "Incompatible shapes|size of tensor|could not be broadcast..."
            with pytest.raises((ValueError, RuntimeError), match=msg):
                x_row @ A
        else:
            xp_assert_close(x_row @ A, x_row @ A_)
            xp_assert_close(A.rdot(x_row), x_row @ A_)
        # b. with column vector (or batch of column vectors)
        with pytest.raises(ValueError, match="Dimension mismatch:..."):
            x_col @ A
        # c. with matrix (or batch of matrices)
        xp_assert_close(xp.conj(x_mat.mT) @ A, xp.conj(x_mat.mT) @ A_)
        xp_assert_close(A.rdot(xp.conj(x_mat.mT)), xp.conj(x_mat.mT) @ A_)
