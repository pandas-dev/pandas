'''Test dtype deprecations.
'''
import numpy as np
import scipy.linalg as linalg
import pytest

from scipy._lib._array_api import SCIPY_ARRAY_API

#
# When the deprecation expires
#     1. Change the test from asserting a DeprecationWarning to assert_raises
#     2. Remove the `_skip_dict` dictionary: it is only used to tell things which
#        already raise errors from things which now emit DeprecationWarnings


# skip these completely
_skip_names = [
    # not functions
    'LinAlgError', 'LinAlgWarning',
    # do not touch the low-level `scipy.linalg.{blas,lapack}`
    'get_lapack_funcs', 'get_blas_funcs', 'find_best_blas_type',
    # qr_update skip for now
    'qr_update', 'qr_insert', 'qr_delete',
    # special matrix constructors do not use LAPACK
    'block_diag',
    'toeplitz', 'circulant', 'hankel', 'hadamard', 'leslie', 'companion',
    'helmert', 'hilbert', 'invhilbert', 'pascal', 'invpascal', 'dft',
    'fiedler', 'fiedler_companion', 'convolution_matrix',
    # _sketches does not use LAPACK
    'clarkson_woodruff_transform',
    # these did not accepted weird dtypes already
    'expm', 'tanm', 'tanhm', 'sinm', 'cosm', 'sinhm', 'coshm', 'expm_cond',
    'solve_continuous_are',
    'rsf2csf',
]


# use this for dtype-specific skips
# (e.g. if something only accepts real arrays, add complex typecodes here.)
_skip_dict = {
    # TypeError: Only real arrays currently supported
    'eigh_tridiagonal': 'G',
    'eigvalsh_tridiagonal': 'G',
    # TypeError: No matching signature found
    'issymmetric': 'eGmM',
    'ishermitian': 'eGmM',
    'bandwidth': 'egGmM',
    # already raise
    'det': 'gGmM',
    'expm': 'gG',

    # datetime specific + float16 specific
    'cdf2rdf' : 'mM',
    'expm_frechet' : 'mM' + 'e',
    'pinv': 'mM',
    'signm': 'mM',
    'fractional_matrix_power' : 'mM' + 'e',
    'khatri_rao' : 'mM',
    'logm' : 'mM' + 'e',
    'orthogonal_procrustes' : 'mM' + 'e',
    'solve_circulant' : 'mM',
    'solve_continuous_lyapunov': 'M',
    'solve_discrete_lyapunov': 'mM' + 'e',
    'solve_lyapunov' : 'M',
    'solve_sylvester' : 'mM',
    'solve_discrete_are': 'mM',
    'matmul_toeplitz': 'mM',
    'solve_toeplitz': 'mM',
}


# these raise if SCIPY_ARRAY_API env var is defined
_array_api_skip_dict = {
    'cho_factor' : 'mM',
    'cho_solve' : 'mM',
    'cho_solve_banded' : 'mM',
    'cholesky' : 'mM',
    'cholesky_banded' : 'mM',
    'cossin' : 'mM',
    'diagsvd' : 'mM',
    'eig_banded' : 'mM',
    'eigh' : 'mM',
    'eigh_tridiagonal' : 'mM',
    'eigvals_banded' : 'mM',
    'eigvalsh' : 'mM',
    'eigvalsh_tridiagonal' : 'mM',
    'funm' : 'mM',
    'hessenberg' : 'mM',
    'ldl' : 'mM',
    'lu_factor' : 'mM',
    'lu_solve' : 'mM',
    'matrix_balance' : 'mM',
    'null_space' : 'mM',
    'ordqz' : 'mM',
    'orth' : 'mM',
    'pinvh' : 'mM',
    'polar' : 'mM',
    'qr' : 'mM',
    'qr_multiply' : 'mM',
    'qz' : 'mM',
    'rq' : 'mM',
    'schur' : 'mM',
    'solve_banded' : 'mM',
    'solve_continuous_lyapunov' : 'mM',
    'solve_lyapunov' : 'mM',
    'solve_triangular' : 'mM',
    'solveh_banded' : 'mM',
    'subspace_angles' : 'mM',
}


if SCIPY_ARRAY_API:
    #_skip_dict.update(**_array_api_skip_dict)
    for key, val in _array_api_skip_dict.items():
        if key in _skip_dict:
            _skip_dict[key] += val
        else:
            _skip_dict[key] = val


# loop over all names
_names = linalg.__all__
_objs = {name: getattr(linalg, name)
    for name in _names
    if name not in _skip_names
}
_objs = {name: obj for name,obj in _objs.items() if callable(obj)}


# Some gymnastics to prepare a set of arguments which all functions can accept
_two_arg_names = set([
    'cdf2rdf', 'rsf2csf', 'solve', 'solve_triangular', 'lstsq', 'khatri_rao',
    'subspace_angles', 'qz', 'ordqz', 'orthogonal_procrustes', 'solve_circulant',
    'solve_discrete_lyapunov', 'solve_continuous_lyapunov', 'solve_lyapunov',
    'expm_frechet',
    'matmul_toeplitz', 'solve_toeplitz', 'lu_solve',
    'eigh_tridiagonal', 'eigvalsh_tridiagonal',
    'solveh_banded', 'solve_banded', 'solveh_banded', 'qr_multiply',
])

_three_arg_names = set(['solve_sylvester'])

_four_arg_names = set(['solve_continuous_are', 'solve_discrete_are'])

_arr_and_int_names = set([
    'clarkson_woodruff_transform', 'convolution_matrix', 'fractional_matrix_power',
])

_arr_and_two_int_names = set(['diagsvd', 'cossin'])


def _patch_args(func_name, args):
    '''Make sure func(*args) does not raise because *args is wrong.'''
    if func_name == 'cdf2rdf':
        args = (args[0][0], args[1])  # cdf2rd(1d, 2d)
    elif func_name == 'funm':
        args = (args[0], lambda x: x)
    elif func_name == 'lu_solve':
        a = args[0]
        piv = np.arange(a.shape[0])
        args = ((a, piv), args[1])
    elif func_name in ('eigh_tridiagonal', 'eigvalsh_tridiagonal'):       
        d, e = args
        args = (d, e[:, :-1])
    elif func_name == 'cossin':
        a = args[0]
        args = (a, a.shape[0]//2, a.shape[1]//2)
    elif func_name == 'solve_banded':
        from scipy.linalg._basic import _to_banded
        args = ((1, 2), _to_banded(1, 2, args[0]), args[0])
    elif func_name == 'solveh_banded':
        from scipy.linalg._basic import _to_banded
        args = (_to_banded(0, 2, args[0])[:, :-1], args[0][0, :-1])
    elif func_name == 'cholesky_banded':
        from scipy.linalg._basic import _to_banded
        args = (_to_banded(0, 2, args[0])[:, :-1],)
    elif func_name == 'solve_toeplitz':
        c = args[0][1]
        args = (c, np.ones_like(c))
    return args


@pytest.mark.skip_xp_backends(np_only=True)
@pytest.mark.filterwarnings(
    'ignore:attempting to conjugate non-numeric dtype:DeprecationWarning'
)
@pytest.mark.filterwarnings('ignore:overflow encountered in dot:RuntimeWarning')
@pytest.mark.filterwarnings('ignore:invalid value encountered in dot:RuntimeWarning')
@pytest.mark.filterwarnings('ignore:logm result may be inaccurate:RuntimeWarning')
@pytest.mark.filterwarnings('ignore:overflow encountered in multiply:RuntimeWarning')
@pytest.mark.parametrize('func_name', _objs.keys())
@pytest.mark.parametrize('tcode', 'e' + 'mMgG')
def test_deprecations(func_name, tcode):
    func = _objs[func_name]

    if tcode in _skip_dict.get(func_name, ''):
        return

    a = np.arange(16).reshape(4, 4) + 8*np.eye(4)
    a = a.T @ a
    a = a.astype(tcode)

    if func_name in _two_arg_names:
        b = a.copy()
        args = (a, b)
    elif func_name in _three_arg_names:
        args = (a, a, a)
    elif func_name in _four_arg_names:
        args = (a, a, a, a)
    elif func_name in _arr_and_int_names:
        args = (a, a.shape[0])
    elif func_name in _arr_and_two_int_names:
        args = (a, a.shape[0], a.shape[1])
    elif func_name in ['cho_solve', 'cho_solve_banded']:
        args = ((a, True), a)
    else:
        args = (a,)

    args = _patch_args(func_name, args)

    with pytest.warns(DeprecationWarning):
        func(*args)

