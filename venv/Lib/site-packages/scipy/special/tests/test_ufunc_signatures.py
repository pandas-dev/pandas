"""Test that all ufuncs have float32-preserving signatures.

This was once guaranteed through the code generation script for
generating ufuncs, `scipy/special/_generate_pyx.py`. Starting with
gh-20260, SciPy developers have begun moving to generate ufuncs
through direct use of the NumPy C API (through C++). Existence of
float32 preserving signatures must now be tested since it is no
longer guaranteed.
"""

import numpy as np
import pytest
import scipy.special._ufuncs
import scipy.special._gufuncs
import scipy.special as special

from numpy.testing import assert_equal


# Single precision is not implemented for these ufuncs;
# floating point inputs must be float64.
exceptions = ['_gen_harmonic', '_normalized_gen_harmonic']

_ufuncs = []
for funcname in dir(scipy.special._ufuncs):
    if funcname not in exceptions:
        _ufuncs.append(getattr(scipy.special._ufuncs, funcname))
for funcname in dir(scipy.special._gufuncs):
    _ufuncs.append(getattr(scipy.special._gufuncs, funcname))

# Not all module members are actually ufuncs
_ufuncs = [func for func in _ufuncs if isinstance(func, np.ufunc)]

@pytest.mark.parametrize("ufunc", _ufuncs)
def test_ufunc_signatures(ufunc):

    # From _generate_pyx.py
    # "Don't add float32 versions of ufuncs with integer arguments, as this
    # can lead to incorrect dtype selection if the integer arguments are
    # arrays, but float arguments are scalars.
    # This may be a NumPy bug, but we need to work around it.
    # cf. gh-4895, https://github.com/numpy/numpy/issues/5895"
    types = set(sig for sig in ufunc.types
                if not ("l" in sig or "i" in sig or "q" in sig or "p" in sig))

    # Generate the full expanded set of signatures which should exist. There
    # should be matching float and double versions of any existing signature.
    expanded_types = set()
    for sig in types:
        expanded_types.update(
            [sig.replace("d", "f").replace("D", "F"),
             sig.replace("f", "d").replace("F", "D")]
        )
    assert types == expanded_types


def _get_nan_val(typecode):
    return {
        "f": np.asarray(np.nan, dtype=np.float32),
        "d": np.asarray(np.nan, dtype=np.float64),
        "g": np.asarray(np.nan, dtype=np.longdouble),
        "F": np.asarray(np.nan + np.nan*1j, dtype=np.complex64),
        "D": np.asarray(np.nan + np.nan*1j, dtype=np.complex128),
        "G": np.asarray(np.nan + np.nan*1j, dtype=np.clongdouble),
    }[typecode]

skips = {
    # These hit https://github.com/scipy/xsf/issues/94
    "eval_chebyc",
    "eval_chebys",
    "eval_chebyt",
    "eval_chebyu",
    "eval_gegenbauer",
    "eval_genlaguerre",
    "eval_jacobi",
    "eval_laguerre",
    "eval_legendre",
    "eval_sh_chebyt",
    "eval_sh_chebyu",
    "eval_sh_jacobi",
    "eval_sh_legendre",
    "hyp1f1",
    "hyp2f1",
    # warns if called with noninteger values of some args
    "bdtr",
    "bdtrc",
    "bdtri",
}
_multi_arg_public_ufuncs = [
    ufunc for ufunc in _ufuncs if ufunc.nargs - ufunc.nout > 1
    and ufunc.__name__ not in skips and ufunc.__name__ in special.__all__
]


@pytest.mark.parametrize("ufunc", _multi_arg_public_ufuncs)
def test_nep50(ufunc):
    # Test that functions with multiple arguments respect nep50 promotion rules.
    rng = np.random.default_rng(1234)
    # As in test_ufunc_signatures, filter out signatures involving integers.
    types = set(sig for sig in ufunc.types
                if not ("l" in sig or "i" in sig or "q" in sig or "p" in sig))
    for sig in types:
        input_types, output_types = sig.split("->")
        # since we only care about dtypes and not values here, just use an appropriately
        # typed nan for each argument.
        args = [_get_nan_val(typecode) for typecode in input_types]
        # swap out a random one of the nans with an appropriately typed numpy scalar.
        idx = rng.choice(len(args))
        args[idx] = float("nan") if np.isrealobj(args[idx]) else complex("nan")
        result = ufunc(*args)
        result = [result] if len(output_types) == 1 else result
        result = np.asarray(result)
        
        # Test that the output is an appropriately typed nan. This also implicitly
        # tests that ufuncs propagate nans correctly.
        assert_equal(
            result, np.asarray([_get_nan_val(typecode) for typecode in output_types]),
            strict=True
        )
