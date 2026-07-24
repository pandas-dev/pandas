# ---------------------------------------------------------------------------
# AUTO-GENERATED - do not edit by hand.
# Re-generate with:  python tools/generate_cython_abi_tests.py linalg
#
# Generated against SciPy 1.17.1 on 2026-03-29
# Python 3.12.13
# ---------------------------------------------------------------------------

"""
Regression test for the cython_blas / cython_lapack binary interface.

The __pyx_capi__ dict of each module contains PyCapsule objects whose
*name* string encodes the C-level function signature that Cython checks
at import time when another module does ``cimport``.  Any change to
those strings is a binary incompatibility that breaks every downstream
package (e.g., scikit-learn, statsmodels) compiled against the previous release,
without any compile-time warning.

The expected signatures are stored in cython_abi_signatures.json.
The test fails if:
  * an existing name disappears, or
  * an existing name's signature string changes.
Adding new names is explicitly allowed (it is not a break).
"""

import ctypes
import importlib
import importlib.resources
import json

import pytest
import scipy

CYTHON_BLAS_ILP64 = (
    scipy.show_config(mode='dicts')['Build Dependencies']['blas']['cython blas ilp64']
)


# ---------------------------------------------------------------------------
# Helper
# ---------------------------------------------------------------------------

_get_capsule_name = ctypes.pythonapi.PyCapsule_GetName
_get_capsule_name.restype = ctypes.c_char_p
_get_capsule_name.argtypes = [ctypes.py_object]


def _extract_capi(module):
    """Return {name: signature} for every entry in module.__pyx_capi__."""
    result = {}
    for name, capsule in module.__pyx_capi__.items():
        raw = _get_capsule_name(capsule)
        if raw is not None:
            result[name] = raw.decode('utf-8')
    return result


_sig_file = importlib.resources.files('scipy.linalg.tests').joinpath(
    'cython_abi_signatures.json'
)
with importlib.resources.as_file(_sig_file) as _f:
    _EXPECTED = json.loads(_f.read_text())


@pytest.mark.skipif(
    CYTHON_BLAS_ILP64,
    reason='cython_blas/cython_lapack ABI depends on ILP64 BLAS',
)
def test_cython_blas_abi_stability():
    """No existing cython_blas signature may change or disappear."""
    expected = _EXPECTED['scipy.linalg.cython_blas']
    mod = importlib.import_module('scipy.linalg.cython_blas')
    actual = _extract_capi(mod)
    errors = []
    for name, expected_sig in expected.items():
        if name not in actual:
            errors.append(
                f"REMOVED  {name!r} (was {expected_sig!r})"
            )
        elif actual[name] != expected_sig:
            errors.append(
                f"CHANGED  {name!r}\n"
                f"  expected: {expected_sig!r}\n"
                f"  actual:   {actual[name]!r}"
            )
    if errors:
        joined = '\n'.join(errors)
        pytest.fail(
            f"scipy.linalg.cython_blas.__pyx_capi__ has {len(errors)} "
            f"ABI breakage(s):\n{joined}"
        )


@pytest.mark.skipif(
    CYTHON_BLAS_ILP64,
    reason='cython_blas/cython_lapack ABI depends on ILP64 BLAS',
)
def test_cython_lapack_abi_stability():
    """No existing cython_lapack signature may change or disappear."""
    expected = _EXPECTED['scipy.linalg.cython_lapack']
    mod = importlib.import_module('scipy.linalg.cython_lapack')
    actual = _extract_capi(mod)
    errors = []
    for name, expected_sig in expected.items():
        if name not in actual:
            errors.append(
                f"REMOVED  {name!r} (was {expected_sig!r})"
            )
        elif actual[name] != expected_sig:
            errors.append(
                f"CHANGED  {name!r}\n"
                f"  expected: {expected_sig!r}\n"
                f"  actual:   {actual[name]!r}"
            )
    if errors:
        joined = '\n'.join(errors)
        pytest.fail(
            f"scipy.linalg.cython_lapack.__pyx_capi__ has {len(errors)} "
            f"ABI breakage(s):\n{joined}"
        )
