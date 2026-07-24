"""
Low-level BLAS functions (:mod:`scipy.linalg.blas`)
===================================================

This module contains low-level functions from the BLAS library.

.. versionadded:: 0.12.0

.. note::

   The common ``overwrite_<>`` option in many routines, allows the
   input arrays to be overwritten to avoid extra memory allocation.
   However this requires the array to satisfy two conditions
   which are memory order and the data type to match exactly the
   order and the type expected by the routine.

   As an example, if you pass a double precision float array to any
   ``S....`` routine which expects single precision arguments, f2py
   will create an intermediate array to match the argument types and
   overwriting will be performed on that intermediate array.

   Similarly, if a C-contiguous array is passed, f2py will pass a
   FORTRAN-contiguous array internally. Please make sure that these
   details are satisfied. More information can be found in the f2py
   documentation.

.. warning::

   These functions do little to no error checking.
   It is possible to cause crashes by mis-using them,
   so prefer using the higher-level routines in `scipy.linalg`.

.. note::

    Prefer using ``get_blas_funcs`` to importing the bare functions directly.
    If you do, for example, ``from scipy.linalg.blas import dgemm``, the ``dgemm``
    function may be either LP64 or ILP64, depending on how SciPy is built.

    The following is more robust:

    >>> from scipy.linalg.blas import get_blas_funcs
    >>> dgemm = get_blas_funcs('gemm', dtype='float64', ilp64='preferred')
    >>> dgemm.int_dtype
    dtype('int32')    # may vary


Finding functions
-----------------

.. autosummary::
   :toctree: generated/

   get_blas_funcs
   find_best_blas_type

BLAS Level 1 functions
----------------------

.. autosummary::
   :toctree: generated/

   sasum
   saxpy
   scasum
   scnrm2
   scopy
   sdot
   snrm2
   srot
   srotg
   srotm
   srotmg
   sscal
   sswap
   dasum
   daxpy
   dcopy
   ddot
   dnrm2
   drot
   drotg
   drotm
   drotmg
   dscal
   dswap
   dzasum
   dznrm2
   icamax
   idamax
   isamax
   izamax
   caxpy
   ccopy
   cdotc
   cdotu
   crotg
   cscal
   csrot
   csscal
   cswap
   zaxpy
   zcopy
   zdotc
   zdotu
   zdrot
   zdscal
   zrotg
   zscal
   zswap

BLAS Level 2 functions
----------------------

.. autosummary::
   :toctree: generated/

   sgbmv
   sgemv
   sger
   ssbmv
   sspmv
   sspr
   sspr2
   ssymv
   ssyr
   ssyr2
   stbmv
   stbsv
   stpmv
   stpsv
   strmv
   strsv
   dgbmv
   dgemv
   dger
   dsbmv
   dspmv
   dspr
   dspr2
   dsymv
   dsyr
   dsyr2
   dtbmv
   dtbsv
   dtpmv
   dtpsv
   dtrmv
   dtrsv
   cgbmv
   cgemv
   cgerc
   cgeru
   chbmv
   chemv
   cher
   cher2
   chpmv
   chpr
   chpr2
   cspmv
   cspr
   csyr
   ctbmv
   ctbsv
   ctpmv
   ctpsv
   ctrmv
   ctrsv
   zgbmv
   zgemv
   zgerc
   zgeru
   zhbmv
   zhemv
   zher
   zher2
   zhpmv
   zhpr
   zhpr2
   zspmv
   zspr
   zsyr
   ztbmv
   ztbsv
   ztpmv
   ztpsv
   ztrmv
   ztrsv

BLAS Level 3 functions
----------------------

.. autosummary::
   :toctree: generated/

   sgemm
   ssymm
   ssyr2k
   ssyrk
   strmm
   strsm
   dgemm
   dsymm
   dsyr2k
   dsyrk
   dtrmm
   dtrsm
   cgemm
   chemm
   cher2k
   cherk
   csymm
   csyr2k
   csyrk
   ctrmm
   ctrsm
   zgemm
   zhemm
   zher2k
   zherk
   zsymm
   zsyr2k
   zsyrk
   ztrmm
   ztrsm

"""
#
# Author: Pearu Peterson, March 2002
#         refactoring by Fabian Pedregosa, March 2010
#

__all__ = ['get_blas_funcs', 'find_best_blas_type']

import numpy as np
import functools
from scipy.__config__ import CONFIG

# If `_fblas` was built, it means the Cython BLAS ABI is LP64, and we're then also
# keeping `linalg.blas` as LP64.
HAS_LP64 = not bool(CONFIG['Build Dependencies']['blas']['cython blas ilp64'])
HAS_ILP64 = CONFIG['Build Dependencies']['blas']['has ilp64']
del CONFIG

_fblas = None
if HAS_LP64:
    from scipy.linalg import _fblas

_fblas_64 = None
if HAS_ILP64:
    from scipy.linalg import _fblas_64

if not (HAS_LP64 or HAS_ILP64):
    raise RuntimeError("SciPy needs either LP64 or ILP64 BLAS.")

if HAS_LP64:
    from scipy.linalg._fblas import *  # noqa: E402, F403
else:
    from scipy.linalg._fblas_64 import *  # noqa: E402, F403

# all numeric dtypes '?bBhHiIlLqQefdgFDGO' that are safe to be converted to

# single precision float   : '?bBhH!!!!!!ef!!!!!!'
# double precision float   : '?bBhHiIlLqQefdg!!!!'
# single precision complex : '?bBhH!!!!!!ef!!F!!!'
# double precision complex : '?bBhHiIlLqQefdgFDG!'

_type_score = {x: 1 for x in '?bBhHef'}
_type_score.update({x: 2 for x in 'iIlLqQd'})

# Handle float128(g) and complex256(G) separately in case non-Windows systems.
# On Windows, the values will be rewritten to the same key with the same value.
_type_score.update({'F': 3, 'D': 4, 'g': 2, 'G': 4})

# Final mapping to the actual prefixes and dtypes
_type_conv = {1: ('s', np.dtype('float32')),
              2: ('d', np.dtype('float64')),
              3: ('c', np.dtype('complex64')),
              4: ('z', np.dtype('complex128'))}

# some convenience alias for complex functions
_blas_alias = {'cnrm2': 'scnrm2', 'znrm2': 'dznrm2',
               'cdot': 'cdotc', 'zdot': 'zdotc',
               'cger': 'cgerc', 'zger': 'zgerc',
               'sdotc': 'sdot', 'sdotu': 'sdot',
               'ddotc': 'ddot', 'ddotu': 'ddot'}


def find_best_blas_type(arrays=(), dtype=None):
    """Find best-matching BLAS/LAPACK type.

    Arrays are used to determine the optimal prefix of BLAS routines.

    Parameters
    ----------
    arrays : sequence of ndarrays, optional
        Arrays can be given to determine optimal prefix of BLAS
        routines. If not given, double-precision routines will be
        used, otherwise the most generic type in arrays will be used.
    dtype : str or dtype, optional
        Data-type specifier. Not used if `arrays` is non-empty.

    Returns
    -------
    prefix : str
        BLAS/LAPACK prefix character.
    dtype : dtype
        Inferred Numpy data type.
    prefer_fortran : bool
        Whether to prefer Fortran order routines over C order.

    Examples
    --------
    >>> import numpy as np
    >>> import scipy.linalg.blas as bla
    >>> rng = np.random.default_rng()
    >>> a = rng.random((10,15))
    >>> b = np.asfortranarray(a)  # Change the memory layout order
    >>> bla.find_best_blas_type((a,))
    ('d', dtype('float64'), False)
    >>> bla.find_best_blas_type((a*1j,))
    ('z', dtype('complex128'), False)
    >>> bla.find_best_blas_type((b,))
    ('d', dtype('float64'), True)

    """
    dtype = np.dtype(dtype)
    max_score = _type_score.get(dtype.char, 5)
    prefer_fortran = False

    if arrays:
        # In most cases, single element is passed through, quicker route
        if len(arrays) == 1:
            max_score = _type_score.get(arrays[0].dtype.char, 5)
            prefer_fortran = arrays[0].flags['FORTRAN']
        else:
            # use the most generic type in arrays
            scores = [_type_score.get(x.dtype.char, 5) for x in arrays]
            max_score = max(scores)
            ind_max_score = scores.index(max_score)
            # safe upcasting for mix of float64 and complex64 --> prefix 'z'
            if max_score == 3 and (2 in scores):
                max_score = 4

            if arrays[ind_max_score].flags['FORTRAN']:
                # prefer Fortran for leading array with column major order
                prefer_fortran = True

    # Get the LAPACK prefix and the corresponding dtype if not fall back
    # to 'd' and double precision float.
    prefix, dtype = _type_conv.get(max_score, ('d', np.dtype('float64')))

    return prefix, dtype, prefer_fortran


def _get_funcs(
    names, arrays, dtype, lib_name, fmodule, fmodule_name, alias, ilp64="preferred"
):
    """
    Return available BLAS/LAPACK functions.

    Used also in lapack.py. See get_blas_funcs for docstring.
    """

    funcs = []
    unpack = False
    dtype = np.dtype(dtype)

    if isinstance(names, str):
        names = (names,)
        unpack = True

    prefix, dtype, _ = find_best_blas_type(arrays, dtype)

    for name in names:
        func_name = prefix + name
        func_name = alias.get(func_name, func_name)
        func = getattr(fmodule, func_name, None)
        if func is None:
            raise ValueError(
                f'{lib_name} function {func_name} could not be found')
        func.module_name, func.typecode, func.dtype = fmodule_name, prefix, dtype
        if not ilp64:
            func.int_dtype = np.dtype(np.intc)
        else:
            func.int_dtype = np.dtype(np.int64)
        func.prefix = prefix  # Backward compatibility
        funcs.append(func)

    if unpack:
        return funcs[0]
    else:
        return funcs


def _memoize_get_funcs(func):
    """
    Memoized fast path for _get_funcs instances
    """
    memo = {}
    func.memo = memo

    @functools.wraps(func)
    def getter(names, arrays=(), dtype=None, ilp64="preferred"):
        key = (names, dtype, ilp64)
        for array in arrays:
            # cf. find_blas_funcs
            key += (array.dtype.char, array.flags.fortran)

        try:
            value = memo.get(key)
        except TypeError:
            # unhashable key etc.
            key = None
            value = None

        if value is not None:
            return value

        value = func(names, arrays, dtype, ilp64)

        if key is not None:
            memo[key] = value

        return value

    return getter


@_memoize_get_funcs
def get_blas_funcs(names, arrays=(), dtype=None, ilp64="preferred"):
    """Return available BLAS function objects from names.

    Arrays are used to determine the optimal prefix of BLAS routines.

    Parameters
    ----------
    names : str or sequence of str
        Name(s) of BLAS functions without type prefix.

    arrays : sequence of ndarrays, optional
        Arrays can be given to determine optimal prefix of BLAS
        routines. If not given, double-precision routines will be
        used, otherwise the most generic type in arrays will be used.

    dtype : str or dtype, optional
        Data-type specifier. Not used if `arrays` is non-empty.

    ilp64 : {True, False, 'preferred'}, optional
        Whether to return ILP64 routine variant.
        Choosing ``'preferred'`` returns ILP64 routine if available,
        and otherwise the 32-bit (LP64) routine. Default: ``'preferred'``.

    Returns
    -------
    funcs : list
        List containing the found function(s).

    Raises
    ------
    RuntimeError
        If the requested LP64/ILP64 variant is not available.

    See Also
    --------
    scipy.linalg.lapack.get_lapack_funcs
        a similar routine for selecting LAPACK functions.

    Notes
    -----
    In BLAS, the naming convention is that all functions start with a
    type prefix, which depends on the type of the principal
    matrix. These can be one of ``{'s', 'd', 'c', 'z'}`` for the NumPy
    types ``{float32, float64, complex64, complex128}`` respectively.
    The code and the dtype are stored in attributes `typecode` and `dtype`
    of the returned functions.

    Examples
    --------
    >>> import numpy as np
    >>> import scipy.linalg as LA
    >>> rng = np.random.default_rng()
    >>> a = rng.random((3, 2))
    >>> x_gemv = LA.get_blas_funcs('gemv', (a,))

    Note that ``x_gemv`` string representation shows the exact BLAS function with the
    prefix (here, ``d-`` because ``a`` is double precision real):

    >>> x_gemv
    <fortran function dgemv>

    The BLAS variant information is also available from the ``typecode`` attribute:

    >>> x_gemv.typecode
    'd'

    For double precision complex arrays, we select the ``z-`` variant, ``zgemv``

    >>> x_gemv = LA.get_blas_funcs('gemv', (a*1j,))
    >>> x_gemv.typecode
    'z'

    If you want to select a specific BLAS variant instead of relying on array types, use
    the ``dtype=`` argument:

    >>> LA.get_blas_funcs('gemv', dtype=np.float32)
    <fortran function sgemv>

    The ``int_dtype`` attribute stores whether the routine is ILP64 (integer arguments
    and outputs are 64-bit) or LP64 (integer arguments and outputs are 32-bit):

    >>> x_gemv.int_dtype
    dtype('int32')   # may vary
    """
    if isinstance(ilp64, str):
        if ilp64 == 'preferred':
            ilp64 = HAS_ILP64
        else:
            raise ValueError(f"Invalid value for {ilp64 = }.")

    if not ilp64:
        return _get_funcs(
            names, arrays, dtype, "BLAS", _fblas, "fblas", _blas_alias, ilp64=False
        )
    else:
        if not HAS_ILP64:
            raise RuntimeError("BLAS ILP64 routine requested, but Scipy "
                               "compiled only with 32-bit BLAS")
        return _get_funcs(
            names, arrays, dtype, "BLAS", _fblas_64, "fblas_64", _blas_alias, ilp64=True
        )
