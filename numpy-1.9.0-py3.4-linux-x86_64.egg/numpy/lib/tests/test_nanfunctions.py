from __future__ import division, absolute_import, print_function

import warnings

import numpy as np
from numpy.testing import (
    run_module_suite, TestCase, assert_, assert_equal, assert_almost_equal,
    assert_raises, assert_array_equal
    )


# Test data
_ndat = np.array([[0.6244, np.nan, 0.2692, 0.0116, np.nan, 0.1170],
                  [0.5351, -0.9403, np.nan, 0.2100, 0.4759, 0.2833],
                  [np.nan, np.nan, np.nan, 0.1042, np.nan, -0.5954],
                  [0.1610, np.nan, np.nan, 0.1859, 0.3146, np.nan]])


# Rows of _ndat with nans removed
_rdat = [np.array([0.6244, 0.2692, 0.0116, 0.1170]),
         np.array([0.5351, -0.9403, 0.2100, 0.4759, 0.2833]),
         np.array([0.1042, -0.5954]),
         np.array([0.1610, 0.1859, 0.3146])]


class TestNanFunctions_MinMax(TestCase):

    nanfuncs = [np.nanmin, np.nanmax]
    stdfuncs = [np.min, np.max]

    def test_mutation(self):
        # Check that passed array is not modified.
        ndat = _ndat.copy()
        for f in self.nanfuncs:
            f(ndat)
            assert_equal(ndat, _ndat)

    def test_keepdims(self):
        mat = np.eye(3)
        for nf, rf in zip(self.nanfuncs, self.stdfuncs):
            for axis in [None, 0, 1]:
                tgt = rf(mat, axis=axis, keepdims=True)
                res = nf(mat, axis=axis, keepdims=True)
                assert_(res.ndim == tgt.ndim)

    def test_out(self):
        mat = np.eye(3)
        for nf, rf in zip(self.nanfuncs, self.stdfuncs):
            resout = np.zeros(3)
            tgt = rf(mat, axis=1)
            res = nf(mat, axis=1, out=resout)
            assert_almost_equal(res, resout)
            assert_almost_equal(res, tgt)

    def test_dtype_from_input(self):
        codes = 'efdgFDG'
        for nf, rf in zip(self.nanfuncs, self.stdfuncs):
            for c in codes:
                mat = np.eye(3, dtype=c)
                tgt = rf(mat, axis=1).dtype.type
                res = nf(mat, axis=1).dtype.type
                assert_(res is tgt)
                # scalar case
                tgt = rf(mat, axis=None).dtype.type
                res = nf(mat, axis=None).dtype.type
                assert_(res is tgt)

    def test_result_values(self):
        for nf, rf in zip(self.nanfuncs, self.stdfuncs):
            tgt = [rf(d) for d in _rdat]
            res = nf(_ndat, axis=1)
            assert_almost_equal(res, tgt)

    def test_allnans(self):
        mat = np.array([np.nan]*9).reshape(3, 3)
        for f in self.nanfuncs:
            for axis in [None, 0, 1]:
                with warnings.catch_warnings(record=True) as w:
                    warnings.simplefilter('always')
                    assert_(np.isnan(f(mat, axis=axis)).all())
                    assert_(len(w) == 1, 'no warning raised')
                    assert_(issubclass(w[0].category, RuntimeWarning))
            # Check scalars
            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter('always')
                assert_(np.isnan(f(np.nan)))
                assert_(len(w) == 1, 'no warning raised')
                assert_(issubclass(w[0].category, RuntimeWarning))

    def test_masked(self):
        mat = np.ma.fix_invalid(_ndat)
        msk = mat._mask.copy()
        for f in [np.nanmin]:
            res = f(mat, axis=1)
            tgt = f(_ndat, axis=1)
            assert_equal(res, tgt)
            assert_equal(mat._mask, msk)
            assert_(not np.isinf(mat).any())

    def test_scalar(self):
        for f in self.nanfuncs:
            assert_(f(0.) == 0.)

    def test_matrices(self):
        # Check that it works and that type and
        # shape are preserved
        mat = np.matrix(np.eye(3))
        for f in self.nanfuncs:
            res = f(mat, axis=0)
            assert_(isinstance(res, np.matrix))
            assert_(res.shape == (1, 3))
            res = f(mat, axis=1)
            assert_(isinstance(res, np.matrix))
            assert_(res.shape == (3, 1))
            res = f(mat)
            assert_(np.isscalar(res))
        # check that rows of nan are dealt with for subclasses (#4628)
        mat[1] = np.nan
        for f in self.nanfuncs:
            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter('always')
                res = f(mat, axis=0)
                assert_(isinstance(res, np.matrix))
                assert_(not np.any(np.isnan(res)))
                assert_(len(w) == 0)

            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter('always')
                res = f(mat, axis=1)
                assert_(isinstance(res, np.matrix))
                assert_(np.isnan(res[1, 0]) and not np.isnan(res[0, 0])
                        and not np.isnan(res[2, 0]))
                assert_(len(w) == 1, 'no warning raised')
                assert_(issubclass(w[0].category, RuntimeWarning))

            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter('always')
                res = f(mat)
                assert_(np.isscalar(res))
                assert_(res != np.nan)
                assert_(len(w) == 0)


class TestNanFunctions_ArgminArgmax(TestCase):

    nanfuncs = [np.nanargmin, np.nanargmax]

    def test_mutation(self):
        # Check that passed array is not modified.
        ndat = _ndat.copy()
        for f in self.nanfuncs:
            f(ndat)
            assert_equal(ndat, _ndat)

    def test_result_values(self):
        for f, fcmp in zip(self.nanfuncs, [np.greater, np.less]):
            for row in _ndat:
                with warnings.catch_warnings(record=True):
                    warnings.simplefilter('always')
                    ind = f(row)
                    val = row[ind]
                    # comparing with NaN is tricky as the result
                    # is always false except for NaN != NaN
                    assert_(not np.isnan(val))
                    assert_(not fcmp(val, row).any())
                    assert_(not np.equal(val, row[:ind]).any())

    def test_allnans(self):
        mat = np.array([np.nan]*9).reshape(3, 3)
        for f in self.nanfuncs:
            for axis in [None, 0, 1]:
                assert_raises(ValueError, f, mat, axis=axis)
            assert_raises(ValueError, f, np.nan)

    def test_empty(self):
        mat = np.zeros((0, 3))
        for f in self.nanfuncs:
            for axis in [0, None]:
                assert_raises(ValueError, f, mat, axis=axis)
            for axis in [1]:
                res = f(mat, axis=axis)
                assert_equal(res, np.zeros(0))

    def test_scalar(self):
        for f in self.nanfuncs:
            assert_(f(0.) == 0.)

    def test_matrices(self):
        # Check that it works and that type and
        # shape are preserved
        mat = np.matrix(np.eye(3))
        for f in self.nanfuncs:
            res = f(mat, axis=0)
            assert_(isinstance(res, np.matrix))
            assert_(res.shape == (1, 3))
            res = f(mat, axis=1)
            assert_(isinstance(res, np.matrix))
            assert_(res.shape == (3, 1))
            res = f(mat)
            assert_(np.isscalar(res))


class TestNanFunctions_IntTypes(TestCase):

    int_types = (np.int8, np.int16, np.int32, np.int64, np.uint8,
                 np.uint16, np.uint32, np.uint64)

    mat = np.array([127, 39, 93, 87, 46])

    def integer_arrays(self):
        for dtype in self.int_types:
            yield self.mat.astype(dtype)

    def test_nanmin(self):
        tgt = np.min(self.mat)
        for mat in self.integer_arrays():
            assert_equal(np.nanmin(mat), tgt)

    def test_nanmax(self):
        tgt = np.max(self.mat)
        for mat in self.integer_arrays():
            assert_equal(np.nanmax(mat), tgt)

    def test_nanargmin(self):
        tgt = np.argmin(self.mat)
        for mat in self.integer_arrays():
            assert_equal(np.nanargmin(mat), tgt)

    def test_nanargmax(self):
        tgt = np.argmax(self.mat)
        for mat in self.integer_arrays():
            assert_equal(np.nanargmax(mat), tgt)

    def test_nansum(self):
        tgt = np.sum(self.mat)
        for mat in self.integer_arrays():
            assert_equal(np.nansum(mat), tgt)

    def test_nanmean(self):
        tgt = np.mean(self.mat)
        for mat in self.integer_arrays():
            assert_equal(np.nanmean(mat), tgt)

    def test_nanvar(self):
        tgt = np.var(self.mat)
        for mat in self.integer_arrays():
            assert_equal(np.nanvar(mat), tgt)

        tgt = np.var(mat, ddof=1)
        for mat in self.integer_arrays():
            assert_equal(np.nanvar(mat, ddof=1), tgt)

    def test_nanstd(self):
        tgt = np.std(self.mat)
        for mat in self.integer_arrays():
            assert_equal(np.nanstd(mat), tgt)

        tgt = np.std(self.mat, ddof=1)
        for mat in self.integer_arrays():
            assert_equal(np.nanstd(mat, ddof=1), tgt)


class TestNanFunctions_Sum(TestCase):

    def test_mutation(self):
        # Check that passed array is not modified.
        ndat = _ndat.copy()
        np.nansum(ndat)
        assert_equal(ndat, _ndat)

    def test_keepdims(self):
        mat = np.eye(3)
        for axis in [None, 0, 1]:
            tgt = np.sum(mat, axis=axis, keepdims=True)
            res = np.nansum(mat, axis=axis, keepdims=True)
            assert_(res.ndim == tgt.ndim)

    def test_out(self):
        mat = np.eye(3)
        resout = np.zeros(3)
        tgt = np.sum(mat, axis=1)
        res = np.nansum(mat, axis=1, out=resout)
        assert_almost_equal(res, resout)
        assert_almost_equal(res, tgt)

    def test_dtype_from_dtype(self):
        mat = np.eye(3)
        codes = 'efdgFDG'
        for c in codes:
            tgt = np.sum(mat, dtype=np.dtype(c), axis=1).dtype.type
            res = np.nansum(mat, dtype=np.dtype(c), axis=1).dtype.type
            assert_(res is tgt)
            # scalar case
            tgt = np.sum(mat, dtype=np.dtype(c), axis=None).dtype.type
            res = np.nansum(mat, dtype=np.dtype(c), axis=None).dtype.type
            assert_(res is tgt)

    def test_dtype_from_char(self):
        mat = np.eye(3)
        codes = 'efdgFDG'
        for c in codes:
            tgt = np.sum(mat, dtype=c, axis=1).dtype.type
            res = np.nansum(mat, dtype=c, axis=1).dtype.type
            assert_(res is tgt)
            # scalar case
            tgt = np.sum(mat, dtype=c, axis=None).dtype.type
            res = np.nansum(mat, dtype=c, axis=None).dtype.type
            assert_(res is tgt)

    def test_dtype_from_input(self):
        codes = 'efdgFDG'
        for c in codes:
            mat = np.eye(3, dtype=c)
            tgt = np.sum(mat, axis=1).dtype.type
            res = np.nansum(mat, axis=1).dtype.type
            assert_(res is tgt)
            # scalar case
            tgt = np.sum(mat, axis=None).dtype.type
            res = np.nansum(mat, axis=None).dtype.type
            assert_(res is tgt)

    def test_result_values(self):
            tgt = [np.sum(d) for d in _rdat]
            res = np.nansum(_ndat, axis=1)
            assert_almost_equal(res, tgt)

    def test_allnans(self):
        # Check for FutureWarning
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter('always')
            res = np.nansum([np.nan]*3, axis=None)
            assert_(res == 0, 'result is not 0')
            assert_(len(w) == 0, 'warning raised')
            # Check scalar
            res = np.nansum(np.nan)
            assert_(res == 0, 'result is not 0')
            assert_(len(w) == 0, 'warning raised')
            # Check there is no warning for not all-nan
            np.nansum([0]*3, axis=None)
            assert_(len(w) == 0, 'unwanted warning raised')

    def test_empty(self):
        mat = np.zeros((0, 3))
        tgt = [0]*3
        res = np.nansum(mat, axis=0)
        assert_equal(res, tgt)
        tgt = []
        res = np.nansum(mat, axis=1)
        assert_equal(res, tgt)
        tgt = 0
        res = np.nansum(mat, axis=None)
        assert_equal(res, tgt)

    def test_scalar(self):
        assert_(np.nansum(0.) == 0.)

    def test_matrices(self):
        # Check that it works and that type and
        # shape are preserved
        mat = np.matrix(np.eye(3))
        res = np.nansum(mat, axis=0)
        assert_(isinstance(res, np.matrix))
        assert_(res.shape == (1, 3))
        res = np.nansum(mat, axis=1)
        assert_(isinstance(res, np.matrix))
        assert_(res.shape == (3, 1))
        res = np.nansum(mat)
        assert_(np.isscalar(res))


class TestNanFunctions_MeanVarStd(TestCase):

    nanfuncs = [np.nanmean, np.nanvar, np.nanstd]
    stdfuncs = [np.mean, np.var, np.std]

    def test_mutation(self):
        # Check that passed array is not modified.
        ndat = _ndat.copy()
        for f in self.nanfuncs:
            f(ndat)
            assert_equal(ndat, _ndat)

    def test_dtype_error(self):
        for f in self.nanfuncs:
            for dtype in [np.bool_, np.int_, np.object]:
                assert_raises(TypeError, f, _ndat, axis=1, dtype=np.int)

    def test_out_dtype_error(self):
        for f in self.nanfuncs:
            for dtype in [np.bool_, np.int_, np.object]:
                out = np.empty(_ndat.shape[0], dtype=dtype)
                assert_raises(TypeError, f, _ndat, axis=1, out=out)

    def test_keepdims(self):
        mat = np.eye(3)
        for nf, rf in zip(self.nanfuncs, self.stdfuncs):
            for axis in [None, 0, 1]:
                tgt = rf(mat, axis=axis, keepdims=True)
                res = nf(mat, axis=axis, keepdims=True)
                assert_(res.ndim == tgt.ndim)

    def test_out(self):
        mat = np.eye(3)
        for nf, rf in zip(self.nanfuncs, self.stdfuncs):
            resout = np.zeros(3)
            tgt = rf(mat, axis=1)
            res = nf(mat, axis=1, out=resout)
            assert_almost_equal(res, resout)
            assert_almost_equal(res, tgt)

    def test_dtype_from_dtype(self):
        mat = np.eye(3)
        codes = 'efdgFDG'
        for nf, rf in zip(self.nanfuncs, self.stdfuncs):
            for c in codes:
                tgt = rf(mat, dtype=np.dtype(c), axis=1).dtype.type
                res = nf(mat, dtype=np.dtype(c), axis=1).dtype.type
                assert_(res is tgt)
                # scalar case
                tgt = rf(mat, dtype=np.dtype(c), axis=None).dtype.type
                res = nf(mat, dtype=np.dtype(c), axis=None).dtype.type
                assert_(res is tgt)

    def test_dtype_from_char(self):
        mat = np.eye(3)
        codes = 'efdgFDG'
        for nf, rf in zip(self.nanfuncs, self.stdfuncs):
            for c in codes:
                tgt = rf(mat, dtype=c, axis=1).dtype.type
                res = nf(mat, dtype=c, axis=1).dtype.type
                assert_(res is tgt)
                # scalar case
                tgt = rf(mat, dtype=c, axis=None).dtype.type
                res = nf(mat, dtype=c, axis=None).dtype.type
                assert_(res is tgt)

    def test_dtype_from_input(self):
        codes = 'efdgFDG'
        for nf, rf in zip(self.nanfuncs, self.stdfuncs):
            for c in codes:
                mat = np.eye(3, dtype=c)
                tgt = rf(mat, axis=1).dtype.type
                res = nf(mat, axis=1).dtype.type
                assert_(res is tgt, "res %s, tgt %s" % (res, tgt))
                # scalar case
                tgt = rf(mat, axis=None).dtype.type
                res = nf(mat, axis=None).dtype.type
                assert_(res is tgt)

    def test_ddof(self):
        nanfuncs = [np.nanvar, np.nanstd]
        stdfuncs = [np.var, np.std]
        for nf, rf in zip(nanfuncs, stdfuncs):
            for ddof in [0, 1]:
                tgt = [rf(d, ddof=ddof) for d in _rdat]
                res = nf(_ndat, axis=1, ddof=ddof)
                assert_almost_equal(res, tgt)

    def test_ddof_too_big(self):
        nanfuncs = [np.nanvar, np.nanstd]
        stdfuncs = [np.var, np.std]
        dsize = [len(d) for d in _rdat]
        for nf, rf in zip(nanfuncs, stdfuncs):
            for ddof in range(5):
                with warnings.catch_warnings(record=True) as w:
                    warnings.simplefilter('always')
                    tgt = [ddof >= d for d in dsize]
                    res = nf(_ndat, axis=1, ddof=ddof)
                    assert_equal(np.isnan(res), tgt)
                    if any(tgt):
                        assert_(len(w) == 1)
                        assert_(issubclass(w[0].category, RuntimeWarning))
                    else:
                        assert_(len(w) == 0)

    def test_result_values(self):
        for nf, rf in zip(self.nanfuncs, self.stdfuncs):
            tgt = [rf(d) for d in _rdat]
            res = nf(_ndat, axis=1)
            assert_almost_equal(res, tgt)

    def test_allnans(self):
        mat = np.array([np.nan]*9).reshape(3, 3)
        for f in self.nanfuncs:
            for axis in [None, 0, 1]:
                with warnings.catch_warnings(record=True) as w:
                    warnings.simplefilter('always')
                    assert_(np.isnan(f(mat, axis=axis)).all())
                    assert_(len(w) == 1)
                    assert_(issubclass(w[0].category, RuntimeWarning))
                    # Check scalar
                    assert_(np.isnan(f(np.nan)))
                    assert_(len(w) == 2)
                    assert_(issubclass(w[0].category, RuntimeWarning))

    def test_empty(self):
        mat = np.zeros((0, 3))
        for f in self.nanfuncs:
            for axis in [0, None]:
                with warnings.catch_warnings(record=True) as w:
                    warnings.simplefilter('always')
                    assert_(np.isnan(f(mat, axis=axis)).all())
                    assert_(len(w) == 1)
                    assert_(issubclass(w[0].category, RuntimeWarning))
            for axis in [1]:
                with warnings.catch_warnings(record=True) as w:
                    warnings.simplefilter('always')
                    assert_equal(f(mat, axis=axis), np.zeros([]))
                    assert_(len(w) == 0)

    def test_scalar(self):
        for f in self.nanfuncs:
            assert_(f(0.) == 0.)

    def test_matrices(self):
        # Check that it works and that type and
        # shape are preserved
        mat = np.matrix(np.eye(3))
        for f in self.nanfuncs:
            res = f(mat, axis=0)
            assert_(isinstance(res, np.matrix))
            assert_(res.shape == (1, 3))
            res = f(mat, axis=1)
            assert_(isinstance(res, np.matrix))
            assert_(res.shape == (3, 1))
            res = f(mat)
            assert_(np.isscalar(res))


class TestNanFunctions_Median(TestCase):

    def test_mutation(self):
        # Check that passed array is not modified.
        ndat = _ndat.copy()
        np.nanmedian(ndat)
        assert_equal(ndat, _ndat)

    def test_keepdims(self):
        mat = np.eye(3)
        for axis in [None, 0, 1]:
            tgt = np.median(mat, axis=axis, out=None, overwrite_input=False)
            res = np.nanmedian(mat, axis=axis, out=None, overwrite_input=False)
            assert_(res.ndim == tgt.ndim)

        d = np.ones((3, 5, 7, 11))
        # Randomly set some elements to NaN:
        w = np.random.random((4, 200)) * np.array(d.shape)[:, None]
        w = w.astype(np.intp)
        d[tuple(w)] = np.nan
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter('always', RuntimeWarning)
            res = np.nanmedian(d, axis=None, keepdims=True)
            assert_equal(res.shape, (1, 1, 1, 1))
            res = np.nanmedian(d, axis=(0, 1), keepdims=True)
            assert_equal(res.shape, (1, 1, 7, 11))
            res = np.nanmedian(d, axis=(0, 3), keepdims=True)
            assert_equal(res.shape, (1, 5, 7, 1))
            res = np.nanmedian(d, axis=(1,), keepdims=True)
            assert_equal(res.shape, (3, 1, 7, 11))
            res = np.nanmedian(d, axis=(0, 1, 2, 3), keepdims=True)
            assert_equal(res.shape, (1, 1, 1, 1))
            res = np.nanmedian(d, axis=(0, 1, 3), keepdims=True)
            assert_equal(res.shape, (1, 1, 7, 1))

    def test_out(self):
        mat = np.random.rand(3, 3)
        nan_mat = np.insert(mat, [0, 2], np.nan, axis=1)
        resout = np.zeros(3)
        tgt = np.median(mat, axis=1)
        res = np.nanmedian(nan_mat, axis=1, out=resout)
        assert_almost_equal(res, resout)
        assert_almost_equal(res, tgt)
        # 0-d output:
        resout = np.zeros(())
        tgt = np.median(mat, axis=None)
        res = np.nanmedian(nan_mat, axis=None, out=resout)
        assert_almost_equal(res, resout)
        assert_almost_equal(res, tgt)
        res = np.nanmedian(nan_mat, axis=(0, 1), out=resout)
        assert_almost_equal(res, resout)
        assert_almost_equal(res, tgt)

    def test_small_large(self):
        # test the small and large code paths, current cutoff 400 elements
        for s in [5, 20, 51, 200, 1000]:
            d = np.random.randn(4, s)
            # Randomly set some elements to NaN:
            w = np.random.randint(0, d.size, size=d.size // 5)
            d.ravel()[w] = np.nan
            d[:,0] = 1.  # ensure at least one good value
            # use normal median without nans to compare
            tgt = []
            for x in d:
                nonan = np.compress(~np.isnan(x), x)
                tgt.append(np.median(nonan, overwrite_input=True))

            assert_array_equal(np.nanmedian(d, axis=-1), tgt)

    def test_result_values(self):
            tgt = [np.median(d) for d in _rdat]
            res = np.nanmedian(_ndat, axis=1)
            assert_almost_equal(res, tgt)

    def test_allnans(self):
        mat = np.array([np.nan]*9).reshape(3, 3)
        for axis in [None, 0, 1]:
            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter('always')
                assert_(np.isnan(np.nanmedian(mat, axis=axis)).all())
                if axis is None:
                    assert_(len(w) == 1)
                else:
                    assert_(len(w) == 3)
                assert_(issubclass(w[0].category, RuntimeWarning))
                # Check scalar
                assert_(np.isnan(np.nanmedian(np.nan)))
                if axis is None:
                    assert_(len(w) == 2)
                else:
                    assert_(len(w) == 4)
                assert_(issubclass(w[0].category, RuntimeWarning))

    def test_empty(self):
        mat = np.zeros((0, 3))
        for axis in [0, None]:
            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter('always')
                assert_(np.isnan(np.nanmedian(mat, axis=axis)).all())
                assert_(len(w) == 1)
                assert_(issubclass(w[0].category, RuntimeWarning))
        for axis in [1]:
            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter('always')
                assert_equal(np.nanmedian(mat, axis=axis), np.zeros([]))
                assert_(len(w) == 0)

    def test_scalar(self):
        assert_(np.nanmedian(0.) == 0.)

    def test_extended_axis_invalid(self):
        d = np.ones((3, 5, 7, 11))
        assert_raises(IndexError, np.nanmedian, d, axis=-5)
        assert_raises(IndexError, np.nanmedian, d, axis=(0, -5))
        assert_raises(IndexError, np.nanmedian, d, axis=4)
        assert_raises(IndexError, np.nanmedian, d, axis=(0, 4))
        assert_raises(ValueError, np.nanmedian, d, axis=(1, 1))


class TestNanFunctions_Percentile(TestCase):

    def test_mutation(self):
        # Check that passed array is not modified.
        ndat = _ndat.copy()
        np.nanpercentile(ndat, 30)
        assert_equal(ndat, _ndat)

    def test_keepdims(self):
        mat = np.eye(3)
        for axis in [None, 0, 1]:
            tgt = np.percentile(mat, 70, axis=axis, out=None,
                                overwrite_input=False)
            res = np.nanpercentile(mat, 70, axis=axis, out=None,
                                   overwrite_input=False)
            assert_(res.ndim == tgt.ndim)

        d = np.ones((3, 5, 7, 11))
        # Randomly set some elements to NaN:
        w = np.random.random((4, 200)) * np.array(d.shape)[:, None]
        w = w.astype(np.intp)
        d[tuple(w)] = np.nan
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter('always', RuntimeWarning)
            res = np.nanpercentile(d, 90, axis=None, keepdims=True)
            assert_equal(res.shape, (1, 1, 1, 1))
            res = np.nanpercentile(d, 90, axis=(0, 1), keepdims=True)
            assert_equal(res.shape, (1, 1, 7, 11))
            res = np.nanpercentile(d, 90, axis=(0, 3), keepdims=True)
            assert_equal(res.shape, (1, 5, 7, 1))
            res = np.nanpercentile(d, 90, axis=(1,), keepdims=True)
            assert_equal(res.shape, (3, 1, 7, 11))
            res = np.nanpercentile(d, 90, axis=(0, 1, 2, 3), keepdims=True)
            assert_equal(res.shape, (1, 1, 1, 1))
            res = np.nanpercentile(d, 90, axis=(0, 1, 3), keepdims=True)
            assert_equal(res.shape, (1, 1, 7, 1))

    def test_out(self):
        mat = np.random.rand(3, 3)
        nan_mat = np.insert(mat, [0, 2], np.nan, axis=1)
        resout = np.zeros(3)
        tgt = np.percentile(mat, 42, axis=1)
        res = np.nanpercentile(nan_mat, 42, axis=1, out=resout)
        assert_almost_equal(res, resout)
        assert_almost_equal(res, tgt)
        # 0-d output:
        resout = np.zeros(())
        tgt = np.percentile(mat, 42, axis=None)
        res = np.nanpercentile(nan_mat, 42, axis=None, out=resout)
        assert_almost_equal(res, resout)
        assert_almost_equal(res, tgt)
        res = np.nanpercentile(nan_mat, 42, axis=(0, 1), out=resout)
        assert_almost_equal(res, resout)
        assert_almost_equal(res, tgt)

    def test_result_values(self):
        tgt = [np.percentile(d, 28) for d in _rdat]
        res = np.nanpercentile(_ndat, 28, axis=1)
        assert_almost_equal(res, tgt)
        tgt = [np.percentile(d, (28, 98)) for d in _rdat]
        res = np.nanpercentile(_ndat, (28, 98), axis=1)
        assert_almost_equal(res, tgt)

    def test_allnans(self):
        mat = np.array([np.nan]*9).reshape(3, 3)
        for axis in [None, 0, 1]:
            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter('always')
                assert_(np.isnan(np.nanpercentile(mat, 60, axis=axis)).all())
                if axis is None:
                    assert_(len(w) == 1)
                else:
                    assert_(len(w) == 3)
                assert_(issubclass(w[0].category, RuntimeWarning))
                # Check scalar
                assert_(np.isnan(np.nanpercentile(np.nan, 60)))
                if axis is None:
                    assert_(len(w) == 2)
                else:
                    assert_(len(w) == 4)
                assert_(issubclass(w[0].category, RuntimeWarning))

    def test_empty(self):
        mat = np.zeros((0, 3))
        for axis in [0, None]:
            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter('always')
                assert_(np.isnan(np.nanpercentile(mat, 40, axis=axis)).all())
                assert_(len(w) == 1)
                assert_(issubclass(w[0].category, RuntimeWarning))
        for axis in [1]:
            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter('always')
                assert_equal(np.nanpercentile(mat, 40, axis=axis), np.zeros([]))
                assert_(len(w) == 0)

    def test_scalar(self):
        assert_(np.nanpercentile(0., 100) == 0.)

    def test_extended_axis_invalid(self):
        d = np.ones((3, 5, 7, 11))
        assert_raises(IndexError, np.nanpercentile, d, q=5, axis=-5)
        assert_raises(IndexError, np.nanpercentile, d, q=5, axis=(0, -5))
        assert_raises(IndexError, np.nanpercentile, d, q=5, axis=4)
        assert_raises(IndexError, np.nanpercentile, d, q=5, axis=(0, 4))
        assert_raises(ValueError, np.nanpercentile, d, q=5, axis=(1, 1))


if __name__ == "__main__":
    run_module_suite()
