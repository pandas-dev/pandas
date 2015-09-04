# -*- coding: utf-8 -*-
from __future__ import division, print_function

from functools import partial

import numpy as np
from pandas import Series
from pandas.core.common import isnull, is_integer_dtype
import pandas.core.nanops as nanops
import pandas.util.testing as tm

use_bn = nanops._USE_BOTTLENECK

class TestnanopsDataFrame(tm.TestCase):

    def setUp(self):
        np.random.seed(11235)
        nanops._USE_BOTTLENECK = False

        self.arr_shape = (11, 7, 5)

        self.arr_float = np.random.randn(*self.arr_shape)
        self.arr_float1 = np.random.randn(*self.arr_shape)
        self.arr_complex = self.arr_float + self.arr_float1*1j
        self.arr_int = np.random.randint(-10, 10, self.arr_shape)
        self.arr_bool = np.random.randint(0, 2, self.arr_shape) == 0
        self.arr_str = np.abs(self.arr_float).astype('S')
        self.arr_utf = np.abs(self.arr_float).astype('U')
        self.arr_date = np.random.randint(0, 20000,
                                          self.arr_shape).astype('M8[ns]')
        self.arr_tdelta = np.random.randint(0, 20000,
                                            self.arr_shape).astype('m8[ns]')

        self.arr_nan = np.tile(np.nan, self.arr_shape)
        self.arr_float_nan = np.vstack([self.arr_float, self.arr_nan])
        self.arr_float1_nan = np.vstack([self.arr_float1, self.arr_nan])
        self.arr_nan_float1 = np.vstack([self.arr_nan, self.arr_float1])
        self.arr_nan_nan = np.vstack([self.arr_nan, self.arr_nan])

        self.arr_inf = self.arr_float*np.inf
        self.arr_float_inf = np.vstack([self.arr_float, self.arr_inf])
        self.arr_float1_inf = np.vstack([self.arr_float1, self.arr_inf])
        self.arr_inf_float1 = np.vstack([self.arr_inf, self.arr_float1])
        self.arr_inf_inf = np.vstack([self.arr_inf, self.arr_inf])

        self.arr_nan_inf = np.vstack([self.arr_nan, self.arr_inf])
        self.arr_float_nan_inf = np.vstack([self.arr_float,
                                            self.arr_nan,
                                            self.arr_inf])
        self.arr_nan_float1_inf = np.vstack([self.arr_float,
                                             self.arr_inf,
                                             self.arr_nan])
        self.arr_nan_nan_inf = np.vstack([self.arr_nan,
                                          self.arr_nan,
                                          self.arr_inf])
        self.arr_obj = np.vstack([self.arr_float.astype('O'),
                                  self.arr_int.astype('O'),
                                  self.arr_bool.astype('O'),
                                  self.arr_complex.astype('O'),
                                  self.arr_str.astype('O'),
                                  self.arr_utf.astype('O'),
                                  self.arr_date.astype('O'),
                                  self.arr_tdelta.astype('O')])

        self.arr_nan_nanj = self.arr_nan + self.arr_nan*1j
        self.arr_complex_nan = np.vstack([self.arr_complex, self.arr_nan_nanj])

        self.arr_nan_infj = self.arr_inf*1j
        self.arr_complex_nan_infj = np.vstack([self.arr_complex,
                                              self.arr_nan_infj])

        self.arr_float_2d = self.arr_float[:, :, 0]
        self.arr_float1_2d = self.arr_float1[:, :, 0]
        self.arr_complex_2d = self.arr_complex[:, :, 0]
        self.arr_int_2d = self.arr_int[:, :, 0]
        self.arr_bool_2d = self.arr_bool[:, :, 0]
        self.arr_str_2d = self.arr_str[:, :, 0]
        self.arr_utf_2d = self.arr_utf[:, :, 0]
        self.arr_date_2d = self.arr_date[:, :, 0]
        self.arr_tdelta_2d = self.arr_tdelta[:, :, 0]

        self.arr_nan_2d = self.arr_nan[:, :, 0]
        self.arr_float_nan_2d = self.arr_float_nan[:, :, 0]
        self.arr_float1_nan_2d = self.arr_float1_nan[:, :, 0]
        self.arr_nan_float1_2d = self.arr_nan_float1[:, :, 0]
        self.arr_nan_nan_2d = self.arr_nan_nan[:, :, 0]
        self.arr_nan_nanj_2d = self.arr_nan_nanj[:, :, 0]
        self.arr_complex_nan_2d = self.arr_complex_nan[:, :, 0]

        self.arr_inf_2d = self.arr_inf[:, :, 0]
        self.arr_float_inf_2d = self.arr_float_inf[:, :, 0]
        self.arr_nan_inf_2d = self.arr_nan_inf[:, :, 0]
        self.arr_float_nan_inf_2d = self.arr_float_nan_inf[:, :, 0]
        self.arr_nan_nan_inf_2d = self.arr_nan_nan_inf[:, :, 0]

        self.arr_float_1d = self.arr_float[:, 0, 0]
        self.arr_float1_1d = self.arr_float1[:, 0, 0]
        self.arr_complex_1d = self.arr_complex[:, 0, 0]
        self.arr_int_1d = self.arr_int[:, 0, 0]
        self.arr_bool_1d = self.arr_bool[:, 0, 0]
        self.arr_str_1d = self.arr_str[:, 0, 0]
        self.arr_utf_1d = self.arr_utf[:, 0, 0]
        self.arr_date_1d = self.arr_date[:, 0, 0]
        self.arr_tdelta_1d = self.arr_tdelta[:, 0, 0]

        self.arr_nan_1d = self.arr_nan[:, 0, 0]
        self.arr_float_nan_1d = self.arr_float_nan[:, 0, 0]
        self.arr_float1_nan_1d = self.arr_float1_nan[:, 0, 0]
        self.arr_nan_float1_1d = self.arr_nan_float1[:, 0, 0]
        self.arr_nan_nan_1d = self.arr_nan_nan[:, 0, 0]
        self.arr_nan_nanj_1d = self.arr_nan_nanj[:, 0, 0]
        self.arr_complex_nan_1d = self.arr_complex_nan[:, 0, 0]

        self.arr_inf_1d = self.arr_inf.ravel()
        self.arr_float_inf_1d = self.arr_float_inf[:, 0, 0]
        self.arr_nan_inf_1d = self.arr_nan_inf[:, 0, 0]
        self.arr_float_nan_inf_1d = self.arr_float_nan_inf[:, 0, 0]
        self.arr_nan_nan_inf_1d = self.arr_nan_nan_inf[:, 0, 0]

    def tearDown(self):
        nanops._USE_BOTTLENECK = use_bn

    def check_results(self, targ, res, axis):
        res = getattr(res, 'asm8', res)
        res = getattr(res, 'values', res)

        # timedeltas are a beast here
        def _coerce_tds(targ, res):
            if targ.dtype == 'm8[ns]':
                if len(targ) == 1:
                    targ = targ[0].item()
                    res = res.item()
                else:
                    targ = targ.view('i8')
            return targ, res

        try:
            if axis != 0 and hasattr(targ, 'shape') and targ.ndim:
                res = np.split(res, [targ.shape[0]], axis=0)[0]
        except:
            targ, res = _coerce_tds(targ, res)

        try:
            tm.assert_almost_equal(targ, res)
        except:

            if targ.dtype == 'm8[ns]':
                targ, res = _coerce_tds(targ, res)
                tm.assert_almost_equal(targ, res)
                return

            # There are sometimes rounding errors with
            # complex and object dtypes.
            # If it isn't one of those, re-raise the error.
            if not hasattr(res, 'dtype') or res.dtype.kind not in ['c', 'O']:
                raise
            # convert object dtypes to something that can be split into
            # real and imaginary parts
            if res.dtype.kind == 'O':
                if targ.dtype.kind != 'O':
                    res = res.astype(targ.dtype)
                else:
                    try:
                        res = res.astype('c16')
                    except:
                        res = res.astype('f8')
                    try:
                        targ = targ.astype('c16')
                    except:
                        targ = targ.astype('f8')
            # there should never be a case where numpy returns an object
            # but nanops doesn't, so make that an exception
            elif targ.dtype.kind == 'O':
                raise
            tm.assert_almost_equal(targ.real, res.real)
            tm.assert_almost_equal(targ.imag, res.imag)

    def check_fun_data(self, testfunc, targfunc,
                       testarval, targarval, targarnanval, **kwargs):
        for axis in list(range(targarval.ndim))+[None]:
            for skipna in [False, True]:
                targartempval = targarval if skipna else targarnanval
                try:
                    targ = targfunc(targartempval, axis=axis, **kwargs)
                    res = testfunc(testarval, axis=axis, skipna=skipna,
                                   **kwargs)
                    self.check_results(targ, res, axis)
                    if skipna:
                        res = testfunc(testarval, axis=axis, **kwargs)
                        self.check_results(targ, res, axis)
                    if axis is None:
                        res = testfunc(testarval, skipna=skipna, **kwargs)
                        self.check_results(targ, res, axis)
                    if skipna and axis is None:
                        res = testfunc(testarval, **kwargs)
                        self.check_results(targ, res, axis)
                except BaseException as exc:
                    exc.args += ('axis: %s of %s' % (axis, testarval.ndim-1),
                                 'skipna: %s' % skipna,
                                 'kwargs: %s' % kwargs)
                    raise

        if testarval.ndim <= 1:
            return

        try:
            testarval2 = np.take(testarval, 0, axis=-1)
            targarval2 = np.take(targarval, 0, axis=-1)
            targarnanval2 = np.take(targarnanval, 0, axis=-1)
        except ValueError:
            return
        self.check_fun_data(testfunc, targfunc,
                            testarval2, targarval2, targarnanval2,
                            **kwargs)

    def check_fun(self, testfunc, targfunc,
                  testar, targar=None, targarnan=None,
                  **kwargs):
        if targar is None:
            targar = testar
        if targarnan is None:
            targarnan = testar
        testarval = getattr(self, testar)
        targarval = getattr(self, targar)
        targarnanval = getattr(self, targarnan)
        try:
            self.check_fun_data(testfunc, targfunc,
                                testarval, targarval, targarnanval, **kwargs)
        except BaseException as exc:
            exc.args += ('testar: %s' % testar,
                         'targar: %s' % targar,
                         'targarnan: %s' % targarnan)
            raise

    def check_funs(self, testfunc, targfunc,
                   allow_complex=True, allow_all_nan=True, allow_str=True,
                   allow_date=True, allow_tdelta=True, allow_obj=True,
                   **kwargs):
        self.check_fun(testfunc, targfunc, 'arr_float', **kwargs)
        self.check_fun(testfunc, targfunc, 'arr_float_nan', 'arr_float',
                       **kwargs)
        self.check_fun(testfunc, targfunc, 'arr_int', **kwargs)
        self.check_fun(testfunc, targfunc, 'arr_bool', **kwargs)
        objs = [self.arr_float.astype('O'),
                self.arr_int.astype('O'),
                self.arr_bool.astype('O')]

        if allow_all_nan:
            self.check_fun(testfunc, targfunc, 'arr_nan', **kwargs)

        if allow_complex:
            self.check_fun(testfunc, targfunc, 'arr_complex', **kwargs)
            self.check_fun(testfunc, targfunc,
                           'arr_complex_nan', 'arr_complex', **kwargs)
            if allow_all_nan:
                self.check_fun(testfunc, targfunc, 'arr_nan_nanj', **kwargs)
            objs += [self.arr_complex.astype('O')]

        if allow_str:
            self.check_fun(testfunc, targfunc, 'arr_str', **kwargs)
            self.check_fun(testfunc, targfunc, 'arr_utf', **kwargs)
            objs += [self.arr_str.astype('O'),
                     self.arr_utf.astype('O')]

        if allow_date:
            try:
                targfunc(self.arr_date)
            except TypeError:
                pass
            else:
                self.check_fun(testfunc, targfunc, 'arr_date', **kwargs)
                objs += [self.arr_date.astype('O')]

        if allow_tdelta:
            try:
                targfunc(self.arr_tdelta)
            except TypeError:
                pass
            else:
                self.check_fun(testfunc, targfunc, 'arr_tdelta', **kwargs)
                objs += [self.arr_tdelta.astype('O')]

        if allow_obj:
            self.arr_obj = np.vstack(objs)
            # some nanops handle object dtypes better than their numpy
            # counterparts, so the numpy functions need to be given something
            # else
            if allow_obj == 'convert':
                targfunc = partial(self._badobj_wrap,
                                   func=targfunc, allow_complex=allow_complex)
            self.check_fun(testfunc, targfunc, 'arr_obj', **kwargs)

    def check_funs_ddof(self, testfunc, targfunc,
                        allow_complex=True, allow_all_nan=True, allow_str=True,
                        allow_date=False, allow_tdelta=False, allow_obj=True,):
        for ddof in range(3):
            try:
                self.check_funs(testfunc, targfunc,
                                allow_complex, allow_all_nan, allow_str,
                                allow_date, allow_tdelta, allow_obj,
                                ddof=ddof)
            except BaseException as exc:
                exc.args += ('ddof %s' % ddof,)
                raise

    def _badobj_wrap(self, value, func, allow_complex=True, **kwargs):
        if value.dtype.kind == 'O':
            if allow_complex:
                value = value.astype('c16')
            else:
                value = value.astype('f8')
        return func(value, **kwargs)

    def test_nanany(self):
        self.check_funs(nanops.nanany, np.any,
                        allow_all_nan=False, allow_str=False, allow_date=False, allow_tdelta=False)

    def test_nanall(self):
        self.check_funs(nanops.nanall, np.all,
                        allow_all_nan=False, allow_str=False, allow_date=False, allow_tdelta=False)

    def test_nansum(self):
        self.check_funs(nanops.nansum, np.sum,
                        allow_str=False, allow_date=False, allow_tdelta=True)

    def test_nanmean(self):
        self.check_funs(nanops.nanmean, np.mean,
                        allow_complex=False, allow_obj=False,
                        allow_str=False, allow_date=False, allow_tdelta=True)

    def test_nanmean_overflow(self):
        # GH 10155
        # In the previous implementation mean can overflow for int dtypes, it
        # is now consistent with numpy

        # numpy < 1.9.0 is not computing this correctly
        from distutils.version import LooseVersion
        if LooseVersion(np.__version__) >= '1.9.0':
            for a in [2 ** 55, -2 ** 55, 20150515061816532]:
                s = Series(a, index=range(500), dtype=np.int64)
                result = s.mean()
                np_result = s.values.mean()
                self.assertEqual(result, a)
                self.assertEqual(result, np_result)
                self.assertTrue(result.dtype == np.float64)

    def test_returned_dtype(self):

        dtypes = [np.int16, np.int32, np.int64, np.float32, np.float64]
        if hasattr(np,'float128'):
            dtypes.append(np.float128)

        for dtype in dtypes:
            s = Series(range(10), dtype=dtype)
            group_a = ['mean', 'std', 'var', 'skew', 'kurt']
            group_b = ['min', 'max']
            for method in group_a + group_b:
                result = getattr(s, method)()
                if is_integer_dtype(dtype) and method in group_a:
                    self.assertTrue(result.dtype == np.float64,
                                    "return dtype expected from %s is np.float64, got %s instead" % (method, result.dtype))
                else:
                    self.assertTrue(result.dtype == dtype,
                                    "return dtype expected from %s is %s, got %s instead" % (method, dtype, result.dtype))

    def test_nanmedian(self):
        self.check_funs(nanops.nanmedian, np.median,
                        allow_complex=False, allow_str=False, allow_date=False,
                        allow_tdelta=True,
                        allow_obj='convert')

    def test_nanvar(self):
        self.check_funs_ddof(nanops.nanvar, np.var,
                             allow_complex=False,
                             allow_str=False,
                             allow_date=False,
                             allow_tdelta=True,
                             allow_obj='convert')

    def test_nanstd(self):
        self.check_funs_ddof(nanops.nanstd, np.std,
                             allow_complex=False,
                             allow_str=False,
                             allow_date=False,
                             allow_tdelta=True,
                             allow_obj='convert')

    def test_nansem(self):
        tm.skip_if_no_package('scipy.stats')
        from scipy.stats import sem
        self.check_funs_ddof(nanops.nansem, sem,
                             allow_complex=False,
                             allow_str=False,
                             allow_date=False,
                             allow_tdelta=True,
                             allow_obj='convert')

    def _minmax_wrap(self, value, axis=None, func=None):
        res = func(value, axis)
        if res.dtype.kind == 'm':
            res = np.atleast_1d(res)
        return res

    def test_nanmin(self):
        func = partial(self._minmax_wrap, func=np.min)
        self.check_funs(nanops.nanmin, func,
                        allow_str=False, allow_obj=False)

    def test_nanmax(self):
        func = partial(self._minmax_wrap, func=np.max)
        self.check_funs(nanops.nanmax, func,
                        allow_str=False, allow_obj=False)

    def _argminmax_wrap(self, value, axis=None, func=None):
        res = func(value, axis)
        nans = np.min(value, axis)
        nullnan = isnull(nans)
        if res.ndim:
            res[nullnan] = -1
        elif (hasattr(nullnan, 'all') and nullnan.all() or
              not hasattr(nullnan, 'all') and nullnan):
            res = -1
        return res

    def test_nanargmax(self):
        func = partial(self._argminmax_wrap, func=np.argmax)
        self.check_funs(nanops.nanargmax, func,
                        allow_str=False, allow_obj=False,
                        allow_date=True,
                        allow_tdelta=True)

    def test_nanargmin(self):
        func = partial(self._argminmax_wrap, func=np.argmin)
        if tm.sys.version_info[0:2] == (2, 6):
            self.check_funs(nanops.nanargmin, func,
                            allow_date=True,
                            allow_tdelta=True,
                            allow_str=False, allow_obj=False)
        else:
            self.check_funs(nanops.nanargmin, func,
                            allow_str=False, allow_obj=False)

    def _skew_kurt_wrap(self, values, axis=None, func=None):
        if not isinstance(values.dtype.type, np.floating):
            values = values.astype('f8')
        result = func(values, axis=axis, bias=False)
        # fix for handling cases where all elements in an axis are the same
        if isinstance(result, np.ndarray):
            result[np.max(values, axis=axis) == np.min(values, axis=axis)] = 0
            return result
        elif np.max(values) == np.min(values):
            return 0.
        return result

    def test_nanskew(self):
        tm.skip_if_no_package('scipy.stats')
        from scipy.stats import skew
        func = partial(self._skew_kurt_wrap, func=skew)
        self.check_funs(nanops.nanskew, func,
                        allow_complex=False, allow_str=False, allow_date=False, allow_tdelta=False)

    def test_nankurt(self):
        tm.skip_if_no_package('scipy.stats')
        from scipy.stats import kurtosis
        func1 = partial(kurtosis, fisher=True)
        func = partial(self._skew_kurt_wrap, func=func1)
        self.check_funs(nanops.nankurt, func,
                        allow_complex=False, allow_str=False, allow_date=False, allow_tdelta=False)

    def test_nanprod(self):
        self.check_funs(nanops.nanprod, np.prod,
                        allow_str=False, allow_date=False, allow_tdelta=False)

    def check_nancorr_nancov_2d(self, checkfun, targ0, targ1, **kwargs):
        res00 = checkfun(self.arr_float_2d, self.arr_float1_2d,
                         **kwargs)
        res01 = checkfun(self.arr_float_2d, self.arr_float1_2d,
                         min_periods=len(self.arr_float_2d)-1,
                         **kwargs)
        tm.assert_almost_equal(targ0, res00)
        tm.assert_almost_equal(targ0, res01)

        res10 = checkfun(self.arr_float_nan_2d, self.arr_float1_nan_2d,
                         **kwargs)
        res11 = checkfun(self.arr_float_nan_2d, self.arr_float1_nan_2d,
                         min_periods=len(self.arr_float_2d)-1,
                         **kwargs)
        tm.assert_almost_equal(targ1, res10)
        tm.assert_almost_equal(targ1, res11)

        targ2 = np.nan
        res20 = checkfun(self.arr_nan_2d, self.arr_float1_2d,
                         **kwargs)
        res21 = checkfun(self.arr_float_2d, self.arr_nan_2d,
                         **kwargs)
        res22 = checkfun(self.arr_nan_2d, self.arr_nan_2d,
                         **kwargs)
        res23 = checkfun(self.arr_float_nan_2d, self.arr_nan_float1_2d,
                         **kwargs)
        res24 = checkfun(self.arr_float_nan_2d, self.arr_nan_float1_2d,
                         min_periods=len(self.arr_float_2d)-1,
                         **kwargs)
        res25 = checkfun(self.arr_float_2d, self.arr_float1_2d,
                         min_periods=len(self.arr_float_2d)+1,
                         **kwargs)
        tm.assert_almost_equal(targ2, res20)
        tm.assert_almost_equal(targ2, res21)
        tm.assert_almost_equal(targ2, res22)
        tm.assert_almost_equal(targ2, res23)
        tm.assert_almost_equal(targ2, res24)
        tm.assert_almost_equal(targ2, res25)

    def check_nancorr_nancov_1d(self, checkfun, targ0, targ1, **kwargs):
        res00 = checkfun(self.arr_float_1d, self.arr_float1_1d,
                         **kwargs)
        res01 = checkfun(self.arr_float_1d, self.arr_float1_1d,
                         min_periods=len(self.arr_float_1d)-1,
                         **kwargs)
        tm.assert_almost_equal(targ0, res00)
        tm.assert_almost_equal(targ0, res01)

        res10 = checkfun(self.arr_float_nan_1d,
                         self.arr_float1_nan_1d,
                         **kwargs)
        res11 = checkfun(self.arr_float_nan_1d,
                         self.arr_float1_nan_1d,
                         min_periods=len(self.arr_float_1d)-1,
                         **kwargs)
        tm.assert_almost_equal(targ1, res10)
        tm.assert_almost_equal(targ1, res11)

        targ2 = np.nan
        res20 = checkfun(self.arr_nan_1d, self.arr_float1_1d,
                         **kwargs)
        res21 = checkfun(self.arr_float_1d, self.arr_nan_1d,
                         **kwargs)
        res22 = checkfun(self.arr_nan_1d, self.arr_nan_1d,
                         **kwargs)
        res23 = checkfun(self.arr_float_nan_1d,
                         self.arr_nan_float1_1d,
                         **kwargs)
        res24 = checkfun(self.arr_float_nan_1d,
                         self.arr_nan_float1_1d,
                         min_periods=len(self.arr_float_1d)-1,
                         **kwargs)
        res25 = checkfun(self.arr_float_1d,
                         self.arr_float1_1d,
                         min_periods=len(self.arr_float_1d)+1,
                         **kwargs)
        tm.assert_almost_equal(targ2, res20)
        tm.assert_almost_equal(targ2, res21)
        tm.assert_almost_equal(targ2, res22)
        tm.assert_almost_equal(targ2, res23)
        tm.assert_almost_equal(targ2, res24)
        tm.assert_almost_equal(targ2, res25)

    def test_nancorr(self):
        targ0 = np.corrcoef(self.arr_float_2d, self.arr_float1_2d)[0, 1]
        targ1 = np.corrcoef(self.arr_float_2d.flat,
                            self.arr_float1_2d.flat)[0, 1]
        self.check_nancorr_nancov_2d(nanops.nancorr, targ0, targ1)
        targ0 = np.corrcoef(self.arr_float_1d, self.arr_float1_1d)[0, 1]
        targ1 = np.corrcoef(self.arr_float_1d.flat,
                            self.arr_float1_1d.flat)[0, 1]
        self.check_nancorr_nancov_1d(nanops.nancorr, targ0, targ1,
                                     method='pearson')

    def test_nancorr_pearson(self):
        targ0 = np.corrcoef(self.arr_float_2d, self.arr_float1_2d)[0, 1]
        targ1 = np.corrcoef(self.arr_float_2d.flat,
                            self.arr_float1_2d.flat)[0, 1]
        self.check_nancorr_nancov_2d(nanops.nancorr, targ0, targ1,
                                     method='pearson')
        targ0 = np.corrcoef(self.arr_float_1d, self.arr_float1_1d)[0, 1]
        targ1 = np.corrcoef(self.arr_float_1d.flat,
                            self.arr_float1_1d.flat)[0, 1]
        self.check_nancorr_nancov_1d(nanops.nancorr, targ0, targ1,
                                     method='pearson')

    def test_nancorr_kendall(self):
        tm.skip_if_no_package('scipy.stats')
        from scipy.stats import kendalltau
        targ0 = kendalltau(self.arr_float_2d, self.arr_float1_2d)[0]
        targ1 = kendalltau(self.arr_float_2d.flat, self.arr_float1_2d.flat)[0]
        self.check_nancorr_nancov_2d(nanops.nancorr, targ0, targ1,
                                     method='kendall')
        targ0 = kendalltau(self.arr_float_1d, self.arr_float1_1d)[0]
        targ1 = kendalltau(self.arr_float_1d.flat, self.arr_float1_1d.flat)[0]
        self.check_nancorr_nancov_1d(nanops.nancorr, targ0, targ1,
                                     method='kendall')

    def test_nancorr_spearman(self):
        tm.skip_if_no_package('scipy.stats')
        from scipy.stats import spearmanr
        targ0 = spearmanr(self.arr_float_2d, self.arr_float1_2d)[0]
        targ1 = spearmanr(self.arr_float_2d.flat, self.arr_float1_2d.flat)[0]
        self.check_nancorr_nancov_2d(nanops.nancorr, targ0, targ1,
                                     method='spearman')
        targ0 = spearmanr(self.arr_float_1d, self.arr_float1_1d)[0]
        targ1 = spearmanr(self.arr_float_1d.flat, self.arr_float1_1d.flat)[0]
        self.check_nancorr_nancov_1d(nanops.nancorr, targ0, targ1,
                                     method='spearman')

    def test_nancov(self):
        targ0 = np.cov(self.arr_float_2d, self.arr_float1_2d)[0, 1]
        targ1 = np.cov(self.arr_float_2d.flat, self.arr_float1_2d.flat)[0, 1]
        self.check_nancorr_nancov_2d(nanops.nancov, targ0, targ1)
        targ0 = np.cov(self.arr_float_1d, self.arr_float1_1d)[0, 1]
        targ1 = np.cov(self.arr_float_1d.flat, self.arr_float1_1d.flat)[0, 1]
        self.check_nancorr_nancov_1d(nanops.nancov, targ0, targ1)

    def check_nancomp(self, checkfun, targ0):
        arr_float = self.arr_float
        arr_float1 = self.arr_float1
        arr_nan = self.arr_nan
        arr_nan_nan = self.arr_nan_nan
        arr_float_nan = self.arr_float_nan
        arr_float1_nan = self.arr_float1_nan
        arr_nan_float1 = self.arr_nan_float1

        while targ0.ndim:
            try:
                res0 = checkfun(arr_float, arr_float1)
                tm.assert_almost_equal(targ0, res0)

                if targ0.ndim > 1:
                    targ1 = np.vstack([targ0, arr_nan])
                else:
                    targ1 = np.hstack([targ0, arr_nan])
                res1 = checkfun(arr_float_nan, arr_float1_nan)
                tm.assert_almost_equal(targ1, res1)

                targ2 = arr_nan_nan
                res2 = checkfun(arr_float_nan, arr_nan_float1)
                tm.assert_almost_equal(targ2, res2)
            except Exception as exc:
                exc.args += ('ndim: %s' % arr_float.ndim,)
                raise

            try:
                arr_float = np.take(arr_float, 0, axis=-1)
                arr_float1 = np.take(arr_float1, 0, axis=-1)
                arr_nan = np.take(arr_nan, 0, axis=-1)
                arr_nan_nan = np.take(arr_nan_nan, 0, axis=-1)
                arr_float_nan = np.take(arr_float_nan, 0, axis=-1)
                arr_float1_nan = np.take(arr_float1_nan, 0, axis=-1)
                arr_nan_float1 = np.take(arr_nan_float1, 0, axis=-1)
                targ0 = np.take(targ0, 0, axis=-1)
            except ValueError:
                break

    def test_nangt(self):
        targ0 = self.arr_float > self.arr_float1
        self.check_nancomp(nanops.nangt, targ0)

    def test_nange(self):
        targ0 = self.arr_float >= self.arr_float1
        self.check_nancomp(nanops.nange, targ0)

    def test_nanlt(self):
        targ0 = self.arr_float < self.arr_float1
        self.check_nancomp(nanops.nanlt, targ0)

    def test_nanle(self):
        targ0 = self.arr_float <= self.arr_float1
        self.check_nancomp(nanops.nanle, targ0)

    def test_naneq(self):
        targ0 = self.arr_float == self.arr_float1
        self.check_nancomp(nanops.naneq, targ0)

    def test_nanne(self):
        targ0 = self.arr_float != self.arr_float1
        self.check_nancomp(nanops.nanne, targ0)

    def check_bool(self, func, value, correct, *args, **kwargs):
        while getattr(value, 'ndim', True):
            try:
                res0 = func(value, *args, **kwargs)
                if correct:
                    self.assertTrue(res0)
                else:
                    self.assertFalse(res0)
            except BaseException as exc:
                exc.args += ('dim: %s' % getattr(value, 'ndim', value),)
                raise
            if not hasattr(value, 'ndim'):
                break
            try:
                value = np.take(value, 0, axis=-1)
            except ValueError:
                break

    def test__has_infs(self):
        pairs = [('arr_complex', False),
                 ('arr_int', False),
                 ('arr_bool', False),
                 ('arr_str', False),
                 ('arr_utf', False),
                 ('arr_complex', False),
                 ('arr_complex_nan', False),

                 ('arr_nan_nanj', False),
                 ('arr_nan_infj', True),
                 ('arr_complex_nan_infj', True)]
        pairs_float = [('arr_float', False),
                       ('arr_nan', False),
                       ('arr_float_nan', False),
                       ('arr_nan_nan', False),

                       ('arr_float_inf', True),
                       ('arr_inf', True),
                       ('arr_nan_inf', True),
                       ('arr_float_nan_inf', True),
                       ('arr_nan_nan_inf', True)]

        for arr, correct in pairs:
            val = getattr(self, arr)
            try:
                self.check_bool(nanops._has_infs, val, correct)
            except BaseException as exc:
                exc.args += (arr,)
                raise

        for arr, correct in pairs_float:
            val = getattr(self, arr)
            try:
                self.check_bool(nanops._has_infs, val, correct)
                self.check_bool(nanops._has_infs, val.astype('f4'), correct)
                self.check_bool(nanops._has_infs, val.astype('f2'), correct)
            except BaseException as exc:
                exc.args += (arr,)
                raise

    def test__isfinite(self):
        pairs = [('arr_complex', False),
                 ('arr_int', False),
                 ('arr_bool', False),
                 ('arr_str', False),
                 ('arr_utf', False),
                 ('arr_complex', False),
                 ('arr_complex_nan', True),

                 ('arr_nan_nanj', True),
                 ('arr_nan_infj', True),
                 ('arr_complex_nan_infj', True)]
        pairs_float = [('arr_float', False),
                       ('arr_nan', True),
                       ('arr_float_nan', True),
                       ('arr_nan_nan', True),

                       ('arr_float_inf', True),
                       ('arr_inf', True),
                       ('arr_nan_inf', True),
                       ('arr_float_nan_inf', True),
                       ('arr_nan_nan_inf', True)]

        func1 = lambda x: np.any(nanops._isfinite(x).ravel())
        func2 = lambda x: np.any(nanops._isfinite(x).values.ravel())
        for arr, correct in pairs:
            val = getattr(self, arr)
            try:
                self.check_bool(func1, val, correct)
            except BaseException as exc:
                exc.args += (arr,)
                raise

        for arr, correct in pairs_float:
            val = getattr(self, arr)
            try:
                self.check_bool(func1, val, correct)
                self.check_bool(func1, val.astype('f4'), correct)
                self.check_bool(func1, val.astype('f2'), correct)
            except BaseException as exc:
                exc.args += (arr,)
                raise

    def test__bn_ok_dtype(self):
        self.assertTrue(nanops._bn_ok_dtype(self.arr_float.dtype, 'test'))
        self.assertTrue(nanops._bn_ok_dtype(self.arr_complex.dtype, 'test'))
        self.assertTrue(nanops._bn_ok_dtype(self.arr_int.dtype, 'test'))
        self.assertTrue(nanops._bn_ok_dtype(self.arr_bool.dtype, 'test'))
        self.assertTrue(nanops._bn_ok_dtype(self.arr_str.dtype, 'test'))
        self.assertTrue(nanops._bn_ok_dtype(self.arr_utf.dtype, 'test'))
        self.assertFalse(nanops._bn_ok_dtype(self.arr_date.dtype, 'test'))
        self.assertFalse(nanops._bn_ok_dtype(self.arr_tdelta.dtype, 'test'))
        self.assertFalse(nanops._bn_ok_dtype(self.arr_obj.dtype, 'test'))


class TestEnsureNumeric(tm.TestCase):
    def test_numeric_values(self):
        # Test integer
        self.assertEqual(nanops._ensure_numeric(1), 1, 'Failed for int')
        # Test float
        self.assertEqual(nanops._ensure_numeric(1.1), 1.1, 'Failed for float')
        # Test complex
        self.assertEqual(nanops._ensure_numeric(1 + 2j), 1 + 2j,
                         'Failed for complex')

    def test_ndarray(self):
        # Test numeric ndarray
        values = np.array([1, 2, 3])
        self.assertTrue(np.allclose(nanops._ensure_numeric(values), values),
                        'Failed for numeric ndarray')

        # Test object ndarray
        o_values = values.astype(object)
        self.assertTrue(np.allclose(nanops._ensure_numeric(o_values), values),
                        'Failed for object ndarray')

        # Test convertible string ndarray
        s_values = np.array(['1', '2', '3'], dtype=object)
        self.assertTrue(np.allclose(nanops._ensure_numeric(s_values), values),
                        'Failed for convertible string ndarray')

        # Test non-convertible string ndarray
        s_values = np.array(['foo', 'bar', 'baz'], dtype=object)
        self.assertRaises(ValueError,
                          lambda: nanops._ensure_numeric(s_values))

    def test_convertable_values(self):
        self.assertTrue(np.allclose(nanops._ensure_numeric('1'), 1.0),
                        'Failed for convertible integer string')
        self.assertTrue(np.allclose(nanops._ensure_numeric('1.1'), 1.1),
                        'Failed for convertible float string')
        self.assertTrue(np.allclose(nanops._ensure_numeric('1+1j'), 1 + 1j),
                        'Failed for convertible complex string')

    def test_non_convertable_values(self):
        self.assertRaises(TypeError,
                          lambda: nanops._ensure_numeric('foo'))
        self.assertRaises(TypeError,
                          lambda: nanops._ensure_numeric({}))
        self.assertRaises(TypeError,
                          lambda: nanops._ensure_numeric([]))


class TestNanvarFixedValues(tm.TestCase):

    # xref GH10242

    def setUp(self):
        # Samples from a normal distribution.
        self.variance = variance = 3.0
        self.samples = self.prng.normal(scale=variance ** 0.5, size=100000)

    def test_nanvar_all_finite(self):
        samples = self.samples
        actual_variance = nanops.nanvar(samples)
        np.testing.assert_almost_equal(
            actual_variance, self.variance, decimal=2)

    def test_nanvar_nans(self):
        samples = np.nan * np.ones(2 * self.samples.shape[0])
        samples[::2] = self.samples

        actual_variance = nanops.nanvar(samples, skipna=True)
        np.testing.assert_almost_equal(
            actual_variance, self.variance, decimal=2)

        actual_variance = nanops.nanvar(samples, skipna=False)
        np.testing.assert_almost_equal(
            actual_variance, np.nan, decimal=2)

    def test_nanstd_nans(self):
        samples = np.nan * np.ones(2 * self.samples.shape[0])
        samples[::2] = self.samples

        actual_std = nanops.nanstd(samples, skipna=True)
        np.testing.assert_almost_equal(
            actual_std, self.variance ** 0.5, decimal=2)

        actual_std = nanops.nanvar(samples, skipna=False)
        np.testing.assert_almost_equal(
            actual_std, np.nan, decimal=2)

    def test_nanvar_axis(self):
        # Generate some sample data.
        samples_norm = self.samples
        samples_unif = self.prng.uniform(size=samples_norm.shape[0])
        samples = np.vstack([samples_norm, samples_unif])

        actual_variance = nanops.nanvar(samples, axis=1)
        np.testing.assert_array_almost_equal(
            actual_variance, np.array([self.variance, 1.0 / 12]), decimal=2)

    def test_nanvar_ddof(self):
        n = 5
        samples = self.prng.uniform(size=(10000, n+1))
        samples[:, -1] = np.nan  # Force use of our own algorithm.

        variance_0 = nanops.nanvar(samples, axis=1, skipna=True, ddof=0).mean()
        variance_1 = nanops.nanvar(samples, axis=1, skipna=True, ddof=1).mean()
        variance_2 = nanops.nanvar(samples, axis=1, skipna=True, ddof=2).mean()

        # The unbiased estimate.
        var = 1.0 / 12
        np.testing.assert_almost_equal(variance_1, var, decimal=2)
        # The underestimated variance.
        np.testing.assert_almost_equal(
            variance_0,  (n - 1.0) / n * var, decimal=2)
        # The overestimated variance.
        np.testing.assert_almost_equal(
            variance_2,  (n - 1.0) / (n - 2.0) * var, decimal=2)

    def test_ground_truth(self):
        # Test against values that were precomputed with Numpy.
        samples = np.empty((4, 4))
        samples[:3, :3] = np.array([[0.97303362, 0.21869576, 0.55560287],
                                    [0.72980153, 0.03109364, 0.99155171],
                                    [0.09317602, 0.60078248, 0.15871292]])
        samples[3] = samples[:, 3] = np.nan

        # Actual variances along axis=0, 1 for ddof=0, 1, 2
        variance = np.array(
            [[[0.13762259, 0.05619224, 0.11568816],
              [0.20643388, 0.08428837, 0.17353224],
              [0.41286776, 0.16857673, 0.34706449]],
             [[0.09519783, 0.16435395, 0.05082054],
              [0.14279674, 0.24653093, 0.07623082],
              [0.28559348, 0.49306186, 0.15246163]]]
        )

        # Test nanvar.
        for axis in range(2):
            for ddof in range(3):
                var = nanops.nanvar(samples, skipna=True, axis=axis, ddof=ddof)
                np.testing.assert_array_almost_equal(
                    var[:3], variance[axis, ddof]
                )
                np.testing.assert_equal(var[3], np.nan)

        # Test nanstd.
        for axis in range(2):
            for ddof in range(3):
                std = nanops.nanstd(samples, skipna=True, axis=axis, ddof=ddof)
                np.testing.assert_array_almost_equal(
                    std[:3], variance[axis, ddof] ** 0.5
                )
                np.testing.assert_equal(std[3], np.nan)

    def test_nanstd_roundoff(self):
        # Regression test for GH 10242 (test data taken from GH 10489). Ensure
        # that variance is stable.
        data = Series(766897346 * np.ones(10))
        for ddof in range(3):
            result = data.std(ddof=ddof)
            self.assertEqual(result, 0.0)

    @property
    def prng(self):
        return np.random.RandomState(1234)


if __name__ == '__main__':
    import nose
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure',
                         '-s'], exit=False)
