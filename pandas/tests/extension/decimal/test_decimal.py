import operator
import decimal

import numpy as np
import pandas as pd
import pandas.util.testing as tm
import pytest

from pandas.tests.extension import base

from .array import DecimalDtype, DecimalArray, make_data


@pytest.fixture
def dtype():
    return DecimalDtype()


@pytest.fixture
def data():
    return DecimalArray(make_data())


@pytest.fixture
def data_missing():
    return DecimalArray([decimal.Decimal('NaN'), decimal.Decimal(1)])


@pytest.fixture
def data_for_sorting():
    return DecimalArray([decimal.Decimal('1'),
                         decimal.Decimal('2'),
                         decimal.Decimal('0')])


@pytest.fixture
def data_missing_for_sorting():
    return DecimalArray([decimal.Decimal('1'),
                         decimal.Decimal('NaN'),
                         decimal.Decimal('0')])


@pytest.fixture
def na_cmp():
    return lambda x, y: x.is_nan() and y.is_nan()


@pytest.fixture
def na_value():
    return decimal.Decimal("NaN")


@pytest.fixture
def data_for_grouping():
    b = decimal.Decimal('1.0')
    a = decimal.Decimal('0.0')
    c = decimal.Decimal('2.0')
    na = decimal.Decimal('NaN')
    return DecimalArray([b, b, na, na, a, a, b, c])


class BaseDecimal(object):

    def assert_series_equal(self, left, right, *args, **kwargs):

        left_na = left.isna()
        right_na = right.isna()

        tm.assert_series_equal(left_na, right_na)
        return tm.assert_series_equal(left[~left_na],
                                      right[~right_na],
                                      *args, **kwargs)

    def assert_frame_equal(self, left, right, *args, **kwargs):
        # TODO(EA): select_dtypes
        tm.assert_index_equal(
            left.columns, right.columns,
            exact=kwargs.get('check_column_type', 'equiv'),
            check_names=kwargs.get('check_names', True),
            check_exact=kwargs.get('check_exact', False),
            check_categorical=kwargs.get('check_categorical', True),
            obj='{obj}.columns'.format(obj=kwargs.get('obj', 'DataFrame')))

        decimals = (left.dtypes == 'decimal').index

        for col in decimals:
            self.assert_series_equal(left[col], right[col],
                                     *args, **kwargs)

        left = left.drop(columns=decimals)
        right = right.drop(columns=decimals)
        tm.assert_frame_equal(left, right, *args, **kwargs)


class TestDtype(BaseDecimal, base.BaseDtypeTests):
    pass


class TestInterface(BaseDecimal, base.BaseInterfaceTests):
    pass


class TestConstructors(BaseDecimal, base.BaseConstructorsTests):

    @pytest.mark.xfail(reason="not implemented constructor from dtype")
    def test_from_dtype(self, data):
        # construct from our dtype & string dtype
        pass


class TestReshaping(BaseDecimal, base.BaseReshapingTests):
    pass


class TestGetitem(BaseDecimal, base.BaseGetitemTests):

    def test_take_na_value_other_decimal(self):
        arr = DecimalArray([decimal.Decimal('1.0'),
                            decimal.Decimal('2.0')])
        result = arr.take([0, -1], allow_fill=True,
                          fill_value=decimal.Decimal('-1.0'))
        expected = DecimalArray([decimal.Decimal('1.0'),
                                 decimal.Decimal('-1.0')])
        self.assert_extension_array_equal(result, expected)


class TestMissing(BaseDecimal, base.BaseMissingTests):
    pass


class TestMethods(BaseDecimal, base.BaseMethodsTests):
    @pytest.mark.parametrize('dropna', [True, False])
    @pytest.mark.xfail(reason="value_counts not implemented yet.")
    def test_value_counts(self, all_data, dropna):
        all_data = all_data[:10]
        if dropna:
            other = np.array(all_data[~all_data.isna()])
        else:
            other = all_data

        result = pd.Series(all_data).value_counts(dropna=dropna).sort_index()
        expected = pd.Series(other).value_counts(dropna=dropna).sort_index()

        tm.assert_series_equal(result, expected)


class TestCasting(BaseDecimal, base.BaseCastingTests):
    pass


class TestGroupby(BaseDecimal, base.BaseGroupbyTests):
    pass


class TestSetitem(BaseDecimal, base.BaseSetitemTests):
    pass


# TODO(extension)
@pytest.mark.xfail(reason=(
    "raising AssertionError as this is not implemented, "
    "though easy enough to do"))
def test_series_constructor_coerce_data_to_extension_dtype_raises():
    xpr = ("Cannot cast data to extension dtype 'decimal'. Pass the "
           "extension array directly.")
    with tm.assert_raises_regex(ValueError, xpr):
        pd.Series([0, 1, 2], dtype=DecimalDtype())


def test_series_constructor_with_dtype():
    arr = DecimalArray([decimal.Decimal('10.0')])
    result = pd.Series(arr, dtype=DecimalDtype())
    expected = pd.Series(arr)
    tm.assert_series_equal(result, expected)

    result = pd.Series(arr, dtype='int64')
    expected = pd.Series([10])
    tm.assert_series_equal(result, expected)


def test_dataframe_constructor_with_dtype():
    arr = DecimalArray([decimal.Decimal('10.0')])

    result = pd.DataFrame({"A": arr}, dtype=DecimalDtype())
    expected = pd.DataFrame({"A": arr})
    tm.assert_frame_equal(result, expected)

    arr = DecimalArray([decimal.Decimal('10.0')])
    result = pd.DataFrame({"A": arr}, dtype='int64')
    expected = pd.DataFrame({"A": [10]})
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize("frame", [True, False])
def test_astype_dispatches(frame):
    # This is a dtype-specific test that ensures Series[decimal].astype
    # gets all the way through to ExtensionArray.astype
    # Designing a reliable smoke test that works for arbitrary data types
    # is difficult.
    data = pd.Series(DecimalArray([decimal.Decimal(2)]), name='a')
    ctx = decimal.Context()
    ctx.prec = 5

    if frame:
        data = data.to_frame()

    result = data.astype(DecimalDtype(ctx))

    if frame:
        result = result['a']

    assert result.dtype.context.prec == ctx.prec


class TestArithmeticOps(BaseDecimal, base.BaseArithmeticOpsTests):

    def check_opname(self, s, op_name, other, exc=None):
        super(TestArithmeticOps, self).check_opname(s, op_name,
                                                    other, exc=None)

    def test_arith_series_with_array(self, data, all_arithmetic_operators):
        op_name = all_arithmetic_operators
        s = pd.Series(data)

        context = decimal.getcontext()
        divbyzerotrap = context.traps[decimal.DivisionByZero]
        invalidoptrap = context.traps[decimal.InvalidOperation]
        context.traps[decimal.DivisionByZero] = 0
        context.traps[decimal.InvalidOperation] = 0

        # Decimal supports ops with int, but not float
        other = pd.Series([int(d * 100) for d in data])
        self.check_opname(s, op_name, other)

        if "mod" not in op_name:
            self.check_opname(s, op_name, s * 2)

        self.check_opname(s, op_name, 0)
        self.check_opname(s, op_name, 5)
        context.traps[decimal.DivisionByZero] = divbyzerotrap
        context.traps[decimal.InvalidOperation] = invalidoptrap

    @pytest.mark.skip(reason="divmod not appropriate for decimal")
    def test_divmod(self, data):
        pass

    def test_error(self):
        pass


class TestComparisonOps(BaseDecimal, base.BaseComparisonOpsTests):

    def check_opname(self, s, op_name, other, exc=None):
        super(TestComparisonOps, self).check_opname(s, op_name,
                                                    other, exc=None)

    def _compare_other(self, s, data, op_name, other):
        self.check_opname(s, op_name, other)

    def test_compare_scalar(self, data, all_compare_operators):
        op_name = all_compare_operators
        s = pd.Series(data)
        self._compare_other(s, data, op_name, 0.5)

    def test_compare_array(self, data, all_compare_operators):
        op_name = all_compare_operators
        s = pd.Series(data)

        alter = np.random.choice([-1, 0, 1], len(data))
        # Randomly double, halve or keep same value
        other = pd.Series(data) * [decimal.Decimal(pow(2.0, i))
                                   for i in alter]
        self._compare_other(s, data, op_name, other)


class DecimalArrayWithoutFromSequence(DecimalArray):
    """Helper class for testing error handling in _from_sequence."""
    def _from_sequence(cls, scalars, dtype=None, copy=False):
        raise KeyError("For the test")


class DecimalArrayWithoutCoercion(DecimalArrayWithoutFromSequence):
    @classmethod
    def _create_arithmetic_method(cls, op):
        return cls._create_method(op, coerce_to_dtype=False)


DecimalArrayWithoutCoercion._add_arithmetic_ops()


def test_combine_from_sequence_raises():
    # https://github.com/pandas-dev/pandas/issues/22850
    ser = pd.Series(DecimalArrayWithoutFromSequence([
        decimal.Decimal("1.0"),
        decimal.Decimal("2.0")
    ]))
    result = ser.combine(ser, operator.add)

    # note: object dtype
    expected = pd.Series([decimal.Decimal("2.0"),
                          decimal.Decimal("4.0")], dtype="object")
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize("class_", [DecimalArrayWithoutFromSequence,
                                    DecimalArrayWithoutCoercion])
def test_scalar_ops_from_sequence_raises(class_):
    # op(EA, EA) should return an EA, or an ndarray if it's not possible
    # to return an EA with the return values.
    arr = class_([
        decimal.Decimal("1.0"),
        decimal.Decimal("2.0")
    ])
    result = arr + arr
    expected = np.array([decimal.Decimal("2.0"), decimal.Decimal("4.0")],
                        dtype="object")
    tm.assert_numpy_array_equal(result, expected)
