import calendar
import datetime
import decimal
import json
import locale
import math
import re
import time

import dateutil
import numpy as np
import pytest
import pytz

import pandas._libs.json as ujson
from pandas._libs.tslib import Timestamp
import pandas.compat as compat

from pandas import DataFrame, DatetimeIndex, Index, NaT, Series, date_range
import pandas._testing as tm


def _clean_dict(d):
    """
    Sanitize dictionary for JSON by converting all keys to strings.

    Parameters
    ----------
    d : dict
        The dictionary to convert.

    Returns
    -------
    cleaned_dict : dict
    """

    return {str(k): v for k, v in d.items()}


@pytest.fixture(
    params=[None, "split", "records", "values", "index"]  # Column indexed by default.
)
def orient(request):
    return request.param


@pytest.fixture(params=[None, True])
def numpy(request):
    return request.param


class TestUltraJSONTests:
    @pytest.mark.skipif(
        compat.is_platform_32bit(), reason="not compliant on 32-bit, xref #15865"
    )
    def test_encode_decimal(self):
        sut = decimal.Decimal("1337.1337")
        encoded = ujson.encode(sut, double_precision=15)
        decoded = ujson.decode(encoded)
        assert decoded == 1337.1337

        sut = decimal.Decimal("0.95")
        encoded = ujson.encode(sut, double_precision=1)
        assert encoded == "1.0"

        decoded = ujson.decode(encoded)
        assert decoded == 1.0

        sut = decimal.Decimal("0.94")
        encoded = ujson.encode(sut, double_precision=1)
        assert encoded == "0.9"

        decoded = ujson.decode(encoded)
        assert decoded == 0.9

        sut = decimal.Decimal("1.95")
        encoded = ujson.encode(sut, double_precision=1)
        assert encoded == "2.0"

        decoded = ujson.decode(encoded)
        assert decoded == 2.0

        sut = decimal.Decimal("-1.95")
        encoded = ujson.encode(sut, double_precision=1)
        assert encoded == "-2.0"

        decoded = ujson.decode(encoded)
        assert decoded == -2.0

        sut = decimal.Decimal("0.995")
        encoded = ujson.encode(sut, double_precision=2)
        assert encoded == "1.0"

        decoded = ujson.decode(encoded)
        assert decoded == 1.0

        sut = decimal.Decimal("0.9995")
        encoded = ujson.encode(sut, double_precision=3)
        assert encoded == "1.0"

        decoded = ujson.decode(encoded)
        assert decoded == 1.0

        sut = decimal.Decimal("0.99999999999999944")
        encoded = ujson.encode(sut, double_precision=15)
        assert encoded == "1.0"

        decoded = ujson.decode(encoded)
        assert decoded == 1.0

    @pytest.mark.parametrize("ensure_ascii", [True, False])
    def test_encode_string_conversion(self, ensure_ascii):
        string_input = "A string \\ / \b \f \n \r \t </script> &"
        not_html_encoded = '"A string \\\\ \\/ \\b \\f \\n \\r \\t <\\/script> &"'
        html_encoded = (
            '"A string \\\\ \\/ \\b \\f \\n \\r \\t \\u003c\\/script\\u003e \\u0026"'
        )

        def helper(expected_output, **encode_kwargs):
            output = ujson.encode(
                string_input, ensure_ascii=ensure_ascii, **encode_kwargs
            )

            assert output == expected_output
            assert string_input == json.loads(output)
            assert string_input == ujson.decode(output)

        # Default behavior assumes encode_html_chars=False.
        helper(not_html_encoded)

        # Make sure explicit encode_html_chars=False works.
        helper(not_html_encoded, encode_html_chars=False)

        # Make sure explicit encode_html_chars=True does the encoding.
        helper(html_encoded, encode_html_chars=True)

    @pytest.mark.parametrize(
        "long_number", [-4342969734183514, -12345678901234.56789012, -528656961.4399388]
    )
    def test_double_long_numbers(self, long_number):
        sut = {"a": long_number}
        encoded = ujson.encode(sut, double_precision=15)

        decoded = ujson.decode(encoded)
        assert sut == decoded

    def test_encode_non_c_locale(self):
        lc_category = locale.LC_NUMERIC

        # We just need one of these locales to work.
        for new_locale in ("it_IT.UTF-8", "Italian_Italy"):
            if tm.can_set_locale(new_locale, lc_category):
                with tm.set_locale(new_locale, lc_category):
                    assert ujson.loads(ujson.dumps(4.78e60)) == 4.78e60
                    assert ujson.loads("4.78", precise_float=True) == 4.78
                break

    def test_decimal_decode_test_precise(self):
        sut = {"a": 4.56}
        encoded = ujson.encode(sut)
        decoded = ujson.decode(encoded, precise_float=True)
        assert sut == decoded

    def test_encode_double_tiny_exponential(self):
        num = 1e-40
        assert num == ujson.decode(ujson.encode(num))
        num = 1e-100
        assert num == ujson.decode(ujson.encode(num))
        num = -1e-45
        assert num == ujson.decode(ujson.encode(num))
        num = -1e-145
        assert np.allclose(num, ujson.decode(ujson.encode(num)))

    @pytest.mark.parametrize("unicode_key", ["key1", "بن"])
    def test_encode_dict_with_unicode_keys(self, unicode_key):
        unicode_dict = {unicode_key: "value1"}
        assert unicode_dict == ujson.decode(ujson.encode(unicode_dict))

    @pytest.mark.parametrize(
        "double_input", [math.pi, -math.pi]  # Should work with negatives too.
    )
    def test_encode_double_conversion(self, double_input):
        output = ujson.encode(double_input)
        assert round(double_input, 5) == round(json.loads(output), 5)
        assert round(double_input, 5) == round(ujson.decode(output), 5)

    def test_encode_with_decimal(self):
        decimal_input = 1.0
        output = ujson.encode(decimal_input)

        assert output == "1.0"

    def test_encode_array_of_nested_arrays(self):
        nested_input = [[[[]]]] * 20
        output = ujson.encode(nested_input)

        assert nested_input == json.loads(output)
        assert nested_input == ujson.decode(output)

        nested_input = np.array(nested_input)
        tm.assert_numpy_array_equal(
            nested_input, ujson.decode(output, numpy=True, dtype=nested_input.dtype)
        )

    def test_encode_array_of_doubles(self):
        doubles_input = [31337.31337, 31337.31337, 31337.31337, 31337.31337] * 10
        output = ujson.encode(doubles_input)

        assert doubles_input == json.loads(output)
        assert doubles_input == ujson.decode(output)

        tm.assert_numpy_array_equal(
            np.array(doubles_input), ujson.decode(output, numpy=True)
        )

    def test_double_precision(self):
        double_input = 30.012345678901234
        output = ujson.encode(double_input, double_precision=15)

        assert double_input == json.loads(output)
        assert double_input == ujson.decode(output)

        for double_precision in (3, 9):
            output = ujson.encode(double_input, double_precision=double_precision)
            rounded_input = round(double_input, double_precision)

            assert rounded_input == json.loads(output)
            assert rounded_input == ujson.decode(output)

    @pytest.mark.parametrize("invalid_val", [20, -1, "9", None])
    def test_invalid_double_precision(self, invalid_val):
        double_input = 30.12345678901234567890
        expected_exception = ValueError if isinstance(invalid_val, int) else TypeError

        with pytest.raises(expected_exception):
            ujson.encode(double_input, double_precision=invalid_val)

    def test_encode_string_conversion2(self):
        string_input = "A string \\ / \b \f \n \r \t"
        output = ujson.encode(string_input)

        assert string_input == json.loads(output)
        assert string_input == ujson.decode(output)
        assert output == '"A string \\\\ \\/ \\b \\f \\n \\r \\t"'

    @pytest.mark.parametrize(
        "unicode_input",
        ["Räksmörgås اسامة بن محمد بن عوض بن لادن", "\xe6\x97\xa5\xd1\x88"],
    )
    def test_encode_unicode_conversion(self, unicode_input):
        enc = ujson.encode(unicode_input)
        dec = ujson.decode(enc)

        assert enc == json.dumps(unicode_input)
        assert dec == json.loads(enc)

    def test_encode_control_escaping(self):
        escaped_input = "\x19"
        enc = ujson.encode(escaped_input)
        dec = ujson.decode(enc)

        assert escaped_input == dec
        assert enc == json.dumps(escaped_input)

    def test_encode_unicode_surrogate_pair(self):
        surrogate_input = "\xf0\x90\x8d\x86"
        enc = ujson.encode(surrogate_input)
        dec = ujson.decode(enc)

        assert enc == json.dumps(surrogate_input)
        assert dec == json.loads(enc)

    def test_encode_unicode_4bytes_utf8(self):
        four_bytes_input = "\xf0\x91\x80\xb0TRAILINGNORMAL"
        enc = ujson.encode(four_bytes_input)
        dec = ujson.decode(enc)

        assert enc == json.dumps(four_bytes_input)
        assert dec == json.loads(enc)

    def test_encode_unicode_4bytes_utf8highest(self):
        four_bytes_input = "\xf3\xbf\xbf\xbfTRAILINGNORMAL"
        enc = ujson.encode(four_bytes_input)

        dec = ujson.decode(enc)

        assert enc == json.dumps(four_bytes_input)
        assert dec == json.loads(enc)

    def test_encode_array_in_array(self):
        arr_in_arr_input = [[[[]]]]
        output = ujson.encode(arr_in_arr_input)

        assert arr_in_arr_input == json.loads(output)
        assert output == json.dumps(arr_in_arr_input)
        assert arr_in_arr_input == ujson.decode(output)

        tm.assert_numpy_array_equal(
            np.array(arr_in_arr_input), ujson.decode(output, numpy=True)
        )

    @pytest.mark.parametrize(
        "num_input",
        [
            31337,
            -31337,  # Negative number.
            -9223372036854775808,  # Large negative number.
        ],
    )
    def test_encode_num_conversion(self, num_input):
        output = ujson.encode(num_input)
        assert num_input == json.loads(output)
        assert output == json.dumps(num_input)
        assert num_input == ujson.decode(output)

    def test_encode_list_conversion(self):
        list_input = [1, 2, 3, 4]
        output = ujson.encode(list_input)

        assert list_input == json.loads(output)
        assert list_input == ujson.decode(output)

        tm.assert_numpy_array_equal(
            np.array(list_input), ujson.decode(output, numpy=True)
        )

    def test_encode_dict_conversion(self):
        dict_input = {"k1": 1, "k2": 2, "k3": 3, "k4": 4}
        output = ujson.encode(dict_input)

        assert dict_input == json.loads(output)
        assert dict_input == ujson.decode(output)

    @pytest.mark.parametrize("builtin_value", [None, True, False])
    def test_encode_builtin_values_conversion(self, builtin_value):
        output = ujson.encode(builtin_value)
        assert builtin_value == json.loads(output)
        assert output == json.dumps(builtin_value)
        assert builtin_value == ujson.decode(output)

    def test_encode_datetime_conversion(self):
        datetime_input = datetime.datetime.fromtimestamp(time.time())
        output = ujson.encode(datetime_input, date_unit="s")
        expected = calendar.timegm(datetime_input.utctimetuple())

        assert int(expected) == json.loads(output)
        assert int(expected) == ujson.decode(output)

    def test_encode_date_conversion(self):
        date_input = datetime.date.fromtimestamp(time.time())
        output = ujson.encode(date_input, date_unit="s")

        tup = (date_input.year, date_input.month, date_input.day, 0, 0, 0)
        expected = calendar.timegm(tup)

        assert int(expected) == json.loads(output)
        assert int(expected) == ujson.decode(output)

    @pytest.mark.parametrize(
        "test",
        [datetime.time(), datetime.time(1, 2, 3), datetime.time(10, 12, 15, 343243)],
    )
    def test_encode_time_conversion_basic(self, test):
        output = ujson.encode(test)
        expected = f'"{test.isoformat()}"'
        assert expected == output

    def test_encode_time_conversion_pytz(self):
        # see gh-11473: to_json segfaults with timezone-aware datetimes
        test = datetime.time(10, 12, 15, 343243, pytz.utc)
        output = ujson.encode(test)
        expected = f'"{test.isoformat()}"'
        assert expected == output

    def test_encode_time_conversion_dateutil(self):
        # see gh-11473: to_json segfaults with timezone-aware datetimes
        test = datetime.time(10, 12, 15, 343243, dateutil.tz.tzutc())
        output = ujson.encode(test)
        expected = f'"{test.isoformat()}"'
        assert expected == output

    @pytest.mark.parametrize(
        "decoded_input", [NaT, np.datetime64("NaT"), np.nan, np.inf, -np.inf]
    )
    def test_encode_as_null(self, decoded_input):
        assert ujson.encode(decoded_input) == "null", "Expected null"

    def test_datetime_units(self):
        val = datetime.datetime(2013, 8, 17, 21, 17, 12, 215504)
        stamp = Timestamp(val)

        roundtrip = ujson.decode(ujson.encode(val, date_unit="s"))
        assert roundtrip == stamp.value // 10 ** 9

        roundtrip = ujson.decode(ujson.encode(val, date_unit="ms"))
        assert roundtrip == stamp.value // 10 ** 6

        roundtrip = ujson.decode(ujson.encode(val, date_unit="us"))
        assert roundtrip == stamp.value // 10 ** 3

        roundtrip = ujson.decode(ujson.encode(val, date_unit="ns"))
        assert roundtrip == stamp.value

        msg = "Invalid value 'foo' for option 'date_unit'"
        with pytest.raises(ValueError, match=msg):
            ujson.encode(val, date_unit="foo")

    def test_encode_to_utf8(self):
        unencoded = "\xe6\x97\xa5\xd1\x88"

        enc = ujson.encode(unencoded, ensure_ascii=False)
        dec = ujson.decode(enc)

        assert enc == json.dumps(unencoded, ensure_ascii=False)
        assert dec == json.loads(enc)

    def test_decode_from_unicode(self):
        unicode_input = '{"obj": 31337}'

        dec1 = ujson.decode(unicode_input)
        dec2 = ujson.decode(str(unicode_input))

        assert dec1 == dec2

    def test_encode_recursion_max(self):
        # 8 is the max recursion depth

        class O2:
            member = 0
            pass

        class O1:
            member = 0
            pass

        decoded_input = O1()
        decoded_input.member = O2()
        decoded_input.member.member = decoded_input

        with pytest.raises(OverflowError):
            ujson.encode(decoded_input)

    def test_decode_jibberish(self):
        jibberish = "fdsa sda v9sa fdsa"

        with pytest.raises(ValueError):
            ujson.decode(jibberish)

    @pytest.mark.parametrize(
        "broken_json",
        [
            "[",  # Broken array start.
            "{",  # Broken object start.
            "]",  # Broken array end.
            "}",  # Broken object end.
        ],
    )
    def test_decode_broken_json(self, broken_json):
        with pytest.raises(ValueError):
            ujson.decode(broken_json)

    @pytest.mark.parametrize("too_big_char", ["[", "{"])
    def test_decode_depth_too_big(self, too_big_char):
        with pytest.raises(ValueError):
            ujson.decode(too_big_char * (1024 * 1024))

    @pytest.mark.parametrize(
        "bad_string",
        [
            '"TESTING',  # Unterminated.
            '"TESTING\\"',  # Unterminated escape.
            "tru",  # Broken True.
            "fa",  # Broken False.
            "n",  # Broken None.
        ],
    )
    def test_decode_bad_string(self, bad_string):
        with pytest.raises(ValueError):
            ujson.decode(bad_string)

    @pytest.mark.parametrize("broken_json", ['{{1337:""}}', '{{"key":"}', "[[[true"])
    def test_decode_broken_json_leak(self, broken_json):
        for _ in range(1000):
            with pytest.raises(ValueError):
                ujson.decode(broken_json)

    @pytest.mark.parametrize(
        "invalid_dict",
        [
            "{{{{31337}}}}",  # No key.
            '{{{{"key":}}}}',  # No value.
            '{{{{"key"}}}}',  # No colon or value.
        ],
    )
    def test_decode_invalid_dict(self, invalid_dict):
        with pytest.raises(ValueError):
            ujson.decode(invalid_dict)

    @pytest.mark.parametrize(
        "numeric_int_as_str", ["31337", "-31337"]  # Should work with negatives.
    )
    def test_decode_numeric_int(self, numeric_int_as_str):
        assert int(numeric_int_as_str) == ujson.decode(numeric_int_as_str)

    def test_encode_null_character(self):
        wrapped_input = "31337 \x00 1337"
        output = ujson.encode(wrapped_input)

        assert wrapped_input == json.loads(output)
        assert output == json.dumps(wrapped_input)
        assert wrapped_input == ujson.decode(output)

        alone_input = "\x00"
        output = ujson.encode(alone_input)

        assert alone_input == json.loads(output)
        assert output == json.dumps(alone_input)
        assert alone_input == ujson.decode(output)
        assert '"  \\u0000\\r\\n "' == ujson.dumps("  \u0000\r\n ")

    def test_decode_null_character(self):
        wrapped_input = '"31337 \\u0000 31337"'
        assert ujson.decode(wrapped_input) == json.loads(wrapped_input)

    def test_encode_list_long_conversion(self):
        long_input = [
            9223372036854775807,
            9223372036854775807,
            9223372036854775807,
            9223372036854775807,
            9223372036854775807,
            9223372036854775807,
        ]
        output = ujson.encode(long_input)

        assert long_input == json.loads(output)
        assert long_input == ujson.decode(output)

        tm.assert_numpy_array_equal(
            np.array(long_input), ujson.decode(output, numpy=True, dtype=np.int64)
        )

    def test_encode_long_conversion(self):
        long_input = 9223372036854775807
        output = ujson.encode(long_input)

        assert long_input == json.loads(output)
        assert output == json.dumps(long_input)
        assert long_input == ujson.decode(output)

    @pytest.mark.parametrize(
        "int_exp", ["1337E40", "1.337E40", "1337E+9", "1.337e+40", "1.337E-4"]
    )
    def test_decode_numeric_int_exp(self, int_exp):
        assert ujson.decode(int_exp) == json.loads(int_exp)

    def test_loads_non_str_bytes_raises(self):
        msg = "Expected 'str' or 'bytes'"
        with pytest.raises(TypeError, match=msg):
            ujson.loads(None)

    def test_version(self):
        assert re.match(
            r"^\d+\.\d+(\.\d+)?$", ujson.__version__
        ), "ujson.__version__ must be a string like '1.4.0'"

    def test_encode_numeric_overflow(self):
        with pytest.raises(OverflowError):
            ujson.encode(12839128391289382193812939)

    def test_encode_numeric_overflow_nested(self):
        class Nested:
            x = 12839128391289382193812939

        for _ in range(0, 100):
            with pytest.raises(OverflowError):
                ujson.encode(Nested())

    @pytest.mark.parametrize("val", [3590016419, 2 ** 31, 2 ** 32, (2 ** 32) - 1])
    def test_decode_number_with_32bit_sign_bit(self, val):
        # Test that numbers that fit within 32 bits but would have the
        # sign bit set (2**31 <= x < 2**32) are decoded properly.
        doc = f'{{"id": {val}}}'
        assert ujson.decode(doc)["id"] == val

    def test_encode_big_escape(self):
        # Make sure no Exception is raised.
        for _ in range(10):
            base = "\u00e5".encode("utf-8")
            escape_input = base * 1024 * 1024 * 2
            ujson.encode(escape_input)

    def test_decode_big_escape(self):
        # Make sure no Exception is raised.
        for _ in range(10):
            base = "\u00e5".encode("utf-8")
            quote = b'"'

            escape_input = quote + (base * 1024 * 1024 * 2) + quote
            ujson.decode(escape_input)

    def test_to_dict(self):
        d = {"key": 31337}

        class DictTest:
            def toDict(self):
                return d

        o = DictTest()
        output = ujson.encode(o)

        dec = ujson.decode(output)
        assert dec == d

    def test_default_handler(self):
        class _TestObject:
            def __init__(self, val):
                self.val = val

            @property
            def recursive_attr(self):
                return _TestObject("recursive_attr")

            def __str__(self) -> str:
                return str(self.val)

        msg = "Maximum recursion level reached"
        with pytest.raises(OverflowError, match=msg):
            ujson.encode(_TestObject("foo"))
        assert '"foo"' == ujson.encode(_TestObject("foo"), default_handler=str)

        def my_handler(_):
            return "foobar"

        assert '"foobar"' == ujson.encode(
            _TestObject("foo"), default_handler=my_handler
        )

        def my_handler_raises(_):
            raise TypeError("I raise for anything")

        with pytest.raises(TypeError, match="I raise for anything"):
            ujson.encode(_TestObject("foo"), default_handler=my_handler_raises)

        def my_int_handler(_):
            return 42

        assert (
            ujson.decode(
                ujson.encode(_TestObject("foo"), default_handler=my_int_handler)
            )
            == 42
        )

        def my_obj_handler(_):
            return datetime.datetime(2013, 2, 3)

        assert ujson.decode(
            ujson.encode(datetime.datetime(2013, 2, 3))
        ) == ujson.decode(
            ujson.encode(_TestObject("foo"), default_handler=my_obj_handler)
        )

        obj_list = [_TestObject("foo"), _TestObject("bar")]
        assert json.loads(json.dumps(obj_list, default=str)) == ujson.decode(
            ujson.encode(obj_list, default_handler=str)
        )


class TestNumpyJSONTests:
    @pytest.mark.parametrize("bool_input", [True, False])
    def test_bool(self, bool_input):
        b = np.bool(bool_input)
        assert ujson.decode(ujson.encode(b)) == b

    def test_bool_array(self):
        bool_array = np.array(
            [True, False, True, True, False, True, False, False], dtype=np.bool
        )
        output = np.array(ujson.decode(ujson.encode(bool_array)), dtype=np.bool)
        tm.assert_numpy_array_equal(bool_array, output)

    def test_int(self, any_int_dtype):
        klass = np.dtype(any_int_dtype).type
        num = klass(1)

        assert klass(ujson.decode(ujson.encode(num))) == num

    def test_int_array(self, any_int_dtype):
        arr = np.arange(100, dtype=np.int)
        arr_input = arr.astype(any_int_dtype)

        arr_output = np.array(
            ujson.decode(ujson.encode(arr_input)), dtype=any_int_dtype
        )
        tm.assert_numpy_array_equal(arr_input, arr_output)

    def test_int_max(self, any_int_dtype):
        if any_int_dtype in ("int64", "uint64") and compat.is_platform_32bit():
            pytest.skip("Cannot test 64-bit integer on 32-bit platform")

        klass = np.dtype(any_int_dtype).type

        # uint64 max will always overflow,
        # as it's encoded to signed.
        if any_int_dtype == "uint64":
            num = np.iinfo("int64").max
        else:
            num = np.iinfo(any_int_dtype).max

        assert klass(ujson.decode(ujson.encode(num))) == num

    def test_float(self, float_dtype):
        klass = np.dtype(float_dtype).type
        num = klass(256.2013)

        assert klass(ujson.decode(ujson.encode(num))) == num

    def test_float_array(self, float_dtype):
        arr = np.arange(12.5, 185.72, 1.7322, dtype=np.float)
        float_input = arr.astype(float_dtype)

        float_output = np.array(
            ujson.decode(ujson.encode(float_input, double_precision=15)),
            dtype=float_dtype,
        )
        tm.assert_almost_equal(float_input, float_output)

    def test_float_max(self, float_dtype):
        klass = np.dtype(float_dtype).type
        num = klass(np.finfo(float_dtype).max / 10)

        tm.assert_almost_equal(
            klass(ujson.decode(ujson.encode(num, double_precision=15))), num
        )

    def test_array_basic(self):
        arr = np.arange(96)
        arr = arr.reshape((2, 2, 2, 2, 3, 2))

        tm.assert_numpy_array_equal(np.array(ujson.decode(ujson.encode(arr))), arr)
        tm.assert_numpy_array_equal(ujson.decode(ujson.encode(arr), numpy=True), arr)

    @pytest.mark.parametrize("shape", [(10, 10), (5, 5, 4), (100, 1)])
    def test_array_reshaped(self, shape):
        arr = np.arange(100)
        arr = arr.reshape(shape)

        tm.assert_numpy_array_equal(np.array(ujson.decode(ujson.encode(arr))), arr)
        tm.assert_numpy_array_equal(ujson.decode(ujson.encode(arr), numpy=True), arr)

    def test_array_list(self):
        arr_list = [
            "a",
            list(),
            dict(),
            dict(),
            list(),
            42,
            97.8,
            ["a", "b"],
            {"key": "val"},
        ]
        arr = np.array(arr_list, dtype=object)
        result = np.array(ujson.decode(ujson.encode(arr)), dtype=object)
        tm.assert_numpy_array_equal(result, arr)

    def test_array_float(self):
        dtype = np.float32

        arr = np.arange(100.202, 200.202, 1, dtype=dtype)
        arr = arr.reshape((5, 5, 4))

        arr_out = np.array(ujson.decode(ujson.encode(arr)), dtype=dtype)
        tm.assert_almost_equal(arr, arr_out)

        arr_out = ujson.decode(ujson.encode(arr), numpy=True, dtype=dtype)
        tm.assert_almost_equal(arr, arr_out)

    def test_0d_array(self):
        # gh-18878
        msg = re.escape("array(1) (0d array) is not JSON serializable at the moment")
        with pytest.raises(TypeError, match=msg):
            ujson.encode(np.array(1))

    @pytest.mark.parametrize(
        "bad_input,exc_type,kwargs",
        [
            ([{}, []], ValueError, {}),
            ([42, None], TypeError, {}),
            ([["a"], 42], ValueError, {}),
            ([42, {}, "a"], TypeError, {}),
            ([42, ["a"], 42], ValueError, {}),
            (["a", "b", [], "c"], ValueError, {}),
            ([{"a": "b"}], ValueError, dict(labelled=True)),
            ({"a": {"b": {"c": 42}}}, ValueError, dict(labelled=True)),
            ([{"a": 42, "b": 23}, {"c": 17}], ValueError, dict(labelled=True)),
        ],
    )
    def test_array_numpy_except(self, bad_input, exc_type, kwargs):
        with pytest.raises(exc_type):
            ujson.decode(ujson.dumps(bad_input), numpy=True, **kwargs)

    def test_array_numpy_labelled(self):
        labelled_input = {"a": []}
        output = ujson.loads(ujson.dumps(labelled_input), numpy=True, labelled=True)
        assert (np.empty((1, 0)) == output[0]).all()
        assert (np.array(["a"]) == output[1]).all()
        assert output[2] is None

        labelled_input = [{"a": 42}]
        output = ujson.loads(ujson.dumps(labelled_input), numpy=True, labelled=True)
        assert (np.array(["a"]) == output[2]).all()
        assert (np.array([42]) == output[0]).all()
        assert output[1] is None

        # see gh-10837: write out the dump explicitly
        # so there is no dependency on iteration order
        input_dumps = '[{"a": 42, "b":31}, {"a": 24, "c": 99}, {"a": 2.4, "b": 78}]'
        output = ujson.loads(input_dumps, numpy=True, labelled=True)
        expected_vals = np.array([42, 31, 24, 99, 2.4, 78], dtype=int).reshape((3, 2))
        assert (expected_vals == output[0]).all()
        assert output[1] is None
        assert (np.array(["a", "b"]) == output[2]).all()

        input_dumps = (
            '{"1": {"a": 42, "b":31}, "2": {"a": 24, "c": 99}, '
            '"3": {"a": 2.4, "b": 78}}'
        )
        output = ujson.loads(input_dumps, numpy=True, labelled=True)
        expected_vals = np.array([42, 31, 24, 99, 2.4, 78], dtype=int).reshape((3, 2))
        assert (expected_vals == output[0]).all()
        assert (np.array(["1", "2", "3"]) == output[1]).all()
        assert (np.array(["a", "b"]) == output[2]).all()


class TestPandasJSONTests:
    def test_dataframe(self, orient, numpy):
        if orient == "records" and numpy:
            pytest.skip("Not idiomatic pandas")

        df = DataFrame(
            [[1, 2, 3], [4, 5, 6]], index=["a", "b"], columns=["x", "y", "z"]
        )
        encode_kwargs = {} if orient is None else dict(orient=orient)
        decode_kwargs = {} if numpy is None else dict(numpy=numpy)

        output = ujson.decode(ujson.encode(df, **encode_kwargs), **decode_kwargs)

        # Ensure proper DataFrame initialization.
        if orient == "split":
            dec = _clean_dict(output)
            output = DataFrame(**dec)
        else:
            output = DataFrame(output)

        # Corrections to enable DataFrame comparison.
        if orient == "values":
            df.columns = [0, 1, 2]
            df.index = [0, 1]
        elif orient == "records":
            df.index = [0, 1]
        elif orient == "index":
            df = df.transpose()

        tm.assert_frame_equal(output, df, check_dtype=False)

    def test_dataframe_nested(self, orient):
        df = DataFrame(
            [[1, 2, 3], [4, 5, 6]], index=["a", "b"], columns=["x", "y", "z"]
        )

        nested = {"df1": df, "df2": df.copy()}
        kwargs = {} if orient is None else dict(orient=orient)

        exp = {
            "df1": ujson.decode(ujson.encode(df, **kwargs)),
            "df2": ujson.decode(ujson.encode(df, **kwargs)),
        }
        assert ujson.decode(ujson.encode(nested, **kwargs)) == exp

    def test_dataframe_numpy_labelled(self, orient):
        if orient in ("split", "values"):
            pytest.skip("Incompatible with labelled=True")

        df = DataFrame(
            [[1, 2, 3], [4, 5, 6]],
            index=["a", "b"],
            columns=["x", "y", "z"],
            dtype=np.int,
        )
        kwargs = {} if orient is None else dict(orient=orient)

        output = DataFrame(
            *ujson.decode(ujson.encode(df, **kwargs), numpy=True, labelled=True)
        )

        if orient is None:
            df = df.T
        elif orient == "records":
            df.index = [0, 1]

        tm.assert_frame_equal(output, df)

    def test_series(self, orient, numpy):
        s = Series(
            [10, 20, 30, 40, 50, 60], name="series", index=[6, 7, 8, 9, 10, 15]
        ).sort_values()

        encode_kwargs = {} if orient is None else dict(orient=orient)
        decode_kwargs = {} if numpy is None else dict(numpy=numpy)

        output = ujson.decode(ujson.encode(s, **encode_kwargs), **decode_kwargs)

        if orient == "split":
            dec = _clean_dict(output)
            output = Series(**dec)
        else:
            output = Series(output)

        if orient in (None, "index"):
            s.name = None
            output = output.sort_values()
            s.index = ["6", "7", "8", "9", "10", "15"]
        elif orient in ("records", "values"):
            s.name = None
            s.index = [0, 1, 2, 3, 4, 5]

        tm.assert_series_equal(output, s, check_dtype=False)

    def test_series_nested(self, orient):
        s = Series(
            [10, 20, 30, 40, 50, 60], name="series", index=[6, 7, 8, 9, 10, 15]
        ).sort_values()
        nested = {"s1": s, "s2": s.copy()}
        kwargs = {} if orient is None else dict(orient=orient)

        exp = {
            "s1": ujson.decode(ujson.encode(s, **kwargs)),
            "s2": ujson.decode(ujson.encode(s, **kwargs)),
        }
        assert ujson.decode(ujson.encode(nested, **kwargs)) == exp

    def test_index(self):
        i = Index([23, 45, 18, 98, 43, 11], name="index")

        # Column indexed.
        output = Index(ujson.decode(ujson.encode(i)), name="index")
        tm.assert_index_equal(i, output)

        output = Index(ujson.decode(ujson.encode(i), numpy=True), name="index")
        tm.assert_index_equal(i, output)

        dec = _clean_dict(ujson.decode(ujson.encode(i, orient="split")))
        output = Index(**dec)

        tm.assert_index_equal(i, output)
        assert i.name == output.name

        dec = _clean_dict(ujson.decode(ujson.encode(i, orient="split"), numpy=True))
        output = Index(**dec)

        tm.assert_index_equal(i, output)
        assert i.name == output.name

        output = Index(ujson.decode(ujson.encode(i, orient="values")), name="index")
        tm.assert_index_equal(i, output)

        output = Index(
            ujson.decode(ujson.encode(i, orient="values"), numpy=True), name="index"
        )
        tm.assert_index_equal(i, output)

        output = Index(ujson.decode(ujson.encode(i, orient="records")), name="index")
        tm.assert_index_equal(i, output)

        output = Index(
            ujson.decode(ujson.encode(i, orient="records"), numpy=True), name="index"
        )
        tm.assert_index_equal(i, output)

        output = Index(ujson.decode(ujson.encode(i, orient="index")), name="index")
        tm.assert_index_equal(i, output)

        output = Index(
            ujson.decode(ujson.encode(i, orient="index"), numpy=True), name="index"
        )
        tm.assert_index_equal(i, output)

    def test_datetime_index(self):
        date_unit = "ns"

        rng = date_range("1/1/2000", periods=20)
        encoded = ujson.encode(rng, date_unit=date_unit)

        decoded = DatetimeIndex(np.array(ujson.decode(encoded)))
        tm.assert_index_equal(rng, decoded)

        ts = Series(np.random.randn(len(rng)), index=rng)
        decoded = Series(ujson.decode(ujson.encode(ts, date_unit=date_unit)))

        idx_values = decoded.index.values.astype(np.int64)
        decoded.index = DatetimeIndex(idx_values)
        tm.assert_series_equal(ts, decoded)

    @pytest.mark.parametrize(
        "invalid_arr",
        [
            "[31337,]",  # Trailing comma.
            "[,31337]",  # Leading comma.
            "[]]",  # Unmatched bracket.
            "[,]",  # Only comma.
        ],
    )
    def test_decode_invalid_array(self, invalid_arr):
        with pytest.raises(ValueError):
            ujson.decode(invalid_arr)

    @pytest.mark.parametrize("arr", [[], [31337]])
    def test_decode_array(self, arr):
        assert arr == ujson.decode(str(arr))

    @pytest.mark.parametrize("extreme_num", [9223372036854775807, -9223372036854775808])
    def test_decode_extreme_numbers(self, extreme_num):
        assert extreme_num == ujson.decode(str(extreme_num))

    @pytest.mark.parametrize(
        "too_extreme_num", ["9223372036854775808", "-90223372036854775809"]
    )
    def test_decode_too_extreme_numbers(self, too_extreme_num):
        with pytest.raises(ValueError):
            ujson.decode(too_extreme_num)

    def test_decode_with_trailing_whitespaces(self):
        assert {} == ujson.decode("{}\n\t ")

    def test_decode_with_trailing_non_whitespaces(self):
        with pytest.raises(ValueError):
            ujson.decode("{}\n\t a")

    def test_decode_array_with_big_int(self):
        with pytest.raises(ValueError):
            ujson.loads("[18446098363113800555]")

    @pytest.mark.parametrize(
        "float_number",
        [
            1.1234567893,
            1.234567893,
            1.34567893,
            1.4567893,
            1.567893,
            1.67893,
            1.7893,
            1.893,
            1.3,
        ],
    )
    @pytest.mark.parametrize("sign", [-1, 1])
    def test_decode_floating_point(self, sign, float_number):
        float_number *= sign
        tm.assert_almost_equal(
            float_number, ujson.loads(str(float_number)), check_less_precise=15
        )

    def test_encode_big_set(self):
        s = set()

        for x in range(0, 100000):
            s.add(x)

        # Make sure no Exception is raised.
        ujson.encode(s)

    def test_encode_empty_set(self):
        assert "[]" == ujson.encode(set())

    def test_encode_set(self):
        s = {1, 2, 3, 4, 5, 6, 7, 8, 9}
        enc = ujson.encode(s)
        dec = ujson.decode(enc)

        for v in dec:
            assert v in s
