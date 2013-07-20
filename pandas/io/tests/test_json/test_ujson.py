import unittest
from unittest import TestCase

import pandas.json as ujson
try:
    import json
except ImportError:
    import simplejson as json
import math
import nose
import platform
import sys
import time
import datetime
import calendar
import StringIO
import re
import random
import decimal
from functools import partial
import pandas.util.py3compat as py3compat

import numpy as np
from pandas.util.testing import assert_almost_equal
from numpy.testing import (assert_array_equal,
                           assert_array_almost_equal_nulp,
                           assert_approx_equal)
from pandas import DataFrame, Series, Index
import pandas.util.testing as tm


def _skip_if_python_ver(skip_major, skip_minor=None):
    major, minor = sys.version_info[:2]
    if major == skip_major and (skip_minor is None or minor == skip_minor):
        raise nose.SkipTest

json_unicode = (json.dumps if sys.version_info[0] >= 3
                else partial(json.dumps, encoding="utf-8"))

class UltraJSONTests(TestCase):

    def test_encodeDecimal(self):
        sut = decimal.Decimal("1337.1337")
        encoded = ujson.encode(sut, double_precision=15)
        decoded = ujson.decode(encoded)
        self.assertEquals(decoded, 1337.1337)

    def test_encodeStringConversion(self):
        input = "A string \\ / \b \f \n \r \t </script> &"
        not_html_encoded = '"A string \\\\ \\/ \\b \\f \\n \\r \\t <\\/script> &"'
        html_encoded = '"A string \\\\ \\/ \\b \\f \\n \\r \\t \\u003c\\/script\\u003e \\u0026"'

        def helper(expected_output, **encode_kwargs):
            output = ujson.encode(input, **encode_kwargs)
            self.assertEquals(input, json.loads(output))
            self.assertEquals(output, expected_output)
            self.assertEquals(input, ujson.decode(output))

        # Default behavior assumes encode_html_chars=False.
        helper(not_html_encoded, ensure_ascii=True)
        helper(not_html_encoded, ensure_ascii=False)

        # Make sure explicit encode_html_chars=False works.
        helper(not_html_encoded, ensure_ascii=True, encode_html_chars=False)
        helper(not_html_encoded, ensure_ascii=False, encode_html_chars=False)

        # Make sure explicit encode_html_chars=True does the encoding.
        helper(html_encoded, ensure_ascii=True, encode_html_chars=True)
        helper(html_encoded, ensure_ascii=False, encode_html_chars=True)

    def test_doubleLongIssue(self):
        sut = {u'a': -4342969734183514}
        encoded = json.dumps(sut)
        decoded = json.loads(encoded)
        self.assertEqual(sut, decoded)
        encoded = ujson.encode(sut, double_precision=15)
        decoded = ujson.decode(encoded)
        self.assertEqual(sut, decoded)

    def test_doubleLongDecimalIssue(self):
        sut = {u'a': -12345678901234.56789012}
        encoded = json.dumps(sut)
        decoded = json.loads(encoded)
        self.assertEqual(sut, decoded)
        encoded = ujson.encode(sut, double_precision=15)
        decoded = ujson.decode(encoded)
        self.assertEqual(sut, decoded)


    def test_encodeDecodeLongDecimal(self):
        sut = {u'a': -528656961.4399388}
        encoded = ujson.dumps(sut, double_precision=15)
        ujson.decode(encoded)

    def test_decimalDecodeTestPrecise(self):
        sut = {u'a': 4.56}
        encoded = ujson.encode(sut)
        decoded = ujson.decode(encoded, precise_float=True)
        self.assertEqual(sut, decoded)

    def test_encodeDoubleTinyExponential(self):
        num = 1e-40
        self.assertEqual(num, ujson.decode(ujson.encode(num)))
        num = 1e-100
        self.assertEqual(num, ujson.decode(ujson.encode(num)))
        num = -1e-45
        self.assertEqual(num, ujson.decode(ujson.encode(num)))
        num = -1e-145
        self.assert_(np.allclose(num, ujson.decode(ujson.encode(num))))

    def test_encodeDictWithUnicodeKeys(self):
        input = { u"key1": u"value1", u"key1": u"value1", u"key1": u"value1", u"key1": u"value1", u"key1": u"value1", u"key1": u"value1" }
        output = ujson.encode(input)

        input = { u"بن": u"value1", u"بن": u"value1", u"بن": u"value1", u"بن": u"value1", u"بن": u"value1", u"بن": u"value1", u"بن": u"value1" }
        output = ujson.encode(input)

        pass

    def test_encodeDoubleConversion(self):
        input = math.pi
        output = ujson.encode(input)
        self.assertEquals(round(input, 5), round(json.loads(output), 5))
        self.assertEquals(round(input, 5), round(ujson.decode(output), 5))

    def test_encodeWithDecimal(self):
        input = 1.0
        output = ujson.encode(input)
        self.assertEquals(output, "1.0")

    def test_encodeDoubleNegConversion(self):
        input = -math.pi
        output = ujson.encode(input)

        self.assertEquals(round(input, 5), round(json.loads(output), 5))
        self.assertEquals(round(input, 5), round(ujson.decode(output), 5))

    def test_encodeArrayOfNestedArrays(self):
        input = [[[[]]]] * 20
        output = ujson.encode(input)
        self.assertEquals(input, json.loads(output))
        #self.assertEquals(output, json.dumps(input))
        self.assertEquals(input, ujson.decode(output))
        input = np.array(input)
        assert_array_equal(input, ujson.decode(output, numpy=True, dtype=input.dtype))

    def test_encodeArrayOfDoubles(self):
        input = [ 31337.31337, 31337.31337, 31337.31337, 31337.31337] * 10
        output = ujson.encode(input)
        self.assertEquals(input, json.loads(output))
        #self.assertEquals(output, json.dumps(input))
        self.assertEquals(input, ujson.decode(output))
        assert_array_equal(np.array(input), ujson.decode(output, numpy=True))

    def test_doublePrecisionTest(self):
        input = 30.012345678901234
        output = ujson.encode(input, double_precision = 15)
        self.assertEquals(input, json.loads(output))
        self.assertEquals(input, ujson.decode(output))

        output = ujson.encode(input, double_precision = 9)
        self.assertEquals(round(input, 9), json.loads(output))
        self.assertEquals(round(input, 9), ujson.decode(output))

        output = ujson.encode(input, double_precision = 3)
        self.assertEquals(round(input, 3), json.loads(output))
        self.assertEquals(round(input, 3), ujson.decode(output))

    def test_invalidDoublePrecision(self):
        input = 30.12345678901234567890

        self.assertRaises(ValueError, ujson.encode, input, double_precision = 20)
        self.assertRaises(ValueError, ujson.encode, input, double_precision = -1)

        # will throw typeError
        self.assertRaises(TypeError, ujson.encode, input, double_precision = '9')
        # will throw typeError
        self.assertRaises(TypeError, ujson.encode, input, double_precision = None)


    def test_encodeStringConversion(self):
        input = "A string \\ / \b \f \n \r \t"
        output = ujson.encode(input)
        self.assertEquals(input, json.loads(output))
        self.assertEquals(output, '"A string \\\\ \\/ \\b \\f \\n \\r \\t"')
        self.assertEquals(input, ujson.decode(output))
        pass

    def test_decodeUnicodeConversion(self):
        pass

    def test_encodeUnicodeConversion1(self):
        input = "Räksmörgås اسامة بن محمد بن عوض بن لادن"
        enc = ujson.encode(input)
        dec = ujson.decode(enc)
        self.assertEquals(enc, json_unicode(input))
        self.assertEquals(dec, json.loads(enc))

    def test_encodeControlEscaping(self):
        input = "\x19"
        enc = ujson.encode(input)
        dec = ujson.decode(enc)
        self.assertEquals(input, dec)
        self.assertEquals(enc, json_unicode(input))


    def test_encodeUnicodeConversion2(self):
        input = "\xe6\x97\xa5\xd1\x88"
        enc = ujson.encode(input)
        dec = ujson.decode(enc)
        self.assertEquals(enc, json_unicode(input))
        self.assertEquals(dec, json.loads(enc))

    def test_encodeUnicodeSurrogatePair(self):
        _skip_if_python_ver(2, 5)
        _skip_if_python_ver(2, 6)
        input = "\xf0\x90\x8d\x86"
        enc = ujson.encode(input)
        dec = ujson.decode(enc)

        self.assertEquals(enc, json_unicode(input))
        self.assertEquals(dec, json.loads(enc))

    def test_encodeUnicode4BytesUTF8(self):
        _skip_if_python_ver(2, 5)
        _skip_if_python_ver(2, 6)
        input = "\xf0\x91\x80\xb0TRAILINGNORMAL"
        enc = ujson.encode(input)
        dec = ujson.decode(enc)

        self.assertEquals(enc, json_unicode(input))
        self.assertEquals(dec, json.loads(enc))

    def test_encodeUnicode4BytesUTF8Highest(self):
        _skip_if_python_ver(2, 5)
        _skip_if_python_ver(2, 6)
        input = "\xf3\xbf\xbf\xbfTRAILINGNORMAL"
        enc = ujson.encode(input)

        dec = ujson.decode(enc)

        self.assertEquals(enc, json_unicode(input))
        self.assertEquals(dec, json.loads(enc))


    def test_encodeArrayInArray(self):
        input = [[[[]]]]
        output = ujson.encode(input)

        self.assertEquals(input, json.loads(output))
        self.assertEquals(output, json.dumps(input))
        self.assertEquals(input, ujson.decode(output))
        assert_array_equal(np.array(input), ujson.decode(output, numpy=True))
        pass

    def test_encodeIntConversion(self):
        input = 31337
        output = ujson.encode(input)
        self.assertEquals(input, json.loads(output))
        self.assertEquals(output, json.dumps(input))
        self.assertEquals(input, ujson.decode(output))
        pass

    def test_encodeIntNegConversion(self):
        input = -31337
        output = ujson.encode(input)
        self.assertEquals(input, json.loads(output))
        self.assertEquals(output, json.dumps(input))
        self.assertEquals(input, ujson.decode(output))
        pass


    def test_encodeLongNegConversion(self):
        input = -9223372036854775808
        output = ujson.encode(input)

        outputjson = json.loads(output)
        outputujson = ujson.decode(output)

        self.assertEquals(input, json.loads(output))
        self.assertEquals(output, json.dumps(input))
        self.assertEquals(input, ujson.decode(output))
        pass

    def test_encodeListConversion(self):
        input = [ 1, 2, 3, 4 ]
        output = ujson.encode(input)
        self.assertEquals(input, json.loads(output))
        self.assertEquals(input, ujson.decode(output))
        assert_array_equal(np.array(input), ujson.decode(output, numpy=True))
        pass

    def test_encodeDictConversion(self):
        input = { "k1": 1, "k2":  2, "k3": 3, "k4": 4 }
        output = ujson.encode(input)
        self.assertEquals(input, json.loads(output))
        self.assertEquals(input, ujson.decode(output))
        self.assertEquals(input, ujson.decode(output))
        pass

    def test_encodeNoneConversion(self):
        input = None
        output = ujson.encode(input)
        self.assertEquals(input, json.loads(output))
        self.assertEquals(output, json.dumps(input))
        self.assertEquals(input, ujson.decode(output))
        pass

    def test_encodeTrueConversion(self):
        input = True
        output = ujson.encode(input)
        self.assertEquals(input, json.loads(output))
        self.assertEquals(output, json.dumps(input))
        self.assertEquals(input, ujson.decode(output))
        pass

    def test_encodeFalseConversion(self):
        input = False
        output = ujson.encode(input)
        self.assertEquals(input, json.loads(output))
        self.assertEquals(output, json.dumps(input))
        self.assertEquals(input, ujson.decode(output))
        pass

    # def test_encodeDatetimeConversion(self):
    #     ts = time.time()
    #     input = datetime.datetime.fromtimestamp(ts)
    #     output = ujson.encode(input)
    #     expected = calendar.timegm(input.utctimetuple())
    #     self.assertEquals(int(expected), json.loads(output))
    #     self.assertEquals(int(expected), ujson.decode(output))
    #     pass

    # def test_encodeDateConversion(self):
    #     ts = time.time()
    #     input = datetime.date.fromtimestamp(ts)

    #     output = ujson.encode(input)
    #     tup = ( input.year, input.month, input.day, 0, 0, 0 )

    #     expected = calendar.timegm(tup)
    #     self.assertEquals(int(expected), json.loads(output))
    #     self.assertEquals(int(expected), ujson.decode(output))

    def test_datetime_nanosecond_unit(self):
        from datetime import datetime
        from pandas.lib import Timestamp

        val = datetime.now()
        stamp = Timestamp(val)

        roundtrip = ujson.decode(ujson.encode(val))
        self.assert_(roundtrip == stamp.value)

    def test_encodeToUTF8(self):
        _skip_if_python_ver(2, 5)
        input = "\xe6\x97\xa5\xd1\x88"
        enc = ujson.encode(input, ensure_ascii=False)
        dec = ujson.decode(enc)
        self.assertEquals(enc, json_unicode(input, ensure_ascii=False))
        self.assertEquals(dec, json.loads(enc))

    def test_decodeFromUnicode(self):
        input = u"{\"obj\": 31337}"
        dec1 = ujson.decode(input)
        dec2 = ujson.decode(str(input))
        self.assertEquals(dec1, dec2)

    def test_encodeRecursionMax(self):
        # 8 is the max recursion depth

        class O2:
            member = 0
            pass

        class O1:
            member = 0
            pass

        input = O1()
        input.member = O2()
        input.member.member = input

        try:
            output = ujson.encode(input)
            assert False, "Expected overflow exception"
        except(OverflowError):
            pass

    def test_encodeDoubleNan(self):
        input = np.nan
        assert ujson.encode(input) == 'null', "Expected null"

    def test_encodeDoubleInf(self):
        input = np.inf
        assert ujson.encode(input) == 'null', "Expected null"

    def test_encodeDoubleNegInf(self):
        input = -np.inf
        assert ujson.encode(input) == 'null', "Expected null"


    def test_decodeJibberish(self):
        input = "fdsa sda v9sa fdsa"
        try:
            ujson.decode(input)
            assert False, "Expected exception!"
        except(ValueError):
            return
        assert False, "Wrong exception"

    def test_decodeBrokenArrayStart(self):
        input = "["
        try:
            ujson.decode(input)
            assert False, "Expected exception!"
        except(ValueError):
            return
        assert False, "Wrong exception"

    def test_decodeBrokenObjectStart(self):
        input = "{"
        try:
            ujson.decode(input)
            assert False, "Expected exception!"
        except(ValueError):
            return
        assert False, "Wrong exception"

    def test_decodeBrokenArrayEnd(self):
        input = "]"
        try:
            ujson.decode(input)
            assert False, "Expected exception!"
        except(ValueError):
            return
        assert False, "Wrong exception"

    def test_decodeArrayDepthTooBig(self):
        input = '[' * (1024 * 1024)
        try:
            ujson.decode(input)
            assert False, "Expected exception!"
        except(ValueError):
            return
        assert False, "Wrong exception"

    def test_decodeBrokenObjectEnd(self):
        input = "}"
        try:
            ujson.decode(input)
            assert False, "Expected exception!"
        except(ValueError):
            return
        assert False, "Wrong exception"

    def test_decodeObjectDepthTooBig(self):
        input = '{' * (1024 * 1024)
        try:
            ujson.decode(input)
            assert False, "Expected exception!"
        except(ValueError):
            return
        assert False, "Wrong exception"

    def test_decodeStringUnterminated(self):
        input = "\"TESTING"
        try:
            ujson.decode(input)
            assert False, "Expected exception!"
        except(ValueError):
            return
        assert False, "Wrong exception"

    def test_decodeStringUntermEscapeSequence(self):
        input = "\"TESTING\\\""
        try:
            ujson.decode(input)
            assert False, "Expected exception!"
        except(ValueError):
            return
        assert False, "Wrong exception"

    def test_decodeStringBadEscape(self):
        input = "\"TESTING\\\""
        try:
            ujson.decode(input)
            assert False, "Expected exception!"
        except(ValueError):
            return
        assert False, "Wrong exception"

    def test_decodeTrueBroken(self):
        input = "tru"
        try:
            ujson.decode(input)
            assert False, "Expected exception!"
        except(ValueError):
            return
        assert False, "Wrong exception"

    def test_decodeFalseBroken(self):
        input = "fa"
        try:
            ujson.decode(input)
            assert False, "Expected exception!"
        except(ValueError):
            return
        assert False, "Wrong exception"

    def test_decodeNullBroken(self):
        input = "n"
        try:
            ujson.decode(input)
            assert False, "Expected exception!"
        except(ValueError):
            return
        assert False, "Wrong exception"


    def test_decodeBrokenDictKeyTypeLeakTest(self):
        input = '{{1337:""}}'
        for x in xrange(1000):
            try:
                ujson.decode(input)
                assert False, "Expected exception!"
            except(ValueError),e:
                continue

            assert False, "Wrong exception"

    def test_decodeBrokenDictLeakTest(self):
        input = '{{"key":"}'
        for x in xrange(1000):
            try:
                ujson.decode(input)
                assert False, "Expected exception!"
            except(ValueError):
                continue

            assert False, "Wrong exception"

    def test_decodeBrokenListLeakTest(self):
        input = '[[[true'
        for x in xrange(1000):
            try:
                ujson.decode(input)
                assert False, "Expected exception!"
            except(ValueError):
                continue

            assert False, "Wrong exception"

    def test_decodeDictWithNoKey(self):
        input = "{{{{31337}}}}"
        try:
            ujson.decode(input)
            assert False, "Expected exception!"
        except(ValueError):
            return

        assert False, "Wrong exception"

    def test_decodeDictWithNoColonOrValue(self):
        input = "{{{{\"key\"}}}}"
        try:
            ujson.decode(input)
            assert False, "Expected exception!"
        except(ValueError):
            return

        assert False, "Wrong exception"

    def test_decodeDictWithNoValue(self):
        input = "{{{{\"key\":}}}}"
        try:
            ujson.decode(input)
            assert False, "Expected exception!"
        except(ValueError):
            return

        assert False, "Wrong exception"

    def test_decodeNumericIntPos(self):
        input = "31337"
        self.assertEquals (31337, ujson.decode(input))

    def test_decodeNumericIntNeg(self):
        input = "-31337"
        self.assertEquals (-31337, ujson.decode(input))

    def test_encodeUnicode4BytesUTF8Fail(self):
        _skip_if_python_ver(3)
        input = "\xfd\xbf\xbf\xbf\xbf\xbf"
        try:
            enc = ujson.encode(input)
            assert False, "Expected exception"
        except OverflowError:
            pass

    def test_encodeNullCharacter(self):
        input = "31337 \x00 1337"
        output = ujson.encode(input)
        self.assertEquals(input, json.loads(output))
        self.assertEquals(output, json.dumps(input))
        self.assertEquals(input, ujson.decode(output))

        input = "\x00"
        output = ujson.encode(input)
        self.assertEquals(input, json.loads(output))
        self.assertEquals(output, json.dumps(input))
        self.assertEquals(input, ujson.decode(output))

        self.assertEquals('"  \\u0000\\r\\n "', ujson.dumps(u"  \u0000\r\n "))
        pass

    def test_decodeNullCharacter(self):
        input = "\"31337 \\u0000 31337\""
        self.assertEquals(ujson.decode(input), json.loads(input))


    def test_encodeListLongConversion(self):
        input = [9223372036854775807, 9223372036854775807, 9223372036854775807,
                 9223372036854775807, 9223372036854775807, 9223372036854775807 ]
        output = ujson.encode(input)
        self.assertEquals(input, json.loads(output))
        self.assertEquals(input, ujson.decode(output))
        assert_array_equal(np.array(input), ujson.decode(output, numpy=True,
                                                         dtype=np.int64))
        pass

    def test_encodeLongConversion(self):
        input = 9223372036854775807
        output = ujson.encode(input)
        self.assertEquals(input, json.loads(output))
        self.assertEquals(output, json.dumps(input))
        self.assertEquals(input, ujson.decode(output))
        pass

    def test_numericIntExp(self):
        input = "1337E40"
        output = ujson.decode(input)
        self.assertEquals(output, json.loads(input))

    def test_numericIntFrcExp(self):
        input = "1.337E40"
        output = ujson.decode(input)
        self.assertAlmostEqual(output, json.loads(input))

    def test_decodeNumericIntExpEPLUS(self):
        input = "1337E+9"
        output = ujson.decode(input)
        self.assertAlmostEqual(output, json.loads(input))

    def test_decodeNumericIntExpePLUS(self):
        input = "1.337e+40"
        output = ujson.decode(input)
        self.assertAlmostEqual(output, json.loads(input))

    def test_decodeNumericIntExpE(self):
        input = "1337E40"
        output = ujson.decode(input)
        self.assertAlmostEqual(output, json.loads(input))

    def test_decodeNumericIntExpe(self):
        input = "1337e40"
        output = ujson.decode(input)
        self.assertAlmostEqual(output, json.loads(input))

    def test_decodeNumericIntExpEMinus(self):
        input = "1.337E-4"
        output = ujson.decode(input)
        self.assertAlmostEqual(output, json.loads(input))

    def test_decodeNumericIntExpeMinus(self):
        input = "1.337e-4"
        output = ujson.decode(input)
        self.assertAlmostEqual(output, json.loads(input))

    def test_dumpToFile(self):
        f = StringIO.StringIO()
        ujson.dump([1, 2, 3], f)
        self.assertEquals("[1,2,3]", f.getvalue())

    def test_dumpToFileLikeObject(self):
        class filelike:
            def __init__(self):
                self.bytes = ''
            def write(self, bytes):
                self.bytes += bytes
        f = filelike()
        ujson.dump([1, 2, 3], f)
        self.assertEquals("[1,2,3]", f.bytes)

    def test_dumpFileArgsError(self):
        try:
            ujson.dump([], '')
        except TypeError:
            pass
        else:
            assert False, 'expected TypeError'

    def test_loadFile(self):
        f = StringIO.StringIO("[1,2,3,4]")
        self.assertEquals([1, 2, 3, 4], ujson.load(f))
        f = StringIO.StringIO("[1,2,3,4]")
        assert_array_equal(np.array([1, 2, 3, 4]), ujson.load(f, numpy=True))

    def test_loadFileLikeObject(self):
        class filelike:
            def read(self):
                try:
                    self.end
                except AttributeError:
                    self.end = True
                    return "[1,2,3,4]"
        f = filelike()
        self.assertEquals([1, 2, 3, 4], ujson.load(f))
        f = filelike()
        assert_array_equal(np.array([1, 2, 3, 4]), ujson.load(f, numpy=True))

    def test_loadFileArgsError(self):
        try:
            ujson.load("[]")
        except TypeError:
            pass
        else:
            assert False, "expected TypeError"

    def test_version(self):
        assert re.match(r'^\d+\.\d+(\.\d+)?$', ujson.__version__), \
               "ujson.__version__ must be a string like '1.4.0'"

    def test_encodeNumericOverflow(self):
        try:
            ujson.encode(12839128391289382193812939)
        except OverflowError:
            pass
        else:
            assert False, "expected OverflowError"

    def test_encodeNumericOverflowNested(self):
        for n in xrange(0, 100):
            class Nested:
                x = 12839128391289382193812939

            nested = Nested()

            try:
                ujson.encode(nested)
            except OverflowError:
                pass
            else:
                assert False, "expected OverflowError"

    def test_decodeNumberWith32bitSignBit(self):
        #Test that numbers that fit within 32 bits but would have the
        # sign bit set (2**31 <= x < 2**32) are decoded properly.
        boundary1 = 2**31
        boundary2 = 2**32
        docs = (
            '{"id": 3590016419}',
            '{"id": %s}' % 2**31,
            '{"id": %s}' % 2**32,
            '{"id": %s}' % ((2**32)-1),
        )
        results = (3590016419, 2**31, 2**32, 2**32-1)
        for doc,result in zip(docs, results):
            self.assertEqual(ujson.decode(doc)['id'], result)

    def test_encodeBigEscape(self):
        for x in xrange(10):
            if py3compat.PY3:
                base = '\u00e5'.encode('utf-8')
            else:
                base = "\xc3\xa5"
            input = base * 1024 * 1024 * 2
            output = ujson.encode(input)

    def test_decodeBigEscape(self):
        for x in xrange(10):
            if py3compat.PY3:
                base = '\u00e5'.encode('utf-8')
            else:
                base = "\xc3\xa5"
            quote = py3compat.str_to_bytes("\"")
            input = quote + (base * 1024 * 1024 * 2) + quote
            output = ujson.decode(input)

    def test_toDict(self):
        d = {u"key": 31337}

        class DictTest:
            def toDict(self):
                return d

        o = DictTest()
        output = ujson.encode(o)
        dec = ujson.decode(output)
        self.assertEquals(dec, d)


class NumpyJSONTests(TestCase):

    def testBool(self):
        b = np.bool(True)
        self.assertEqual(ujson.decode(ujson.encode(b)), b)

    def testBoolArray(self):
        inpt = np.array([True, False, True, True, False, True, False , False],
                         dtype=np.bool)
        outp = np.array(ujson.decode(ujson.encode(inpt)), dtype=np.bool)
        assert_array_equal(inpt, outp)

    def testInt(self):
        num = np.int(2562010)
        self.assertEqual(np.int(ujson.decode(ujson.encode(num))), num)

        num = np.int8(127)
        self.assertEqual(np.int8(ujson.decode(ujson.encode(num))), num)

        num = np.int16(2562010)
        self.assertEqual(np.int16(ujson.decode(ujson.encode(num))), num)

        num = np.int32(2562010)
        self.assertEqual(np.int32(ujson.decode(ujson.encode(num))), num)

        num = np.int64(2562010)
        self.assertEqual(np.int64(ujson.decode(ujson.encode(num))), num)

        num = np.uint8(255)
        self.assertEqual(np.uint8(ujson.decode(ujson.encode(num))), num)

        num = np.uint16(2562010)
        self.assertEqual(np.uint16(ujson.decode(ujson.encode(num))), num)

        num = np.uint32(2562010)
        self.assertEqual(np.uint32(ujson.decode(ujson.encode(num))), num)

        num = np.uint64(2562010)
        self.assertEqual(np.uint64(ujson.decode(ujson.encode(num))), num)

    def testIntArray(self):
        arr = np.arange(100, dtype=np.int)
        dtypes = (np.int, np.int8, np.int16, np.int32, np.int64,
                  np.uint, np.uint8, np.uint16, np.uint32, np.uint64)
        for dtype in dtypes:
            inpt = arr.astype(dtype)
            outp = np.array(ujson.decode(ujson.encode(inpt)), dtype=dtype)
            assert_array_equal(inpt, outp)

    def testIntMax(self):
        num = np.int(np.iinfo(np.int).max)
        self.assertEqual(np.int(ujson.decode(ujson.encode(num))), num)

        num = np.int8(np.iinfo(np.int8).max)
        self.assertEqual(np.int8(ujson.decode(ujson.encode(num))), num)

        num = np.int16(np.iinfo(np.int16).max)
        self.assertEqual(np.int16(ujson.decode(ujson.encode(num))), num)

        num = np.int32(np.iinfo(np.int32).max)
        self.assertEqual(np.int32(ujson.decode(ujson.encode(num))), num)

        num = np.uint8(np.iinfo(np.uint8).max)
        self.assertEqual(np.uint8(ujson.decode(ujson.encode(num))), num)

        num = np.uint16(np.iinfo(np.uint16).max)
        self.assertEqual(np.uint16(ujson.decode(ujson.encode(num))), num)

        num = np.uint32(np.iinfo(np.uint32).max)
        self.assertEqual(np.uint32(ujson.decode(ujson.encode(num))), num)

        if platform.architecture()[0] != '32bit':
            num = np.int64(np.iinfo(np.int64).max)
            self.assertEqual(np.int64(ujson.decode(ujson.encode(num))), num)

            # uint64 max will always overflow as it's encoded to signed
            num = np.uint64(np.iinfo(np.int64).max)
            self.assertEqual(np.uint64(ujson.decode(ujson.encode(num))), num)

    def testFloat(self):
        num = np.float(256.2013)
        self.assertEqual(np.float(ujson.decode(ujson.encode(num))), num)

        num = np.float32(256.2013)
        self.assertEqual(np.float32(ujson.decode(ujson.encode(num))), num)

        num = np.float64(256.2013)
        self.assertEqual(np.float64(ujson.decode(ujson.encode(num))), num)

    def testFloatArray(self):
        arr = np.arange(12.5, 185.72, 1.7322, dtype=np.float)
        dtypes = (np.float, np.float32, np.float64)

        for dtype in dtypes:
            inpt = arr.astype(dtype)
            outp = np.array(ujson.decode(ujson.encode(inpt, double_precision=15)), dtype=dtype)
            assert_array_almost_equal_nulp(inpt, outp)

    def testFloatMax(self):
        num = np.float(np.finfo(np.float).max/10)
        assert_approx_equal(np.float(ujson.decode(ujson.encode(num, double_precision=15))), num, 15)

        num = np.float32(np.finfo(np.float32).max/10)
        assert_approx_equal(np.float32(ujson.decode(ujson.encode(num, double_precision=15))), num, 15)

        num = np.float64(np.finfo(np.float64).max/10)
        assert_approx_equal(np.float64(ujson.decode(ujson.encode(num, double_precision=15))), num, 15)

    def testArrays(self):
        arr = np.arange(100);

        arr = arr.reshape((10, 10))
        assert_array_equal(np.array(ujson.decode(ujson.encode(arr))), arr)
        assert_array_equal(ujson.decode(ujson.encode(arr), numpy=True), arr)

        arr = arr.reshape((5, 5, 4))
        assert_array_equal(np.array(ujson.decode(ujson.encode(arr))), arr)
        assert_array_equal(ujson.decode(ujson.encode(arr), numpy=True), arr)

        arr = arr.reshape((100, 1))
        assert_array_equal(np.array(ujson.decode(ujson.encode(arr))), arr)
        assert_array_equal(ujson.decode(ujson.encode(arr), numpy=True), arr)

        arr = np.arange(96);
        arr = arr.reshape((2, 2, 2, 2, 3, 2))
        assert_array_equal(np.array(ujson.decode(ujson.encode(arr))), arr)
        assert_array_equal(ujson.decode(ujson.encode(arr), numpy=True), arr)

        l = ['a', list(), dict(), dict(), list(),
             42, 97.8, ['a', 'b'], {'key': 'val'}]
        arr = np.array(l)
        assert_array_equal(np.array(ujson.decode(ujson.encode(arr))), arr)

        arr = np.arange(100.202, 200.202, 1, dtype=np.float32);
        arr = arr.reshape((5, 5, 4))
        outp = np.array(ujson.decode(ujson.encode(arr)), dtype=np.float32)
        assert_array_almost_equal_nulp(arr, outp)
        outp = ujson.decode(ujson.encode(arr), numpy=True, dtype=np.float32)
        assert_array_almost_equal_nulp(arr, outp)

    def testArrayNumpyExcept(self):

        input = ujson.dumps([42, {}, 'a'])
        try:
            ujson.decode(input, numpy=True)
            assert False, "Expected exception!"
        except(TypeError):
            pass
        except:
            assert False, "Wrong exception"

        input = ujson.dumps(['a', 'b', [], 'c'])
        try:
            ujson.decode(input, numpy=True)
            assert False, "Expected exception!"
        except(ValueError):
            pass
        except:
            assert False, "Wrong exception"

        input = ujson.dumps([['a'], 42])
        try:
            ujson.decode(input, numpy=True)
            assert False, "Expected exception!"
        except(ValueError):
            pass
        except:
            assert False, "Wrong exception"

        input = ujson.dumps([42, ['a'], 42])
        try:
            ujson.decode(input, numpy=True)
            assert False, "Expected exception!"
        except(ValueError):
            pass
        except:
            assert False, "Wrong exception"

        input = ujson.dumps([{}, []])
        try:
            ujson.decode(input, numpy=True)
            assert False, "Expected exception!"
        except(ValueError):
            pass
        except:
            assert False, "Wrong exception"

        input = ujson.dumps([42, None])
        try:
            ujson.decode(input, numpy=True)
            assert False, "Expected exception!"
        except(TypeError):
            pass
        except:
            assert False, "Wrong exception"

        input = ujson.dumps([{'a': 'b'}])
        try:
            ujson.decode(input, numpy=True, labelled=True)
            assert False, "Expected exception!"
        except(ValueError):
            pass
        except:
            assert False, "Wrong exception"

        input = ujson.dumps({'a': {'b': {'c': 42}}})
        try:
            ujson.decode(input, numpy=True, labelled=True)
            assert False, "Expected exception!"
        except(ValueError):
            pass
        except:
            assert False, "Wrong exception"

        input = ujson.dumps([{'a': 42, 'b': 23}, {'c': 17}])
        try:
            ujson.decode(input, numpy=True, labelled=True)
            assert False, "Expected exception!"
        except(ValueError):
            pass
        except:
            assert False, "Wrong exception"

    def testArrayNumpyLabelled(self):
        input = {'a': []}
        output = ujson.loads(ujson.dumps(input), numpy=True, labelled=True)
        self.assertTrue((np.empty((1, 0)) == output[0]).all())
        self.assertTrue((np.array(['a']) == output[1]).all())
        self.assertTrue(output[2] is None)

        input = [{'a': 42}]
        output = ujson.loads(ujson.dumps(input), numpy=True, labelled=True)
        self.assertTrue((np.array([42]) == output[0]).all())
        self.assertTrue(output[1] is None)
        self.assertTrue((np.array([u'a']) == output[2]).all())

        # py3 is non-determinstic on the ordering......
        if not py3compat.PY3:
            input = [{'a': 42, 'b':31}, {'a': 24, 'c': 99}, {'a': 2.4, 'b': 78}]
            output = ujson.loads(ujson.dumps(input), numpy=True, labelled=True)
            expectedvals = np.array([42, 31, 24, 99, 2.4, 78], dtype=int).reshape((3,2))
            self.assertTrue((expectedvals == output[0]).all())
            self.assertTrue(output[1] is None)
            self.assertTrue((np.array([u'a', 'b']) == output[2]).all())


            input = {1: {'a': 42, 'b':31}, 2: {'a': 24, 'c': 99}, 3: {'a': 2.4, 'b': 78}}
            output = ujson.loads(ujson.dumps(input), numpy=True, labelled=True)
            expectedvals = np.array([42, 31, 24, 99, 2.4, 78], dtype=int).reshape((3,2))
            self.assertTrue((expectedvals == output[0]).all())
            self.assertTrue((np.array(['1','2','3']) == output[1]).all())
            self.assertTrue((np.array(['a', 'b']) == output[2]).all())

class PandasJSONTests(TestCase):

    def testDataFrame(self):
        df = DataFrame([[1,2,3], [4,5,6]], index=['a', 'b'], columns=['x', 'y', 'z'])

        # column indexed
        outp = DataFrame(ujson.decode(ujson.encode(df)))
        self.assertTrue((df == outp).values.all())
        assert_array_equal(df.columns, outp.columns)
        assert_array_equal(df.index, outp.index)

        dec = _clean_dict(ujson.decode(ujson.encode(df, orient="split")))
        outp = DataFrame(**dec)
        self.assertTrue((df == outp).values.all())
        assert_array_equal(df.columns, outp.columns)
        assert_array_equal(df.index, outp.index)

        outp = DataFrame(ujson.decode(ujson.encode(df, orient="records")))
        outp.index = df.index
        self.assertTrue((df == outp).values.all())
        assert_array_equal(df.columns, outp.columns)

        outp = DataFrame(ujson.decode(ujson.encode(df, orient="values")))
        outp.index = df.index
        self.assertTrue((df.values == outp.values).all())

        outp = DataFrame(ujson.decode(ujson.encode(df, orient="index")))
        self.assertTrue((df.transpose() == outp).values.all())
        assert_array_equal(df.transpose().columns, outp.columns)
        assert_array_equal(df.transpose().index, outp.index)


    def testDataFrameNumpy(self):
        df = DataFrame([[1,2,3], [4,5,6]], index=['a', 'b'], columns=['x', 'y', 'z'])

        # column indexed
        outp = DataFrame(ujson.decode(ujson.encode(df), numpy=True))
        self.assertTrue((df == outp).values.all())
        assert_array_equal(df.columns, outp.columns)
        assert_array_equal(df.index, outp.index)

        dec = _clean_dict(ujson.decode(ujson.encode(df, orient="split"),
                          numpy=True))
        outp = DataFrame(**dec)
        self.assertTrue((df == outp).values.all())
        assert_array_equal(df.columns, outp.columns)
        assert_array_equal(df.index, outp.index)

        outp = DataFrame(ujson.decode(ujson.encode(df, orient="index"), numpy=True))
        self.assertTrue((df.transpose() == outp).values.all())
        assert_array_equal(df.transpose().columns, outp.columns)
        assert_array_equal(df.transpose().index, outp.index)

    def testDataFrameNested(self):
        df = DataFrame([[1,2,3], [4,5,6]], index=['a', 'b'], columns=['x', 'y', 'z'])

        nested = {'df1': df, 'df2': df.copy()}

        exp = {'df1': ujson.decode(ujson.encode(df)),
               'df2': ujson.decode(ujson.encode(df))}
        self.assertTrue(ujson.decode(ujson.encode(nested)) == exp)

        exp = {'df1': ujson.decode(ujson.encode(df, orient="index")),
               'df2': ujson.decode(ujson.encode(df, orient="index"))}
        self.assertTrue(ujson.decode(ujson.encode(nested, orient="index")) == exp)

        exp = {'df1': ujson.decode(ujson.encode(df, orient="records")),
               'df2': ujson.decode(ujson.encode(df, orient="records"))}
        self.assertTrue(ujson.decode(ujson.encode(nested, orient="records")) == exp)

        exp = {'df1': ujson.decode(ujson.encode(df, orient="values")),
               'df2': ujson.decode(ujson.encode(df, orient="values"))}
        self.assertTrue(ujson.decode(ujson.encode(nested, orient="values")) == exp)

        exp = {'df1': ujson.decode(ujson.encode(df, orient="split")),
               'df2': ujson.decode(ujson.encode(df, orient="split"))}
        self.assertTrue(ujson.decode(ujson.encode(nested, orient="split")) == exp)

    def testDataFrameNumpyLabelled(self):
        df = DataFrame([[1,2,3], [4,5,6]], index=['a', 'b'], columns=['x', 'y', 'z'])

        # column indexed
        outp = DataFrame(*ujson.decode(ujson.encode(df), numpy=True, labelled=True))
        self.assertTrue((df.T == outp).values.all())
        assert_array_equal(df.T.columns, outp.columns)
        assert_array_equal(df.T.index, outp.index)

        outp = DataFrame(*ujson.decode(ujson.encode(df, orient="records"), numpy=True, labelled=True))
        outp.index = df.index
        self.assertTrue((df == outp).values.all())
        assert_array_equal(df.columns, outp.columns)

        outp = DataFrame(*ujson.decode(ujson.encode(df, orient="index"), numpy=True, labelled=True))
        self.assertTrue((df == outp).values.all())
        assert_array_equal(df.columns, outp.columns)
        assert_array_equal(df.index, outp.index)

    def testSeries(self):
        s = Series([10, 20, 30, 40, 50, 60], name="series", index=[6,7,8,9,10,15])
        s.sort()

        # column indexed
        outp = Series(ujson.decode(ujson.encode(s)))
        outp.sort()
        self.assertTrue((s == outp).values.all())

        outp = Series(ujson.decode(ujson.encode(s), numpy=True))
        outp.sort()
        self.assertTrue((s == outp).values.all())

        dec = _clean_dict(ujson.decode(ujson.encode(s, orient="split")))
        outp = Series(**dec)
        self.assertTrue((s == outp).values.all())
        self.assertTrue(s.name == outp.name)

        dec = _clean_dict(ujson.decode(ujson.encode(s, orient="split"),
                          numpy=True))
        outp = Series(**dec)
        self.assertTrue((s == outp).values.all())
        self.assertTrue(s.name == outp.name)

        outp = Series(ujson.decode(ujson.encode(s, orient="records"), numpy=True))
        self.assertTrue((s == outp).values.all())

        outp = Series(ujson.decode(ujson.encode(s, orient="records")))
        self.assertTrue((s == outp).values.all())

        outp = Series(ujson.decode(ujson.encode(s, orient="values"), numpy=True))
        self.assertTrue((s == outp).values.all())

        outp = Series(ujson.decode(ujson.encode(s, orient="values")))
        self.assertTrue((s == outp).values.all())

        outp = Series(ujson.decode(ujson.encode(s, orient="index")))
        outp.sort()
        self.assertTrue((s == outp).values.all())

        outp = Series(ujson.decode(ujson.encode(s, orient="index"), numpy=True))
        outp.sort()
        self.assertTrue((s == outp).values.all())

    def testSeriesNested(self):
        s = Series([10, 20, 30, 40, 50, 60], name="series", index=[6,7,8,9,10,15])
        s.sort()

        nested = {'s1': s, 's2': s.copy()}

        exp = {'s1': ujson.decode(ujson.encode(s)),
               's2': ujson.decode(ujson.encode(s))}
        self.assertTrue(ujson.decode(ujson.encode(nested)) == exp)

        exp = {'s1': ujson.decode(ujson.encode(s, orient="split")),
               's2': ujson.decode(ujson.encode(s, orient="split"))}
        self.assertTrue(ujson.decode(ujson.encode(nested, orient="split")) == exp)

        exp = {'s1': ujson.decode(ujson.encode(s, orient="records")),
               's2': ujson.decode(ujson.encode(s, orient="records"))}
        self.assertTrue(ujson.decode(ujson.encode(nested, orient="records")) == exp)

        exp = {'s1': ujson.decode(ujson.encode(s, orient="values")),
               's2': ujson.decode(ujson.encode(s, orient="values"))}
        self.assertTrue(ujson.decode(ujson.encode(nested, orient="values")) == exp)

        exp = {'s1': ujson.decode(ujson.encode(s, orient="index")),
               's2': ujson.decode(ujson.encode(s, orient="index"))}
        self.assertTrue(ujson.decode(ujson.encode(nested, orient="index")) == exp)

    def testIndex(self):
        i = Index([23, 45, 18, 98, 43, 11], name="index")

        # column indexed
        outp = Index(ujson.decode(ujson.encode(i)))
        self.assert_(i.equals(outp))

        outp = Index(ujson.decode(ujson.encode(i), numpy=True))
        self.assert_(i.equals(outp))

        dec = _clean_dict(ujson.decode(ujson.encode(i, orient="split")))
        outp = Index(**dec)
        self.assert_(i.equals(outp))
        self.assertTrue(i.name == outp.name)

        dec = _clean_dict(ujson.decode(ujson.encode(i, orient="split"),
                          numpy=True))
        outp = Index(**dec)
        self.assert_(i.equals(outp))
        self.assertTrue(i.name == outp.name)

        outp = Index(ujson.decode(ujson.encode(i, orient="values")))
        self.assert_(i.equals(outp))

        outp = Index(ujson.decode(ujson.encode(i, orient="values"), numpy=True))
        self.assert_(i.equals(outp))

        outp = Index(ujson.decode(ujson.encode(i, orient="records")))
        self.assert_(i.equals(outp))

        outp = Index(ujson.decode(ujson.encode(i, orient="records"), numpy=True))
        self.assert_(i.equals(outp))

        outp = Index(ujson.decode(ujson.encode(i, orient="index")))
        self.assert_(i.equals(outp))

        outp = Index(ujson.decode(ujson.encode(i, orient="index"), numpy=True))
        self.assert_(i.equals(outp))

    def test_datetimeindex(self):
        from pandas.tseries.index import date_range, DatetimeIndex

        rng = date_range('1/1/2000', periods=20)

        encoded = ujson.encode(rng)
        decoded = DatetimeIndex(np.array(ujson.decode(encoded)))

        self.assert_(rng.equals(decoded))

        ts = Series(np.random.randn(len(rng)), index=rng)
        decoded = Series(ujson.decode(ujson.encode(ts)))
        idx_values = decoded.index.values.astype(np.int64)
        decoded.index = DatetimeIndex(idx_values)
        tm.assert_series_equal(ts, decoded)

    def test_decodeArrayTrailingCommaFail(self):
        input = "[31337,]"
        try:
            ujson.decode(input)
        except ValueError:
            pass
        else:
            assert False, "expected ValueError"

    def test_decodeArrayLeadingCommaFail(self):
        input = "[,31337]"
        try:
            ujson.decode(input)
        except ValueError:
            pass
        else:
            assert False, "expected ValueError"

    def test_decodeArrayOnlyCommaFail(self):
        input = "[,]"
        try:
            ujson.decode(input)
        except ValueError:
            pass
        else:
            assert False, "expected ValueError"

    def test_decodeArrayUnmatchedBracketFail(self):
        input = "[]]"
        try:
            ujson.decode(input)
        except ValueError:
            pass
        else:
            assert False, "expected ValueError"

    def test_decodeArrayEmpty(self):
        input = "[]"
        ujson.decode(input)

    def test_decodeArrayOneItem(self):
        input = "[31337]"
        ujson.decode(input)

    def test_decodeBigValue(self):
        input = "9223372036854775807"
        ujson.decode(input)

    def test_decodeSmallValue(self):
        input = "-9223372036854775808"
        ujson.decode(input)

    def test_decodeTooBigValue(self):
        try:
            input = "9223372036854775808"
            ujson.decode(input)
        except ValueError, e:
            pass
        else:
            assert False, "expected ValueError"

    def test_decodeTooSmallValue(self):
        try:
            input = "-90223372036854775809"
            ujson.decode(input)
        except ValueError,e:
            pass
        else:
            assert False, "expected ValueError"

    def test_decodeVeryTooBigValue(self):
        try:
            input = "9223372036854775808"
            ujson.decode(input)
        except ValueError:
            pass
        else:
            assert False, "expected ValueError"

    def test_decodeVeryTooSmallValue(self):
        try:
            input = "-90223372036854775809"
            ujson.decode(input)
        except ValueError:
            pass
        else:
            assert False, "expected ValueError"

    def test_decodeWithTrailingWhitespaces(self):
        input = "{}\n\t "
        ujson.decode(input)

    def test_decodeWithTrailingNonWhitespaces(self):
        try:
            input = "{}\n\t a"
            ujson.decode(input)
        except ValueError:
            pass
        else:
            assert False, "expected ValueError"

    def test_decodeArrayWithBigInt(self):
        try:
            ujson.loads('[18446098363113800555]')
        except ValueError:
            pass
        else:
            assert False, "expected ValueError"

    def test_decodeArrayFaultyUnicode(self):
        try:
            ujson.loads('[18446098363113800555]')
        except ValueError:
            pass
        else:
            assert False, "expected ValueError"


    def test_decodeFloatingPointAdditionalTests(self):
        places = 15

        self.assertAlmostEquals(-1.1234567893, ujson.loads("-1.1234567893"), places=places)
        self.assertAlmostEquals(-1.234567893, ujson.loads("-1.234567893"), places=places)
        self.assertAlmostEquals(-1.34567893, ujson.loads("-1.34567893"), places=places)
        self.assertAlmostEquals(-1.4567893, ujson.loads("-1.4567893"), places=places)
        self.assertAlmostEquals(-1.567893, ujson.loads("-1.567893"), places=places)
        self.assertAlmostEquals(-1.67893, ujson.loads("-1.67893"), places=places)
        self.assertAlmostEquals(-1.7893, ujson.loads("-1.7893"), places=places)
        self.assertAlmostEquals(-1.893, ujson.loads("-1.893"), places=places)
        self.assertAlmostEquals(-1.3, ujson.loads("-1.3"), places=places)

        self.assertAlmostEquals(1.1234567893, ujson.loads("1.1234567893"), places=places)
        self.assertAlmostEquals(1.234567893, ujson.loads("1.234567893"), places=places)
        self.assertAlmostEquals(1.34567893, ujson.loads("1.34567893"), places=places)
        self.assertAlmostEquals(1.4567893, ujson.loads("1.4567893"), places=places)
        self.assertAlmostEquals(1.567893, ujson.loads("1.567893"), places=places)
        self.assertAlmostEquals(1.67893, ujson.loads("1.67893"), places=places)
        self.assertAlmostEquals(1.7893, ujson.loads("1.7893"), places=places)
        self.assertAlmostEquals(1.893, ujson.loads("1.893"), places=places)
        self.assertAlmostEquals(1.3, ujson.loads("1.3"), places=places)

    def test_encodeBigSet(self):
        s = set()
        for x in xrange(0, 100000):
            s.add(x)
        ujson.encode(s)

    def test_encodeEmptySet(self):
        s = set()
        self.assertEquals("[]", ujson.encode(s))

    def test_encodeSet(self):
        s = set([1,2,3,4,5,6,7,8,9])
        enc = ujson.encode(s)
        dec = ujson.decode(enc)

        for v in dec:
            self.assertTrue(v in s)


"""
def test_decodeNumericIntFrcOverflow(self):
input = "X.Y"
raise NotImplementedError("Implement this test!")


def test_decodeStringUnicodeEscape(self):
input = "\u3131"
raise NotImplementedError("Implement this test!")

def test_decodeStringUnicodeBrokenEscape(self):
input = "\u3131"
raise NotImplementedError("Implement this test!")

def test_decodeStringUnicodeInvalidEscape(self):
input = "\u3131"
raise NotImplementedError("Implement this test!")

def test_decodeStringUTF8(self):
input = "someutfcharacters"
raise NotImplementedError("Implement this test!")



"""

def _clean_dict(d):
    return dict((str(k), v) for k, v in d.iteritems())

if __name__ == '__main__':
    # unittest.main()
    import nose
    # nose.runmodule(argv=[__file__,'-vvs','-x', '--ipdb-failure'],
    #                exit=False)
    nose.runmodule(argv=[__file__,'-vvs','-x','--pdb', '--pdb-failure'],
                   exit=False)
