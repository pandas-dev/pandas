import unittest
import nose
import pandas.core.encoding as en
from pandas.util import py3compat

try:
    next
except NameError:  # pragma: no cover
    # Python < 2.6
    def next(x):
        return x.next()

class TestEncoding(unittest.TestCase):
    def setUp(self):
        self.u = u"\u03c3"
        self.u_e = self.u.encode('utf-8')
        self.seq =self.assertEqual

    def tearDown(self):
        pass

    def test_decode(self):
        if py3compat.PY3:
            raise nose.SkipTest()

        self.seq(en.decode_catch_errors([]),[])
        self.seq(en.decode_catch_errors([1,2.5,True]),[1,2.5,True])
        self.seq(en.decode_catch_errors([u"1","2",3]), [u'1','2', 3])
        self.seq(en.decode_catch_errors([self.u,2,3]),[self.u, 2, 3])
        self.seq(en.decode_catch_errors([self.u_e,"a"]),[self.u, 'a'])

    def test_decode_with_nested(self):
        if py3compat.PY3:
            raise nose.SkipTest()

        self.seq(en.decode_catch_errors([self.u_e,["a"],u"a"]),[self.u, ['a'], u'a']) # ascii left alone


    def test_decode_with_immutable_seq(self):
        if py3compat.PY3:
            raise nose.SkipTest()

        # mutables are not altered
        self.assertTrue(en.decode_catch_errors((self.u_e,))==(self.u_e,))
        self.assertTrue(en.decode_catch_errors(["abc",(u"abc",)])==["abc", (u"abc",)]) # mutables not converted

    def test_decode_with_nested_and_dicts(self):
        if py3compat.PY3:
            raise nose.SkipTest()

        self.seq(en.decode_catch_errors({"a":"b"}), {u'a': u'b'})

        r=[u'a',u'b', 1, 2.5, True, {u'a': u'b'},
           [u'a',  u'b',  1,  2.5,  True,  {u'a': u'b'},
            [u'a', u'b', 1, 2.5, True, {u'a': u'b'}],
            [u'a', u'b', 1, 2.5, True, {u'a': u'b'}]]]

        self.seq(en.decode_catch_errors(["a",u"b",1,2.5,True,{"a":"b"},
                                 ["a",u"b",1,2.5,True,{"a":"b"},
                                  ["a",u"b",1,2.5,True,{"a":"b"}],
                                  ["a",u"b",1,2.5,True,{"a":"b"}]]]),r)

        r= [{"k": [self.u, [self.u,1], u'b']}]
        self.seq(en.decode_catch_errors([{"k":[self.u_e,[self.u_e,1],u"b"]}]),r)

    def test_decode_non_seq(self):
        self.seq(en.decode_catch_errors("abcd"),"abcd")
        self.seq(en.decode_catch_errors(u"abcd"),u"abcd")

    def test_extract_text(self):
        if py3compat.PY3:
            raise nose.SkipTest()

        # test with self.seq, pure str, pure unicode
        g=en._extract_txt_from_obj(u"abcd")

        try:
            next(g)
        except StopIteration:
            pass
        else:
            self.fail("erroneous yield")
        #        self.assertRaises(StopIteration,next(g))

        g=en._extract_txt_from_obj("abcd")
        self.seq(next(g),"abcd")

        g=en._extract_txt_from_obj("\xcc")
        self.seq(next(g),"\xcc")

        g=en._extract_txt_from_obj(["abcd","\xcc"])
        self.seq(next(g),"abcd")
        self.seq(next(g),"\xcc")

    def test_recursion_limit_safe(self):
        "Test against recursion limit"
        import sys

        a=["a"]
        for i in range(sys.getrecursionlimit()+1):
            a=["a",a]

        try:
            en.decode_catch_errors(a)
        except RuntimeError:
            self.fail("en.decode_self.seq() Implementation cannot handle deeply-nested self.sequences")

    def test_ordered_dict_key_ordering(self):
        "Test That OrderedDicts keep their key ordering"
        import string,random,sys

        if sys.version_info[:2]<(2,7):
            raise nose.SkipTest

        from collections import OrderedDict
        self.seq=self.assertEqual

        for i in range(100):
            keys=[string.ascii_letters[random.randint(1,20)] for x in range(20)]
            d=OrderedDict.fromkeys(keys)
            # after decoding, is the order of keys is maintained?
            self.seq( en.decode_catch_errors([d])[0].keys(),map(unicode,d.keys()))

    def test_detect_text_enc(self):
        import string
        if en._can_import("chardet"):
            res=en._detect_encoding(string.ascii_letters,min_cnt=10)
            self.assertTrue(isinstance(res,dict))
            self.assertTrue('confidence' in res and 'encoding' in res) # keys in result dict
            res=en._detect_encoding("a") # not enough confidence, return empty
            self.assertTrue(res=={})

    def test_detector_detects_enc(self):
        s='\xf9\xec\xe5\xed \xf8\xe1 \xf9\xe5\xe1\xea'+\
          '\xf6\xe9\xf4\xe5\xf8\xe4 \xf0\xe7\xee\xe3\xfa'

        if en._can_import("chardet"):
            res=en._detect_encoding(s,min_cnt=0)
            self.assertTrue(isinstance(res,dict))
            self.assertTrue('confidence' in res and 'encoding' in res) # keys in result dict
            self.assertEqual(res['encoding'],"windows-1255") # keys in result dict


    def test_text_extract_limit_iter(self):
        if en._can_import("chardet"):
            seq=["a","a","b"]
            for x in en._extract_txt_from_obj(seq,2):
                self.assertNotEqual(x,"b")
