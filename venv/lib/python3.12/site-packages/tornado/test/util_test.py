import re
import sys
import datetime
import textwrap
import unittest

import tornado
from tornado.escape import utf8
from tornado.util import (
    raise_exc_info,
    Configurable,
    exec_in,
    ArgReplacer,
    timedelta_to_seconds,
    import_object,
    re_unescape,
)

from typing import cast, Dict, Any


class RaiseExcInfoTest(unittest.TestCase):
    def test_two_arg_exception(self):
        # This test would fail on python 3 if raise_exc_info were simply
        # a three-argument raise statement, because TwoArgException
        # doesn't have a "copy constructor"
        class TwoArgException(Exception):
            def __init__(self, a, b):
                super().__init__()
                self.a, self.b = a, b

        try:
            raise TwoArgException(1, 2)
        except TwoArgException:
            exc_info = sys.exc_info()
        try:
            raise_exc_info(exc_info)
            self.fail("didn't get expected exception")
        except TwoArgException as e:
            self.assertIs(e, exc_info[1])


class TestConfigurable(Configurable):
    @classmethod
    def configurable_base(cls):
        return TestConfigurable

    @classmethod
    def configurable_default(cls):
        return TestConfig1


class TestConfig1(TestConfigurable):
    def initialize(self, pos_arg=None, a=None):
        self.a = a
        self.pos_arg = pos_arg


class TestConfig2(TestConfigurable):
    def initialize(self, pos_arg=None, b=None):
        self.b = b
        self.pos_arg = pos_arg


class TestConfig3(TestConfigurable):
    # TestConfig3 is a configuration option that is itself configurable.
    @classmethod
    def configurable_base(cls):
        return TestConfig3

    @classmethod
    def configurable_default(cls):
        return TestConfig3A


class TestConfig3A(TestConfig3):
    def initialize(self, a=None):
        self.a = a


class TestConfig3B(TestConfig3):
    def initialize(self, b=None):
        self.b = b


class ConfigurableTest(unittest.TestCase):
    def setUp(self):
        self.saved = TestConfigurable._save_configuration()
        self.saved3 = TestConfig3._save_configuration()

    def tearDown(self):
        TestConfigurable._restore_configuration(self.saved)
        TestConfig3._restore_configuration(self.saved3)

    def checkSubclasses(self):
        # no matter how the class is configured, it should always be
        # possible to instantiate the subclasses directly
        self.assertIsInstance(TestConfig1(), TestConfig1)
        self.assertIsInstance(TestConfig2(), TestConfig2)

        obj = TestConfig1(a=1)
        self.assertEqual(obj.a, 1)
        obj2 = TestConfig2(b=2)
        self.assertEqual(obj2.b, 2)

    def test_default(self):
        # In these tests we combine a typing.cast to satisfy mypy with
        # a runtime type-assertion. Without the cast, mypy would only
        # let us access attributes of the base class.
        obj = cast(TestConfig1, TestConfigurable())
        self.assertIsInstance(obj, TestConfig1)
        self.assertIsNone(obj.a)

        obj = cast(TestConfig1, TestConfigurable(a=1))
        self.assertIsInstance(obj, TestConfig1)
        self.assertEqual(obj.a, 1)

        self.checkSubclasses()

    def test_config_class(self):
        TestConfigurable.configure(TestConfig2)
        obj = cast(TestConfig2, TestConfigurable())
        self.assertIsInstance(obj, TestConfig2)
        self.assertIsNone(obj.b)

        obj = cast(TestConfig2, TestConfigurable(b=2))
        self.assertIsInstance(obj, TestConfig2)
        self.assertEqual(obj.b, 2)

        self.checkSubclasses()

    def test_config_str(self):
        TestConfigurable.configure("tornado.test.util_test.TestConfig2")
        obj = cast(TestConfig2, TestConfigurable())
        self.assertIsInstance(obj, TestConfig2)
        self.assertIsNone(obj.b)

        obj = cast(TestConfig2, TestConfigurable(b=2))
        self.assertIsInstance(obj, TestConfig2)
        self.assertEqual(obj.b, 2)

        self.checkSubclasses()

    def test_config_args(self):
        TestConfigurable.configure(None, a=3)
        obj = cast(TestConfig1, TestConfigurable())
        self.assertIsInstance(obj, TestConfig1)
        self.assertEqual(obj.a, 3)

        obj = cast(TestConfig1, TestConfigurable(42, a=4))
        self.assertIsInstance(obj, TestConfig1)
        self.assertEqual(obj.a, 4)
        self.assertEqual(obj.pos_arg, 42)

        self.checkSubclasses()
        # args bound in configure don't apply when using the subclass directly
        obj = TestConfig1()
        self.assertIsNone(obj.a)

    def test_config_class_args(self):
        TestConfigurable.configure(TestConfig2, b=5)
        obj = cast(TestConfig2, TestConfigurable())
        self.assertIsInstance(obj, TestConfig2)
        self.assertEqual(obj.b, 5)

        obj = cast(TestConfig2, TestConfigurable(42, b=6))
        self.assertIsInstance(obj, TestConfig2)
        self.assertEqual(obj.b, 6)
        self.assertEqual(obj.pos_arg, 42)

        self.checkSubclasses()
        # args bound in configure don't apply when using the subclass directly
        obj = TestConfig2()
        self.assertIsNone(obj.b)

    def test_config_multi_level(self):
        TestConfigurable.configure(TestConfig3, a=1)
        obj = cast(TestConfig3A, TestConfigurable())
        self.assertIsInstance(obj, TestConfig3A)
        self.assertEqual(obj.a, 1)

        TestConfigurable.configure(TestConfig3)
        TestConfig3.configure(TestConfig3B, b=2)
        obj2 = cast(TestConfig3B, TestConfigurable())
        self.assertIsInstance(obj2, TestConfig3B)
        self.assertEqual(obj2.b, 2)

    def test_config_inner_level(self):
        # The inner level can be used even when the outer level
        # doesn't point to it.
        obj = TestConfig3()
        self.assertIsInstance(obj, TestConfig3A)

        TestConfig3.configure(TestConfig3B)
        obj = TestConfig3()
        self.assertIsInstance(obj, TestConfig3B)

        # Configuring the base doesn't configure the inner.
        obj2 = TestConfigurable()
        self.assertIsInstance(obj2, TestConfig1)
        TestConfigurable.configure(TestConfig2)

        obj3 = TestConfigurable()
        self.assertIsInstance(obj3, TestConfig2)

        obj = TestConfig3()
        self.assertIsInstance(obj, TestConfig3B)


class UnicodeLiteralTest(unittest.TestCase):
    def test_unicode_escapes(self):
        self.assertEqual(utf8("\u00e9"), b"\xc3\xa9")


class ExecInTest(unittest.TestCase):
    def test_no_inherit_future(self):
        # Two files: the first has "from __future__ import annotations", and it executes the second
        # which doesn't. The second file should not be affected by the first's __future__ imports.
        #
        # The annotations future became available in python 3.7 but has been replaced by PEP 649, so
        # it should remain supported but off-by-default for the foreseeable future.
        code1 = textwrap.dedent(
            """
            from __future__ import annotations
            from tornado.util import exec_in

            exec_in(code2, globals())
            """
        )

        code2 = textwrap.dedent(
            """
            def f(x: int) -> int:
                return x + 1
            output[0] = f.__annotations__
            """
        )

        # Make a mutable container to pass the result back to the caller
        output = [None]
        exec_in(code1, dict(code2=code2, output=output))
        # If the annotations future were in effect, these would be strings instead of the int type
        # object.
        self.assertEqual(output[0], {"x": int, "return": int})


class ArgReplacerTest(unittest.TestCase):
    def setUp(self):
        def function(x, y, callback=None, z=None):
            pass

        self.replacer = ArgReplacer(function, "callback")

    def test_omitted(self):
        args = (1, 2)
        kwargs: Dict[str, Any] = dict()
        self.assertIsNone(self.replacer.get_old_value(args, kwargs))
        self.assertEqual(
            self.replacer.replace("new", args, kwargs),
            (None, (1, 2), dict(callback="new")),
        )

    def test_position(self):
        args = (1, 2, "old", 3)
        kwargs: Dict[str, Any] = dict()
        self.assertEqual(self.replacer.get_old_value(args, kwargs), "old")
        self.assertEqual(
            self.replacer.replace("new", args, kwargs),
            ("old", [1, 2, "new", 3], dict()),
        )

    def test_keyword(self):
        args = (1,)
        kwargs = dict(y=2, callback="old", z=3)
        self.assertEqual(self.replacer.get_old_value(args, kwargs), "old")
        self.assertEqual(
            self.replacer.replace("new", args, kwargs),
            ("old", (1,), dict(y=2, callback="new", z=3)),
        )


class TimedeltaToSecondsTest(unittest.TestCase):
    def test_timedelta_to_seconds(self):
        time_delta = datetime.timedelta(hours=1)
        self.assertEqual(timedelta_to_seconds(time_delta), 3600.0)


class ImportObjectTest(unittest.TestCase):
    def test_import_member(self):
        self.assertIs(import_object("tornado.escape.utf8"), utf8)

    def test_import_member_unicode(self):
        self.assertIs(import_object("tornado.escape.utf8"), utf8)

    def test_import_module(self):
        self.assertIs(import_object("tornado.escape"), tornado.escape)

    def test_import_module_unicode(self):
        # The internal implementation of __import__ differs depending on
        # whether the thing being imported is a module or not.
        # This variant requires a byte string in python 2.
        self.assertIs(import_object("tornado.escape"), tornado.escape)


class ReUnescapeTest(unittest.TestCase):
    def test_re_unescape(self):
        test_strings = ("/favicon.ico", "index.html", "Hello, World!", "!$@#%;")
        for string in test_strings:
            self.assertEqual(string, re_unescape(re.escape(string)))

    def test_re_unescape_raises_error_on_invalid_input(self):
        with self.assertRaises(ValueError):
            re_unescape("\\d")
        with self.assertRaises(ValueError):
            re_unescape("\\b")
        with self.assertRaises(ValueError):
            re_unescape("\\Z")


class VersionInfoTest(unittest.TestCase):
    def assert_version_info_compatible(self, version, version_info):
        # We map our version identifier string (a subset of
        # https://packaging.python.org/en/latest/specifications/version-specifiers/#public-version-identifiers)
        # to a 4-tuple of integers for easy comparisons. The last component is
        # 0 for a final release, negative for a pre-release, and would be positive for a
        # post-release if we did any of those. This test is not a promise that these are the
        # only formats we will ever use, but it does catch accidents like
        # https://github.com/tornadoweb/tornado/issues/3406.
        major = minor = patch = "0"
        is_pre = False
        if m := re.fullmatch(r"(\d+)\.(\d+)\.(\d+)", version):
            # Regular 3-component version number
            major, minor, patch = m.groups()
        elif m := re.fullmatch(r"(\d+)\.(\d+)", version):
            # Two-component version number, equivalent to major.minor.0
            major, minor = m.groups()
        elif m := re.fullmatch(r"(\d+)\.(\d+)\.(\d+)(?:\.dev|a|b|rc)\d+", version):
            # Pre-release 3-component version number.
            major, minor, patch = m.groups()
            is_pre = True
        elif m := re.fullmatch(r"(\d+)\.(\d+)(?:\.dev|a|b|rc)\d+", version):
            # Pre-release 2-component version number.
            major, minor = m.groups()
            is_pre = True
        else:
            self.fail(f"Unrecognized version format: {version}")

        self.assertEqual(version_info[:3], (int(major), int(minor), int(patch)))
        if is_pre:
            self.assertLess(int(version_info[3]), 0)
        else:
            self.assertEqual(int(version_info[3]), 0)

    def test_version_info_compatible(self):
        self.assert_version_info_compatible("6.5.0", (6, 5, 0, 0))
        self.assert_version_info_compatible("6.5", (6, 5, 0, 0))
        self.assert_version_info_compatible("6.5.1", (6, 5, 1, 0))
        self.assert_version_info_compatible("6.6.dev1", (6, 6, 0, -100))
        self.assert_version_info_compatible("6.6a1", (6, 6, 0, -100))
        self.assert_version_info_compatible("6.6b1", (6, 6, 0, -100))
        self.assert_version_info_compatible("6.6rc1", (6, 6, 0, -100))
        self.assertRaises(
            AssertionError, self.assert_version_info_compatible, "6.5.0", (6, 5, 0, 1)
        )
        self.assertRaises(
            AssertionError, self.assert_version_info_compatible, "6.5.0", (6, 4, 0, 0)
        )
        self.assertRaises(
            AssertionError, self.assert_version_info_compatible, "6.5.1", (6, 5, 0, 1)
        )

    def test_current_version(self):
        self.assert_version_info_compatible(tornado.version, tornado.version_info)
