#!/usr/bin/python
# -*- coding: utf-8 -*-
import pandas as pd
import unittest
import warnings


class TestConfig(unittest.TestCase):
    _multiprocess_can_split_ = True

    def __init__(self, *args):
        super(TestConfig, self).__init__(*args)

        from copy import deepcopy
        self.cf = pd.core.config
        self.gc = deepcopy(getattr(self.cf, '_global_config'))
        self.do = deepcopy(getattr(self.cf, '_deprecated_options'))
        self.ro = deepcopy(getattr(self.cf, '_registered_options'))

    def setUp(self):
        setattr(self.cf, '_global_config', {})
        setattr(
            self.cf, 'options', self.cf.DictWrapper(self.cf._global_config))
        setattr(self.cf, '_deprecated_options', {})
        setattr(self.cf, '_registered_options', {})

    def tearDown(self):
        setattr(self.cf, '_global_config', self.gc)
        setattr(self.cf, '_deprecated_options', self.do)
        setattr(self.cf, '_registered_options', self.ro)

    def test_api(self):

        # the pandas object exposes the user API
        self.assertTrue(hasattr(pd, 'get_option'))
        self.assertTrue(hasattr(pd, 'set_option'))
        self.assertTrue(hasattr(pd, 'reset_option'))
        self.assertTrue(hasattr(pd, 'describe_option'))

    def test_is_one_of_factory(self):
        v = self.cf.is_one_of_factory([None, 12])

        v(12)
        v(None)
        self.assertRaises(ValueError, v, 1.1)

    def test_register_option(self):
        self.cf.register_option('a', 1, 'doc')

        # can't register an already registered option
        self.assertRaises(KeyError, self.cf.register_option, 'a', 1, 'doc')

        # can't register an already registered option
        self.assertRaises(KeyError, self.cf.register_option, 'a.b.c.d1', 1,
                          'doc')
        self.assertRaises(KeyError, self.cf.register_option, 'a.b.c.d2', 1,
                          'doc')

        # no python keywords
        self.assertRaises(ValueError, self.cf.register_option, 'for', 0)
        self.assertRaises(ValueError, self.cf.register_option, 'a.for.b', 0)
        # must be valid identifier (ensure attribute access works)
        self.assertRaises(ValueError, self.cf.register_option,
                          'Oh my Goddess!', 0)

        # we can register options several levels deep
        # without predefining the intermediate steps
        # and we can define differently named options
        # in the same namespace
        self.cf.register_option('k.b.c.d1', 1, 'doc')
        self.cf.register_option('k.b.c.d2', 1, 'doc')

    def test_describe_option(self):
        self.cf.register_option('a', 1, 'doc')
        self.cf.register_option('b', 1, 'doc2')
        self.cf.deprecate_option('b')

        self.cf.register_option('c.d.e1', 1, 'doc3')
        self.cf.register_option('c.d.e2', 1, 'doc4')
        self.cf.register_option('f', 1)
        self.cf.register_option('g.h', 1)
        self.cf.register_option('k', 2)
        self.cf.deprecate_option('g.h', rkey="k")
        self.cf.register_option('l', "foo")

        # non-existent keys raise KeyError
        self.assertRaises(KeyError, self.cf.describe_option, 'no.such.key')

        # we can get the description for any key we registered
        self.assertTrue(
            'doc' in self.cf.describe_option('a', _print_desc=False))
        self.assertTrue(
            'doc2' in self.cf.describe_option('b', _print_desc=False))
        self.assertTrue(
            'precated' in self.cf.describe_option('b', _print_desc=False))

        self.assertTrue(
            'doc3' in self.cf.describe_option('c.d.e1', _print_desc=False))
        self.assertTrue(
            'doc4' in self.cf.describe_option('c.d.e2', _print_desc=False))

        # if no doc is specified we get a default message
        # saying "description not available"
        self.assertTrue(
            'vailable' in self.cf.describe_option('f', _print_desc=False))
        self.assertTrue(
            'vailable' in self.cf.describe_option('g.h', _print_desc=False))
        self.assertTrue(
            'precated' in self.cf.describe_option('g.h', _print_desc=False))
        self.assertTrue(
            'k' in self.cf.describe_option('g.h', _print_desc=False))

        # default is reported
        self.assertTrue(
            'foo' in self.cf.describe_option('l', _print_desc=False))
        # current value is reported
        self.assertFalse(
            'bar' in self.cf.describe_option('l', _print_desc=False))
        self.cf.set_option("l", "bar")
        self.assertTrue(
            'bar' in self.cf.describe_option('l', _print_desc=False))

    def test_case_insensitive(self):
        self.cf.register_option('KanBAN', 1, 'doc')

        self.assertTrue(
            'doc' in self.cf.describe_option('kanbaN', _print_desc=False))
        self.assertEqual(self.cf.get_option('kanBaN'), 1)
        self.cf.set_option('KanBan', 2)
        self.assertEqual(self.cf.get_option('kAnBaN'), 2)

        # gets of non-existent keys fail
        self.assertRaises(KeyError, self.cf.get_option, 'no_such_option')
        self.cf.deprecate_option('KanBan')

        self.assertTrue(self.cf._is_deprecated('kAnBaN'))

    def test_get_option(self):
        self.cf.register_option('a', 1, 'doc')
        self.cf.register_option('b.c', 'hullo', 'doc2')
        self.cf.register_option('b.b', None, 'doc2')

        # gets of existing keys succeed
        self.assertEqual(self.cf.get_option('a'), 1)
        self.assertEqual(self.cf.get_option('b.c'), 'hullo')
        self.assertTrue(self.cf.get_option('b.b') is None)

        # gets of non-existent keys fail
        self.assertRaises(KeyError, self.cf.get_option, 'no_such_option')

    def test_set_option(self):
        self.cf.register_option('a', 1, 'doc')
        self.cf.register_option('b.c', 'hullo', 'doc2')
        self.cf.register_option('b.b', None, 'doc2')

        self.assertEqual(self.cf.get_option('a'), 1)
        self.assertEqual(self.cf.get_option('b.c'), 'hullo')
        self.assertTrue(self.cf.get_option('b.b') is None)

        self.cf.set_option('a', 2)
        self.cf.set_option('b.c', 'wurld')
        self.cf.set_option('b.b', 1.1)

        self.assertEqual(self.cf.get_option('a'), 2)
        self.assertEqual(self.cf.get_option('b.c'), 'wurld')
        self.assertEqual(self.cf.get_option('b.b'), 1.1)

        self.assertRaises(KeyError, self.cf.set_option, 'no.such.key', None)

    def test_set_option_empty_args(self):
        self.assertRaises(ValueError, self.cf.set_option)

    def test_set_option_uneven_args(self):
        self.assertRaises(ValueError, self.cf.set_option, 'a.b', 2, 'b.c')

    def test_set_option_invalid_single_argument_type(self):
        self.assertRaises(ValueError, self.cf.set_option, 2)

    def test_set_option_multiple(self):
        self.cf.register_option('a', 1, 'doc')
        self.cf.register_option('b.c', 'hullo', 'doc2')
        self.cf.register_option('b.b', None, 'doc2')

        self.assertEqual(self.cf.get_option('a'), 1)
        self.assertEqual(self.cf.get_option('b.c'), 'hullo')
        self.assertTrue(self.cf.get_option('b.b') is None)

        self.cf.set_option('a', '2', 'b.c', None, 'b.b', 10.0)

        self.assertEqual(self.cf.get_option('a'), '2')
        self.assertTrue(self.cf.get_option('b.c') is None)
        self.assertEqual(self.cf.get_option('b.b'), 10.0)

    def test_validation(self):
        self.cf.register_option('a', 1, 'doc', validator=self.cf.is_int)
        self.cf.register_option('b.c', 'hullo', 'doc2',
                                validator=self.cf.is_text)
        self.assertRaises(ValueError, self.cf.register_option, 'a.b.c.d2',
                          'NO', 'doc', validator=self.cf.is_int)

        self.cf.set_option('a', 2)  # int is_int
        self.cf.set_option('b.c', 'wurld')  # str is_str

        self.assertRaises(
            ValueError, self.cf.set_option, 'a', None)  # None not is_int
        self.assertRaises(ValueError, self.cf.set_option, 'a', 'ab')
        self.assertRaises(ValueError, self.cf.set_option, 'b.c', 1)

        validator = self.cf.is_one_of_factory([None, self.cf.is_callable])
        self.cf.register_option('b', lambda: None, 'doc',
                                validator=validator)
        self.cf.set_option('b', '%.1f'.format)  # Formatter is callable
        self.cf.set_option('b', None)  # Formatter is none (default)
        self.assertRaises(ValueError, self.cf.set_option, 'b', '%.1f')

    def test_reset_option(self):
        self.cf.register_option('a', 1, 'doc', validator=self.cf.is_int)
        self.cf.register_option('b.c', 'hullo', 'doc2',
                                validator=self.cf.is_str)
        self.assertEqual(self.cf.get_option('a'), 1)
        self.assertEqual(self.cf.get_option('b.c'), 'hullo')

        self.cf.set_option('a', 2)
        self.cf.set_option('b.c', 'wurld')
        self.assertEqual(self.cf.get_option('a'), 2)
        self.assertEqual(self.cf.get_option('b.c'), 'wurld')

        self.cf.reset_option('a')
        self.assertEqual(self.cf.get_option('a'), 1)
        self.assertEqual(self.cf.get_option('b.c'), 'wurld')
        self.cf.reset_option('b.c')
        self.assertEqual(self.cf.get_option('a'), 1)
        self.assertEqual(self.cf.get_option('b.c'), 'hullo')

    def test_reset_option_all(self):
        self.cf.register_option('a', 1, 'doc', validator=self.cf.is_int)
        self.cf.register_option('b.c', 'hullo', 'doc2',
                                validator=self.cf.is_str)
        self.assertEqual(self.cf.get_option('a'), 1)
        self.assertEqual(self.cf.get_option('b.c'), 'hullo')

        self.cf.set_option('a', 2)
        self.cf.set_option('b.c', 'wurld')
        self.assertEqual(self.cf.get_option('a'), 2)
        self.assertEqual(self.cf.get_option('b.c'), 'wurld')

        self.cf.reset_option("all")
        self.assertEqual(self.cf.get_option('a'), 1)
        self.assertEqual(self.cf.get_option('b.c'), 'hullo')

    def test_deprecate_option(self):
        # we can deprecate non-existent options
        self.cf.deprecate_option('foo')

        self.assertTrue(self.cf._is_deprecated('foo'))
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter('always')
            try:
                self.cf.get_option('foo')
            except KeyError:
                pass
            else:
                self.fail("Nonexistent option didn't raise KeyError")

            self.assertEqual(len(w), 1)  # should have raised one warning
            self.assertTrue(
                'deprecated' in str(w[-1]))  # we get the default message

        self.cf.register_option('a', 1, 'doc', validator=self.cf.is_int)
        self.cf.register_option('b.c', 'hullo', 'doc2')
        self.cf.register_option('foo', 'hullo', 'doc2')

        self.cf.deprecate_option('a', removal_ver='nifty_ver')
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter('always')
            self.cf.get_option('a')

            self.assertEqual(len(w), 1)  # should have raised one warning
            self.assertTrue(
                'eprecated' in str(w[-1]))  # we get the default message
            self.assertTrue(
                'nifty_ver' in str(w[-1]))  # with the removal_ver quoted

            self.assertRaises(
                KeyError, self.cf.deprecate_option, 'a')  # can't depr. twice

        self.cf.deprecate_option('b.c', 'zounds!')
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter('always')
            self.cf.get_option('b.c')

            self.assertEqual(len(w), 1)  # should have raised one warning
            self.assertTrue(
                'zounds!' in str(w[-1]))  # we get the custom message

        # test rerouting keys
        self.cf.register_option('d.a', 'foo', 'doc2')
        self.cf.register_option('d.dep', 'bar', 'doc2')
        self.assertEqual(self.cf.get_option('d.a'), 'foo')
        self.assertEqual(self.cf.get_option('d.dep'), 'bar')

        self.cf.deprecate_option('d.dep', rkey='d.a')  # reroute d.dep to d.a
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter('always')
            self.assertEqual(self.cf.get_option('d.dep'), 'foo')

            self.assertEqual(len(w), 1)  # should have raised one warning
            self.assertTrue(
                'eprecated' in str(w[-1]))  # we get the custom message

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter('always')
            self.cf.set_option('d.dep', 'baz')  # should overwrite "d.a"

            self.assertEqual(len(w), 1)  # should have raised one warning
            self.assertTrue(
                'eprecated' in str(w[-1]))  # we get the custom message

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter('always')
            self.assertEqual(self.cf.get_option('d.dep'), 'baz')

            self.assertEqual(len(w), 1)  # should have raised one warning
            self.assertTrue(
                'eprecated' in str(w[-1]))  # we get the custom message

    def test_config_prefix(self):
        with self.cf.config_prefix("base"):
            self.cf.register_option('a', 1, "doc1")
            self.cf.register_option('b', 2, "doc2")
            self.assertEqual(self.cf.get_option('a'), 1)
            self.assertEqual(self.cf.get_option('b'), 2)

            self.cf.set_option('a', 3)
            self.cf.set_option('b', 4)
            self.assertEqual(self.cf.get_option('a'), 3)
            self.assertEqual(self.cf.get_option('b'), 4)

        self.assertEqual(self.cf.get_option('base.a'), 3)
        self.assertEqual(self.cf.get_option('base.b'), 4)
        self.assertTrue(
            'doc1' in self.cf.describe_option('base.a', _print_desc=False))
        self.assertTrue(
            'doc2' in self.cf.describe_option('base.b', _print_desc=False))

        self.cf.reset_option('base.a')
        self.cf.reset_option('base.b')

        with self.cf.config_prefix("base"):
            self.assertEqual(self.cf.get_option('a'), 1)
            self.assertEqual(self.cf.get_option('b'), 2)

    def test_callback(self):
        k = [None]
        v = [None]

        def callback(key):
            k.append(key)
            v.append(self.cf.get_option(key))

        self.cf.register_option('d.a', 'foo', cb=callback)
        self.cf.register_option('d.b', 'foo', cb=callback)

        del k[-1], v[-1]
        self.cf.set_option("d.a", "fooz")
        self.assertEqual(k[-1], "d.a")
        self.assertEqual(v[-1], "fooz")

        del k[-1], v[-1]
        self.cf.set_option("d.b", "boo")
        self.assertEqual(k[-1], "d.b")
        self.assertEqual(v[-1], "boo")

        del k[-1], v[-1]
        self.cf.reset_option("d.b")
        self.assertEqual(k[-1], "d.b")

    def test_set_ContextManager(self):
        def eq(val):
            self.assertEqual(self.cf.get_option("a"), val)

        self.cf.register_option('a', 0)
        eq(0)
        with self.cf.option_context("a", 15):
            eq(15)
            with self.cf.option_context("a", 25):
                eq(25)
            eq(15)
        eq(0)

        self.cf.set_option("a", 17)
        eq(17)

    def test_attribute_access(self):
        holder = []

        def f():
            options.b = 1

        def f2():
            options.display = 1

        def f3(key):
            holder.append(True)

        self.cf.register_option('a', 0)
        self.cf.register_option('c', 0, cb=f3)
        options = self.cf.options

        self.assertEqual(options.a, 0)
        with self.cf.option_context("a", 15):
            self.assertEqual(options.a, 15)

        options.a = 500
        self.assertEqual(self.cf.get_option("a"), 500)

        self.cf.reset_option("a")
        self.assertEqual(options.a, self.cf.get_option("a", 0))

        self.assertRaises(KeyError, f)
        self.assertRaises(KeyError, f2)

        # make sure callback kicks when using this form of setting
        options.c = 1
        self.assertEqual(len(holder), 1)

    def test_option_context_scope(self):
        # Ensure that creating a context does not affect the existing
        # environment as it is supposed to be used with the `with` statement.
        # See https://github.com/pydata/pandas/issues/8514

        original_value = 60
        context_value = 10
        option_name = 'a'

        self.cf.register_option(option_name, original_value)

        # Ensure creating contexts didn't affect the current context.
        ctx = self.cf.option_context(option_name, context_value)
        self.assertEqual(self.cf.get_option(option_name), original_value)

        # Ensure the correct value is available inside the context.
        with ctx:
            self.assertEqual(self.cf.get_option(option_name), context_value)

        # Ensure the current context is reset
        self.assertEqual(self.cf.get_option(option_name), original_value)
