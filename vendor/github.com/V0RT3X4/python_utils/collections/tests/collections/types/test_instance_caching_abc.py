#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 19 14:02:03 2018
@author: richard
"""
import unittest

from vortexa_utils.collections.types.instance_caching_abc import (
        InstanceCachingABC,
        instance_caching)


class InstanceCachingABCTests(unittest.TestCase):

    def register_class(self, klass):
        setattr(self, klass.__name__, klass)
        return klass

    def setUp(self):
        @self.register_class
        class Foo(object, metaclass=InstanceCachingABC):
            pass

        @self.register_class
        class Bar(object):
            pass

    def test_signiture(self):
        self.assertEqual(repr(self.Foo), repr(self.Bar).replace('Bar', 'Foo'))

    def test_instance_cache(self):
        # no instances
        self.assertFalse(list(self.Foo))

        # one instance
        foo = self.Foo()
        foos = list(self.Foo)
        self.assertEqual(len(foos), 1)
        klass_name, instance = foos[0]
        self.assertEqual(instance, foo)
        self.assertEqual(klass_name, 'Foo')

        # more instances
        foo2 = self.Foo()
        foos = list(self.Foo)
        self.assertEqual(len(foos), 2)
        klass_name, instance = foos[-1]
        self.assertEqual(instance, foo2)
        self.assertEqual(klass_name, 'Foo')


class InstanceCachingDecoratorTests(InstanceCachingABCTests):

    def setUp(self):
        register = self.register_class

        @register
        class Foo(object):
            pass

        self._Foo = Foo
        self.Foo = Foo = instance_caching(Foo)

        @register
        class Bar(Foo):
            pass

        @register
        class Baz(Bar):
            pass

        @register
        class Bo(Foo):
            pass

        @register
        class Bill(Bo):
            pass

    def test_signiture(self):
        self.assertEqual(repr(self.Foo), repr(self._Foo))

    def test_list_subclasses(self):
        self.assertEqual(
                set(self.Foo._allsubclasses()),
                set((self.Foo, self.Bar, self.Baz, self.Bo, self.Bill))
            )
        self.assertEqual(
                set(self.Bar._allsubclasses()),
                set((self.Bar, self.Baz))
            )
        self.assertEqual(
                set(self.Bo._allsubclasses()),
                set((self.Bill, self.Bo))
            )

    def test_instance_cache(self):
        super().test_instance_cache()
        # no instances in subclasses
        for klass in self.Bar._allsubclasses():
            self.assertFalse(list(klass))

        for klass in self.Bo._allsubclasses():
            self.assertFalse(list(klass))

        self.assertEqual(len(list(self.Foo)), 2)
        # one instance
        bar = self.Bar()
        foos = list(self.Foo)
        bars = list(self.Bar)
        self.assertEqual(len(foos), 3)
        self.assertEqual(len(bars), 1)
        klass_name, instance = bars[0]
        self.assertEqual(instance, bar)
        self.assertEqual(klass_name, 'Bar')

        baz = self.Baz()
        foos = list(self.Foo)
        bars = list(self.Bar)
        bazs = list(self.Baz)
        self.assertEqual(len(foos), 4)
        self.assertEqual(len(bars), 2)
        self.assertEqual(len(bazs), 1)
        klass_name, instance = bazs[0]
        self.assertEqual(instance, baz)
        self.assertEqual(klass_name, 'Baz')

        for klass in self.Bo._allsubclasses():
            self.assertFalse(list(klass))
