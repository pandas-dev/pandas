#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 19 09:57:05 2018
@author: richard
"""
from abc import ABCMeta


class InstanceCachingABC(ABCMeta):
    """Metaclass for defining Instance Caching Abstract Base Classs (ICABC)
    Use this metaclass to create an ICABC. An ICABC will remember the instances
    created from it and can be iterated over to return all instances and sub
    class instances
    """

    def __init__(cls, name, bases, namespace):
        super().__init__(name, bases, namespace)
        cls._instances = list()

    def __call__(cls, *args, **kwargs):
        instance = super().__call__(*args, **kwargs)
        cls._instances.append(instance)
        return instance

    def _allsubclasses(cls):
        yield cls
        for subclass in cls.__subclasses__():
            yield from subclass._allsubclasses()

    # Metamethods, called on class objects:
    def __iter__(cls):
        return ((klass.__name__, instance)
                for klass in cls._allsubclasses()
                for instance in klass._instances)


def instance_caching(klass):
    class Decorated(klass, metaclass=InstanceCachingABC):
        pass

    Decorated.__name__ = klass.__name__
    Decorated.__qualname__ = klass.__qualname__
    Decorated.__module__ = klass.__module__
    return Decorated
