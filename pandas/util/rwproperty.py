# Read & write properties
#
# Copyright (c) 2006 by Philipp "philiKON" von Weitershausen
#                       philikon@philikon.de
#
# Freely distributable under the terms of the Zope Public License, v2.1.
#
# See rwproperty.txt for detailed explanations
#
import sys

__all__ = ['getproperty', 'setproperty', 'delproperty']

class rwproperty(object):

    def __new__(cls, func):
        name = func.__name__

        # ugly, but common hack
        frame = sys._getframe(1)
        locals = frame.f_locals

        if name not in locals:
            return cls.createProperty(func)

        oldprop = locals[name]
        if isinstance(oldprop, property):
            return cls.enhanceProperty(oldprop, func)

        raise TypeError("read & write properties cannot be mixed with "
                        "other attributes except regular property objects.")

    # this might not be particularly elegant, but it's easy on the eyes

    @staticmethod
    def createProperty(func):
        raise NotImplementedError

    @staticmethod
    def enhanceProperty(oldprop, func):
        raise NotImplementedError

class getproperty(rwproperty):

    @staticmethod
    def createProperty(func):
        return property(func)

    @staticmethod
    def enhanceProperty(oldprop, func):
        return property(func, oldprop.fset, oldprop.fdel)

class setproperty(rwproperty):

    @staticmethod
    def createProperty(func):
        return property(None, func)

    @staticmethod
    def enhanceProperty(oldprop, func):
        return property(oldprop.fget, func, oldprop.fdel)

class delproperty(rwproperty):

    @staticmethod
    def createProperty(func):
        return property(None, None, func)

    @staticmethod
    def enhanceProperty(oldprop, func):
        return property(oldprop.fget, oldprop.fset, func)

if __name__ == "__main__":
    import doctest
    doctest.testfile('rwproperty.txt')
