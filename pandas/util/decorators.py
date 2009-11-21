"""
Pierre G-M's caching decorators
"""

import warnings

__all__ = ['resettable_cache','cache_readonly', 'cache_writable']

#-------------------------------------------------------------------------------
# Pierre G-M's caching decorators

class CacheWriteWarning(UserWarning):
    pass


class ResettableCache(dict):
    """
    Dictionary whose elements mey depend one from another.

    If entry `B` depends on entry `A`, changing the values of entry `A` will
    reset the value of entry `B` to a default (None); deleteing entry `A` will
    delete entry `B`.  The connections between entries are stored in a
    `_resetdict` private attribute.

    Parameters
    ----------
    reset : dictionary, optional
        An optional dictionary, associated a sequence of entries to any key
        of the object.
    items : var, optional
        An optional dictionary used to initialize the dictionary

    Examples
    --------
    >>> reset = dict(a=('b',), b=('c',))
    >>> cache = resettable_cache(a=0, b=1, c=2, reset=reset)
    >>> assert_equal(cache, dict(a=0, b=1, c=2))

    >>> print "Try resetting a"
    >>> cache['a'] = 1
    >>> assert_equal(cache, dict(a=1, b=None, c=None))
    >>> cache['c'] = 2
    >>> assert_equal(cache, dict(a=1, b=None, c=2))
    >>> cache['b'] = 0
    >>> assert_equal(cache, dict(a=1, b=0, c=None))

    >>> print "Try deleting b"
    >>> del(cache['a'])
    >>> assert_equal(cache, {})
    """

    def __init__(self, reset=None, **items):
        self._resetdict = reset or {}
        dict.__init__(self, **items)

    def __setitem__(self, key, value):
        dict.__setitem__(self, key, value)
        for mustreset in self._resetdict.get(key, []):
            self[mustreset] = None

    def __delitem__(self, key):
        dict.__delitem__(self, key)
        for mustreset in self._resetdict.get(key, []):
            del(self[mustreset])

resettable_cache = ResettableCache

class CachedAttribute(object):

    def __init__(self, func, cachename=None, resetlist=None):
        self.fget = func
        self.name = func.__name__
        self.cachename = cachename or '_cache'
        self.resetlist = resetlist or ()

    def __get__(self, obj, type=None):
        if obj is None:
            return self.fget
        # Get the cache or set a default one if needed
        _cachename = self.cachename
        _cache = getattr(obj, _cachename, None)
        if _cache is None:
            setattr(obj, _cachename, resettable_cache())
            _cache = getattr(obj, _cachename)
        # Get the name of the attribute to set and cache
        name = self.name
        _cachedval = _cache.get(name, None)
#        print "[_cachedval=%s]" % _cachedval
        if _cachedval is None:
            # Call the "fget" function
            _cachedval = self.fget(obj)
            # Set the attribute in obj
#            print "Setting %s in cache to %s" % (name, _cachedval)
            try:
                _cache[name] = _cachedval
            except KeyError:
                setattr(_cache, name, _cachedval)
            # Update the reset list if needed (and possible)
            resetlist = self.resetlist
            if resetlist is not ():
                try:
                    _cache._resetdict[name] = self.resetlist
                except AttributeError:
                    pass
#        else:
#            print "Reading %s from cache (%s)" % (name, _cachedval)
        return _cachedval

    def __set__(self, obj, value):
        errmsg = "The attribute '%s' cannot be overwritten" % self.name
        warnings.warn(errmsg, CacheWriteWarning)

class CachedWritableAttribute(CachedAttribute):
    #
    def __set__(self, obj, value):
        _cache = getattr(obj, self.cachename)
        name = self.name
        try:
            _cache[name] = value
        except KeyError:
            setattr(_cache, name, value)

class _cache_readonly(object):
    """
    Decorator for CachedAttribute
    """

    def __init__(self, cachename=None, resetlist=None):
        self.func = None
        self.cachename = cachename
        self.resetlist = resetlist or None

    def __call__(self, func):
        return CachedAttribute(func,
                               cachename=self.cachename,
                               resetlist=self.resetlist)
cache_readonly = _cache_readonly()

class cache_writable(_cache_readonly):
    """
    Decorator for CachedWritableAttribute
    """
    def __call__(self, func):
        return CachedWritableAttribute(func,
                                       cachename=self.cachename,
                                       resetlist=self.resetlist)


if __name__ == "__main__":
### Tests resettable_cache ----------------------------------------------------

    from numpy.testing import *

    reset = dict(a=('b',), b=('c',))
    cache = resettable_cache(a=0, b=1, c=2, reset=reset)
    assert_equal(cache, dict(a=0, b=1, c=2))
    #
    print "Try resetting a"
    cache['a'] = 1
    assert_equal(cache, dict(a=1, b=None, c=None))
    cache['c'] = 2
    assert_equal(cache, dict(a=1, b=None, c=2))
    cache['b'] = 0
    assert_equal(cache, dict(a=1, b=0, c=None))
    #
    print "Try deleting b"
    del(cache['a'])
    assert_equal(cache, {})
### ---------------------------------------------------------------------------


    class Example(object):
        #
        def __init__(self):
            self._cache = resettable_cache()
            self.a = 0
        #
        @cache_readonly
        def b(self):
            return 1
        @cache_writable(resetlist='d')
        def c(self):
            return 2
        @cache_writable(resetlist=('e', 'f'))
        def d(self):
            return self.c + 1
        #
        @cache_readonly
        def e(self):
            return 4
        @cache_readonly
        def f(self):
            return self.e + 1

    ex = Example()
    print "(attrs  : %s)" % str(ex.__dict__)
    print "(cached : %s)" % str(ex._cache)
    print "Try a   :", ex.a
    print "Try accessing/setting a readonly attribute"
    assert_equal(ex.__dict__, dict(a=0, _cache={}))
    print "Try b #1:", ex.b
    b = ex.b
    assert_equal(b, 1)
    assert_equal(ex.__dict__, dict(a=0, _cache=dict(b=1,)))
#   assert_equal(ex.__dict__, dict(a=0, b=1, _cache=dict(b=1)))
    ex.b = -1
    print "Try dict", ex.__dict__
    assert_equal(ex._cache, dict(b=1,))
    #
    print "Try accessing/resetting a cachewritable attribute"
    c = ex.c
    assert_equal(c, 2)
    assert_equal(ex._cache, dict(b=1, c=2))
    d = ex.d
    assert_equal(d, 3)
    assert_equal(ex._cache, dict(b=1, c=2, d=3))
    ex.c = 0
    assert_equal(ex._cache, dict(b=1, c=0, d=None, e=None, f=None))
    d = ex.d
    assert_equal(ex._cache, dict(b=1, c=0, d=1, e=None, f=None))
    ex.d = 5
    assert_equal(ex._cache, dict(b=1, c=0, d=5, e=None, f=None))
