"""Method decorator helpers."""

import functools
import weakref


def warn_cache_none():
    from warnings import warn

    warn(
        "returning `None` from `cache(self)` is deprecated",
        DeprecationWarning,
        stacklevel=3,
    )


def _condition(method, cache, key, lock, cond):
    pending = weakref.WeakKeyDictionary()

    def wrapper(self, *args, **kwargs):
        c = cache(self)
        if c is None:
            warn_cache_none()
            return method(self, *args, **kwargs)
        k = key(self, *args, **kwargs)
        with lock(self):
            p = pending.setdefault(self, set())
            cond(self).wait_for(lambda: k not in p)
            try:
                return c[k]
            except KeyError:
                p.add(k)
        try:
            v = method(self, *args, **kwargs)
            with lock(self):
                try:
                    c[k] = v
                except ValueError:
                    pass  # value too large
                return v
        finally:
            with lock(self):
                pending[self].remove(k)
                cond(self).notify_all()

    def cache_clear(self):
        c = cache(self)
        if c is not None:
            with lock(self):
                c.clear()

    wrapper.cache_clear = cache_clear
    return wrapper


def _locked(method, cache, key, lock):
    def wrapper(self, *args, **kwargs):
        c = cache(self)
        if c is None:
            warn_cache_none()
            return method(self, *args, **kwargs)
        k = key(self, *args, **kwargs)
        with lock(self):
            try:
                return c[k]
            except KeyError:
                pass  # key not found
        v = method(self, *args, **kwargs)
        # in case of a race, prefer the item already in the cache
        with lock(self):
            try:
                return c.setdefault(k, v)
            except ValueError:
                return v  # value too large

    def cache_clear(self):
        c = cache(self)
        if c is not None:
            with lock(self):
                c.clear()

    wrapper.cache_clear = cache_clear
    return wrapper


def _unlocked(method, cache, key):
    def wrapper(self, *args, **kwargs):
        c = cache(self)
        if c is None:
            warn_cache_none()
            return method(self, *args, **kwargs)
        k = key(self, *args, **kwargs)
        try:
            return c[k]
        except KeyError:
            pass  # key not found
        v = method(self, *args, **kwargs)
        try:
            c[k] = v
        except ValueError:
            pass  # value too large
        return v

    def cache_clear(self):
        c = cache(self)
        if c is not None:
            c.clear()

    wrapper.cache_clear = cache_clear
    return wrapper


def _wrapper(method, cache, key, lock=None, cond=None):
    if cond is not None and lock is not None:
        wrapper = _condition(method, cache, key, lock, cond)
    elif cond is not None:
        wrapper = _condition(method, cache, key, cond, cond)
    elif lock is not None:
        wrapper = _locked(method, cache, key, lock)
    else:
        wrapper = _unlocked(method, cache, key)

    wrapper.cache = cache
    wrapper.cache_key = key
    wrapper.cache_lock = lock if lock is not None else cond
    wrapper.cache_condition = cond

    return functools.update_wrapper(wrapper, method)
