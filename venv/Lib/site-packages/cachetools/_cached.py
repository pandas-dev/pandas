"""Function decorator helpers."""

import functools


def _condition_info(func, cache, key, lock, cond, info):
    hits = misses = 0
    pending = set()

    def wrapper(*args, **kwargs):
        nonlocal hits, misses
        k = key(*args, **kwargs)
        with lock:
            cond.wait_for(lambda: k not in pending)
            try:
                result = cache[k]
                hits += 1
                return result
            except KeyError:
                pending.add(k)
                misses += 1
        try:
            v = func(*args, **kwargs)
            with lock:
                try:
                    cache[k] = v
                except ValueError:
                    pass  # value too large
                return v
        finally:
            with lock:
                pending.remove(k)
                cond.notify_all()

    def cache_clear():
        nonlocal hits, misses
        with lock:
            cache.clear()
            hits = misses = 0

    def cache_info():
        with lock:
            return info(hits, misses)

    wrapper.cache_clear = cache_clear
    wrapper.cache_info = cache_info
    return wrapper


def _locked_info(func, cache, key, lock, info):
    hits = misses = 0

    def wrapper(*args, **kwargs):
        nonlocal hits, misses
        k = key(*args, **kwargs)
        with lock:
            try:
                result = cache[k]
                hits += 1
                return result
            except KeyError:
                misses += 1
        v = func(*args, **kwargs)
        with lock:
            try:
                # In case of a race condition, i.e. if another thread
                # stored a value for this key while we were calling
                # func(), prefer the cached value.
                return cache.setdefault(k, v)
            except ValueError:
                return v  # value too large

    def cache_clear():
        nonlocal hits, misses
        with lock:
            cache.clear()
            hits = misses = 0

    def cache_info():
        with lock:
            return info(hits, misses)

    wrapper.cache_clear = cache_clear
    wrapper.cache_info = cache_info
    return wrapper


def _unlocked_info(func, cache, key, info):
    hits = misses = 0

    def wrapper(*args, **kwargs):
        nonlocal hits, misses
        k = key(*args, **kwargs)
        try:
            result = cache[k]
            hits += 1
            return result
        except KeyError:
            misses += 1
        v = func(*args, **kwargs)
        try:
            cache[k] = v
        except ValueError:
            pass  # value too large
        return v

    def cache_clear():
        nonlocal hits, misses
        cache.clear()
        hits = misses = 0

    def cache_info():
        return info(hits, misses)

    wrapper.cache_clear = cache_clear
    wrapper.cache_info = cache_info
    return wrapper


def _uncached_info(func, info):
    misses = 0

    def wrapper(*args, **kwargs):
        nonlocal misses
        misses += 1
        return func(*args, **kwargs)

    def cache_clear():
        nonlocal misses
        misses = 0

    wrapper.cache_clear = cache_clear
    wrapper.cache_info = lambda: info(0, misses)
    return wrapper


def _condition(func, cache, key, lock, cond):
    pending = set()

    def wrapper(*args, **kwargs):
        k = key(*args, **kwargs)
        with lock:
            cond.wait_for(lambda: k not in pending)
            try:
                result = cache[k]
                return result
            except KeyError:
                pending.add(k)
        try:
            v = func(*args, **kwargs)
            with lock:
                try:
                    cache[k] = v
                except ValueError:
                    pass  # value too large
                return v
        finally:
            with lock:
                pending.remove(k)
                cond.notify_all()

    def cache_clear():
        with lock:
            cache.clear()

    wrapper.cache_clear = cache_clear
    return wrapper


def _locked(func, cache, key, lock):
    def wrapper(*args, **kwargs):
        k = key(*args, **kwargs)
        with lock:
            try:
                return cache[k]
            except KeyError:
                pass  # key not found
        v = func(*args, **kwargs)
        with lock:
            try:
                # possible race condition: see above
                return cache.setdefault(k, v)
            except ValueError:
                return v  # value too large

    def cache_clear():
        with lock:
            cache.clear()

    wrapper.cache_clear = cache_clear
    return wrapper


def _unlocked(func, cache, key):
    def wrapper(*args, **kwargs):
        k = key(*args, **kwargs)
        try:
            return cache[k]
        except KeyError:
            pass  # key not found
        v = func(*args, **kwargs)
        try:
            cache[k] = v
        except ValueError:
            pass  # value too large
        return v

    wrapper.cache_clear = lambda: cache.clear()
    return wrapper


def _uncached(func):
    def wrapper(*args, **kwargs):
        return func(*args, **kwargs)

    wrapper.cache_clear = lambda: None
    return wrapper


def _wrapper(func, cache, key, lock=None, cond=None, info=None):
    if info is not None:
        if cache is None:
            wrapper = _uncached_info(func, info)
        elif cond is not None and lock is not None:
            wrapper = _condition_info(func, cache, key, lock, cond, info)
        elif cond is not None:
            wrapper = _condition_info(func, cache, key, cond, cond, info)
        elif lock is not None:
            wrapper = _locked_info(func, cache, key, lock, info)
        else:
            wrapper = _unlocked_info(func, cache, key, info)
    else:
        if cache is None:
            wrapper = _uncached(func)
        elif cond is not None and lock is not None:
            wrapper = _condition(func, cache, key, lock, cond)
        elif cond is not None:
            wrapper = _condition(func, cache, key, cond, cond)
        elif lock is not None:
            wrapper = _locked(func, cache, key, lock)
        else:
            wrapper = _unlocked(func, cache, key)
        wrapper.cache_info = None

    wrapper.cache = cache
    wrapper.cache_key = key
    wrapper.cache_lock = lock if lock is not None else cond
    wrapper.cache_condition = cond

    return functools.update_wrapper(wrapper, func)
