# This file is part of Hypothesis, which may be found at
# https://github.com/HypothesisWorks/hypothesis/
#
# Copyright the Hypothesis Authors.
# Individual contributors are listed in AUTHORS.rst and the git log.
#
# This Source Code Form is subject to the terms of the Mozilla Public License,
# v. 2.0. If a copy of the MPL was not distributed with this file, You can
# obtain one at https://mozilla.org/MPL/2.0/.

import threading

import attr


@attr.s(slots=True)
class Entry:
    key = attr.ib()
    value = attr.ib()
    score = attr.ib()
    pins = attr.ib(default=0)

    @property
    def sort_key(self):
        if self.pins == 0:
            # Unpinned entries are sorted by score.
            return (0, self.score)
        else:
            # Pinned entries sort after unpinned ones. Beyond that, we don't
            # worry about their relative order.
            return (1,)


class GenericCache:
    """Generic supertype for cache implementations.

    Defines a dict-like mapping with a maximum size, where as well as mapping
    to a value, each key also maps to a score. When a write would cause the
    dict to exceed its maximum size, it first evicts the existing key with
    the smallest score, then adds the new key to the map.

    A key has the following lifecycle:

    1. key is written for the first time, the key is given the score
       self.new_entry(key, value)
    2. whenever an existing key is read or written, self.on_access(key, value,
       score) is called. This returns a new score for the key.
    3. When a key is evicted, self.on_evict(key, value, score) is called.

    The cache will be in a valid state in all of these cases.

    Implementations are expected to implement new_entry and optionally
    on_access and on_evict to implement a specific scoring strategy.
    """

    __slots__ = ("max_size", "_threadlocal")

    def __init__(self, max_size):
        self.max_size = max_size

        # Implementation: We store a binary heap of Entry objects in self.data,
        # with the heap property requiring that a parent's score is <= that of
        # its children. keys_to_index then maps keys to their index in the
        # heap. We keep these two in sync automatically - the heap is never
        # reordered without updating the index.
        self._threadlocal = threading.local()

    @property
    def keys_to_indices(self):
        try:
            return self._threadlocal.keys_to_indices
        except AttributeError:
            self._threadlocal.keys_to_indices = {}
            return self._threadlocal.keys_to_indices

    @property
    def data(self):
        try:
            return self._threadlocal.data
        except AttributeError:
            self._threadlocal.data = []
            return self._threadlocal.data

    @property
    def __pinned_entry_count(self):
        return getattr(self._threadlocal, "_pinned_entry_count", 0)

    @__pinned_entry_count.setter
    def __pinned_entry_count(self, value):
        self._threadlocal._pinned_entry_count = value

    def __len__(self):
        assert len(self.keys_to_indices) == len(self.data)
        return len(self.data)

    def __contains__(self, key):
        return key in self.keys_to_indices

    def __getitem__(self, key):
        i = self.keys_to_indices[key]
        result = self.data[i]
        self.on_access(result.key, result.value, result.score)
        self.__balance(i)
        return result.value

    def __setitem__(self, key, value):
        if self.max_size == 0:
            return
        evicted = None
        try:
            i = self.keys_to_indices[key]
        except KeyError:
            if self.max_size == self.__pinned_entry_count:
                raise ValueError(
                    "Cannot increase size of cache where all keys have been pinned."
                ) from None
            entry = Entry(key, value, self.new_entry(key, value))
            if len(self.data) >= self.max_size:
                evicted = self.data[0]
                assert evicted.pins == 0
                del self.keys_to_indices[evicted.key]
                i = 0
                self.data[0] = entry
            else:
                i = len(self.data)
                self.data.append(entry)
            self.keys_to_indices[key] = i
        else:
            entry = self.data[i]
            assert entry.key == key
            entry.value = value
            entry.score = self.on_access(entry.key, entry.value, entry.score)

        self.__balance(i)

        if evicted is not None:
            if self.data[0] is not entry:
                assert evicted.score <= self.data[0].score
            self.on_evict(evicted.key, evicted.value, evicted.score)

    def __iter__(self):
        return iter(self.keys_to_indices)

    def pin(self, key):
        """Mark ``key`` as pinned. That is, it may not be evicted until
        ``unpin(key)`` has been called. The same key may be pinned multiple
        times and will not be unpinned until the same number of calls to
        unpin have been made."""
        i = self.keys_to_indices[key]
        entry = self.data[i]
        entry.pins += 1
        if entry.pins == 1:
            self.__pinned_entry_count += 1
            assert self.__pinned_entry_count <= self.max_size
            self.__balance(i)

    def unpin(self, key):
        """Undo one previous call to ``pin(key)``. Once all calls are
        undone this key may be evicted as normal."""
        i = self.keys_to_indices[key]
        entry = self.data[i]
        if entry.pins == 0:
            raise ValueError(f"Key {key!r} has not been pinned")
        entry.pins -= 1
        if entry.pins == 0:
            self.__pinned_entry_count -= 1
            self.__balance(i)

    def is_pinned(self, key):
        """Returns True if the key is currently pinned."""
        i = self.keys_to_indices[key]
        return self.data[i].pins > 0

    def clear(self):
        """Remove all keys, clearing their pinned status."""
        del self.data[:]
        self.keys_to_indices.clear()
        self.__pinned_entry_count = 0

    def __repr__(self):
        return "{" + ", ".join(f"{e.key!r}: {e.value!r}" for e in self.data) + "}"

    def new_entry(self, key, value):
        """Called when a key is written that does not currently appear in the
        map.

        Returns the score to associate with the key.
        """
        raise NotImplementedError

    def on_access(self, key, value, score):
        """Called every time a key that is already in the map is read or
        written.

        Returns the new score for the key.
        """
        return score

    def on_evict(self, key, value, score):
        """Called after a key has been evicted, with the score it had had at
        the point of eviction."""

    def check_valid(self):
        """Debugging method for use in tests.

        Asserts that all of the cache's invariants hold. When everything
        is working correctly this should be an expensive no-op.
        """
        for i, e in enumerate(self.data):
            assert self.keys_to_indices[e.key] == i
            for j in [i * 2 + 1, i * 2 + 2]:
                if j < len(self.data):
                    assert e.score <= self.data[j].score, self.data

    def __swap(self, i, j):
        assert i < j
        assert self.data[j].sort_key < self.data[i].sort_key
        self.data[i], self.data[j] = self.data[j], self.data[i]
        self.keys_to_indices[self.data[i].key] = i
        self.keys_to_indices[self.data[j].key] = j

    def __balance(self, i):
        """When we have made a modification to the heap such that means that
        the heap property has been violated locally around i but previously
        held for all other indexes (and no other values have been modified),
        this fixes the heap so that the heap property holds everywhere."""
        while i > 0:
            parent = (i - 1) // 2
            if self.__out_of_order(parent, i):
                self.__swap(parent, i)
                i = parent
            else:
                # This branch is never taken on versions of Python where dicts
                # preserve their insertion order (pypy or cpython >= 3.7)
                break  # pragma: no cover
        while True:
            children = [j for j in (2 * i + 1, 2 * i + 2) if j < len(self.data)]
            if len(children) == 2:
                children.sort(key=lambda j: self.data[j].score)
            for j in children:
                if self.__out_of_order(i, j):
                    self.__swap(i, j)
                    i = j
                    break
            else:
                break

    def __out_of_order(self, i, j):
        """Returns True if the indices i, j are in the wrong order.

        i must be the parent of j.
        """
        assert i == (j - 1) // 2
        return self.data[j].sort_key < self.data[i].sort_key


class LRUReusedCache(GenericCache):
    """The only concrete implementation of GenericCache we use outside of tests
    currently.

    Adopts a modified least-frequently used eviction policy: It evicts the key
    that has been used least recently, but it will always preferentially evict
    keys that have only ever been accessed once. Among keys that have been
    accessed more than once, it ignores the number of accesses.

    This retains most of the benefits of an LRU cache, but adds an element of
    scan-resistance to the process: If we end up scanning through a large
    number of keys without reusing them, this does not evict the existing
    entries in preference for the new ones.
    """

    __slots__ = ("__tick",)

    def __init__(self, max_size):
        super().__init__(max_size)
        self.__tick = 0

    def tick(self):
        self.__tick += 1
        return self.__tick

    def new_entry(self, key, value):
        return [1, self.tick()]

    def on_access(self, key, value, score):
        score[0] = 2
        score[1] = self.tick()
        return score

    def pin(self, key):
        try:
            super().pin(key)
        except KeyError:
            # The whole point of an LRU cache is that it might drop things for you
            assert key not in self.keys_to_indices

    def unpin(self, key):
        try:
            super().unpin(key)
        except KeyError:
            assert key not in self.keys_to_indices
