# Copyright 2026 Google LLC
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""A Least Recently Used (LRU) cache implementation."""

from collections import OrderedDict
from typing import Generic, Optional, TypeVar

K = TypeVar("K")
V = TypeVar("V")


class LRUCache(Generic[K, V]):
    """A Least Recently Used (LRU) cache implementation using OrderedDict.

    :type capacity: int
    :param capacity: The maximum number of items the cache can hold.
    """

    def __init__(self, capacity: int) -> None:
        if capacity <= 0:
            raise ValueError("Capacity must be greater than 0")
        self._capacity = capacity
        self._cache: OrderedDict[K, V] = OrderedDict()

    @property
    def capacity(self) -> int:
        """Return the capacity of the cache."""
        return self._capacity

    def get(self, key: K, default: Optional[V] = None) -> Optional[V]:
        """Retrieve an item from the cache.

        If the key exists, it is moved to the end (marked as most recently used).

        :type key: Any
        :param key: The key to look up in the cache.

        :type default: Any
        :param default: Default value to return if key is not found.
        """
        if key not in self._cache:
            return default
        self._cache.move_to_end(key)
        return self._cache[key]

    def put(self, key: K, value: V) -> None:
        """Add or update an item in the cache.

        If the key already exists, it is updated and moved to the end.
        If adding the item exceeds capacity, the least recently used item (at the beginning)
        is evicted.

        :type key: Any
        :param key: The key to store.

        :type value: Any
        :param value: The value to store.
        """
        if key in self._cache:
            self._cache.move_to_end(key)
        self._cache[key] = value
        if len(self._cache) > self._capacity:
            self._cache.popitem(last=False)

    def __len__(self) -> int:
        return len(self._cache)

    def __contains__(self, key: K) -> bool:
        return key in self._cache

    def clear(self) -> None:
        """Clear all items from the cache."""
        self._cache.clear()

    def delete(self, key: K) -> None:
        """Remove an item from the cache if it exists."""
        self._cache.pop(key, None)
