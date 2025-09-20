# Licensed to the Apache Software Foundation (ASF) under one
# or more contributor license agreements.  See the NOTICE file
# distributed with this work for additional information
# regarding copyright ownership.  The ASF licenses this file
# to you under the Apache License, Version 2.0 (the
# "License"); you may not use this file except in compliance
# with the License.  You may obtain a copy of the License at
#
#   http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing,
# software distributed under the License is distributed on an
# "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
# KIND, either express or implied.  See the License for the
# specific language governing permissions and limitations
# under the License.
from __future__ import annotations

from typing import (
    Callable,
    Generic,
    Iterable,
    List,
    Optional,
    TypeVar,
)

T = TypeVar("T")


class Bin(Generic[T]):
    def __init__(self, target_weight: int) -> None:
        self.bin_weight = 0
        self.target_weight = target_weight
        self.items: List[T] = []

    def weight(self) -> int:
        return self.bin_weight

    def can_add(self, weight: int) -> bool:
        return self.bin_weight + weight <= self.target_weight

    def add(self, item: T, weight: int) -> None:
        self.bin_weight += weight
        self.items.append(item)


class PackingIterator(Generic[T]):
    bins: List[Bin[T]]

    def __init__(
        self,
        items: Iterable[T],
        target_weight: int,
        lookback: int,
        weight_func: Callable[[T], int],
        largest_bin_first: bool = False,
    ) -> None:
        self.items = iter(items)
        self.target_weight = target_weight
        self.lookback = lookback
        self.weight_func = weight_func
        self.largest_bin_first = largest_bin_first
        self.bins = []

    def __iter__(self) -> PackingIterator[T]:
        """Return an iterator for the PackingIterator class."""
        return self

    def __next__(self) -> List[T]:
        """Return the next item when iterating over the PackingIterator class."""
        while True:
            try:
                item = next(self.items)
                weight = self.weight_func(item)
                bin_ = self.find_bin(weight)
                if bin_ is not None:
                    bin_.add(item, weight)
                else:
                    bin_ = Bin(self.target_weight)
                    bin_.add(item, weight)
                    self.bins.append(bin_)

                    if len(self.bins) > self.lookback:
                        return self.remove_bin().items
            except StopIteration:
                break

        if len(self.bins) == 0:
            raise StopIteration()

        return self.remove_bin().items

    def find_bin(self, weight: int) -> Optional[Bin[T]]:
        for bin_ in self.bins:
            if bin_.can_add(weight):
                return bin_
        return None

    def remove_bin(self) -> Bin[T]:
        if self.largest_bin_first:
            bin_ = max(self.bins, key=lambda b: b.weight())
            self.bins.remove(bin_)
            return bin_
        else:
            return self.bins.pop(0)


class ListPacker(Generic[T]):
    _target_weight: int
    _lookback: int
    _largest_bin_first: bool

    def __init__(self, target_weight: int, lookback: int, largest_bin_first: bool) -> None:
        self._target_weight = target_weight
        self._lookback = lookback
        self._largest_bin_first = largest_bin_first

    def pack(self, items: List[T], weight_func: Callable[[T], int]) -> List[List[T]]:
        return list(
            PackingIterator(
                items=items,
                target_weight=self._target_weight,
                lookback=self._lookback,
                weight_func=weight_func,
                largest_bin_first=self._largest_bin_first,
            )
        )

    def pack_end(self, items: List[T], weight_func: Callable[[T], int]) -> List[List[T]]:
        packed = self.pack(items=list(reversed(items)), weight_func=weight_func)
        return [list(reversed(bin_items)) for bin_items in reversed(packed)]
