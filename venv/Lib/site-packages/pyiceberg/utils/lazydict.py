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

from collections.abc import Iterator, Mapping, Sequence
from typing import (
    TypeVar,
    cast,
)

K = TypeVar("K")
V = TypeVar("V")


class LazyDict(Mapping[K, V]):
    """Lazily build a dictionary from an array of items."""

    __slots__ = ("_contents", "_dict")

    # Since Python's type system is not powerful enough to express the type of the
    # contents of the dictionary, we use specify the type as a sequence of either K or V
    # values.
    #
    # Rather than spending the runtime cost of checking the type of each item, we presume
    # that the developer has correctly used the class and that the contents are valid.
    def __init__(self, contents: Sequence[Sequence[K | V]]):
        self._contents = contents
        self._dict: dict[K, V] | None = None

    def _build_dict(self) -> dict[K, V]:
        self._dict = {}
        for item in self._contents:
            self._dict.update(dict(zip(cast(Sequence[K], item[::2]), cast(Sequence[V], item[1::2]), strict=True)))

        return self._dict

    def __getitem__(self, key: K, /) -> V:
        """Return the value for the given key."""
        source = self._dict or self._build_dict()
        return source[key]

    def __iter__(self) -> Iterator[K]:
        """Return an iterator over the keys of the dictionary."""
        source = self._dict or self._build_dict()
        return iter(source)

    def __len__(self) -> int:
        """Return the number of items in the dictionary."""
        source = self._dict or self._build_dict()
        return len(source)

    def __dict__(self) -> dict[K, V]:  # type: ignore
        """Convert the lazy dict in a dict."""
        return self._dict or self._build_dict()
