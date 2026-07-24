# #############################################################################
# Copyright 2018 Hoffmann-La Roche
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#   http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
# #############################################################################

"""
Functions to work with multiprocessing
"""

from os import PathLike
from typing import TYPE_CHECKING, Any, TypeAlias

if TYPE_CHECKING:
    from .pyreadstat import PyreadstatReadFunction, DataFrame, DictOutput

Input: TypeAlias = "tuple[PyreadstatReadFunction, str | bytes | PathLike, int, int, dict[str, Any]]"


def worker(inpt: Input) -> "DataFrame | DictOutput":
    read_function, path, row_offset, row_limit, kwargs = inpt
    df, meta = read_function(path, row_offset=row_offset, row_limit=row_limit, **kwargs)
    return df
