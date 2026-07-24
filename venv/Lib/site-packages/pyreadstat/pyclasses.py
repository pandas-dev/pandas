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

# Typing

from dataclasses import dataclass, field
from datetime import datetime
from typing import Literal, TypedDict


class MissingRange(TypedDict):
    """A dictionary to hold the definition of a missing range"""

    lo: float
    hi: float


class MRSet(TypedDict):
    """A dictionary to hold the definition of a multiple-response (MR) set"""

    type: Literal["D", "C"]
    is_dichotomy: bool
    counted_value: int | None
    label: str
    variable_list: list[str]


# Classes


@dataclass
class metadata_container:
    """
    This class holds metadata we want to give back to python
    """

    column_names: list[str] = field(default_factory=list)
    column_labels: list[str] = field(default_factory=list)
    column_names_to_labels: dict[str, str] = field(default_factory=dict)
    file_encoding: str | None = None
    file_label: str | None = None
    number_columns: int | None = None
    number_rows: int | None = None
    variable_value_labels: dict[str, dict[float | int, str]] = field(default_factory=dict)
    value_labels: dict[str, dict[float | int, str]] = field(default_factory=dict)
    variable_to_label: dict[str, str] = field(default_factory=dict)
    notes: list[str] = field(default_factory=list)
    original_variable_types: dict[str, str] = field(default_factory=dict)
    readstat_variable_types: dict[str, str] = field(default_factory=dict)
    table_name: str | None = None
    missing_ranges: dict[str, list[int | float | str | MissingRange]] = field(default_factory=dict)
    missing_user_values: dict[str, list[int | float | str | MissingRange]] = field(default_factory=dict)
    variable_storage_width: dict[str, int] = field(default_factory=dict)
    variable_display_width: dict[str, int] = field(default_factory=dict)
    variable_alignment: dict[str, str] = field(default_factory=dict)
    variable_measure: dict[str, Literal["nominal", "ordinal", "scale", "unknown"]] = field(default_factory=dict)
    creation_time: datetime | None = None
    modification_time: datetime | None = None
    mr_sets: dict[str, MRSet] = field(default_factory=dict)
