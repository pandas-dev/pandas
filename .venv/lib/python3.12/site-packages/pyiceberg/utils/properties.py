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

from typing import (
    Any,
)

from pyiceberg.typedef import Properties
from pyiceberg.types import strtobool

HEADER_PREFIX = "header."


def property_as_int(
    properties: dict[str, str],
    property_name: str,
    default: int | None = None,
) -> int | None:
    if value := properties.get(property_name):
        try:
            return int(value)
        except ValueError as e:
            raise ValueError(f"Could not parse table property {property_name} to an integer: {value}") from e
    else:
        return default


def property_as_float(
    properties: dict[str, str],
    property_name: str,
    default: float | None = None,
) -> float | None:
    if value := properties.get(property_name):
        try:
            return float(value)
        except ValueError as e:
            raise ValueError(f"Could not parse table property {property_name} to a float: {value}") from e
    else:
        return default


def property_as_bool(
    properties: dict[str, str],
    property_name: str,
    default: bool,
) -> bool:
    if value := properties.get(property_name):
        try:
            return strtobool(value)
        except ValueError as e:
            raise ValueError(f"Could not parse table property {property_name} to a boolean: {value}") from e
    return default


def get_first_property_value(
    properties: Properties,
    *property_names: str,
) -> Any | None:
    for property_name in property_names:
        if property_value := properties.get(property_name):
            return property_value
    return None


def get_header_properties(
    properties: Properties,
) -> Properties:
    header_prefix_len = len(HEADER_PREFIX)
    return {key[header_prefix_len:]: value for key, value in properties.items() if key.startswith(HEADER_PREFIX)}
