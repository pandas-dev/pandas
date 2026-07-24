#  Licensed to the Apache Software Foundation (ASF) under one
#  or more contributor license agreements.  See the NOTICE file
#  distributed with this work for additional information
#  regarding copyright ownership.  The ASF licenses this file
#  to you under the Apache License, Version 2.0 (the
#  "License"); you may not use this file except in compliance
#  with the License.  You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
#  Unless required by applicable law or agreed to in writing,
#  software distributed under the License is distributed on an
#  "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
#  KIND, either express or implied.  See the License for the
#  specific language governing permissions and limitations
#  under the License.
import re
from re import Pattern

from pyiceberg.exceptions import ValidationError


class ParseNumberFromBrackets:
    """Extracts the size from a string in the form of prefix[22]."""

    regex: Pattern  # type: ignore
    prefix: str

    def __init__(self, prefix: str):
        self.prefix = prefix
        self.regex = re.compile(rf"{prefix}\[(\d+)\]")

    def match(self, str_repr: str) -> int:
        matches = self.regex.search(str_repr)
        if matches:
            return int(matches.group(1))
        raise ValidationError(f"Could not match {str_repr}, expected format {self.prefix}[22]")
