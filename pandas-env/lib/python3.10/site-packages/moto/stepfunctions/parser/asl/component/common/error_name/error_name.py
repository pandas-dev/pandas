from __future__ import annotations

import abc
from typing import Optional

from moto.stepfunctions.parser.asl.component.component import Component


class ErrorName(Component, abc.ABC):
    error_name: Optional[str]

    def __init__(self, error_name: Optional[str]):
        self.error_name = error_name

    def matches(self, error_name: Optional[str]) -> bool:
        return self.error_name == error_name

    def __eq__(self, other):
        if isinstance(other, ErrorName):
            return self.matches(other.error_name)
        return False
