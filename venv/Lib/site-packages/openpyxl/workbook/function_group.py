# Copyright (c) 2010-2024 openpyxl

from openpyxl.descriptors import (
    Integer,
    Sequence,
    String,
)
from openpyxl.descriptors.serialisable import Serialisable


class FunctionGroup(Serialisable):

    tagname = "functionGroup"

    name = String()

    def __init__(self,
                 name=None,
                ):
        self.name = name


class FunctionGroupList(Serialisable):

    tagname = "functionGroups"

    builtInGroupCount = Integer(allow_none=True)
    functionGroup = Sequence(expected_type=FunctionGroup, allow_none=True)

    __elements__ = ('functionGroup',)

    def __init__(self,
                 builtInGroupCount=16,
                 functionGroup=(),
                ):
        self.builtInGroupCount = builtInGroupCount
        self.functionGroup = functionGroup
