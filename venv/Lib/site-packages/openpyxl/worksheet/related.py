# Copyright (c) 2010-2024 openpyxl

from openpyxl.descriptors.excel import Relation
from openpyxl.descriptors.serialisable import Serialisable


class Related(Serialisable):

    id = Relation()


    def __init__(self, id=None):
        self.id = id


    def to_tree(self, tagname, idx=None):
        return super().to_tree(tagname)
