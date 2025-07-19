# Copyright (c) 2010-2024 openpyxl

from openpyxl.descriptors.excel import Relation
from openpyxl.descriptors.serialisable import Serialisable


class Drawing(Serialisable):

    tagname = "drawing"

    id = Relation()

    def __init__(self, id=None):
        self.id = id
