# Copyright (c) 2010-2024 openpyxl

from openpyxl.descriptors.excel import Relation
from openpyxl.descriptors.serialisable import Serialisable
from openpyxl.xml.constants import CHART_NS


class ChartRelation(Serialisable):

    tagname = "chart"
    namespace = CHART_NS

    id = Relation()

    def __init__(self, id):
        self.id = id
