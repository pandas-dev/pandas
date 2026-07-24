from typing import ClassVar

from openpyxl.descriptors.serialisable import Serialisable

class SheetBackgroundPicture(Serialisable):
    tagname: ClassVar[str]
