from typing import Final

from pyasn1.type.constraint import ConstraintsIntersection, ValueRangeConstraint
from pyasn1.type.namedtype import NamedTypes
from pyasn1.type.univ import Integer, OctetString, Sequence

MAXINT: Final[Integer]
rangeInt0ToMaxConstraint: ValueRangeConstraint

class Integer0ToMax(Integer):
    subtypeSpec: ConstraintsIntersection

class Size(Integer0ToMax): ...
class Cookie(OctetString): ...

class RealSearchControlValue(Sequence):
    componentType: NamedTypes

def paged_search_control(criticality: bool = False, size: int = 10, cookie=None): ...
