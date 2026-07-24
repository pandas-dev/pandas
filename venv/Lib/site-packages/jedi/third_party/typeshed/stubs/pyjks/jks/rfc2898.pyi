from pyasn1.type.namedtype import NamedTypes
from pyasn1.type.univ import Sequence

class PBEParameter(Sequence):
    componentType: NamedTypes
