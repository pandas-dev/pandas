from pyasn1.type.namedtype import NamedTypes
from pyasn1.type.namedval import NamedValues
from pyasn1.type.univ import Enumerated, Sequence

class PersistentSearchControl(Sequence):
    componentType: NamedTypes

class ChangeType(Enumerated):
    namedValues: NamedValues

class EntryChangeNotificationControl(Sequence):
    componentType: NamedTypes

def persistent_search_control(change_types, changes_only: bool = True, return_ecs: bool = True, criticality: bool = False): ...
