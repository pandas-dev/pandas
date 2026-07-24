from pyasn1.type.namedtype import NamedTypes
from pyasn1.type.tag import TagSet
from pyasn1.type.univ import Sequence

class SicilyBindResponse(Sequence):
    tagSet: TagSet
    componentType: NamedTypes

class DirSyncControlRequestValue(Sequence):
    componentType: NamedTypes

class DirSyncControlResponseValue(Sequence):
    componentType: NamedTypes

class SdFlags(Sequence):
    componentType: NamedTypes

class ExtendedDN(Sequence):
    componentType: NamedTypes

def dir_sync_control(criticality, object_security, ancestors_first, public_data_only, incremental_values, max_length, cookie): ...
def extended_dn_control(criticality: bool = False, hex_format: bool = False): ...
def show_deleted_control(criticality: bool = False): ...
def security_descriptor_control(criticality: bool = False, sdflags: int = 15): ...
def persistent_search_control(criticality: bool = False): ...
