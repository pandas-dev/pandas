from pyasn1.type.namedtype import NamedTypes
from pyasn1.type.tag import TagSet
from pyasn1.type.univ import Integer, OctetString, Sequence, SequenceOf

NMAS_LDAP_EXT_VERSION: int

class Identity(OctetString):
    encoding: str

class LDAPDN(OctetString):
    tagSet: TagSet
    encoding: str

class Password(OctetString):
    tagSet: TagSet
    encoding: str

class LDAPOID(OctetString):
    tagSet: TagSet
    encoding: str

class GroupCookie(Integer):
    tagSet: TagSet

class NmasVer(Integer):
    tagSet: TagSet

class Error(Integer):
    tagSet: TagSet

class NmasGetUniversalPasswordRequestValue(Sequence):
    componentType: NamedTypes

class NmasGetUniversalPasswordResponseValue(Sequence):
    componentType: NamedTypes

class NmasSetUniversalPasswordRequestValue(Sequence):
    componentType: NamedTypes

class NmasSetUniversalPasswordResponseValue(Sequence):
    componentType: NamedTypes

class ReplicaList(SequenceOf):
    componentType: OctetString

class ReplicaInfoRequestValue(Sequence):
    tagSet: TagSet
    componentType: NamedTypes

class ReplicaInfoResponseValue(Sequence):
    tagSet: TagSet
    componentType: NamedTypes

class CreateGroupTypeRequestValue(Sequence):
    componentType: NamedTypes

class CreateGroupTypeResponseValue(Sequence):
    componentType: NamedTypes

class EndGroupTypeRequestValue(Sequence):
    componentType: NamedTypes

class EndGroupTypeResponseValue(Sequence):
    componentType: NamedTypes

class GroupingControlValue(Sequence):
    componentType: NamedTypes
