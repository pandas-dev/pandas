from pyasn1.type.namedtype import NamedTypes
from pyasn1.type.tag import TagSet
from pyasn1.type.univ import OctetString, Sequence

class UserIdentity(OctetString):
    tagSet: TagSet
    encoding: str

class OldPasswd(OctetString):
    tagSet: TagSet
    encoding: str

class NewPasswd(OctetString):
    tagSet: TagSet
    encoding: str

class GenPasswd(OctetString):
    tagSet: TagSet
    encoding: str

class PasswdModifyRequestValue(Sequence):
    componentType: NamedTypes

class PasswdModifyResponseValue(Sequence):
    componentType: NamedTypes
