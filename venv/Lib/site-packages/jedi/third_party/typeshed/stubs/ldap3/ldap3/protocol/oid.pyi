from _typeshed import Incomplete

OID_CONTROL: str
OID_EXTENSION: str
OID_FEATURE: str
OID_UNSOLICITED_NOTICE: str
OID_ATTRIBUTE_TYPE: str
OID_DIT_CONTENT_RULE: str
OID_LDAP_URL_EXTENSION: str
OID_FAMILY: str
OID_MATCHING_RULE: str
OID_NAME_FORM: str
OID_OBJECT_CLASS: str
OID_ADMINISTRATIVE_ROLE: str
OID_LDAP_SYNTAX: str
CLASS_STRUCTURAL: str
CLASS_ABSTRACT: str
CLASS_AUXILIARY: str
ATTRIBUTE_USER_APPLICATION: str
ATTRIBUTE_DIRECTORY_OPERATION: str
ATTRIBUTE_DISTRIBUTED_OPERATION: str
ATTRIBUTE_DSA_OPERATION: str

def constant_to_oid_kind(oid_kind): ...
def decode_oids(sequence): ...
def decode_syntax(syntax): ...
def oid_to_string(oid): ...

Oids: Incomplete
