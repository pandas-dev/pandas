from pyasn1.type import namedtype, univ

RSA_ENCRYPTION_ASN1_OID: str

class RsaAlgorithmIdentifier(univ.Sequence):
    componentType: namedtype.NamedTypes

class PKCS8PrivateKey(univ.Sequence):
    componentType: namedtype.NamedTypes

class PublicKeyInfo(univ.Sequence):
    componentType: namedtype.NamedTypes

def rsa_private_key_pkcs8_to_pkcs1(pkcs8_key) -> bytes: ...
def rsa_private_key_pkcs1_to_pkcs8(pkcs1_key) -> bytes: ...
def rsa_public_key_pkcs1_to_pkcs8(pkcs1_key) -> bytes: ...
def rsa_public_key_pkcs8_to_pkcs1(pkcs8_key) -> bytes: ...
