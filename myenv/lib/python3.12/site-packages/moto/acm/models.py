import base64
import datetime
import re
from typing import Any, Dict, Iterable, List, Optional, Set, Tuple

import cryptography.hazmat.primitives.asymmetric.rsa
import cryptography.x509
from cryptography.hazmat.backends import default_backend
from cryptography.hazmat.primitives import hashes, serialization
from cryptography.x509 import OID_COMMON_NAME, DNSName, NameOID

from moto import settings
from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.common_models import BaseModel
from moto.core.utils import utcnow

from .exceptions import (
    AWSTooManyTagsException,
    AWSValidationException,
    CertificateNotFound,
)
from .utils import make_arn_for_certificate

AWS_ROOT_CA = b"""-----BEGIN CERTIFICATE-----
MIIESTCCAzGgAwIBAgITBntQXCplJ7wevi2i0ZmY7bibLDANBgkqhkiG9w0BAQsF
ADA5MQswCQYDVQQGEwJVUzEPMA0GA1UEChMGQW1hem9uMRkwFwYDVQQDExBBbWF6
b24gUm9vdCBDQSAxMB4XDTE1MTAyMTIyMjQzNFoXDTQwMTAyMTIyMjQzNFowRjEL
MAkGA1UEBhMCVVMxDzANBgNVBAoTBkFtYXpvbjEVMBMGA1UECxMMU2VydmVyIENB
IDFCMQ8wDQYDVQQDEwZBbWF6b24wggEiMA0GCSqGSIb3DQEBAQUAA4IBDwAwggEK
AoIBAQDCThZn3c68asg3Wuw6MLAd5tES6BIoSMzoKcG5blPVo+sDORrMd4f2AbnZ
cMzPa43j4wNxhplty6aUKk4T1qe9BOwKFjwK6zmxxLVYo7bHViXsPlJ6qOMpFge5
blDP+18x+B26A0piiQOuPkfyDyeR4xQghfj66Yo19V+emU3nazfvpFA+ROz6WoVm
B5x+F2pV8xeKNR7u6azDdU5YVX1TawprmxRC1+WsAYmz6qP+z8ArDITC2FMVy2fw
0IjKOtEXc/VfmtTFch5+AfGYMGMqqvJ6LcXiAhqG5TI+Dr0RtM88k+8XUBCeQ8IG
KuANaL7TiItKZYxK1MMuTJtV9IblAgMBAAGjggE7MIIBNzASBgNVHRMBAf8ECDAG
AQH/AgEAMA4GA1UdDwEB/wQEAwIBhjAdBgNVHQ4EFgQUWaRmBlKge5WSPKOUByeW
dFv5PdAwHwYDVR0jBBgwFoAUhBjMhTTsvAyUlC4IWZzHshBOCggwewYIKwYBBQUH
AQEEbzBtMC8GCCsGAQUFBzABhiNodHRwOi8vb2NzcC5yb290Y2ExLmFtYXpvbnRy
dXN0LmNvbTA6BggrBgEFBQcwAoYuaHR0cDovL2NybC5yb290Y2ExLmFtYXpvbnRy
dXN0LmNvbS9yb290Y2ExLmNlcjA/BgNVHR8EODA2MDSgMqAwhi5odHRwOi8vY3Js
LnJvb3RjYTEuYW1hem9udHJ1c3QuY29tL3Jvb3RjYTEuY3JsMBMGA1UdIAQMMAow
CAYGZ4EMAQIBMA0GCSqGSIb3DQEBCwUAA4IBAQAfsaEKwn17DjAbi/Die0etn+PE
gfY/I6s8NLWkxGAOUfW2o+vVowNARRVjaIGdrhAfeWHkZI6q2pI0x/IJYmymmcWa
ZaW/2R7DvQDtxCkFkVaxUeHvENm6IyqVhf6Q5oN12kDSrJozzx7I7tHjhBK7V5Xo
TyS4NU4EhSyzGgj2x6axDd1hHRjblEpJ80LoiXlmUDzputBXyO5mkcrplcVvlIJi
WmKjrDn2zzKxDX5nwvkskpIjYlJcrQu4iCX1/YwZ1yNqF9LryjlilphHCACiHbhI
RnGfN8j8KLDVmWyTYMk8V+6j0LI4+4zFh2upqGMQHL3VFVFWBek6vCDWhB/b
 -----END CERTIFICATE-----"""
# Added aws root CA as AWS returns chain you gave it + root CA (provided or not)
# so for now a cheap response is just give any old root CA


def datetime_to_epoch(date: datetime.datetime) -> float:
    return date.timestamp()


class TagHolder(Dict[str, Optional[str]]):
    MAX_TAG_COUNT = 50
    MAX_KEY_LENGTH = 128
    MAX_VALUE_LENGTH = 256

    def _validate_kv(self, key: str, value: Optional[str], index: int) -> None:
        if len(key) > self.MAX_KEY_LENGTH:
            raise AWSValidationException(
                "Value '%s' at 'tags.%d.member.key' failed to satisfy constraint: Member must have length less than or equal to %s"
                % (key, index, self.MAX_KEY_LENGTH)
            )
        if value and len(value) > self.MAX_VALUE_LENGTH:
            raise AWSValidationException(
                "Value '%s' at 'tags.%d.member.value' failed to satisfy constraint: Member must have length less than or equal to %s"
                % (value, index, self.MAX_VALUE_LENGTH)
            )
        if key.startswith("aws:"):
            raise AWSValidationException(
                'Invalid Tag Key: "%s". AWS internal tags cannot be changed with this API'
                % key
            )

    def add(self, tags: List[Dict[str, str]]) -> None:
        tags_copy = self.copy()
        for i, tag in enumerate(tags):
            key = tag["Key"]
            value = tag.get("Value")
            self._validate_kv(key, value, i + 1)

            tags_copy[key] = value
        if len(tags_copy) > self.MAX_TAG_COUNT:
            raise AWSTooManyTagsException(
                "the TagSet: '{%s}' contains too many Tags"
                % ", ".join(k + "=" + str(v or "") for k, v in tags_copy.items())
            )

        self.update(tags_copy)

    def remove(self, tags: List[Dict[str, str]]) -> None:
        for i, tag in enumerate(tags):
            key = tag["Key"]
            value = tag.get("Value")
            self._validate_kv(key, value, i + 1)
            try:
                # If value isnt provided, just delete key
                if value is None:
                    del self[key]
                # If value is provided, only delete if it matches what already exists
                elif self[key] == value:
                    del self[key]
            except KeyError:
                pass

    def equals(self, tags: List[Dict[str, str]]) -> bool:
        flat_tags = {t["Key"]: t.get("Value") for t in tags} if tags else {}
        return self == flat_tags


class CertBundle(BaseModel):
    def __init__(
        self,
        account_id: str,
        certificate: bytes,
        private_key: bytes,
        chain: Optional[bytes] = None,
        region: str = "us-east-1",
        arn: Optional[str] = None,
        cert_type: str = "IMPORTED",
        cert_status: str = "ISSUED",
    ):
        self.created_at = utcnow()
        self.cert = certificate
        self.key = private_key
        # AWS always returns your chain + root CA
        self.chain = chain + b"\n" + AWS_ROOT_CA if chain else AWS_ROOT_CA
        self.tags = TagHolder()
        self.type = cert_type  # Should really be an enum
        self.status = cert_status  # Should really be an enum
        self.in_use_by: List[str] = []

        # Takes care of PEM checking
        self._key = self.validate_pk()
        self._cert = self.validate_certificate()
        # Extracting some common fields for ease of use
        # Have to search through cert.subject for OIDs
        self.common_name: Any = self._cert.subject.get_attributes_for_oid(
            OID_COMMON_NAME
        )[0].value
        if chain is not None:
            self.validate_chain()

        # TODO check cert is valid, or if self-signed then a chain is provided, otherwise
        # raise AWSValidationException('Provided certificate is not a valid self signed. Please provide either a valid self-signed certificate or certificate chain.')

        # Used for when one wants to overwrite an arn
        self.arn = arn or make_arn_for_certificate(account_id, region)

    @classmethod
    def generate_cert(
        cls,
        domain_name: str,
        account_id: str,
        region: str,
        sans: Optional[List[str]] = None,
    ) -> "CertBundle":
        unique_sans: Set[str] = set(sans) if sans else set()

        unique_sans.add(domain_name)
        unique_dns_names = [DNSName(item) for item in unique_sans]

        key = cryptography.hazmat.primitives.asymmetric.rsa.generate_private_key(
            public_exponent=65537, key_size=2048, backend=default_backend()
        )
        subject = cryptography.x509.Name(
            [
                cryptography.x509.NameAttribute(NameOID.COUNTRY_NAME, "US"),
                cryptography.x509.NameAttribute(NameOID.STATE_OR_PROVINCE_NAME, "CA"),
                cryptography.x509.NameAttribute(NameOID.LOCALITY_NAME, "San Francisco"),
                cryptography.x509.NameAttribute(
                    NameOID.ORGANIZATION_NAME, "My Company"
                ),
                cryptography.x509.NameAttribute(NameOID.COMMON_NAME, domain_name),
            ]
        )
        issuer = cryptography.x509.Name(
            [  # C = US, O = Amazon, OU = Server CA 1B, CN = Amazon
                cryptography.x509.NameAttribute(NameOID.COUNTRY_NAME, "US"),
                cryptography.x509.NameAttribute(NameOID.ORGANIZATION_NAME, "Amazon"),
                cryptography.x509.NameAttribute(
                    NameOID.ORGANIZATIONAL_UNIT_NAME, "Server CA 1B"
                ),
                cryptography.x509.NameAttribute(NameOID.COMMON_NAME, "Amazon"),
            ]
        )
        cert = (
            cryptography.x509.CertificateBuilder()
            .subject_name(subject)
            .issuer_name(issuer)
            .public_key(key.public_key())
            .serial_number(cryptography.x509.random_serial_number())
            .not_valid_before(utcnow())
            .not_valid_after(utcnow() + datetime.timedelta(days=365))
            .add_extension(
                cryptography.x509.SubjectAlternativeName(unique_dns_names),
                critical=False,
            )
            .sign(key, hashes.SHA512(), default_backend())
        )

        cert_armored = cert.public_bytes(serialization.Encoding.PEM)
        private_key = key.private_bytes(
            encoding=serialization.Encoding.PEM,
            format=serialization.PrivateFormat.TraditionalOpenSSL,
            encryption_algorithm=serialization.NoEncryption(),
        )

        return cls(
            certificate=cert_armored,
            private_key=private_key,
            cert_type="AMAZON_ISSUED",
            cert_status="PENDING_VALIDATION",
            account_id=account_id,
            region=region,
        )

    def validate_pk(self) -> Any:
        try:
            return serialization.load_pem_private_key(
                self.key, password=None, backend=default_backend()
            )
        except Exception as err:
            if isinstance(err, AWSValidationException):
                raise
            raise AWSValidationException(
                "The private key is not PEM-encoded or is not valid."
            )

    def validate_certificate(self) -> cryptography.x509.base.Certificate:
        try:
            _cert = cryptography.x509.load_pem_x509_certificate(
                self.cert, default_backend()
            )

            now = utcnow()
            if self._not_valid_after(_cert) < now:
                raise AWSValidationException(
                    "The certificate has expired, is not valid."
                )

            if self._not_valid_before(_cert) > now:
                raise AWSValidationException(
                    "The certificate is not in effect yet, is not valid."
                )

        except Exception as err:
            if isinstance(err, AWSValidationException):
                raise
            raise AWSValidationException(
                "The certificate is not PEM-encoded or is not valid."
            )
        return _cert

    def _not_valid_after(
        self, _cert: cryptography.x509.base.Certificate
    ) -> datetime.datetime:
        try:
            return _cert.not_valid_after_utc.replace(tzinfo=None)
        except AttributeError:
            return _cert.not_valid_after

    def _not_valid_before(
        self, _cert: cryptography.x509.base.Certificate
    ) -> datetime.datetime:
        try:
            return _cert.not_valid_before_utc.replace(tzinfo=None)
        except AttributeError:
            return _cert.not_valid_before

    def validate_chain(self) -> None:
        try:
            for cert_armored in self.chain.split(b"-\n-"):
                # Fix missing -'s on split
                cert_armored = re.sub(b"^----B", b"-----B", cert_armored)
                cert_armored = re.sub(b"E----$", b"E-----", cert_armored)
                cryptography.x509.load_pem_x509_certificate(
                    cert_armored, default_backend()
                )

                now = utcnow()
                if self._not_valid_after(self._cert) < now:
                    raise AWSValidationException(
                        "The certificate chain has expired, is not valid."
                    )

                if self._not_valid_before(self._cert) > now:
                    raise AWSValidationException(
                        "The certificate chain is not in effect yet, is not valid."
                    )

        except Exception as err:
            if isinstance(err, AWSValidationException):
                raise
            raise AWSValidationException(
                "The certificate is not PEM-encoded or is not valid."
            )

    def check(self) -> None:
        # Basically, if the certificate is pending, and then checked again after a
        # while, it will appear as if its been validated. The default wait time is 60
        # seconds but you can set an environment to change it.
        waited_seconds = (utcnow() - self.created_at).total_seconds()
        if (
            self.type == "AMAZON_ISSUED"
            and self.status == "PENDING_VALIDATION"
            and waited_seconds > settings.ACM_VALIDATION_WAIT
        ):
            self.status = "ISSUED"

    def describe(self) -> Dict[str, Any]:
        # 'RenewalSummary': {},  # Only when cert is amazon issued
        if self._key.key_size == 1024:
            key_algo = "RSA_1024"
        elif self._key.key_size == 2048:
            key_algo = "RSA_2048"
        else:
            key_algo = "EC_prime256v1"

        # Look for SANs
        try:
            san_obj: Any = self._cert.extensions.get_extension_for_oid(
                cryptography.x509.OID_SUBJECT_ALTERNATIVE_NAME
            )
        except cryptography.x509.ExtensionNotFound:
            san_obj = None
        sans = []
        if san_obj is not None:
            sans = [item.value for item in san_obj.value]

        result: Dict[str, Any] = {
            "Certificate": {
                "CertificateArn": self.arn,
                "DomainName": self.common_name,
                "InUseBy": self.in_use_by,
                "Issuer": self._cert.issuer.get_attributes_for_oid(OID_COMMON_NAME)[
                    0
                ].value,
                "KeyAlgorithm": key_algo,
                "NotAfter": datetime_to_epoch(self._not_valid_after(self._cert)),
                "NotBefore": datetime_to_epoch(self._not_valid_before(self._cert)),
                "Serial": str(self._cert.serial_number),
                "SignatureAlgorithm": self._cert.signature_algorithm_oid._name.upper().replace(
                    "ENCRYPTION", ""
                ),
                "Status": self.status,  # One of PENDING_VALIDATION, ISSUED, INACTIVE, EXPIRED, VALIDATION_TIMED_OUT, REVOKED, FAILED.
                "Subject": f"CN={self.common_name}",
                "SubjectAlternativeNames": sans,
                "Type": self.type,  # One of IMPORTED, AMAZON_ISSUED,
                "ExtendedKeyUsages": [],
                "RenewalEligibility": "INELIGIBLE",
                "Options": {"CertificateTransparencyLoggingPreference": "ENABLED"},
            }
        }

        domain_names = set(sans + [self.common_name])
        validation_options = []

        for san in domain_names:
            resource_record = {
                "Name": f"_d930b28be6c5927595552b219965053e.{san}.",
                "Type": "CNAME",
                "Value": "_c9edd76ee4a0e2a74388032f3861cc50.ykybfrwcxw.acm-validations.aws.",
            }
            validation_options.append(
                {
                    "DomainName": san,
                    "ValidationDomain": san,
                    "ValidationStatus": self.status,
                    "ValidationMethod": "DNS",
                    "ResourceRecord": resource_record,
                }
            )

        result["Certificate"]["DomainValidationOptions"] = validation_options

        if self.type == "IMPORTED":
            result["Certificate"]["ImportedAt"] = datetime_to_epoch(self.created_at)
        else:
            result["Certificate"]["CreatedAt"] = datetime_to_epoch(self.created_at)
            result["Certificate"]["IssuedAt"] = datetime_to_epoch(self.created_at)

        return result

    def serialize_pk(self, passphrase_bytes: bytes) -> str:
        pk_bytes = self._key.private_bytes(
            encoding=serialization.Encoding.PEM,
            format=serialization.PrivateFormat.PKCS8,
            encryption_algorithm=serialization.BestAvailableEncryption(
                passphrase_bytes
            ),
        )
        return pk_bytes.decode("utf-8")

    def __str__(self) -> str:
        return self.arn

    def __repr__(self) -> str:
        return "<Certificate>"


class AWSCertificateManagerBackend(BaseBackend):
    def __init__(self, region_name: str, account_id: str):
        super().__init__(region_name, account_id)
        self._certificates: Dict[str, CertBundle] = {}
        self._idempotency_tokens: Dict[str, Any] = {}

    def set_certificate_in_use_by(self, arn: str, load_balancer_name: str) -> None:
        if arn not in self._certificates:
            raise CertificateNotFound(arn=arn, account_id=self.account_id)

        cert_bundle = self._certificates[arn]
        cert_bundle.in_use_by.append(load_balancer_name)

    def _get_arn_from_idempotency_token(self, token: str) -> Optional[str]:
        """
        If token doesnt exist, return None, later it will be
        set with an expiry and arn.

        If token expiry has passed, delete entry and return None

        Else return ARN

        :param token: String token
        :return: None or ARN
        """
        now = utcnow()
        if token in self._idempotency_tokens:
            if self._idempotency_tokens[token]["expires"] < now:
                # Token has expired, new request
                del self._idempotency_tokens[token]
                return None
            else:
                return self._idempotency_tokens[token]["arn"]

        return None

    def _set_idempotency_token_arn(self, token: str, arn: str) -> None:
        self._idempotency_tokens[token] = {
            "arn": arn,
            "expires": utcnow() + datetime.timedelta(hours=1),
        }

    def import_certificate(
        self,
        certificate: bytes,
        private_key: bytes,
        chain: Optional[bytes],
        arn: Optional[str],
        tags: List[Dict[str, str]],
    ) -> str:
        if arn is not None:
            if arn not in self._certificates:
                raise CertificateNotFound(arn=arn, account_id=self.account_id)
            else:
                # Will reuse provided ARN
                bundle = CertBundle(
                    self.account_id,
                    certificate,
                    private_key,
                    chain=chain,
                    region=self.region_name,
                    arn=arn,
                )
        else:
            # Will generate a random ARN
            bundle = CertBundle(
                self.account_id,
                certificate,
                private_key,
                chain=chain,
                region=self.region_name,
            )

        self._certificates[bundle.arn] = bundle

        if tags:
            self.add_tags_to_certificate(bundle.arn, tags)

        return bundle.arn

    def list_certificates(self, statuses: List[str]) -> Iterable[CertBundle]:
        for arn in self._certificates.keys():
            cert = self.get_certificate(arn)
            if not statuses or cert.status in statuses:
                yield cert

    def get_certificate(self, arn: str) -> CertBundle:
        if arn not in self._certificates:
            raise CertificateNotFound(arn=arn, account_id=self.account_id)

        cert_bundle = self._certificates[arn]
        cert_bundle.check()
        return cert_bundle

    def describe_certificate(self, arn: str) -> CertBundle:
        return self.get_certificate(arn)

    def delete_certificate(self, arn: str) -> None:
        if arn not in self._certificates:
            raise CertificateNotFound(arn=arn, account_id=self.account_id)

        del self._certificates[arn]

    def request_certificate(
        self,
        domain_name: str,
        idempotency_token: str,
        subject_alt_names: List[str],
        tags: List[Dict[str, str]],
    ) -> str:
        """
        The parameter DomainValidationOptions has not yet been implemented
        """
        if idempotency_token is not None:
            arn = self._get_arn_from_idempotency_token(idempotency_token)
            if arn and self._certificates[arn].tags.equals(tags):
                return arn

        cert = CertBundle.generate_cert(
            domain_name,
            account_id=self.account_id,
            region=self.region_name,
            sans=subject_alt_names,
        )
        if idempotency_token is not None:
            self._set_idempotency_token_arn(idempotency_token, cert.arn)
        self._certificates[cert.arn] = cert

        if tags:
            cert.tags.add(tags)

        return cert.arn

    def add_tags_to_certificate(self, arn: str, tags: List[Dict[str, str]]) -> None:
        # get_cert does arn check
        cert_bundle = self.get_certificate(arn)
        cert_bundle.tags.add(tags)

    def remove_tags_from_certificate(
        self, arn: str, tags: List[Dict[str, str]]
    ) -> None:
        # get_cert does arn check
        cert_bundle = self.get_certificate(arn)
        cert_bundle.tags.remove(tags)

    def export_certificate(
        self, certificate_arn: str, passphrase: str
    ) -> Tuple[str, str, str]:
        passphrase_bytes = base64.standard_b64decode(passphrase)
        cert_bundle = self.get_certificate(certificate_arn)

        certificate = cert_bundle.cert.decode()
        certificate_chain = cert_bundle.chain.decode()
        private_key = cert_bundle.serialize_pk(passphrase_bytes)

        return certificate, certificate_chain, private_key


acm_backends = BackendDict(AWSCertificateManagerBackend, "acm")
