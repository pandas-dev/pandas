"""ACMPCABackend class with methods for supported APIs."""

import base64
import datetime
from typing import Any, Dict, List, Optional, Tuple, cast

from cryptography import x509
from cryptography.hazmat.backends import default_backend
from cryptography.hazmat.primitives import hashes, serialization
from cryptography.hazmat.primitives.asymmetric import rsa

from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.common_models import BaseModel
from moto.core.utils import unix_time, utcnow
from moto.moto_api._internal import mock_random
from moto.utilities.tagging_service import TaggingService
from moto.utilities.utils import get_partition

from .exceptions import (
    InvalidS3ObjectAclInCrlConfiguration,
    InvalidStateException,
    MalformedCertificateAuthorityException,
    ResourceNotFoundException,
)


class CertificateAuthority(BaseModel):
    def __init__(
        self,
        region: str,
        account_id: str,
        certificate_authority_configuration: Dict[str, Any],
        certificate_authority_type: str,
        revocation_configuration: Dict[str, Any],
        security_standard: Optional[str],
    ):
        self.id = mock_random.uuid4()
        self.arn = f"arn:{get_partition(region)}:acm-pca:{region}:{account_id}:certificate-authority/{self.id}"
        self.account_id = account_id
        self.region_name = region
        self.certificate_authority_configuration = certificate_authority_configuration
        self.certificate_authority_type = certificate_authority_type
        self.revocation_configuration: Dict[str, Any] = {
            "CrlConfiguration": {"Enabled": False}
        }
        self.set_revocation_configuration(revocation_configuration)
        self.created_at = unix_time()
        self.updated_at: Optional[float] = None
        self.status = "PENDING_CERTIFICATE"
        self.usage_mode = "SHORT_LIVED_CERTIFICATE"
        self.security_standard = security_standard or "FIPS_140_2_LEVEL_3_OR_HIGHER"

        self.password = str(mock_random.uuid4()).encode("utf-8")

        private_key = rsa.generate_private_key(public_exponent=65537, key_size=2048)
        self.private_bytes = private_key.private_bytes(
            encoding=serialization.Encoding.PEM,
            format=serialization.PrivateFormat.TraditionalOpenSSL,
            encryption_algorithm=serialization.BestAvailableEncryption(self.password),
        )

        self.certificate_bytes: bytes = b""
        self.certificate_chain: Optional[bytes] = None
        self.issued_certificates: Dict[str, bytes] = dict()

        self.subject = self.certificate_authority_configuration.get("Subject", {})

    def generate_cert(
        self,
        subject: x509.Name,
        public_key: rsa.RSAPublicKey,
        extensions: List[Tuple[x509.ExtensionType, bool]],
    ) -> bytes:
        builder = (
            x509.CertificateBuilder()
            .subject_name(subject)
            .issuer_name(self.issuer)
            .public_key(public_key)
            .serial_number(x509.random_serial_number())
            .not_valid_before(utcnow())
            .not_valid_after(utcnow() + datetime.timedelta(days=365))
        )

        for extension, critical in extensions:
            builder = builder.add_extension(extension, critical)

        cert = builder.sign(self.key, hashes.SHA512(), default_backend())

        return cert.public_bytes(serialization.Encoding.PEM)

    @property
    def key(self) -> rsa.RSAPrivateKey:
        private_key = serialization.load_pem_private_key(
            self.private_bytes,
            password=self.password,
        )
        return cast(rsa.RSAPrivateKey, private_key)

    @property
    def certificate(self) -> Optional[x509.Certificate]:
        if self.certificate_bytes:
            return x509.load_pem_x509_certificate(self.certificate_bytes)
        return None

    @property
    def issuer(self) -> x509.Name:
        name_attributes = []
        if "Country" in self.subject:
            name_attributes.append(
                x509.NameAttribute(x509.NameOID.COUNTRY_NAME, self.subject["Country"])
            )
        if "State" in self.subject:
            name_attributes.append(
                x509.NameAttribute(
                    x509.NameOID.STATE_OR_PROVINCE_NAME, self.subject["State"]
                )
            )
        if "Organization" in self.subject:
            name_attributes.append(
                x509.NameAttribute(
                    x509.NameOID.ORGANIZATION_NAME, self.subject["Organization"]
                )
            )
        if "OrganizationalUnit" in self.subject:
            name_attributes.append(
                x509.NameAttribute(
                    x509.NameOID.ORGANIZATIONAL_UNIT_NAME,
                    self.subject["OrganizationalUnit"],
                )
            )
        if "CommonName" in self.subject:
            name_attributes.append(
                x509.NameAttribute(x509.NameOID.COMMON_NAME, self.subject["CommonName"])
            )
        return x509.Name(name_attributes)

    @property
    def csr(self) -> bytes:
        csr = (
            x509.CertificateSigningRequestBuilder()
            .subject_name(self.issuer)
            .add_extension(
                x509.BasicConstraints(ca=True, path_length=None),
                critical=True,
            )
            .sign(self.key, hashes.SHA256())
        )
        return csr.public_bytes(serialization.Encoding.PEM)

    def issue_certificate(self, csr_bytes: bytes, template_arn: Optional[str]) -> str:
        csr = x509.load_pem_x509_csr(base64.b64decode(csr_bytes))
        extensions = self._x509_extensions(csr, template_arn)
        new_cert = self.generate_cert(
            subject=csr.subject,
            public_key=csr.public_key(),  # type: ignore[arg-type]
            extensions=extensions,
        )

        cert_id = str(mock_random.uuid4()).replace("-", "")
        cert_arn = f"arn:{get_partition(self.region_name)}:acm-pca:{self.region_name}:{self.account_id}:certificate-authority/{self.id}/certificate/{cert_id}"
        self.issued_certificates[cert_arn] = new_cert
        return cert_arn

    def _x509_extensions(
        self, csr: x509.CertificateSigningRequest, template_arn: Optional[str]
    ) -> List[Tuple[x509.ExtensionType, bool]]:
        """
        Uses a PCA certificate template ARN to return a list of X.509 extensions.
        These extensions are part of the constructed certificate.

        See https://docs.aws.amazon.com/privateca/latest/userguide/UsingTemplates.html
        """
        extensions = []

        if template_arn == "arn:aws:acm-pca:::template/RootCACertificate/V1":
            extensions.extend(
                [
                    (
                        x509.BasicConstraints(ca=True, path_length=None),
                        True,
                    ),
                    (
                        x509.KeyUsage(
                            crl_sign=True,
                            key_cert_sign=True,
                            digital_signature=True,
                            content_commitment=False,
                            key_encipherment=False,
                            data_encipherment=False,
                            key_agreement=False,
                            encipher_only=False,
                            decipher_only=False,
                        ),
                        True,
                    ),
                    (
                        x509.SubjectKeyIdentifier.from_public_key(csr.public_key()),
                        False,
                    ),
                ]
            )

        elif template_arn in (
            "arn:aws:acm-pca:::template/EndEntityCertificate/V1",
            None,
        ):
            extensions.extend(
                [
                    (
                        x509.BasicConstraints(ca=False, path_length=None),
                        True,
                    ),
                    (
                        x509.AuthorityKeyIdentifier.from_issuer_public_key(
                            self.key.public_key()
                        ),
                        False,
                    ),
                    (
                        x509.SubjectKeyIdentifier.from_public_key(csr.public_key()),
                        False,
                    ),
                    (
                        x509.KeyUsage(
                            crl_sign=False,
                            key_cert_sign=False,
                            digital_signature=True,
                            content_commitment=False,
                            key_encipherment=True,
                            data_encipherment=False,
                            key_agreement=False,
                            encipher_only=False,
                            decipher_only=False,
                        ),
                        True,
                    ),
                    (
                        x509.ExtendedKeyUsage(
                            [
                                x509.ExtendedKeyUsageOID.SERVER_AUTH,
                                x509.ExtendedKeyUsageOID.CLIENT_AUTH,
                            ]
                        ),
                        False,
                    ),
                ]
            )

        cn = csr.subject.get_attributes_for_oid(x509.NameOID.COMMON_NAME)
        if cn:
            extensions.append(
                (
                    x509.SubjectAlternativeName([x509.DNSName(cn[0].value)]),  # type: ignore[arg-type]
                    False,
                ),
            )

        return extensions

    def get_certificate(self, certificate_arn: str) -> bytes:
        return self.issued_certificates[certificate_arn]

    def set_revocation_configuration(
        self, revocation_configuration: Optional[Dict[str, Any]]
    ) -> None:
        if revocation_configuration is not None:
            self.revocation_configuration = revocation_configuration
            if "CrlConfiguration" in self.revocation_configuration:
                acl = self.revocation_configuration["CrlConfiguration"].get(
                    "S3ObjectAcl", None
                )
                if acl is None:
                    self.revocation_configuration["CrlConfiguration"]["S3ObjectAcl"] = (
                        "PUBLIC_READ"
                    )
                else:
                    if acl not in ["PUBLIC_READ", "BUCKET_OWNER_FULL_CONTROL"]:
                        raise InvalidS3ObjectAclInCrlConfiguration(acl)

    @property
    def not_valid_after(self) -> Optional[float]:
        if self.certificate is None:
            return None
        try:
            return unix_time(self.certificate.not_valid_after_utc.replace(tzinfo=None))
        except AttributeError:
            return unix_time(self.certificate.not_valid_after)

    @property
    def not_valid_before(self) -> Optional[float]:
        if self.certificate is None:
            return None
        try:
            return unix_time(self.certificate.not_valid_before_utc.replace(tzinfo=None))
        except AttributeError:
            return unix_time(self.certificate.not_valid_before)

    def import_certificate_authority_certificate(
        self, certificate: bytes, certificate_chain: Optional[bytes]
    ) -> None:
        try:
            x509.load_pem_x509_certificate(certificate)
        except ValueError:
            raise MalformedCertificateAuthorityException()

        self.certificate_bytes = certificate
        self.certificate_chain = certificate_chain
        self.status = "ACTIVE"
        self.updated_at = unix_time()

    def to_json(self) -> Dict[str, Any]:
        dct = {
            "Arn": self.arn,
            "OwnerAccount": self.account_id,
            "CertificateAuthorityConfiguration": self.certificate_authority_configuration,
            "Type": self.certificate_authority_type,
            "RevocationConfiguration": self.revocation_configuration,
            "CreatedAt": self.created_at,
            "Status": self.status,
            "UsageMode": self.usage_mode,
            "KeyStorageSecurityStandard": self.security_standard,
        }
        if self.updated_at:
            dct["LastStateChangeAt"] = self.updated_at
        if self.certificate:
            dct.update(
                {
                    "NotBefore": self.not_valid_before,
                    "NotAfter": self.not_valid_after,
                }
            )
        return dct


class ACMPCABackend(BaseBackend):
    """Implementation of ACMPCA APIs."""

    def __init__(self, region_name: str, account_id: str):
        super().__init__(region_name, account_id)
        self.certificate_authorities: Dict[str, CertificateAuthority] = dict()
        self.tagger = TaggingService()

    def create_certificate_authority(
        self,
        certificate_authority_configuration: Dict[str, Any],
        revocation_configuration: Dict[str, Any],
        certificate_authority_type: str,
        security_standard: Optional[str],
        tags: List[Dict[str, str]],
    ) -> str:
        """
        The following parameters are not yet implemented: IdempotencyToken, KeyStorageSecurityStandard, UsageMode
        """
        authority = CertificateAuthority(
            region=self.region_name,
            account_id=self.account_id,
            certificate_authority_configuration=certificate_authority_configuration,
            certificate_authority_type=certificate_authority_type,
            revocation_configuration=revocation_configuration,
            security_standard=security_standard,
        )
        self.certificate_authorities[authority.arn] = authority
        if tags:
            self.tagger.tag_resource(authority.arn, tags)
        return authority.arn

    def describe_certificate_authority(
        self, certificate_authority_arn: str
    ) -> CertificateAuthority:
        if certificate_authority_arn not in self.certificate_authorities:
            raise ResourceNotFoundException(certificate_authority_arn)
        return self.certificate_authorities[certificate_authority_arn]

    def get_certificate_authority_certificate(
        self, certificate_authority_arn: str
    ) -> Tuple[bytes, Optional[bytes]]:
        ca = self.describe_certificate_authority(certificate_authority_arn)
        if ca.status != "ACTIVE":
            raise InvalidStateException(certificate_authority_arn)
        return ca.certificate_bytes, ca.certificate_chain

    def get_certificate_authority_csr(self, certificate_authority_arn: str) -> bytes:
        ca = self.describe_certificate_authority(certificate_authority_arn)
        return ca.csr

    def list_tags(
        self, certificate_authority_arn: str
    ) -> Dict[str, List[Dict[str, str]]]:
        """
        Pagination is not yet implemented
        """
        return self.tagger.list_tags_for_resource(certificate_authority_arn)

    def update_certificate_authority(
        self,
        certificate_authority_arn: str,
        revocation_configuration: Dict[str, Any],
        status: str,
    ) -> None:
        ca = self.describe_certificate_authority(certificate_authority_arn)
        if status is not None:
            ca.status = status
        ca.set_revocation_configuration(revocation_configuration)
        ca.updated_at = unix_time()

    def delete_certificate_authority(self, certificate_authority_arn: str) -> None:
        ca = self.describe_certificate_authority(certificate_authority_arn)
        ca.status = "DELETED"

    def issue_certificate(
        self, certificate_authority_arn: str, csr: bytes, template_arn: Optional[str]
    ) -> str:
        """
        The following parameters are not yet implemented: ApiPassthrough, SigningAlgorithm, Validity, ValidityNotBefore, IdempotencyToken
        Some fields of the resulting certificate will have default values, instead of using the CSR
        """
        ca = self.describe_certificate_authority(certificate_authority_arn)
        certificate_arn = ca.issue_certificate(csr, template_arn)
        return certificate_arn

    def get_certificate(
        self, certificate_authority_arn: str, certificate_arn: str
    ) -> Tuple[bytes, Optional[str]]:
        """
        The CertificateChain will always return None for now
        """
        ca = self.describe_certificate_authority(certificate_authority_arn)
        certificate = ca.get_certificate(certificate_arn)
        certificate_chain = None
        return certificate, certificate_chain

    def import_certificate_authority_certificate(
        self,
        certificate_authority_arn: str,
        certificate: bytes,
        certificate_chain: Optional[bytes],
    ) -> None:
        ca = self.describe_certificate_authority(certificate_authority_arn)
        ca.import_certificate_authority_certificate(certificate, certificate_chain)

    def revoke_certificate(
        self,
        certificate_authority_arn: str,
        certificate_serial: str,
        revocation_reason: str,
    ) -> None:
        """
        This is currently a NO-OP
        """

    def tag_certificate_authority(
        self, certificate_authority_arn: str, tags: List[Dict[str, str]]
    ) -> None:
        self.tagger.tag_resource(certificate_authority_arn, tags)

    def untag_certificate_authority(
        self, certificate_authority_arn: str, tags: List[Dict[str, str]]
    ) -> None:
        self.tagger.untag_resource_using_tags(certificate_authority_arn, tags)


acmpca_backends = BackendDict(ACMPCABackend, "acm-pca")
