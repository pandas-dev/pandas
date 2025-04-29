"""Handles incoming acmpca requests, invokes methods, returns responses."""

import base64
import json

from moto.core.responses import BaseResponse

from .models import ACMPCABackend, acmpca_backends


class ACMPCAResponse(BaseResponse):
    """Handler for ACMPCA requests and responses."""

    def __init__(self) -> None:
        super().__init__(service_name="acmpca")

    @property
    def acmpca_backend(self) -> ACMPCABackend:
        """Return backend instance specific for this region."""
        return acmpca_backends[self.current_account][self.region]

    def create_certificate_authority(self) -> str:
        params = json.loads(self.body)
        certificate_authority_configuration = params.get(
            "CertificateAuthorityConfiguration"
        )
        revocation_configuration = params.get("RevocationConfiguration")
        certificate_authority_type = params.get("CertificateAuthorityType")
        security_standard = params.get("KeyStorageSecurityStandard")
        tags = params.get("Tags")
        certificate_authority_arn = self.acmpca_backend.create_certificate_authority(
            certificate_authority_configuration=certificate_authority_configuration,
            revocation_configuration=revocation_configuration,
            certificate_authority_type=certificate_authority_type,
            security_standard=security_standard,
            tags=tags,
        )
        return json.dumps(dict(CertificateAuthorityArn=certificate_authority_arn))

    def describe_certificate_authority(self) -> str:
        params = json.loads(self.body)
        certificate_authority_arn = params.get("CertificateAuthorityArn")
        certificate_authority = self.acmpca_backend.describe_certificate_authority(
            certificate_authority_arn=certificate_authority_arn,
        )
        return json.dumps(dict(CertificateAuthority=certificate_authority.to_json()))

    def get_certificate_authority_certificate(self) -> str:
        params = json.loads(self.body)
        certificate_authority_arn = params.get("CertificateAuthorityArn")
        (
            certificate,
            certificate_chain,
        ) = self.acmpca_backend.get_certificate_authority_certificate(
            certificate_authority_arn=certificate_authority_arn,
        )
        return json.dumps(
            dict(
                Certificate=certificate.decode("utf-8"),
                CertificateChain=certificate_chain,
            )
        )

    def get_certificate_authority_csr(self) -> str:
        params = json.loads(self.body)
        certificate_authority_arn = params.get("CertificateAuthorityArn")
        csr = self.acmpca_backend.get_certificate_authority_csr(
            certificate_authority_arn=certificate_authority_arn,
        )
        return json.dumps(dict(Csr=csr.decode("utf-8")))

    def list_tags(self) -> str:
        params = json.loads(self.body)
        certificate_authority_arn = params.get("CertificateAuthorityArn")
        tags = self.acmpca_backend.list_tags(
            certificate_authority_arn=certificate_authority_arn
        )
        return json.dumps(tags)

    def update_certificate_authority(self) -> str:
        params = json.loads(self.body)
        certificate_authority_arn = params.get("CertificateAuthorityArn")
        revocation_configuration = params.get("RevocationConfiguration")
        status = params.get("Status")
        self.acmpca_backend.update_certificate_authority(
            certificate_authority_arn=certificate_authority_arn,
            revocation_configuration=revocation_configuration,
            status=status,
        )
        return "{}"

    def delete_certificate_authority(self) -> str:
        params = json.loads(self.body)
        certificate_authority_arn = params.get("CertificateAuthorityArn")
        self.acmpca_backend.delete_certificate_authority(
            certificate_authority_arn=certificate_authority_arn
        )
        return "{}"

    def issue_certificate(self) -> str:
        params = json.loads(self.body)
        certificate_authority_arn = params.get("CertificateAuthorityArn")
        template_arn = params.get("TemplateArn")
        csr = params.get("Csr").encode("utf-8")
        certificate_arn = self.acmpca_backend.issue_certificate(
            certificate_authority_arn=certificate_authority_arn,
            csr=csr,
            template_arn=template_arn,
        )
        return json.dumps(dict(CertificateArn=certificate_arn))

    def get_certificate(self) -> str:
        params = json.loads(self.body)
        certificate_authority_arn = params.get("CertificateAuthorityArn")
        certificate_arn = params.get("CertificateArn")
        certificate, certificate_chain = self.acmpca_backend.get_certificate(
            certificate_authority_arn=certificate_authority_arn,
            certificate_arn=certificate_arn,
        )
        return json.dumps(
            dict(
                Certificate=certificate.decode("utf-8"),
                CertificateChain=certificate_chain,
            )
        )

    def import_certificate_authority_certificate(self) -> str:
        params = json.loads(self.body)
        certificate_authority_arn = params.get("CertificateAuthorityArn")
        certificate = params.get("Certificate")
        certificate_bytes = base64.b64decode(certificate)
        certificate_chain = params.get("CertificateChain")
        self.acmpca_backend.import_certificate_authority_certificate(
            certificate_authority_arn=certificate_authority_arn,
            certificate=certificate_bytes,
            certificate_chain=certificate_chain,
        )
        return "{}"

    def revoke_certificate(self) -> str:
        params = json.loads(self.body)
        certificate_authority_arn = params.get("CertificateAuthorityArn")
        certificate_serial = params.get("CertificateSerial")
        revocation_reason = params.get("RevocationReason")
        self.acmpca_backend.revoke_certificate(
            certificate_authority_arn=certificate_authority_arn,
            certificate_serial=certificate_serial,
            revocation_reason=revocation_reason,
        )
        return "{}"

    def tag_certificate_authority(self) -> str:
        params = json.loads(self.body)
        certificate_authority_arn = params.get("CertificateAuthorityArn")
        tags = params.get("Tags")
        self.acmpca_backend.tag_certificate_authority(
            certificate_authority_arn=certificate_authority_arn,
            tags=tags,
        )
        return "{}"

    def untag_certificate_authority(self) -> str:
        params = json.loads(self.body)
        certificate_authority_arn = params.get("CertificateAuthorityArn")
        tags = params.get("Tags")
        self.acmpca_backend.untag_certificate_authority(
            certificate_authority_arn=certificate_authority_arn,
            tags=tags,
        )
        return "{}"
