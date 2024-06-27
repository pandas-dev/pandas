import base64
import json
from typing import Dict, List, Tuple, Union

from moto.core.responses import BaseResponse

from .exceptions import AWSValidationException
from .models import AWSCertificateManagerBackend, acm_backends

GENERIC_RESPONSE_TYPE = Union[str, Tuple[str, Dict[str, int]]]


class AWSCertificateManagerResponse(BaseResponse):
    def __init__(self) -> None:
        super().__init__(service_name="acm")

    @property
    def acm_backend(self) -> AWSCertificateManagerBackend:
        return acm_backends[self.current_account][self.region]

    def add_tags_to_certificate(self) -> GENERIC_RESPONSE_TYPE:
        arn = self._get_param("CertificateArn")
        tags = self._get_param("Tags")

        if arn is None:
            msg = "A required parameter for the specified action is not supplied."
            return (
                json.dumps({"__type": "MissingParameter", "message": msg}),
                dict(status=400),
            )

        self.acm_backend.add_tags_to_certificate(arn, tags)

        return ""

    def delete_certificate(self) -> GENERIC_RESPONSE_TYPE:
        arn = self._get_param("CertificateArn")

        if arn is None:
            msg = "A required parameter for the specified action is not supplied."
            return (
                json.dumps({"__type": "MissingParameter", "message": msg}),
                dict(status=400),
            )

        self.acm_backend.delete_certificate(arn)

        return ""

    def describe_certificate(self) -> GENERIC_RESPONSE_TYPE:
        arn = self._get_param("CertificateArn")

        if arn is None:
            msg = "A required parameter for the specified action is not supplied."
            return (
                json.dumps({"__type": "MissingParameter", "message": msg}),
                dict(status=400),
            )

        cert_bundle = self.acm_backend.describe_certificate(arn)

        return json.dumps(cert_bundle.describe())

    def get_certificate(self) -> GENERIC_RESPONSE_TYPE:
        arn = self._get_param("CertificateArn")

        if arn is None:
            msg = "A required parameter for the specified action is not supplied."
            return (
                json.dumps({"__type": "MissingParameter", "message": msg}),
                dict(status=400),
            )

        cert_bundle = self.acm_backend.get_certificate(arn)

        result = {
            "Certificate": cert_bundle.cert.decode(),
            "CertificateChain": cert_bundle.chain.decode(),
        }
        return json.dumps(result)

    def import_certificate(self) -> str:
        """
        Returns errors on:
        Certificate, PrivateKey or Chain not being properly formatted
        Arn not existing if its provided
        PrivateKey size > 2048
        Certificate expired or is not yet in effect

        Does not return errors on:
        Checking Certificate is legit, or a selfsigned chain is provided

        :return: str(JSON) for response
        """
        certificate = self._get_param("Certificate")
        private_key = self._get_param("PrivateKey")
        chain = self._get_param("CertificateChain")  # Optional
        current_arn = self._get_param("CertificateArn")  # Optional
        tags = self._get_param("Tags")  # Optional

        # Simple parameter decoding. Rather do it here as its a data transport decision not part of the
        # actual data
        try:
            certificate = base64.standard_b64decode(certificate)
        except Exception:
            raise AWSValidationException(
                "The certificate is not PEM-encoded or is not valid."
            )
        try:
            private_key = base64.standard_b64decode(private_key)
        except Exception:
            raise AWSValidationException(
                "The private key is not PEM-encoded or is not valid."
            )
        if chain is not None:
            try:
                chain = base64.standard_b64decode(chain)
            except Exception:
                raise AWSValidationException(
                    "The certificate chain is not PEM-encoded or is not valid."
                )

        arn = self.acm_backend.import_certificate(
            certificate, private_key, chain=chain, arn=current_arn, tags=tags
        )

        return json.dumps({"CertificateArn": arn})

    def list_certificates(self) -> str:
        certs = []
        statuses = self._get_param("CertificateStatuses")
        for cert_bundle in self.acm_backend.list_certificates(statuses):
            certs.append(cert_bundle.describe()["Certificate"])

        result = {"CertificateSummaryList": certs}
        return json.dumps(result)

    def list_tags_for_certificate(self) -> GENERIC_RESPONSE_TYPE:
        arn = self._get_param("CertificateArn")

        if arn is None:
            msg = "A required parameter for the specified action is not supplied."
            return json.dumps({"__type": "MissingParameter", "message": msg}), dict(
                status=400
            )

        cert_bundle = self.acm_backend.get_certificate(arn)

        result: Dict[str, List[Dict[str, str]]] = {"Tags": []}
        # Tag "objects" can not contain the Value part
        for key, value in cert_bundle.tags.items():
            tag_dict = {"Key": key}
            if value is not None:
                tag_dict["Value"] = value
            result["Tags"].append(tag_dict)

        return json.dumps(result)

    def remove_tags_from_certificate(self) -> GENERIC_RESPONSE_TYPE:
        arn = self._get_param("CertificateArn")
        tags = self._get_param("Tags")

        if arn is None:
            msg = "A required parameter for the specified action is not supplied."
            return (
                json.dumps({"__type": "MissingParameter", "message": msg}),
                dict(status=400),
            )

        self.acm_backend.remove_tags_from_certificate(arn, tags)

        return ""

    def request_certificate(self) -> GENERIC_RESPONSE_TYPE:
        domain_name = self._get_param("DomainName")
        idempotency_token = self._get_param("IdempotencyToken")
        subject_alt_names = self._get_param("SubjectAlternativeNames")
        tags = self._get_param("Tags")  # Optional

        if subject_alt_names is not None and len(subject_alt_names) > 10:
            # There is initial AWS limit of 10
            msg = (
                "An ACM limit has been exceeded. Need to request SAN limit to be raised"
            )
            return (
                json.dumps({"__type": "LimitExceededException", "message": msg}),
                dict(status=400),
            )

        arn = self.acm_backend.request_certificate(
            domain_name,
            idempotency_token,
            subject_alt_names,
            tags,
        )

        return json.dumps({"CertificateArn": arn})

    def resend_validation_email(self) -> GENERIC_RESPONSE_TYPE:
        arn = self._get_param("CertificateArn")
        domain = self._get_param("Domain")
        # ValidationDomain not used yet.
        # Contains domain which is equal to or a subset of Domain
        # that AWS will send validation emails to
        # https://docs.aws.amazon.com/acm/latest/APIReference/API_ResendValidationEmail.html
        # validation_domain = self._get_param('ValidationDomain')

        if arn is None:
            msg = "A required parameter for the specified action is not supplied."
            return (
                json.dumps({"__type": "MissingParameter", "message": msg}),
                dict(status=400),
            )

        cert_bundle = self.acm_backend.get_certificate(arn)

        if cert_bundle.common_name != domain:
            msg = "Parameter Domain does not match certificate domain"
            _type = "InvalidDomainValidationOptionsException"
            return json.dumps({"__type": _type, "message": msg}), dict(status=400)

        return ""

    def export_certificate(self) -> GENERIC_RESPONSE_TYPE:
        certificate_arn = self._get_param("CertificateArn")
        passphrase = self._get_param("Passphrase")

        if certificate_arn is None:
            msg = "A required parameter for the specified action is not supplied."
            return (
                json.dumps({"__type": "MissingParameter", "message": msg}),
                dict(status=400),
            )

        (
            certificate,
            certificate_chain,
            private_key,
        ) = self.acm_backend.export_certificate(
            certificate_arn=certificate_arn, passphrase=passphrase
        )
        return json.dumps(
            dict(
                Certificate=certificate,
                CertificateChain=certificate_chain,
                PrivateKey=private_key,
            )
        )
