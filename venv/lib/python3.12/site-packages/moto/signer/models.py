from typing import Any, Dict, List, Optional

from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.common_models import BaseModel
from moto.moto_api._internal import mock_random
from moto.utilities.tagging_service import TaggingService
from moto.utilities.utils import get_partition


class SigningProfile(BaseModel):
    def __init__(
        self,
        backend: "SignerBackend",
        name: str,
        platform_id: str,
        signing_material: Dict[str, str],
        signature_validity_period: Optional[Dict[str, Any]],
    ):
        self.name = name
        self.platform_id = platform_id
        self.signature_validity_period = signature_validity_period or {
            "value": 135,
            "type": "MONTHS",
        }
        self.backend = backend

        self.status = "Active"
        self.arn = f"arn:{get_partition(backend.region_name)}:signer:{backend.region_name}:{backend.account_id}:/signing-profiles/{name}"
        self.profile_version = mock_random.get_random_hex(10)
        self.profile_version_arn = f"{self.arn}/{self.profile_version}"
        self.signing_material = signing_material

    def cancel(self) -> None:
        self.status = "Canceled"

    def to_dict(self, full: bool = True) -> Dict[str, Any]:
        small: Dict[str, Any] = {
            "arn": self.arn,
            "profileVersion": self.profile_version,
            "profileVersionArn": self.profile_version_arn,
        }
        if full:
            small.update(
                {
                    "status": self.status,
                    "profileName": self.name,
                    "platformId": self.platform_id,
                    "signatureValidityPeriod": self.signature_validity_period,
                    "signingMaterial": self.signing_material,
                    "platformDisplayName": next(
                        (
                            p["displayName"]
                            for p in SignerBackend.platforms
                            if p["platformId"] == self.platform_id
                        ),
                        None,
                    ),
                }
            )
            tags = self.backend.list_tags_for_resource(self.arn)
            if tags:
                small.update({"tags": tags})
        return small


class SignerBackend(BaseBackend):
    """Implementation of signer APIs."""

    platforms = [
        {
            "platformId": "AWSIoTDeviceManagement-SHA256-ECDSA",
            "displayName": "AWS IoT Device Management SHA256-ECDSA ",
            "partner": "AWSIoTDeviceManagement",
            "target": "SHA256-ECDSA",
            "category": "AWS",
            "signingConfiguration": {
                "encryptionAlgorithmOptions": {
                    "allowedValues": ["ECDSA"],
                    "defaultValue": "ECDSA",
                },
                "hashAlgorithmOptions": {
                    "allowedValues": ["SHA256"],
                    "defaultValue": "SHA256",
                },
            },
            "signingImageFormat": {
                "supportedFormats": ["JSONDetached"],
                "defaultFormat": "JSONDetached",
            },
            "maxSizeInMB": 2048,
            "revocationSupported": False,
        },
        {
            "platformId": "AWSLambda-SHA384-ECDSA",
            "displayName": "AWS Lambda",
            "partner": "AWSLambda",
            "target": "SHA384-ECDSA",
            "category": "AWS",
            "signingConfiguration": {
                "encryptionAlgorithmOptions": {
                    "allowedValues": ["ECDSA"],
                    "defaultValue": "ECDSA",
                },
                "hashAlgorithmOptions": {
                    "allowedValues": ["SHA384"],
                    "defaultValue": "SHA384",
                },
            },
            "signingImageFormat": {
                "supportedFormats": ["JSONDetached"],
                "defaultFormat": "JSONDetached",
            },
            "maxSizeInMB": 250,
            "revocationSupported": True,
        },
        {
            "platformId": "AmazonFreeRTOS-TI-CC3220SF",
            "displayName": "Amazon FreeRTOS SHA1-RSA CC3220SF-Format",
            "partner": "AmazonFreeRTOS",
            "target": "SHA1-RSA-TISHA1",
            "category": "AWS",
            "signingConfiguration": {
                "encryptionAlgorithmOptions": {
                    "allowedValues": ["RSA"],
                    "defaultValue": "RSA",
                },
                "hashAlgorithmOptions": {
                    "allowedValues": ["SHA1"],
                    "defaultValue": "SHA1",
                },
            },
            "signingImageFormat": {
                "supportedFormats": ["JSONEmbedded"],
                "defaultFormat": "JSONEmbedded",
            },
            "maxSizeInMB": 16,
            "revocationSupported": False,
        },
        {
            "platformId": "AmazonFreeRTOS-Default",
            "displayName": "Amazon FreeRTOS SHA256-ECDSA",
            "partner": "AmazonFreeRTOS",
            "target": "SHA256-ECDSA",
            "category": "AWS",
            "signingConfiguration": {
                "encryptionAlgorithmOptions": {
                    "allowedValues": ["ECDSA", "RSA"],
                    "defaultValue": "ECDSA",
                },
                "hashAlgorithmOptions": {
                    "allowedValues": ["SHA256"],
                    "defaultValue": "SHA256",
                },
            },
            "signingImageFormat": {
                "supportedFormats": ["JSONEmbedded"],
                "defaultFormat": "JSONEmbedded",
            },
            "maxSizeInMB": 16,
            "revocationSupported": False,
        },
    ]

    def __init__(self, region_name: str, account_id: str):
        super().__init__(region_name, account_id)
        self.signing_profiles: Dict[str, SigningProfile] = dict()
        self.tagger = TaggingService()

    def cancel_signing_profile(self, profile_name: str) -> None:
        self.signing_profiles[profile_name].cancel()

    def get_signing_profile(self, profile_name: str) -> SigningProfile:
        return self.signing_profiles[profile_name]

    def put_signing_profile(
        self,
        profile_name: str,
        signature_validity_period: Optional[Dict[str, Any]],
        platform_id: str,
        signing_material: Dict[str, str],
        tags: Dict[str, str],
    ) -> SigningProfile:
        """
        The following parameters are not yet implemented: Overrides, SigningParameters
        """
        profile = SigningProfile(
            backend=self,
            name=profile_name,
            platform_id=platform_id,
            signing_material=signing_material,
            signature_validity_period=signature_validity_period,
        )
        self.signing_profiles[profile_name] = profile
        self.tag_resource(profile.arn, tags)
        return profile

    def list_signing_platforms(self) -> List[Dict[str, Any]]:
        """
        Pagination is not yet implemented. The parameters category, partner, target are not yet implemented
        """
        return SignerBackend.platforms

    def list_tags_for_resource(self, resource_arn: str) -> Dict[str, str]:
        return self.tagger.get_tag_dict_for_resource(resource_arn)

    def tag_resource(self, resource_arn: str, tags: Dict[str, str]) -> None:
        self.tagger.tag_resource(
            resource_arn, TaggingService.convert_dict_to_tags_input(tags)
        )

    def untag_resource(self, resource_arn: str, tag_keys: List[str]) -> None:
        self.tagger.untag_resource_using_names(resource_arn, tag_keys)


signer_backends = BackendDict(SignerBackend, "signer")
