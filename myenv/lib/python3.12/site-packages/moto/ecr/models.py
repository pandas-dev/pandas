import hashlib
import json
import re
import threading
from collections import namedtuple
from datetime import datetime, timezone
from typing import Any, Dict, Iterable, List, Optional, Tuple

from botocore.exceptions import ParamValidationError

from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.common_models import BaseModel, CloudFormationModel
from moto.core.utils import iso_8601_datetime_without_milliseconds, utcnow
from moto.ecr.exceptions import (
    ImageAlreadyExistsException,
    ImageNotFoundException,
    InvalidParameterException,
    LifecyclePolicyNotFoundException,
    LimitExceededException,
    RegistryPolicyNotFoundException,
    RepositoryAlreadyExistsException,
    RepositoryNotEmptyException,
    RepositoryNotFoundException,
    RepositoryPolicyNotFoundException,
    ScanNotFoundException,
    ValidationException,
)
from moto.ecr.policy_validation import EcrLifecyclePolicyValidator
from moto.iam.exceptions import MalformedPolicyDocument
from moto.iam.policy_validation import IAMPolicyDocumentValidator
from moto.moto_api._internal import mock_random as random
from moto.utilities.tagging_service import TaggingService
from moto.utilities.utils import get_partition

ECR_REPOSITORY_ARN_PATTERN = "^arn:(?P<partition>[^:]+):ecr:(?P<region>[^:]+):(?P<account_id>[^:]+):repository/(?P<repo_name>.*)$"
ECR_REPOSITORY_NAME_PATTERN = (
    "(?:[a-z0-9]+(?:[._-][a-z0-9]+)*/)*[a-z0-9]+(?:[._-][a-z0-9]+)*"
)

EcrRepositoryArn = namedtuple(
    "EcrRepositoryArn", ["partition", "region", "account_id", "repo_name"]
)


class BaseObject(BaseModel):
    def camelCase(self, key: str) -> str:
        words = []
        for i, word in enumerate(key.split("_")):
            if i > 0:
                words.append(word.title())
            else:
                words.append(word)
        return "".join(words)

    def gen_response_object(self) -> Dict[str, Any]:
        response_object = dict()
        for key, value in self.__dict__.items():
            if "_" in key:
                response_object[self.camelCase(key)] = value
            else:
                response_object[key] = value
        return response_object

    @property
    def response_object(self) -> Dict[str, Any]:  # type: ignore[misc]
        return self.gen_response_object()


class Repository(BaseObject, CloudFormationModel):
    def __init__(
        self,
        account_id: str,
        region_name: str,
        repository_name: str,
        registry_id: Optional[str],
        encryption_config: Optional[Dict[str, str]],
        image_scan_config: str,
        image_tag_mutablility: str,
    ):
        self.account_id = account_id
        self.region_name = region_name
        self.registry_id = registry_id or account_id
        self.arn = f"arn:{get_partition(region_name)}:ecr:{region_name}:{self.registry_id}:repository/{repository_name}"
        self.name = repository_name
        self.created_at = utcnow()
        self.uri = (
            f"{self.registry_id}.dkr.ecr.{region_name}.amazonaws.com/{repository_name}"
        )
        self.image_tag_mutability = image_tag_mutablility or "MUTABLE"
        self.image_scanning_configuration = image_scan_config or {"scanOnPush": False}
        self.encryption_configuration = self._determine_encryption_config(
            encryption_config
        )
        self.policy: Optional[str] = None
        self.lifecycle_policy: Optional[str] = None
        self.images: List[Image] = []
        self.scanning_config = {
            "repositoryArn": self.arn,
            "repositoryName": self.name,
            "scanOnPush": False,
            "scanFrequency": "MANUAL",
            "appliedScanFilters": [],
        }

    def _determine_encryption_config(
        self, encryption_config: Optional[Dict[str, str]]
    ) -> Dict[str, str]:
        if not encryption_config:
            return {"encryptionType": "AES256"}
        if encryption_config == {"encryptionType": "KMS"}:
            encryption_config["kmsKey"] = (
                f"arn:{get_partition(self.region_name)}:kms:{self.region_name}:{self.account_id}:key/{random.uuid4()}"
            )
        return encryption_config

    def _get_image(
        self, image_tag: Optional[str], image_digest: Optional[str]
    ) -> "Image":
        # you can either search for one or both
        image = next(
            (
                i
                for i in self.images
                if (not image_tag or image_tag in i.image_tags)
                and (not image_digest or image_digest == i.get_image_digest())
            ),
            None,
        )

        if not image:
            idigest = image_digest or "null"
            itag = image_tag or "null"
            image_id_rep = f"{{imageDigest:'{idigest}', imageTag:'{itag}'}}"

            raise ImageNotFoundException(
                image_id=image_id_rep,
                repository_name=self.name,
                registry_id=self.registry_id,
            )

        return image

    @property
    def physical_resource_id(self) -> str:
        return self.name

    @property
    def response_object(self) -> Dict[str, Any]:  # type: ignore[misc]
        response_object = self.gen_response_object()

        response_object["registryId"] = self.registry_id
        response_object["repositoryArn"] = self.arn
        response_object["repositoryName"] = self.name
        response_object["repositoryUri"] = self.uri
        response_object["createdAt"] = iso_8601_datetime_without_milliseconds(
            self.created_at
        )
        del response_object["arn"], response_object["name"], response_object["images"]
        return response_object

    def update(
        self,
        image_scan_config: Optional[Dict[str, Any]] = None,
        image_tag_mutability: Optional[str] = None,
    ) -> None:
        if image_scan_config:
            self.image_scanning_configuration = image_scan_config
        if image_tag_mutability:
            self.image_tag_mutability = image_tag_mutability

    def delete(self, account_id: str, region_name: str) -> None:
        ecr_backend = ecr_backends[account_id][region_name]
        ecr_backend.delete_repository(self.name)

    @classmethod
    def has_cfn_attr(cls, attr: str) -> bool:
        return attr in ["Arn", "RepositoryUri"]

    def get_cfn_attribute(self, attribute_name: str) -> str:
        from moto.cloudformation.exceptions import UnformattedGetAttTemplateException

        if attribute_name == "Arn":
            return self.arn
        elif attribute_name == "RepositoryUri":
            return self.uri

        raise UnformattedGetAttTemplateException()

    @staticmethod
    def cloudformation_name_type() -> str:
        return "RepositoryName"

    @staticmethod
    def cloudformation_type() -> str:
        # https://docs.aws.amazon.com/AWSCloudFormation/latest/UserGuide/aws-resource-ecr-repository.html
        return "AWS::ECR::Repository"

    @classmethod
    def create_from_cloudformation_json(  # type: ignore[misc]
        cls,
        resource_name: str,
        cloudformation_json: Any,
        account_id: str,
        region_name: str,
        **kwargs: Any,
    ) -> "Repository":
        ecr_backend = ecr_backends[account_id][region_name]
        properties = cloudformation_json["Properties"]

        encryption_config = properties.get("EncryptionConfiguration")
        image_scan_config = properties.get("ImageScanningConfiguration")
        image_tag_mutablility = properties.get("ImageTagMutability")
        tags = properties.get("Tags", [])

        return ecr_backend.create_repository(
            # RepositoryName is optional in CloudFormation, thus create a random
            # name if necessary
            repository_name=resource_name,
            registry_id=None,
            encryption_config=encryption_config,
            image_scan_config=image_scan_config,
            image_tag_mutablility=image_tag_mutablility,
            tags=tags,
        )

    @classmethod
    def update_from_cloudformation_json(  # type: ignore[misc]
        cls,
        original_resource: Any,
        new_resource_name: str,
        cloudformation_json: Any,
        account_id: str,
        region_name: str,
    ) -> "Repository":
        ecr_backend = ecr_backends[account_id][region_name]
        properties = cloudformation_json["Properties"]
        encryption_configuration = properties.get(
            "EncryptionConfiguration", {"encryptionType": "AES256"}
        )

        if (
            new_resource_name == original_resource.name
            and encryption_configuration == original_resource.encryption_configuration
        ):
            original_resource.update(
                properties.get("ImageScanningConfiguration"),
                properties.get("ImageTagMutability"),
            )

            ecr_backend.tagger.tag_resource(
                original_resource.arn, properties.get("Tags", [])
            )

            return original_resource
        else:
            original_resource.delete(account_id, region_name)
            return cls.create_from_cloudformation_json(
                new_resource_name, cloudformation_json, account_id, region_name
            )


class Image(BaseObject):
    def __init__(
        self,
        account_id: str,
        tag: str,
        manifest: str,
        repository: str,
        image_manifest_mediatype: Optional[str] = None,
        digest: Optional[str] = None,
        registry_id: Optional[str] = None,
    ):
        self.image_tag = tag
        self.image_tags = [tag] if tag is not None else []
        self.image_manifest = manifest
        self.image_manifest_mediatype = image_manifest_mediatype
        self.repository = repository
        self.registry_id = registry_id or account_id
        self.image_digest = digest
        self.image_pushed_at = int(datetime.now(timezone.utc).timestamp())
        self.last_scan: Optional[datetime] = None

    def _create_digest(self) -> None:
        image_manifest = json.loads(self.image_manifest)
        if "layers" in image_manifest:
            layer_digests = [layer["digest"] for layer in image_manifest["layers"]]
            self.image_digest = (
                "sha256:"
                + hashlib.sha256("".join(layer_digests).encode("utf-8")).hexdigest()
            )
        else:
            random_sha = hashlib.sha256(
                f"{random.randint(0,100)}".encode("utf-8")
            ).hexdigest()
            self.image_digest = f"sha256:{random_sha}"

    def get_image_digest(self) -> str:
        if not self.image_digest:
            self._create_digest()
        return self.image_digest  # type: ignore[return-value]

    def get_image_size_in_bytes(self) -> Optional[int]:
        image_manifest = json.loads(self.image_manifest)
        if "layers" in image_manifest:
            try:
                return image_manifest["config"]["size"]
            except KeyError:
                return 50 * 1024 * 1024
        else:
            return None

    def get_image_manifest(self) -> str:
        return self.image_manifest

    def remove_tag(self, tag: str) -> None:
        if tag is not None and tag in self.image_tags:
            self.image_tags.remove(tag)
            if self.image_tags:
                self.image_tag = self.image_tags[-1]

    def update_tag(self, tag: str) -> None:
        self.image_tag = tag
        if tag not in self.image_tags and tag is not None:
            self.image_tags.append(tag)

    @property
    def response_object(self) -> Dict[str, Any]:  # type: ignore[misc]
        response_object = self.gen_response_object()
        response_object["imageId"] = {}
        response_object["imageId"]["imageTag"] = self.image_tag
        response_object["imageId"]["imageDigest"] = self.get_image_digest()
        response_object["imageManifest"] = self.image_manifest
        response_object["imageManifestMediaType"] = self.image_manifest_mediatype
        response_object["repositoryName"] = self.repository
        response_object["registryId"] = self.registry_id
        return {
            k: v for k, v in response_object.items() if v is not None and v != [None]
        }

    @property
    def response_list_object(self) -> Dict[str, Any]:  # type: ignore[misc]
        response_object = self.gen_response_object()
        response_object["imageTag"] = self.image_tag
        response_object["imageDigest"] = self.get_image_digest()
        return {
            k: v for k, v in response_object.items() if v is not None and v != [None]
        }

    @property
    def response_describe_object(self) -> Dict[str, Any]:  # type: ignore[misc]
        response_object = self.gen_response_object()
        response_object["imageTags"] = self.image_tags
        response_object["imageDigest"] = self.get_image_digest()
        response_object["imageManifest"] = self.image_manifest
        response_object["imageManifestMediaType"] = self.image_manifest_mediatype
        response_object["repositoryName"] = self.repository
        response_object["registryId"] = self.registry_id
        response_object["imageSizeInBytes"] = self.get_image_size_in_bytes()
        response_object["imagePushedAt"] = self.image_pushed_at
        return {k: v for k, v in response_object.items() if v is not None and v != []}

    @property
    def response_batch_get_image(self) -> Dict[str, Any]:  # type: ignore[misc]
        response_object = {
            "imageId": {
                "imageTag": self.image_tag,
                "imageDigest": self.get_image_digest(),
            },
            "imageManifest": self.image_manifest,
            "repositoryName": self.repository,
            "registryId": self.registry_id,
        }
        return {
            k: v
            for k, v in response_object.items()
            if v is not None and v != [None]  # type: ignore
        }

    @property
    def response_batch_delete_image(self) -> Dict[str, Any]:  # type: ignore[misc]
        response_object = {}
        response_object["imageDigest"] = self.get_image_digest()
        response_object["imageTag"] = self.image_tag
        return {
            k: v
            for k, v in response_object.items()
            if v is not None and v != [None]  # type: ignore
        }


class ECRBackend(BaseBackend):
    def __init__(self, region_name: str, account_id: str):
        super().__init__(region_name, account_id)
        self.registry_policy: Optional[str] = None
        self.replication_config: Dict[str, Any] = {"rules": []}
        self.repositories: Dict[str, Repository] = {}
        self.registry_scanning_configuration: Dict[str, Any] = {
            "scanType": "BASIC",
            "rules": [],
        }
        self.registry_scanning_configuration_update_lock = threading.RLock()
        self.tagger = TaggingService(tag_name="tags")

    @staticmethod
    def default_vpc_endpoint_service(  # type: ignore[misc]
        service_region: str, zones: List[str]
    ) -> List[Dict[str, Any]]:
        """Default VPC endpoint service."""
        docker_endpoint = {
            "AcceptanceRequired": False,
            "AvailabilityZones": zones,
            "BaseEndpointDnsNames": [f"dkr.ecr.{service_region}.vpce.amazonaws.com"],
            "ManagesVpcEndpoints": False,
            "Owner": "amazon",
            "PrivateDnsName": f"*.dkr.ecr.{service_region}.amazonaws.com",
            "PrivateDnsNameVerificationState": "verified",
            "PrivateDnsNames": [
                {"PrivateDnsName": f"*.dkr.ecr.{service_region}.amazonaws.com"}
            ],
            "ServiceId": f"vpce-svc-{BaseBackend.vpce_random_number()}",
            "ServiceName": f"com.amazonaws.{service_region}.ecr.dkr",
            "ServiceType": [{"ServiceType": "Interface"}],
            "Tags": [],
            "VpcEndpointPolicySupported": True,
        }
        return BaseBackend.default_vpc_endpoint_service_factory(
            service_region, zones, "api.ecr", special_service_name="ecr.api"
        ) + [docker_endpoint]

    def _get_repository(
        self, name: str, registry_id: Optional[str] = None
    ) -> Repository:
        repo = self.repositories.get(name)
        reg_id = registry_id or self.account_id

        if not repo or repo.registry_id != reg_id:
            raise RepositoryNotFoundException(name, reg_id)
        return repo

    @staticmethod
    def _parse_resource_arn(resource_arn: str) -> EcrRepositoryArn:  # type: ignore[misc]
        match = re.match(ECR_REPOSITORY_ARN_PATTERN, resource_arn)
        if not match:
            raise InvalidParameterException(
                "Invalid parameter at 'resourceArn' failed to satisfy constraint: "
                "'Invalid ARN'"
            )
        return EcrRepositoryArn(**match.groupdict())

    def describe_repositories(
        self,
        registry_id: Optional[str] = None,
        repository_names: Optional[List[str]] = None,
    ) -> List[Dict[str, Any]]:
        """
        maxResults and nextToken not implemented
        """
        if repository_names:
            for repository_name in repository_names:
                if repository_name not in self.repositories:
                    raise RepositoryNotFoundException(
                        repository_name, registry_id or self.account_id
                    )

        repositories = []
        for repository in self.repositories.values():
            # If a registry_id was supplied, ensure this repository matches
            if registry_id:
                if repository.registry_id != registry_id:
                    continue
            # If a list of repository names was supplied, esure this repository
            # is in that list
            if repository_names:
                if repository.name not in repository_names:
                    continue
            repositories.append(repository.response_object)
        return repositories

    def create_repository(
        self,
        repository_name: str,
        registry_id: Optional[str],
        encryption_config: Dict[str, str],
        image_scan_config: Any,
        image_tag_mutablility: str,
        tags: List[Dict[str, str]],
    ) -> Repository:
        if self.repositories.get(repository_name):
            raise RepositoryAlreadyExistsException(repository_name, self.account_id)

        match = re.fullmatch(ECR_REPOSITORY_NAME_PATTERN, repository_name)
        if not match:
            raise InvalidParameterException(
                f"Invalid parameter at 'repositoryName' failed to satisfy constraint: 'must satisfy regular expression '{ECR_REPOSITORY_NAME_PATTERN}'"
            )

        repository = Repository(
            account_id=self.account_id,
            region_name=self.region_name,
            repository_name=repository_name,
            registry_id=registry_id,
            encryption_config=encryption_config,
            image_scan_config=image_scan_config,
            image_tag_mutablility=image_tag_mutablility,
        )
        self.repositories[repository_name] = repository
        self.tagger.tag_resource(repository.arn, tags)

        # check if any of the registry scanning policies applies to the repository
        with self.registry_scanning_configuration_update_lock:
            for rule in self.registry_scanning_configuration["rules"]:
                for repo_filter in rule["repositoryFilters"]:
                    if self._match_repository_filter(
                        repo_filter["filter"], repository_name
                    ):
                        repository.scanning_config["scanFrequency"] = rule[
                            "scanFrequency"
                        ]
                        # AWS testing seems to indicate that this is always overwritten
                        repository.scanning_config["appliedScanFilters"] = [repo_filter]

        return repository

    def delete_repository(
        self,
        repository_name: str,
        registry_id: Optional[str] = None,
        force: bool = False,
    ) -> Repository:
        repo = self._get_repository(repository_name, registry_id)

        if repo.images and not force:
            raise RepositoryNotEmptyException(
                repository_name, registry_id or self.account_id
            )

        self.tagger.delete_all_tags_for_resource(repo.arn)
        return self.repositories.pop(repository_name)

    def list_images(
        self, repository_name: str, registry_id: Optional[str] = None
    ) -> List[Image]:
        """
        maxResults and filtering not implemented
        """
        repository = None
        found = False
        if repository_name in self.repositories:
            repository = self.repositories[repository_name]
            if registry_id:
                if repository.registry_id == registry_id:
                    found = True
            else:
                found = True

        if not found:
            raise RepositoryNotFoundException(
                repository_name, registry_id or self.account_id
            )

        return list(repository.images)  # type: ignore[union-attr]

    def describe_images(
        self,
        repository_name: str,
        registry_id: Optional[str] = None,
        image_ids: Optional[List[Dict[str, str]]] = None,
    ) -> Iterable[Image]:
        repository = self._get_repository(repository_name, registry_id)

        if image_ids:
            return set(
                repository._get_image(
                    image_id.get("imageTag"), image_id.get("imageDigest")
                )
                for image_id in image_ids
            )

        else:
            return list(repository.images)

    def put_image(
        self,
        repository_name: str,
        image_manifest: str,
        image_tag: str,
        image_manifest_mediatype: Optional[str] = None,
        digest: Optional[str] = None,
    ) -> Image:
        if repository_name in self.repositories:
            repository = self.repositories[repository_name]
        else:
            raise Exception(f"{repository_name} is not a repository")

        try:
            parsed_image_manifest = json.loads(image_manifest)
        except json.JSONDecodeError:
            raise Exception(
                "Invalid parameter at 'ImageManifest' failed to satisfy constraint: 'Invalid JSON syntax'"
            )

        if image_manifest_mediatype:
            parsed_image_manifest["imageManifest"] = image_manifest_mediatype
        else:
            if "mediaType" not in parsed_image_manifest:
                raise InvalidParameterException(
                    message="image manifest mediatype not provided in manifest or parameter"
                )
            else:
                image_manifest_mediatype = parsed_image_manifest["mediaType"]

        existing_images_with_matching_manifest = list(
            filter(
                lambda x: x.response_object["imageManifest"] == image_manifest,
                repository.images,
            )
        )

        # if an image with a matching manifest exists and it is tagged,
        # trying to put the same image with the same tag will result in an
        # ImageAlreadyExistsException

        try:
            existing_images_with_matching_tag = list(
                filter(
                    lambda x: image_tag in x.image_tags,
                    repository.images,
                )
            )
        except KeyError:
            existing_images_with_matching_tag = []

        if not existing_images_with_matching_manifest:
            # this image is not in ECR yet
            image = Image(
                self.account_id,
                image_tag,
                image_manifest,
                repository_name,
                image_manifest_mediatype,
                digest,
            )
            repository.images.append(image)
            if existing_images_with_matching_tag:
                # Tags are unique, so delete any existing image with this tag first
                # (or remove the tag if the image has more than one tag)
                self.batch_delete_image(
                    repository_name=repository_name, image_ids=[{"imageTag": image_tag}]
                )
            return image
        else:
            # this image is in ECR
            image = existing_images_with_matching_manifest[0]
            if image.image_tag == image_tag:
                raise ImageAlreadyExistsException(
                    registry_id=repository.registry_id,
                    image_tag=image_tag,
                    digest=image.get_image_digest(),
                    repository_name=repository_name,
                )
            else:
                # Tags are unique, so delete any existing image with this tag first
                # (or remove the tag if the image has more than one tag)
                self.batch_delete_image(
                    repository_name=repository_name, image_ids=[{"imageTag": image_tag}]
                )
                # update existing image
                image.update_tag(image_tag)
                return image

    def batch_get_image(
        self,
        repository_name: str,
        registry_id: Optional[str] = None,
        image_ids: Optional[List[Dict[str, Any]]] = None,
    ) -> Dict[str, Any]:
        """
        The parameter AcceptedMediaTypes has not yet been implemented
        """
        if repository_name in self.repositories:
            repository = self.repositories[repository_name]
        else:
            raise RepositoryNotFoundException(
                repository_name, registry_id or self.account_id
            )

        if not image_ids:
            raise ParamValidationError(
                msg='Missing required parameter in input: "imageIds"'
            )

        response: Dict[str, Any] = {"images": [], "failures": []}

        for image_id in image_ids:
            found = False
            for image in repository.images:
                if (
                    "imageDigest" in image_id
                    and image.get_image_digest() == image_id["imageDigest"]
                ) or (
                    "imageTag" in image_id and image_id["imageTag"] in image.image_tags
                ):
                    found = True
                    response["images"].append(image.response_batch_get_image)

            if not found:
                response["failures"].append(
                    {
                        "imageId": {"imageTag": image_id.get("imageTag", "null")},
                        "failureCode": "ImageNotFound",
                        "failureReason": "Requested image not found",
                    }
                )

        return response

    def batch_delete_image(
        self,
        repository_name: str,
        registry_id: Optional[str] = None,
        image_ids: Optional[List[Dict[str, str]]] = None,
    ) -> Dict[str, Any]:
        if repository_name in self.repositories:
            repository = self.repositories[repository_name]
        else:
            raise RepositoryNotFoundException(
                repository_name, registry_id or self.account_id
            )

        if not image_ids:
            raise ParamValidationError(
                msg='Missing required parameter in input: "imageIds"'
            )

        response: Dict[str, Any] = {"imageIds": [], "failures": []}

        for image_id in image_ids:
            image_found = False

            # Is request missing both digest and tag?
            if "imageDigest" not in image_id and "imageTag" not in image_id:
                response["failures"].append(
                    {
                        "imageId": {},
                        "failureCode": "MissingDigestAndTag",
                        "failureReason": "Invalid request parameters: both tag and digest cannot be null",
                    }
                )
                continue

            # If we have a digest, is it valid?
            if "imageDigest" in image_id:
                pattern = re.compile(r"^[0-9a-zA-Z_+\.-]+:[0-9a-fA-F]{64}")
                if not pattern.match(image_id["imageDigest"]):
                    response["failures"].append(
                        {
                            "imageId": {"imageDigest": image_id["imageDigest"]},
                            "failureCode": "InvalidImageDigest",
                            "failureReason": "Invalid request parameters: image digest should satisfy the regex '[a-zA-Z0-9-_+.]+:[a-fA-F0-9]+'",
                        }
                    )
                    continue

            for num, image in enumerate(repository.images):
                # Search by matching both digest and tag
                if "imageDigest" in image_id and "imageTag" in image_id:
                    if (
                        image_id["imageDigest"] == image.get_image_digest()
                        and image_id["imageTag"] in image.image_tags
                    ):
                        image_found = True
                        for image_tag in reversed(image.image_tags):
                            repository.images[num].image_tag = image_tag
                            response["imageIds"].append(
                                image.response_batch_delete_image
                            )
                            repository.images[num].remove_tag(image_tag)
                        del repository.images[num]

                # Search by matching digest
                elif (
                    "imageDigest" in image_id
                    and image.get_image_digest() == image_id["imageDigest"]
                ):
                    image_found = True
                    for image_tag in reversed(image.image_tags):
                        repository.images[num].image_tag = image_tag
                        response["imageIds"].append(image.response_batch_delete_image)
                        repository.images[num].remove_tag(image_tag)
                    del repository.images[num]

                # Search by matching tag
                elif (
                    "imageTag" in image_id and image_id["imageTag"] in image.image_tags
                ):
                    image_found = True
                    repository.images[num].image_tag = image_id["imageTag"]
                    response["imageIds"].append(image.response_batch_delete_image)
                    if len(image.image_tags) > 1:
                        repository.images[num].remove_tag(image_id["imageTag"])
                    else:
                        repository.images.remove(image)

            if not image_found:
                failure_response: Dict[str, Any] = {
                    "imageId": {},
                    "failureCode": "ImageNotFound",
                    "failureReason": "Requested image not found",
                }

                if "imageDigest" in image_id:
                    failure_response["imageId"]["imageDigest"] = image_id.get(
                        "imageDigest", "null"
                    )

                if "imageTag" in image_id:
                    failure_response["imageId"]["imageTag"] = image_id.get(
                        "imageTag", "null"
                    )

                response["failures"].append(failure_response)

        return response

    def batch_get_repository_scanning_configuration(
        self, names: List[str]
    ) -> Tuple[List[Dict[str, Any]], List[str]]:
        configs = []
        failing = []
        for name in names:
            try:
                configs.append(
                    self._get_repository(name=name, registry_id=None).scanning_config
                )
            except RepositoryNotFoundException:
                failing.append(name)
        return configs, failing

    def list_tags_for_resource(self, arn: str) -> Dict[str, List[Dict[str, str]]]:
        resource = self._parse_resource_arn(arn)
        repo = self._get_repository(resource.repo_name, resource.account_id)

        return self.tagger.list_tags_for_resource(repo.arn)

    def tag_resource(self, arn: str, tags: List[Dict[str, str]]) -> None:
        resource = self._parse_resource_arn(arn)
        repo = self._get_repository(resource.repo_name, resource.account_id)
        self.tagger.tag_resource(repo.arn, tags)

    def untag_resource(self, arn: str, tag_keys: List[str]) -> None:
        resource = self._parse_resource_arn(arn)
        repo = self._get_repository(resource.repo_name, resource.account_id)
        self.tagger.untag_resource_using_names(repo.arn, tag_keys)

    def put_image_tag_mutability(
        self, registry_id: str, repository_name: str, image_tag_mutability: str
    ) -> Dict[str, str]:
        if image_tag_mutability not in ["IMMUTABLE", "MUTABLE"]:
            raise InvalidParameterException(
                "Invalid parameter at 'imageTagMutability' failed to satisfy constraint: "
                "'Member must satisfy enum value set: [IMMUTABLE, MUTABLE]'"
            )

        repo = self._get_repository(repository_name, registry_id)
        repo.update(image_tag_mutability=image_tag_mutability)

        return {
            "registryId": repo.registry_id,
            "repositoryName": repository_name,
            "imageTagMutability": repo.image_tag_mutability,
        }

    def put_image_scanning_configuration(
        self, registry_id: str, repository_name: str, image_scan_config: Dict[str, Any]
    ) -> Dict[str, Any]:
        repo = self._get_repository(repository_name, registry_id)
        repo.update(image_scan_config=image_scan_config)

        return {
            "registryId": repo.registry_id,
            "repositoryName": repository_name,
            "imageScanningConfiguration": repo.image_scanning_configuration,
        }

    def set_repository_policy(
        self, registry_id: str, repository_name: str, policy_text: str
    ) -> Dict[str, Any]:
        repo = self._get_repository(repository_name, registry_id)

        try:
            iam_policy_document_validator = IAMPolicyDocumentValidator(policy_text)
            # the repository policy can be defined without a resource field
            iam_policy_document_validator._validate_resource_exist = lambda: None  # type: ignore
            # the repository policy can have the old version 2008-10-17
            iam_policy_document_validator._validate_version = lambda: None  # type: ignore
            iam_policy_document_validator.validate()
        except MalformedPolicyDocument:
            raise InvalidParameterException(
                "Invalid parameter at 'PolicyText' failed to satisfy constraint: "
                "'Invalid repository policy provided'"
            )

        repo.policy = policy_text

        return {
            "registryId": repo.registry_id,
            "repositoryName": repository_name,
            "policyText": repo.policy,
        }

    def get_repository_policy(
        self, registry_id: str, repository_name: str
    ) -> Dict[str, Any]:
        repo = self._get_repository(repository_name, registry_id)

        if not repo.policy:
            raise RepositoryPolicyNotFoundException(repository_name, repo.registry_id)

        return {
            "registryId": repo.registry_id,
            "repositoryName": repository_name,
            "policyText": repo.policy,
        }

    def delete_repository_policy(
        self, registry_id: str, repository_name: str
    ) -> Dict[str, Any]:
        repo = self._get_repository(repository_name, registry_id)
        policy = repo.policy

        if not policy:
            raise RepositoryPolicyNotFoundException(repository_name, repo.registry_id)

        repo.policy = None

        return {
            "registryId": repo.registry_id,
            "repositoryName": repository_name,
            "policyText": policy,
        }

    def put_lifecycle_policy(
        self, registry_id: str, repository_name: str, lifecycle_policy_text: str
    ) -> Dict[str, Any]:
        repo = self._get_repository(repository_name, registry_id)

        validator = EcrLifecyclePolicyValidator(lifecycle_policy_text)
        validator.validate()

        repo.lifecycle_policy = lifecycle_policy_text

        return {
            "registryId": repo.registry_id,
            "repositoryName": repository_name,
            "lifecyclePolicyText": repo.lifecycle_policy,
        }

    def get_lifecycle_policy(
        self, registry_id: str, repository_name: str
    ) -> Dict[str, Any]:
        repo = self._get_repository(repository_name, registry_id)

        if not repo.lifecycle_policy:
            raise LifecyclePolicyNotFoundException(repository_name, repo.registry_id)

        return {
            "registryId": repo.registry_id,
            "repositoryName": repository_name,
            "lifecyclePolicyText": repo.lifecycle_policy,
            "lastEvaluatedAt": iso_8601_datetime_without_milliseconds(utcnow()),
        }

    def delete_lifecycle_policy(
        self, registry_id: str, repository_name: str
    ) -> Dict[str, Any]:
        repo = self._get_repository(repository_name, registry_id)
        policy = repo.lifecycle_policy

        if not policy:
            raise LifecyclePolicyNotFoundException(repository_name, repo.registry_id)

        repo.lifecycle_policy = None

        return {
            "registryId": repo.registry_id,
            "repositoryName": repository_name,
            "lifecyclePolicyText": policy,
            "lastEvaluatedAt": iso_8601_datetime_without_milliseconds(utcnow()),
        }

    def _validate_registry_policy_action(self, policy_text: str) -> None:
        # only CreateRepository & ReplicateImage actions are allowed
        VALID_ACTIONS = {"ecr:CreateRepository", "ecr:ReplicateImage"}

        policy = json.loads(policy_text)
        for statement in policy["Statement"]:
            action = statement["Action"]
            if isinstance(action, str):
                action = [action]
            if set(action) - VALID_ACTIONS:
                raise MalformedPolicyDocument()

    def put_registry_policy(self, policy_text: str) -> Dict[str, Any]:
        try:
            iam_policy_document_validator = IAMPolicyDocumentValidator(policy_text)
            iam_policy_document_validator.validate()

            self._validate_registry_policy_action(policy_text)
        except MalformedPolicyDocument:
            raise InvalidParameterException(
                "Invalid parameter at 'PolicyText' failed to satisfy constraint: "
                "'Invalid registry policy provided'"
            )

        self.registry_policy = policy_text

        return {
            "registryId": self.account_id,
            "policyText": policy_text,
        }

    def get_registry_policy(self) -> Dict[str, Any]:
        if not self.registry_policy:
            raise RegistryPolicyNotFoundException(self.account_id)

        return {
            "registryId": self.account_id,
            "policyText": self.registry_policy,
        }

    def delete_registry_policy(self) -> Dict[str, Any]:
        policy = self.registry_policy
        if not policy:
            raise RegistryPolicyNotFoundException(self.account_id)

        self.registry_policy = None

        return {
            "registryId": self.account_id,
            "policyText": policy,
        }

    def start_image_scan(
        self, registry_id: str, repository_name: str, image_id: Dict[str, str]
    ) -> Dict[str, Any]:
        repo = self._get_repository(repository_name, registry_id)

        image = repo._get_image(image_id.get("imageTag"), image_id.get("imageDigest"))

        # scanning an image is only allowed once per day
        if image.last_scan and image.last_scan.date() == datetime.today().date():
            raise LimitExceededException()

        image.last_scan = datetime.today()

        return {
            "registryId": repo.registry_id,
            "repositoryName": repository_name,
            "imageId": {
                "imageDigest": image.image_digest,
                "imageTag": image.image_tag,
            },
            "imageScanStatus": {"status": "IN_PROGRESS"},
        }

    def describe_image_scan_findings(
        self, registry_id: str, repository_name: str, image_id: Dict[str, Any]
    ) -> Dict[str, Any]:
        repo = self._get_repository(repository_name, registry_id)

        image = repo._get_image(image_id.get("imageTag"), image_id.get("imageDigest"))

        if not image.last_scan:
            idigest = image_id.get("imageDigest") or "null"
            itag = image_id.get("imageTag") or "null"
            image_id_rep = f"{{imageDigest:'{idigest}', imageTag:'{itag}'}}"
            raise ScanNotFoundException(
                image_id=image_id_rep,
                repository_name=repository_name,
                registry_id=repo.registry_id,
            )

        return {
            "registryId": repo.registry_id,
            "repositoryName": repository_name,
            "imageId": {
                "imageDigest": image.image_digest,
                "imageTag": image.image_tag,
            },
            "imageScanStatus": {
                "status": "COMPLETE",
                "description": "The scan was completed successfully.",
            },
            "imageScanFindings": {
                "imageScanCompletedAt": iso_8601_datetime_without_milliseconds(
                    image.last_scan
                ),
                "vulnerabilitySourceUpdatedAt": iso_8601_datetime_without_milliseconds(
                    utcnow()
                ),
                "findings": [
                    {
                        "name": "CVE-9999-9999",
                        "uri": "https://cve.mitre.org/cgi-bin/cvename.cgi?name=CVE-9999-9999",
                        "severity": "HIGH",
                        "attributes": [
                            {"key": "package_version", "value": "9.9.9"},
                            {"key": "package_name", "value": "moto_fake"},
                            {
                                "key": "CVSS2_VECTOR",
                                "value": "AV:N/AC:L/Au:N/C:P/I:P/A:P",
                            },
                            {"key": "CVSS2_SCORE", "value": "7.5"},
                        ],
                    }
                ],
                "findingSeverityCounts": {"HIGH": 1},
            },
        }

    def put_replication_configuration(
        self, replication_config: Dict[str, Any]
    ) -> Dict[str, Any]:
        rules = replication_config["rules"]
        if len(rules) > 1:
            raise ValidationException("This feature is disabled")

        if len(rules) == 1:
            for dest in rules[0]["destinations"]:
                if (
                    dest["region"] == self.region_name
                    and dest["registryId"] == self.account_id
                ):
                    raise InvalidParameterException(
                        "Invalid parameter at 'replicationConfiguration' failed to satisfy constraint: "
                        "'Replication destination cannot be the same as the source registry'"
                    )

        self.replication_config = replication_config

        return {"replicationConfiguration": replication_config}

    def _match_repository_filter(self, filter: str, repository_name: str) -> bool:
        filter_regex = filter.replace("*", ".*")
        return filter in repository_name or bool(
            re.match(filter_regex, repository_name)
        )

    def get_registry_scanning_configuration(self) -> Dict[str, Any]:
        return self.registry_scanning_configuration

    def put_registry_scanning_configuration(
        self, scan_type: str, rules: List[Dict[str, Any]]
    ) -> None:
        # locking here to avoid simultaneous updates which leads to inconsistent state
        with self.registry_scanning_configuration_update_lock:
            self.registry_scanning_configuration = {
                "scanType": scan_type,
                "rules": rules,
            }

            # reset all rules first
            for repo in self.repositories.values():
                repo.scanning_config["scanFrequency"] = "MANUAL"
                repo.scanning_config["appliedScanFilters"] = []

            for rule in rules:
                for repo_filter in rule["repositoryFilters"]:
                    for repo in self.repositories.values():
                        if self._match_repository_filter(
                            repo_filter["filter"], repo.name
                        ):
                            repo.scanning_config["scanFrequency"] = rule[
                                "scanFrequency"
                            ]
                            # AWS testing seems to indicate that this is always overwritten
                            repo.scanning_config["appliedScanFilters"] = [repo_filter]

    def describe_registry(self) -> Dict[str, Any]:
        return {
            "registryId": self.account_id,
            "replicationConfiguration": self.replication_config,
        }


ecr_backends = BackendDict(ECRBackend, "ecr")
