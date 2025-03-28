import json
import time
from base64 import b64encode
from datetime import datetime

from moto.core.responses import BaseResponse

from .models import ECRBackend, ecr_backends


class ECRResponse(BaseResponse):
    def __init__(self) -> None:
        super().__init__(service_name="ecr")

    @property
    def ecr_backend(self) -> ECRBackend:
        return ecr_backends[self.current_account][self.region]

    def create_repository(self) -> str:
        repository_name = self._get_param("repositoryName")
        registry_id = self._get_param("registryId")
        encryption_config = self._get_param("encryptionConfiguration")
        image_scan_config = self._get_param("imageScanningConfiguration")
        image_tag_mutablility = self._get_param("imageTagMutability")
        tags = self._get_param("tags", [])

        repository = self.ecr_backend.create_repository(
            repository_name=repository_name,
            registry_id=registry_id,
            encryption_config=encryption_config,
            image_scan_config=image_scan_config,
            image_tag_mutablility=image_tag_mutablility,
            tags=tags,
        )
        return json.dumps({"repository": repository.response_object})

    def describe_repositories(self) -> str:
        describe_repositories_name = self._get_param("repositoryNames")
        registry_id = self._get_param("registryId")

        repositories = self.ecr_backend.describe_repositories(
            repository_names=describe_repositories_name, registry_id=registry_id
        )
        return json.dumps({"repositories": repositories, "failures": []})

    def delete_repository(self) -> str:
        repository_str = self._get_param("repositoryName")
        registry_id = self._get_param("registryId")
        force = self._get_param("force")

        repository = self.ecr_backend.delete_repository(
            repository_str, registry_id, force
        )
        return json.dumps({"repository": repository.response_object})

    def put_image(self) -> str:
        repository_str = self._get_param("repositoryName")
        image_manifest = self._get_param("imageManifest")
        image_tag = self._get_param("imageTag")
        image_manifest_media_type = self._get_param("imageManifestMediaType")
        digest = self._get_param("imageDigest")
        image = self.ecr_backend.put_image(
            repository_str, image_manifest, image_tag, image_manifest_media_type, digest
        )

        return json.dumps({"image": image.response_object})

    def list_images(self) -> str:
        repository_str = self._get_param("repositoryName")
        registry_id = self._get_param("registryId")
        images = self.ecr_backend.list_images(repository_str, registry_id)
        resp = []
        for image in images:
            for tag in image.image_tags or [None]:  # type: ignore
                resp.append({"imageDigest": image.image_digest, "imageTag": tag})
        return json.dumps({"imageIds": resp})

    def describe_images(self) -> str:
        repository_str = self._get_param("repositoryName")
        registry_id = self._get_param("registryId")
        image_ids = self._get_param("imageIds")
        images = self.ecr_backend.describe_images(
            repository_str, registry_id, image_ids
        )
        return json.dumps(
            {"imageDetails": [image.response_describe_object for image in images]}
        )

    def batch_check_layer_availability(self) -> None:
        self.error_on_dryrun()
        raise NotImplementedError(
            "ECR.batch_check_layer_availability is not yet implemented"
        )

    def batch_delete_image(self) -> str:
        repository_str = self._get_param("repositoryName")
        registry_id = self._get_param("registryId")
        image_ids = self._get_param("imageIds")

        response = self.ecr_backend.batch_delete_image(
            repository_str, registry_id, image_ids
        )
        return json.dumps(response)

    def batch_get_image(self) -> str:
        repository_str = self._get_param("repositoryName")
        registry_id = self._get_param("registryId")
        image_ids = self._get_param("imageIds")

        response = self.ecr_backend.batch_get_image(
            repository_str, registry_id, image_ids
        )
        return json.dumps(response)

    def batch_get_repository_scanning_configuration(self) -> str:
        names = self._get_param("repositoryNames")
        configs, missing = self.ecr_backend.batch_get_repository_scanning_configuration(
            names
        )
        return json.dumps(
            {
                "scanningConfigurations": configs,
                "failures": [
                    {
                        "repositoryName": m,
                        "failureCode": "REPOSITORY_NOT_FOUND",
                        "failureReason": "REPOSITORY_NOT_FOUND",
                    }
                    for m in missing
                ],
            }
        )

    def complete_layer_upload(self) -> None:
        self.error_on_dryrun()
        raise NotImplementedError("ECR.complete_layer_upload is not yet implemented")

    def delete_repository_policy(self) -> str:
        registry_id = self._get_param("registryId")
        repository_name = self._get_param("repositoryName")

        return json.dumps(
            self.ecr_backend.delete_repository_policy(
                registry_id=registry_id, repository_name=repository_name
            )
        )

    def get_authorization_token(self) -> str:
        registry_ids = self._get_param("registryIds")
        if not registry_ids:
            registry_ids = [self.current_account]
        auth_data = []
        for registry_id in registry_ids:
            password = f"{registry_id}-auth-token"
            auth_token = b64encode(f"AWS:{password}".encode("ascii")).decode()
            auth_data.append(
                {
                    "authorizationToken": auth_token,
                    "expiresAt": time.mktime(datetime(2015, 1, 1).timetuple()),
                    "proxyEndpoint": f"https://{registry_id}.dkr.ecr.{self.region}.amazonaws.com",
                }
            )
        return json.dumps({"authorizationData": auth_data})

    def get_download_url_for_layer(self) -> None:
        self.error_on_dryrun()
        raise NotImplementedError(
            "ECR.get_download_url_for_layer is not yet implemented"
        )

    def get_repository_policy(self) -> str:
        registry_id = self._get_param("registryId")
        repository_name = self._get_param("repositoryName")

        return json.dumps(
            self.ecr_backend.get_repository_policy(
                registry_id=registry_id, repository_name=repository_name
            )
        )

    def initiate_layer_upload(self) -> None:
        self.error_on_dryrun()
        raise NotImplementedError("ECR.initiate_layer_upload is not yet implemented")

    def set_repository_policy(self) -> str:
        registry_id = self._get_param("registryId")
        repository_name = self._get_param("repositoryName")
        policy_text = self._get_param("policyText")
        # this is usually a safety flag to prevent accidental repository lock outs
        # but this would need a much deeper validation of the provided policy
        # force = self._get_param("force")

        return json.dumps(
            self.ecr_backend.set_repository_policy(
                registry_id=registry_id,
                repository_name=repository_name,
                policy_text=policy_text,
            )
        )

    def upload_layer_part(self) -> None:
        self.error_on_dryrun()
        raise NotImplementedError("ECR.upload_layer_part is not yet implemented")

    def list_tags_for_resource(self) -> str:
        arn = self._get_param("resourceArn")

        return json.dumps(self.ecr_backend.list_tags_for_resource(arn))

    def tag_resource(self) -> str:
        arn = self._get_param("resourceArn")
        tags = self._get_param("tags", [])

        self.ecr_backend.tag_resource(arn, tags)
        return "{}"

    def untag_resource(self) -> str:
        arn = self._get_param("resourceArn")
        tag_keys = self._get_param("tagKeys", [])

        self.ecr_backend.untag_resource(arn, tag_keys)
        return "{}"

    def put_image_tag_mutability(self) -> str:
        registry_id = self._get_param("registryId")
        repository_name = self._get_param("repositoryName")
        image_tag_mutability = self._get_param("imageTagMutability")

        return json.dumps(
            self.ecr_backend.put_image_tag_mutability(
                registry_id=registry_id,
                repository_name=repository_name,
                image_tag_mutability=image_tag_mutability,
            )
        )

    def put_image_scanning_configuration(self) -> str:
        registry_id = self._get_param("registryId")
        repository_name = self._get_param("repositoryName")
        image_scan_config = self._get_param("imageScanningConfiguration")

        return json.dumps(
            self.ecr_backend.put_image_scanning_configuration(
                registry_id=registry_id,
                repository_name=repository_name,
                image_scan_config=image_scan_config,
            )
        )

    def put_lifecycle_policy(self) -> str:
        registry_id = self._get_param("registryId")
        repository_name = self._get_param("repositoryName")
        lifecycle_policy_text = self._get_param("lifecyclePolicyText")

        return json.dumps(
            self.ecr_backend.put_lifecycle_policy(
                registry_id=registry_id,
                repository_name=repository_name,
                lifecycle_policy_text=lifecycle_policy_text,
            )
        )

    def get_lifecycle_policy(self) -> str:
        registry_id = self._get_param("registryId")
        repository_name = self._get_param("repositoryName")

        return json.dumps(
            self.ecr_backend.get_lifecycle_policy(
                registry_id=registry_id, repository_name=repository_name
            )
        )

    def delete_lifecycle_policy(self) -> str:
        registry_id = self._get_param("registryId")
        repository_name = self._get_param("repositoryName")

        return json.dumps(
            self.ecr_backend.delete_lifecycle_policy(
                registry_id=registry_id, repository_name=repository_name
            )
        )

    def put_registry_policy(self) -> str:
        policy_text = self._get_param("policyText")

        return json.dumps(self.ecr_backend.put_registry_policy(policy_text=policy_text))

    def get_registry_policy(self) -> str:
        return json.dumps(self.ecr_backend.get_registry_policy())

    def delete_registry_policy(self) -> str:
        return json.dumps(self.ecr_backend.delete_registry_policy())

    def start_image_scan(self) -> str:
        registry_id = self._get_param("registryId")
        repository_name = self._get_param("repositoryName")
        image_id = self._get_param("imageId")

        return json.dumps(
            self.ecr_backend.start_image_scan(
                registry_id=registry_id,
                repository_name=repository_name,
                image_id=image_id,
            )
        )

    def describe_image_scan_findings(self) -> str:
        registry_id = self._get_param("registryId")
        repository_name = self._get_param("repositoryName")
        image_id = self._get_param("imageId")

        return json.dumps(
            self.ecr_backend.describe_image_scan_findings(
                registry_id=registry_id,
                repository_name=repository_name,
                image_id=image_id,
            )
        )

    def put_replication_configuration(self) -> str:
        replication_config = self._get_param("replicationConfiguration")

        return json.dumps(
            self.ecr_backend.put_replication_configuration(
                replication_config=replication_config
            )
        )

    def get_registry_scanning_configuration(self) -> str:
        registry_scanning_config = (
            self.ecr_backend.get_registry_scanning_configuration()
        )
        return json.dumps(
            {
                "registryId": self.current_account,
                "scanningConfiguration": registry_scanning_config,
            }
        )

    def put_registry_scanning_configuration(self) -> str:
        scan_type = self._get_param("scanType")
        rules = self._get_param("rules")
        self.ecr_backend.put_registry_scanning_configuration(scan_type, rules)
        return json.dumps(
            {"registryScanningConfiguration": {"scanType": scan_type, "rules": rules}}
        )

    def describe_registry(self) -> str:
        return json.dumps(self.ecr_backend.describe_registry())
