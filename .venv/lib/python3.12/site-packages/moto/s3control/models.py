import uuid
from collections import defaultdict
from datetime import datetime, timezone
from typing import Any, Optional

from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.common_models import BaseModel
from moto.moto_api._internal import mock_random
from moto.s3 import s3_backends
from moto.s3.exceptions import (
    InvalidPublicAccessBlockConfiguration,
    WrongPublicAccessBlockAccountIdError,
)
from moto.s3.models import PublicAccessBlock, S3Backend
from moto.utilities.paginator import paginate
from moto.utilities.tagging_service import TaggingService
from moto.utilities.utils import PARTITION_NAMES, get_partition

from .exceptions import (
    AccessPointNotFound,
    AccessPointPolicyNotFound,
    InvalidRequestException,
    MultiRegionAccessPointNotFound,
    MultiRegionAccessPointOperationNotFound,
    MultiRegionAccessPointPolicyNotFound,
    NoSuchPublicAccessBlockConfiguration,
    StorageLensConfigurationNotFound,
)

PAGINATION_MODEL = {
    "list_storage_lens_configurations": {
        "input_token": "next_token",
        "limit_default": 100,
        "unique_attribute": "id",
    },
    "list_access_points": {
        "input_token": "next_token",
        "limit_key": "max_results",
        "limit_default": 1000,
        "unique_attribute": "name",
    },
    "list_multi_region_access_points": {
        "input_token": "next_token",
        "limit_key": "max_results",
        "limit_default": 100,
        "unique_attribute": "name",
    },
}


class AccessPoint(BaseModel):
    def __init__(
        self,
        account_id: str,
        region_name: str,
        name: str,
        bucket: str,
        vpc_configuration: dict[str, Any],
        public_access_block_configuration: dict[str, Any],
    ):
        self.name = name
        self.alias = f"{name}-{mock_random.get_random_hex(34)}-s3alias"
        self.bucket = bucket
        self.created = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%S.%f")
        self.arn = f"arn:{get_partition(region_name)}:s3:us-east-1:{account_id}:accesspoint/{name}"
        self.policy: Optional[str] = None
        self.network_origin = "VPC" if vpc_configuration else "Internet"
        self.vpc_id = (vpc_configuration or {}).get("VpcId")
        pubc = public_access_block_configuration or {}
        self.pubc = {
            "BlockPublicAcls": pubc.get("BlockPublicAcls", "true"),
            "IgnorePublicAcls": pubc.get("IgnorePublicAcls", "true"),
            "BlockPublicPolicy": pubc.get("BlockPublicPolicy", "true"),
            "RestrictPublicBuckets": pubc.get("RestrictPublicBuckets", "true"),
        }

    def delete_policy(self) -> None:
        self.policy = None

    def set_policy(self, policy: str) -> None:
        self.policy = policy

    def has_policy(self) -> bool:
        return self.policy is not None


class MultiRegionAccessPoint(BaseModel):
    def __init__(
        self,
        name: str,
        public_access_block: dict[str, Any],
        regions: list[dict[str, str]],
    ):
        self.name = name
        self.alias = f"{name}-{mock_random.get_random_hex(10)}.mrap"
        self.created_at = datetime.now(timezone.utc)
        self.public_access_block = public_access_block or {
            "BlockPublicAcls": True,
            "IgnorePublicAcls": True,
            "BlockPublicPolicy": True,
            "RestrictPublicBuckets": True,
        }
        self.regions = regions
        self.status = "READY"
        self.policy: Optional[str] = None

    def set_policy(self, policy: str) -> None:
        self.policy = policy

    def has_policy(self) -> bool:
        return self.policy is not None

    def to_dict(self) -> dict[str, Any]:
        return {
            "Name": self.name or "",
            "Alias": self.alias,
            "CreatedAt": self.created_at.strftime("%Y-%m-%dT%H:%M:%S.%fZ"),
            "PublicAccessBlock": self.public_access_block,
            "Status": self.status,
            "Regions": self.regions,
        }


class MultiRegionAccessPointOperation(BaseModel):
    def __init__(
        self,
        account_id: str,
        operation: str,
        region_name: str,
        name: Optional[str] = None,
        details: Optional[dict[str, Any]] = None,
    ):
        self.request_token_arn = f"arn:aws:s3:{region_name}:{account_id}:async-request/mrap/{operation.lower()}/{uuid.uuid4().hex[:24]}"
        self.request_status = "SUCCEEDED"
        self.name = name
        self.operation = operation
        self.created_at = datetime.now(timezone.utc)
        self.details = details or {}

    def to_dict(self) -> dict[str, Any]:
        response: dict[str, Any] = {
            "RequestTokenARN": self.request_token_arn,
            "RequestStatus": self.request_status,
            "CreationTime": self.created_at.strftime("%Y-%m-%dT%H:%M:%S.%fZ"),
            "Operation": self.operation,
        }

        if self.operation == "CreateMultiRegionAccessPoint" and self.details:
            response["RequestParameters"] = {
                "CreateMultiRegionAccessPointRequest": {
                    "Name": self.name or "",
                    "Regions": self.details.get("Regions", []),
                    "PublicAccessBlock": self.details.get("PublicAccessBlock", {}),
                }
            }
            regions_response = [
                {"Name": r.get("Region", ""), "RequestStatus": self.request_status}
                for r in self.details.get("Regions", [])
            ]
            response["ResponseDetails"] = {
                "MultiRegionAccessPointDetails": {
                    "Regions": regions_response,
                }
            }
        elif self.operation == "DeleteMultiRegionAccessPoint":
            response["RequestParameters"] = {
                "DeleteMultiRegionAccessPointRequest": {
                    "Name": self.name or "",
                }
            }
            response["ResponseDetails"] = {
                "MultiRegionAccessPointDetails": {
                    "Regions": [],
                }
            }
        elif self.operation == "PutMultiRegionAccessPointPolicy":
            response["RequestParameters"] = {
                "PutMultiRegionAccessPointPolicyRequest": {
                    "Name": self.name or "",
                    "Policy": self.details.get("Policy", ""),
                }
            }
            response["ResponseDetails"] = {}

        return response


class StorageLensConfiguration(BaseModel):
    def __init__(
        self,
        account_id: str,
        config_id: str,
        storage_lens_configuration: dict[str, Any],
        tags: Optional[dict[str, str]] = None,
    ):
        self.account_id = account_id
        self.config_id = config_id
        self.config = storage_lens_configuration
        self.tags = tags or {}
        self.arn = f"arn:{get_partition('us-east-1')}:s3:us-east-1:{account_id}:storage-lens/{config_id}"


class S3ControlBackend(BaseBackend):
    def __init__(self, region_name: str, account_id: str):
        super().__init__(region_name, account_id)
        self.public_access_block: Optional[PublicAccessBlock] = None
        self.access_points: dict[str, dict[str, AccessPoint]] = defaultdict(dict)
        self.multi_region_access_points: dict[
            str, dict[str, MultiRegionAccessPoint]
        ] = defaultdict(dict)
        self.mrap_operations: dict[str, dict[str, MultiRegionAccessPointOperation]] = (
            defaultdict(dict)
        )
        self.storage_lens_configs: dict[str, StorageLensConfiguration] = {}
        self.tagger = TaggingService()

    def get_public_access_block(self, account_id: str) -> PublicAccessBlock:
        if account_id != self.account_id:
            raise WrongPublicAccessBlockAccountIdError()

        if not self.public_access_block:
            raise NoSuchPublicAccessBlockConfiguration()

        return self.public_access_block

    def delete_public_access_block(self, account_id: str) -> None:
        if account_id != self.account_id:
            raise WrongPublicAccessBlockAccountIdError()

        self.public_access_block = None

    def put_public_access_block(
        self, account_id: str, pub_block_config: dict[str, Any]
    ) -> None:
        if account_id != self.account_id:
            raise WrongPublicAccessBlockAccountIdError()

        if not pub_block_config:
            raise InvalidPublicAccessBlockConfiguration()

        self.public_access_block = PublicAccessBlock(
            pub_block_config.get("BlockPublicAcls"),
            pub_block_config.get("IgnorePublicAcls"),
            pub_block_config.get("BlockPublicPolicy"),
            pub_block_config.get("RestrictPublicBuckets"),
        )

    def create_access_point(
        self,
        account_id: str,
        name: str,
        bucket: str,
        vpc_configuration: dict[str, Any],
        public_access_block_configuration: dict[str, Any],
    ) -> AccessPoint:
        access_point = AccessPoint(
            account_id,
            region_name=self.region_name,
            name=name,
            bucket=bucket,
            vpc_configuration=vpc_configuration,
            public_access_block_configuration=public_access_block_configuration,
        )
        self.access_points[account_id][name] = access_point
        return access_point

    def delete_access_point(self, account_id: str, name: str) -> None:
        self.access_points[account_id].pop(name, None)

    def get_access_point(self, account_id: str, name: str) -> AccessPoint:
        if name not in self.access_points[account_id]:
            raise AccessPointNotFound(name)
        return self.access_points[account_id][name]

    def put_access_point_policy(self, account_id: str, name: str, policy: str) -> None:
        access_point = self.get_access_point(account_id, name)
        access_point.set_policy(policy)

    def get_access_point_policy(self, account_id: str, name: str) -> str:
        access_point = self.get_access_point(account_id, name)
        if access_point.has_policy():
            return access_point.policy  # type: ignore[return-value]
        raise AccessPointPolicyNotFound(name)

    def delete_access_point_policy(self, account_id: str, name: str) -> None:
        access_point = self.get_access_point(account_id, name)
        access_point.delete_policy()

    def get_access_point_policy_status(self, account_id: str, name: str) -> bool:
        self.get_access_point_policy(account_id, name)
        return True

    def create_multi_region_access_point(
        self,
        account_id: str,
        name: str,
        public_access_block: dict[str, Any],
        regions: list[dict[str, str]],
        region_name: str,
    ) -> MultiRegionAccessPointOperation:
        if name in self.multi_region_access_points[account_id]:
            raise InvalidRequestException(
                f"Multi-Region Access Point {name} already exists"
            )

        processed_regions: list[dict[str, str]] = []
        for region_item in regions:
            bucket_name = region_item.get("Bucket", "")
            found_region = "us-east-1"

            found = False
            for account_id_key in s3_backends:
                if found:
                    break
                for region_key, backend in s3_backends[account_id_key].items():
                    if region_key == "aws":
                        continue

                    if bucket_name in backend.buckets:
                        found_region = region_key
                        found = True
                        break

            processed_regions.append({"Bucket": bucket_name, "Region": found_region})

        mrap = MultiRegionAccessPoint(
            name=name,
            public_access_block=public_access_block,
            regions=processed_regions,
        )
        self.multi_region_access_points[account_id][name] = mrap

        operation = MultiRegionAccessPointOperation(
            account_id=account_id,
            operation="CreateMultiRegionAccessPoint",
            region_name=region_name,
            name=name,
            details={
                "Regions": processed_regions,
                "PublicAccessBlock": public_access_block,
            },
        )
        self.mrap_operations[account_id][operation.request_token_arn] = operation

        return operation

    def delete_multi_region_access_point(
        self,
        account_id: str,
        name: str,
        region_name: str,
    ) -> MultiRegionAccessPointOperation:
        if name not in self.multi_region_access_points[account_id]:
            raise MultiRegionAccessPointNotFound(name)

        del self.multi_region_access_points[account_id][name]

        operation = MultiRegionAccessPointOperation(
            account_id=account_id,
            operation="DeleteMultiRegionAccessPoint",
            region_name=region_name,
            name=name,
        )
        self.mrap_operations[account_id][operation.request_token_arn] = operation

        return operation

    def describe_multi_region_access_point_operation(
        self,
        account_id: str,
        request_token_arn: str,
    ) -> MultiRegionAccessPointOperation:
        if request_token_arn not in self.mrap_operations[account_id]:
            raise MultiRegionAccessPointOperationNotFound(request_token_arn)

        return self.mrap_operations[account_id][request_token_arn]

    def get_multi_region_access_point(
        self,
        account_id: str,
        name: str,
    ) -> MultiRegionAccessPoint:
        if name not in self.multi_region_access_points[account_id]:
            raise MultiRegionAccessPointNotFound(name)

        return self.multi_region_access_points[account_id][name]

    def get_multi_region_access_point_policy(
        self,
        account_id: str,
        name: str,
    ) -> str:
        mrap = self.get_multi_region_access_point(account_id, name)
        if not mrap.has_policy():
            raise MultiRegionAccessPointPolicyNotFound(name)
        return mrap.policy  # type: ignore[return-value]

    def _is_policy_public(
        self, policy: str, public_access_block: dict[str, Any]
    ) -> bool:
        block_public = public_access_block.get("BlockPublicPolicy")
        if block_public is True or block_public == "true":
            return False

        policy_no_spaces = policy.replace(" ", "")
        if '"Principal":"*"' in policy_no_spaces:
            return True
        if '"Principal":{"AWS":"*"}' in policy_no_spaces:
            return True

        return False

    def get_multi_region_access_point_policy_status(
        self,
        account_id: str,
        name: str,
    ) -> dict[str, Any]:
        mrap = self.get_multi_region_access_point(account_id, name)
        if not mrap.has_policy():
            raise MultiRegionAccessPointPolicyNotFound(name)
        is_public = self._is_policy_public(mrap.policy, mrap.public_access_block)  # type: ignore[arg-type]
        return {"IsPublic": is_public}

    @paginate(pagination_model=PAGINATION_MODEL)
    def list_multi_region_access_points(
        self,
        account_id: str,
        max_results: Optional[int] = None,
        next_token: Optional[str] = None,
    ) -> list[MultiRegionAccessPoint]:
        return list(self.multi_region_access_points[account_id].values())

    def put_multi_region_access_point_policy(
        self,
        account_id: str,
        name: str,
        policy: str,
        region_name: str,
    ) -> MultiRegionAccessPointOperation:
        mrap = self.get_multi_region_access_point(account_id, name)
        mrap.set_policy(policy)

        operation = MultiRegionAccessPointOperation(
            account_id=account_id,
            operation="PutMultiRegionAccessPointPolicy",
            region_name=region_name,
            name=name,
            details={"Policy": policy},
        )
        self.mrap_operations[account_id][operation.request_token_arn] = operation

        return operation

    def put_storage_lens_configuration(
        self,
        config_id: str,
        account_id: str,
        storage_lens_configuration: dict[str, Any],
        tags: Optional[dict[str, str]] = None,
    ) -> None:
        if account_id != self.account_id:
            raise WrongPublicAccessBlockAccountIdError()

        storage_lens = StorageLensConfiguration(
            account_id=account_id,
            config_id=config_id,
            storage_lens_configuration=storage_lens_configuration,
            tags=tags,
        )
        self.storage_lens_configs[config_id] = storage_lens

    def get_storage_lens_configuration(
        self, config_id: str, account_id: str
    ) -> StorageLensConfiguration:
        if config_id not in self.storage_lens_configs:
            raise StorageLensConfigurationNotFound(config_id)
        storage_lens_configuration = self.storage_lens_configs[config_id]
        return storage_lens_configuration

    def delete_storage_lens_configuration(
        self, config_id: str, account_id: str
    ) -> None:
        if account_id != self.account_id:
            raise WrongPublicAccessBlockAccountIdError()
        if config_id not in self.storage_lens_configs:
            raise StorageLensConfigurationNotFound(config_id)
        del self.storage_lens_configs[config_id]

    @paginate(pagination_model=PAGINATION_MODEL)
    def list_storage_lens_configurations(
        self, account_id: str
    ) -> list[StorageLensConfiguration]:
        storage_lens_configuration_list = list(self.storage_lens_configs.values())
        return storage_lens_configuration_list

    @paginate(pagination_model=PAGINATION_MODEL)
    def list_access_points(
        self,
        account_id: str,
        bucket: Optional[str] = None,
        max_results: Optional[int] = None,
        next_token: Optional[str] = None,
    ) -> list[AccessPoint]:
        account_access_points = self.access_points.get(account_id, {})
        all_access_points = list(account_access_points.values())

        if bucket:
            return [ap for ap in all_access_points if ap.bucket == bucket]
        return all_access_points

    def put_storage_lens_configuration_tagging(
        self, config_id: str, account_id: str, tags: dict[str, str]
    ) -> None:
        if account_id != self.account_id:
            raise WrongPublicAccessBlockAccountIdError()

        if config_id not in self.storage_lens_configs:
            raise AccessPointNotFound(config_id)

        self.storage_lens_configs[config_id].tags = tags

    def get_storage_lens_configuration_tagging(
        self, config_id: str, account_id: str
    ) -> dict[str, str]:
        if account_id != self.account_id:
            raise WrongPublicAccessBlockAccountIdError()
        if config_id not in self.storage_lens_configs:
            raise AccessPointNotFound(config_id)

        return self.storage_lens_configs[config_id].tags

    def list_tags_for_resource(self, resource_arn: str) -> list[dict[str, str]]:
        backend: S3Backend = s3_backends[self.account_id][self.partition]
        return backend.tagger.list_tags_for_resource(resource_arn)["Tags"]

    def tag_resource(self, resource_arn: str, tags: list[dict[str, str]]) -> None:
        backend: S3Backend = s3_backends[self.account_id][self.partition]
        backend.tagger.tag_resource(resource_arn, tags=tags)

    def untag_resource(self, resource_arn: str, tag_keys: list[str]) -> None:
        backend: S3Backend = s3_backends[self.account_id][self.partition]
        backend.tagger.untag_resource_using_names(resource_arn, tag_names=tag_keys)


s3control_backends = BackendDict(
    S3ControlBackend,
    "s3control",
    use_boto3_regions=False,
    additional_regions=PARTITION_NAMES,
)
