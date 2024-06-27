from collections import defaultdict
from datetime import datetime
from typing import Any, Dict, Optional

from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.common_models import BaseModel
from moto.moto_api._internal import mock_random
from moto.s3.exceptions import (
    InvalidPublicAccessBlockConfiguration,
    NoSuchPublicAccessBlockConfiguration,
    WrongPublicAccessBlockAccountIdError,
)
from moto.s3.models import PublicAccessBlock
from moto.utilities.utils import PARTITION_NAMES, get_partition

from .exceptions import AccessPointNotFound, AccessPointPolicyNotFound


class AccessPoint(BaseModel):
    def __init__(
        self,
        account_id: str,
        region_name: str,
        name: str,
        bucket: str,
        vpc_configuration: Dict[str, Any],
        public_access_block_configuration: Dict[str, Any],
    ):
        self.name = name
        self.alias = f"{name}-{mock_random.get_random_hex(34)}-s3alias"
        self.bucket = bucket
        self.created = datetime.now().strftime("%Y-%m-%dT%H:%M:%S.%f")
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


class S3ControlBackend(BaseBackend):
    def __init__(self, region_name: str, account_id: str):
        super().__init__(region_name, account_id)
        self.public_access_block: Optional[PublicAccessBlock] = None
        self.access_points: Dict[str, Dict[str, AccessPoint]] = defaultdict(dict)

    def get_public_access_block(self, account_id: str) -> PublicAccessBlock:
        # The account ID should equal the account id that is set for Moto:
        if account_id != self.account_id:
            raise WrongPublicAccessBlockAccountIdError()

        if not self.public_access_block:
            raise NoSuchPublicAccessBlockConfiguration()

        return self.public_access_block

    def delete_public_access_block(self, account_id: str) -> None:
        # The account ID should equal the account id that is set for Moto:
        if account_id != self.account_id:
            raise WrongPublicAccessBlockAccountIdError()

        self.public_access_block = None

    def put_public_access_block(
        self, account_id: str, pub_block_config: Dict[str, Any]
    ) -> None:
        # The account ID should equal the account id that is set for Moto:
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
        vpc_configuration: Dict[str, Any],
        public_access_block_configuration: Dict[str, Any],
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
        """
        We assume the policy status is always public
        """
        self.get_access_point_policy(account_id, name)
        return True


s3control_backends = BackendDict(
    S3ControlBackend,
    "s3control",
    use_boto3_regions=False,
    additional_regions=PARTITION_NAMES,
)
