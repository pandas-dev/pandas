"""AuroraDSQLBackend class with methods for supported APIs."""

from collections import OrderedDict
from typing import Any, Dict, Optional

from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.common_models import BaseModel
from moto.core.utils import iso_8601_datetime_with_milliseconds, utcnow
from moto.moto_api._internal import mock_random
from moto.moto_api._internal.managed_state_model import ManagedState
from moto.utilities.utils import get_partition

from .exceptions import ValidationException


class Cluster(BaseModel, ManagedState):
    """Model for an AuroraDSQL cluster."""

    def to_dict(self) -> Dict[str, Any]:
        dct = {
            "identifier": self.identifier,
            "arn": self.arn,
            "status": self.status,
            "creationTime": iso_8601_datetime_with_milliseconds(self.creation_time),
            "deletionProtectionEnabled": self.deletion_protection_enabled,
        }
        return {k: v for k, v in dct.items() if v is not None}

    def __init__(
        self,
        region_name: str,
        account_id: str,
        deletion_protection_enabled: Optional[bool],
        tags: Optional[Dict[str, str]],
        client_token: Optional[str],
    ):
        ManagedState.__init__(
            self, "dsql::cluster", transitions=[("CREATING", "ACTIVE")]
        )
        self.region_name = region_name
        self.account_id = account_id
        self.identifier = mock_random.get_random_hex(26)
        self.arn = f"arn:{get_partition(self.region_name)}:dsql:{self.region_name}:{self.account_id}:cluster/{self.identifier}"
        self.creation_time = utcnow()
        self.deletion_protection_enabled = deletion_protection_enabled
        self.tags = tags
        self.client_token = client_token


class AuroraDSQLBackend(BaseBackend):
    """Implementation of AuroraDSQL APIs."""

    def __init__(self, region_name: str, account_id: str):
        super().__init__(region_name, account_id)
        self.region_name = region_name
        self.account_id = account_id
        self.partition = get_partition(region_name)
        self.clusters: Dict[str, Cluster] = OrderedDict()

    def create_cluster(
        self,
        deletion_protection_enabled: bool,
        tags: Optional[Dict[str, str]],
        client_token: Optional[str],
    ) -> Cluster:
        cluster = Cluster(
            self.region_name,
            self.account_id,
            deletion_protection_enabled,
            tags,
            client_token,
        )
        self.clusters[cluster.identifier] = cluster
        return cluster

    def get_cluster(
        self,
        identifier: str,
    ) -> Cluster:
        cluster = self.clusters.get(identifier)
        if cluster is None:
            raise ValidationException("invalid Cluster Id")

        return cluster


dsql_backends = BackendDict(
    AuroraDSQLBackend,
    "dsql",
    # currently botocore does not provide a dsql endpoint
    # https://github.com/boto/botocore/blob/e07cddc333fe4fb90efcd5d04324dd83f9cc3a57/botocore/data/endpoints.json
    use_boto3_regions=False,
    additional_regions=["us-east-1", "us-east-2"],
)
