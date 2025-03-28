"""DAXBackend class with methods for supported APIs."""

from typing import Any, Dict, Iterable, List

from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.common_models import BaseModel
from moto.core.utils import unix_time
from moto.moto_api._internal import mock_random as random
from moto.moto_api._internal.managed_state_model import ManagedState
from moto.utilities.paginator import paginate
from moto.utilities.tagging_service import TaggingService
from moto.utilities.utils import get_partition

from .exceptions import ClusterNotFoundFault
from .utils import PAGINATION_MODEL


class DaxParameterGroup(BaseModel):
    def __init__(self) -> None:
        self.name = "default.dax1.0"
        self.status = "in-sync"

    def to_json(self) -> Dict[str, Any]:
        return {
            "ParameterGroupName": self.name,
            "ParameterApplyStatus": self.status,
            "NodeIdsToReboot": [],
        }


class DaxNode:
    def __init__(self, endpoint: "DaxEndpoint", name: str, index: int):
        self.node_id = f"{name}-{chr(ord('a')+index)}"  # name-a, name-b, etc
        self.node_endpoint = {
            "Address": f"{self.node_id}.{endpoint.cluster_hex}.nodes.dax-clusters.{endpoint.region}.amazonaws.com",
            "Port": endpoint.port,
        }
        self.create_time = unix_time()
        # AWS spreads nodes across zones, i.e. three nodes will probably end up in us-east-1a, us-east-1b, us-east-1c
        # For simplicity, we'll 'deploy' everything to us-east-1a
        self.availability_zone = f"{endpoint.region}a"
        self.status = "available"
        self.parameter_status = "in-sync"

    def to_json(self) -> Dict[str, Any]:
        return {
            "NodeId": self.node_id,
            "Endpoint": self.node_endpoint,
            "NodeCreateTime": self.create_time,
            "AvailabilityZone": self.availability_zone,
            "NodeStatus": self.status,
            "ParameterGroupStatus": self.parameter_status,
        }


class DaxEndpoint:
    def __init__(self, name: str, cluster_hex: str, region: str):
        self.name = name
        self.cluster_hex = cluster_hex
        self.region = region
        self.port = 8111

    def to_json(self, full: bool = False) -> Dict[str, Any]:
        dct: Dict[str, Any] = {"Port": self.port}
        if full:
            dct["Address"] = (
                f"{self.name}.{self.cluster_hex}.dax-clusters.{self.region}.amazonaws.com"
            )
            dct["URL"] = f"dax://{dct['Address']}"
        return dct


class DaxCluster(BaseModel, ManagedState):
    def __init__(
        self,
        account_id: str,
        region: str,
        name: str,
        description: str,
        node_type: str,
        replication_factor: int,
        iam_role_arn: str,
        sse_specification: Dict[str, Any],
        encryption_type: str,
    ):
        # Configure ManagedState
        super().__init__(
            model_name="dax::cluster",
            transitions=[("creating", "available"), ("deleting", "deleted")],
        )
        # Set internal properties
        self.name = name
        self.description = description
        self.arn = (
            f"arn:{get_partition(region)}:dax:{region}:{account_id}:cache/{self.name}"
        )
        self.node_type = node_type
        self.replication_factor = replication_factor
        self.cluster_hex = random.get_random_hex(6)
        self.endpoint = DaxEndpoint(
            name=name, cluster_hex=self.cluster_hex, region=region
        )
        self.nodes = [self._create_new_node(i) for i in range(0, replication_factor)]
        self.preferred_maintenance_window = "thu:23:30-fri:00:30"
        self.subnet_group = "default"
        self.iam_role_arn = iam_role_arn
        self.parameter_group = DaxParameterGroup()
        self.security_groups = [
            {
                "SecurityGroupIdentifier": f"sg-{random.get_random_hex(10)}",
                "Status": "active",
            }
        ]
        self.sse_specification = sse_specification
        self.encryption_type = encryption_type

    def _create_new_node(self, idx: int) -> DaxNode:
        return DaxNode(endpoint=self.endpoint, name=self.name, index=idx)

    def increase_replication_factor(self, new_replication_factor: int) -> None:
        for idx in range(self.replication_factor, new_replication_factor):
            self.nodes.append(self._create_new_node(idx))
        self.replication_factor = new_replication_factor

    def decrease_replication_factor(
        self, new_replication_factor: int, node_ids_to_remove: List[str]
    ) -> None:
        if node_ids_to_remove:
            self.nodes = [n for n in self.nodes if n.node_id not in node_ids_to_remove]
        else:
            self.nodes = self.nodes[0:new_replication_factor]
        self.replication_factor = new_replication_factor

    def delete(self) -> None:
        self.status = "deleting"

    def is_deleted(self) -> bool:
        return self.status == "deleted"

    def to_json(self) -> Dict[str, Any]:
        use_full_repr = self.status == "available"
        dct = {
            "ClusterName": self.name,
            "Description": self.description,
            "ClusterArn": self.arn,
            "TotalNodes": self.replication_factor,
            "ActiveNodes": 0,
            "NodeType": self.node_type,
            "Status": self.status,
            "ClusterDiscoveryEndpoint": self.endpoint.to_json(use_full_repr),
            "PreferredMaintenanceWindow": self.preferred_maintenance_window,
            "SubnetGroup": self.subnet_group,
            "IamRoleArn": self.iam_role_arn,
            "ParameterGroup": self.parameter_group.to_json(),
            "SSEDescription": {
                "Status": "ENABLED"
                if self.sse_specification.get("Enabled") is True
                else "DISABLED"
            },
            "ClusterEndpointEncryptionType": self.encryption_type,
            "SecurityGroups": self.security_groups,
        }
        if use_full_repr:
            dct["Nodes"] = [n.to_json() for n in self.nodes]
        return dct


class DAXBackend(BaseBackend):
    def __init__(self, region_name: str, account_id: str):
        super().__init__(region_name, account_id)
        self._clusters: Dict[str, DaxCluster] = dict()
        self._tagger = TaggingService()

    @property
    def clusters(self) -> Dict[str, DaxCluster]:
        self._clusters = {
            name: cluster
            for name, cluster in self._clusters.items()
            if cluster.status != "deleted"
        }
        return self._clusters

    def create_cluster(
        self,
        cluster_name: str,
        node_type: str,
        description: str,
        replication_factor: int,
        iam_role_arn: str,
        tags: List[Dict[str, str]],
        sse_specification: Dict[str, Any],
        encryption_type: str,
    ) -> DaxCluster:
        """
        The following parameters are not yet processed:
        AvailabilityZones, SubnetGroupNames, SecurityGroups, PreferredMaintenanceWindow, NotificationTopicArn, ParameterGroupName
        """
        cluster = DaxCluster(
            account_id=self.account_id,
            region=self.region_name,
            name=cluster_name,
            description=description,
            node_type=node_type,
            replication_factor=replication_factor,
            iam_role_arn=iam_role_arn,
            sse_specification=sse_specification,
            encryption_type=encryption_type,
        )
        self.clusters[cluster_name] = cluster
        self._tagger.tag_resource(cluster.arn, tags)
        return cluster

    def delete_cluster(self, cluster_name: str) -> DaxCluster:
        if cluster_name not in self.clusters:
            raise ClusterNotFoundFault()
        self.clusters[cluster_name].delete()
        return self.clusters[cluster_name]

    @paginate(PAGINATION_MODEL)
    def describe_clusters(self, cluster_names: Iterable[str]) -> List[DaxCluster]:
        clusters = self.clusters
        if not cluster_names:
            cluster_names = clusters.keys()

        for name in cluster_names:
            if name in self.clusters:
                self.clusters[name].advance()

        # Clusters may have been deleted while advancing the states
        clusters = self.clusters
        for name in cluster_names:
            if name not in self.clusters:
                raise ClusterNotFoundFault(name)
        return [cluster for name, cluster in clusters.items() if name in cluster_names]

    def list_tags(self, resource_name: str) -> Dict[str, List[Dict[str, str]]]:
        """
        Pagination is not yet implemented
        """
        # resource_name can be the name, or the full ARN
        name = resource_name.split("/")[-1]
        if name not in self.clusters:
            raise ClusterNotFoundFault()
        return self._tagger.list_tags_for_resource(self.clusters[name].arn)

    def increase_replication_factor(
        self, cluster_name: str, new_replication_factor: int
    ) -> DaxCluster:
        """
        The AvailabilityZones-parameter is not yet implemented
        """
        if cluster_name not in self.clusters:
            raise ClusterNotFoundFault()
        self.clusters[cluster_name].increase_replication_factor(new_replication_factor)
        return self.clusters[cluster_name]

    def decrease_replication_factor(
        self,
        cluster_name: str,
        new_replication_factor: int,
        node_ids_to_remove: List[str],
    ) -> DaxCluster:
        """
        The AvailabilityZones-parameter is not yet implemented
        """
        if cluster_name not in self.clusters:
            raise ClusterNotFoundFault()
        self.clusters[cluster_name].decrease_replication_factor(
            new_replication_factor, node_ids_to_remove
        )
        return self.clusters[cluster_name]


dax_backends = BackendDict(DAXBackend, "dax")
