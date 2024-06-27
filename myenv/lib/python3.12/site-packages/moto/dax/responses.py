import json
import re

from moto.core.responses import BaseResponse

from .exceptions import InvalidParameterValueException
from .models import DAXBackend, dax_backends


class DAXResponse(BaseResponse):
    def __init__(self) -> None:
        super().__init__(service_name="dax")

    @property
    def dax_backend(self) -> DAXBackend:
        return dax_backends[self.current_account][self.region]

    def create_cluster(self) -> str:
        params = json.loads(self.body)
        cluster_name = params.get("ClusterName")
        node_type = params.get("NodeType")
        description = params.get("Description")
        replication_factor = params.get("ReplicationFactor")
        iam_role_arn = params.get("IamRoleArn")
        tags = params.get("Tags", [])
        sse_specification = params.get("SSESpecification", {})
        encryption_type = params.get("ClusterEndpointEncryptionType", "NONE")

        self._validate_arn(iam_role_arn)
        self._validate_name(cluster_name)

        cluster = self.dax_backend.create_cluster(
            cluster_name=cluster_name,
            node_type=node_type,
            description=description,
            replication_factor=replication_factor,
            iam_role_arn=iam_role_arn,
            tags=tags,
            sse_specification=sse_specification,
            encryption_type=encryption_type,
        )
        return json.dumps(dict(Cluster=cluster.to_json()))

    def delete_cluster(self) -> str:
        cluster_name = json.loads(self.body).get("ClusterName")
        cluster = self.dax_backend.delete_cluster(cluster_name)
        return json.dumps(dict(Cluster=cluster.to_json()))

    def describe_clusters(self) -> str:
        params = json.loads(self.body)
        cluster_names = params.get("ClusterNames", [])
        max_results = params.get("MaxResults")
        next_token = params.get("NextToken")

        for name in cluster_names:
            self._validate_name(name)

        clusters, next_token = self.dax_backend.describe_clusters(
            cluster_names=cluster_names, max_results=max_results, next_token=next_token
        )
        return json.dumps(
            {"Clusters": [c.to_json() for c in clusters], "NextToken": next_token}
        )

    def _validate_arn(self, arn: str) -> None:
        if not arn.startswith("arn:"):
            raise InvalidParameterValueException(f"ARNs must start with 'arn:': {arn}")
        sections = arn.split(":")
        if len(sections) < 3:
            raise InvalidParameterValueException(
                f"Second colon partition not found: {arn}"
            )
        if len(sections) < 4:
            raise InvalidParameterValueException(f"Third colon vendor not found: {arn}")
        if len(sections) < 5:
            raise InvalidParameterValueException(
                f"Fourth colon (region/namespace delimiter) not found: {arn}"
            )
        if len(sections) < 6:
            raise InvalidParameterValueException(
                f"Fifth colon (namespace/relative-id delimiter) not found: {arn}"
            )

    def _validate_name(self, name: str) -> None:
        msg = "Cluster ID specified is not a valid identifier. Identifiers must begin with a letter; must contain only ASCII letters, digits, and hyphens; and must not end with a hyphen or contain two consecutive hyphens."
        if not re.match("^[a-z][a-z0-9-]+[a-z0-9]$", name):
            raise InvalidParameterValueException(msg)
        if "--" in name:
            raise InvalidParameterValueException(msg)

    def list_tags(self) -> str:
        params = json.loads(self.body)
        resource_name = params.get("ResourceName")
        tags = self.dax_backend.list_tags(resource_name=resource_name)
        return json.dumps(tags)

    def increase_replication_factor(self) -> str:
        params = json.loads(self.body)
        cluster_name = params.get("ClusterName")
        new_replication_factor = params.get("NewReplicationFactor")
        cluster = self.dax_backend.increase_replication_factor(
            cluster_name=cluster_name, new_replication_factor=new_replication_factor
        )
        return json.dumps({"Cluster": cluster.to_json()})

    def decrease_replication_factor(self) -> str:
        params = json.loads(self.body)
        cluster_name = params.get("ClusterName")
        new_replication_factor = params.get("NewReplicationFactor")
        node_ids_to_remove = params.get("NodeIdsToRemove")
        cluster = self.dax_backend.decrease_replication_factor(
            cluster_name=cluster_name,
            new_replication_factor=new_replication_factor,
            node_ids_to_remove=node_ids_to_remove,
        )
        return json.dumps({"Cluster": cluster.to_json()})
