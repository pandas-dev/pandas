"""KafkaBackend class with methods for supported APIs."""

import uuid
from datetime import datetime
from typing import Any, Dict, List, Optional, Tuple

from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.common_models import BaseModel
from moto.utilities.utils import get_partition

from ..utilities.tagging_service import TaggingService


class FakeKafkaCluster(BaseModel):
    def __init__(
        self,
        cluster_name: str,
        account_id: str,
        region_name: str,
        cluster_type: str,
        tags: Optional[Dict[str, str]] = None,
        broker_node_group_info: Optional[Dict[str, Any]] = None,
        kafka_version: Optional[str] = None,
        number_of_broker_nodes: Optional[int] = None,
        configuration_info: Optional[Dict[str, Any]] = None,
        serverless_config: Optional[Dict[str, Any]] = None,
        encryption_info: Optional[Dict[str, Any]] = None,
        enhanced_monitoring: str = "DEFAULT",
        open_monitoring: Optional[Dict[str, Any]] = None,
        logging_info: Optional[Dict[str, Any]] = None,
        storage_mode: str = "LOCAL",
        current_version: str = "1.0",
        client_authentication: Optional[Dict[str, Any]] = None,
        state: str = "CREATING",
        active_operation_arn: Optional[str] = None,
        zookeeper_connect_string: Optional[str] = None,
        zookeeper_connect_string_tls: Optional[str] = None,
    ):
        # General attributes
        self.cluster_id = str(uuid.uuid4())
        self.cluster_name = cluster_name
        self.account_id = account_id
        self.region_name = region_name
        self.cluster_type = cluster_type
        self.tags = tags or {}
        self.state = state
        self.creation_time = datetime.now().isoformat()
        self.current_version = current_version
        self.active_operation_arn = active_operation_arn
        self.arn = self._generate_arn()

        # Attributes specific to PROVISIONED clusters
        self.broker_node_group_info = broker_node_group_info
        self.kafka_version = kafka_version
        self.number_of_broker_nodes = number_of_broker_nodes
        self.configuration_info = configuration_info
        self.encryption_info = encryption_info
        self.enhanced_monitoring = enhanced_monitoring
        self.open_monitoring = open_monitoring
        self.logging_info = logging_info
        self.storage_mode = storage_mode
        self.client_authentication = client_authentication
        self.zookeeper_connect_string = zookeeper_connect_string
        self.zookeeper_connect_string_tls = zookeeper_connect_string_tls

        # Attributes specific to SERVERLESS clusters
        self.serverless_config = serverless_config

    def _generate_arn(self) -> str:
        resource_type = (
            "cluster" if self.cluster_type == "PROVISIONED" else "serverless-cluster"
        )
        partition = get_partition(self.region_name)
        return f"arn:{partition}:kafka:{self.region_name}:{self.account_id}:{resource_type}/{self.cluster_id}"


class KafkaBackend(BaseBackend):
    """Implementation of Kafka APIs."""

    def __init__(self, region_name: str, account_id: str):
        super().__init__(region_name, account_id)
        self.clusters: Dict[str, FakeKafkaCluster] = {}
        self.tagger = TaggingService()

    def create_cluster_v2(
        self,
        cluster_name: str,
        tags: Optional[Dict[str, str]],
        provisioned: Optional[Dict[str, Any]],
        serverless: Optional[Dict[str, Any]],
    ) -> Tuple[str, str, str, str]:
        if provisioned:
            cluster_type = "PROVISIONED"
            broker_node_group_info = provisioned.get("brokerNodeGroupInfo")
            kafka_version = provisioned.get("kafkaVersion", "default-kafka-version")
            number_of_broker_nodes = int(provisioned.get("numberOfBrokerNodes", 1))
            storage_mode = provisioned.get("storageMode", "LOCAL")
            serverless_config = None
        elif serverless:
            cluster_type = "SERVERLESS"
            broker_node_group_info = None
            kafka_version = None
            number_of_broker_nodes = None
            storage_mode = None
            serverless_config = serverless

        new_cluster = FakeKafkaCluster(
            cluster_name=cluster_name,
            account_id=self.account_id,
            region_name=self.region_name,
            cluster_type=cluster_type,
            broker_node_group_info=broker_node_group_info,
            kafka_version=kafka_version,
            number_of_broker_nodes=number_of_broker_nodes,
            serverless_config=serverless_config,
            tags=tags,
            state="CREATING",
            storage_mode=storage_mode if storage_mode else "LOCAL",
            current_version="1.0",
        )

        self.clusters[new_cluster.arn] = new_cluster

        if tags:
            self.tag_resource(new_cluster.arn, tags)

        return (
            new_cluster.arn,
            new_cluster.cluster_name,
            new_cluster.state,
            new_cluster.cluster_type,
        )

    def describe_cluster_v2(self, cluster_arn: str) -> Dict[str, Any]:
        cluster = self.clusters[cluster_arn]

        cluster_info: Dict[str, Any] = {
            "activeOperationArn": "arn:aws:kafka:region:account-id:operation/active-operation",
            "clusterArn": cluster.arn,
            "clusterName": cluster.cluster_name,
            "clusterType": cluster.cluster_type,
            "creationTime": cluster.creation_time,
            "currentVersion": cluster.current_version,
            "state": cluster.state,
            "stateInfo": {
                "code": "string",
                "message": "Cluster state details.",
            },
            "tags": self.list_tags_for_resource(cluster.arn),
        }

        if cluster.cluster_type == "PROVISIONED":
            cluster_info.update(
                {
                    "provisioned": {
                        "brokerNodeGroupInfo": cluster.broker_node_group_info or {},
                        "clientAuthentication": cluster.client_authentication or {},
                        "currentBrokerSoftwareInfo": {
                            "configurationArn": (cluster.configuration_info or {}).get(
                                "arn", "string"
                            ),
                            "configurationRevision": (
                                cluster.configuration_info or {}
                            ).get("revision", 1),
                            "kafkaVersion": cluster.kafka_version,
                        },
                        "encryptionInfo": cluster.encryption_info or {},
                        "enhancedMonitoring": cluster.enhanced_monitoring,
                        "openMonitoring": cluster.open_monitoring or {},
                        "loggingInfo": cluster.logging_info or {},
                        "numberOfBrokerNodes": cluster.number_of_broker_nodes or 0,
                        "zookeeperConnectString": cluster.zookeeper_connect_string
                        or "zookeeper.example.com:2181",
                        "zookeeperConnectStringTls": cluster.zookeeper_connect_string_tls
                        or "zookeeper.example.com:2181",
                        "storageMode": cluster.storage_mode,
                        "customerActionStatus": "NONE",
                    }
                }
            )

        elif cluster.cluster_type == "SERVERLESS":
            cluster_info.update(
                {
                    "serverless": {
                        "vpcConfigs": cluster.serverless_config.get("vpcConfigs", [])
                        if cluster.serverless_config
                        else [],
                        "clientAuthentication": cluster.serverless_config.get(
                            "clientAuthentication", {}
                        )
                        if cluster.serverless_config
                        else {},
                    }
                }
            )

        return cluster_info

    def list_clusters_v2(
        self,
        cluster_name_filter: Optional[str],
        cluster_type_filter: Optional[str],
        max_results: Optional[int],
        next_token: Optional[str],
    ) -> Tuple[List[Dict[str, Any]], Optional[str]]:
        cluster_info_list = []
        for cluster_arn in self.clusters.keys():
            cluster_info = self.describe_cluster_v2(cluster_arn)
            cluster_info_list.append(cluster_info)

        return cluster_info_list, None

    def create_cluster(
        self,
        broker_node_group_info: Dict[str, Any],
        client_authentication: Optional[Dict[str, Any]],
        cluster_name: str,
        configuration_info: Optional[Dict[str, Any]] = None,
        encryption_info: Optional[Dict[str, Any]] = None,
        enhanced_monitoring: str = "DEFAULT",
        open_monitoring: Optional[Dict[str, Any]] = None,
        kafka_version: str = "2.8.1",
        logging_info: Optional[Dict[str, Any]] = None,
        number_of_broker_nodes: int = 1,
        tags: Optional[Dict[str, str]] = None,
        storage_mode: str = "LOCAL",
    ) -> Tuple[str, str, str]:
        new_cluster = FakeKafkaCluster(
            cluster_name=cluster_name,
            account_id=self.account_id,
            region_name=self.region_name,
            cluster_type="PROVISIONED",
            broker_node_group_info=broker_node_group_info,
            client_authentication=client_authentication,
            kafka_version=kafka_version,
            number_of_broker_nodes=number_of_broker_nodes,
            configuration_info=configuration_info,
            encryption_info=encryption_info,
            enhanced_monitoring=enhanced_monitoring,
            open_monitoring=open_monitoring,
            logging_info=logging_info,
            storage_mode=storage_mode,
        )

        self.clusters[new_cluster.arn] = new_cluster

        if tags:
            self.tag_resource(new_cluster.arn, tags)

        return new_cluster.arn, new_cluster.cluster_name, new_cluster.state

    def describe_cluster(self, cluster_arn: str) -> Dict[str, Any]:
        cluster = self.clusters[cluster_arn]

        return {
            "activeOperationArn": "arn:aws:kafka:region:account-id:operation/active-operation",
            "brokerNodeGroupInfo": cluster.broker_node_group_info or {},
            "clientAuthentication": cluster.client_authentication or {},
            "clusterArn": cluster.arn,
            "clusterName": cluster.cluster_name,
            "creationTime": cluster.creation_time,
            "currentBrokerSoftwareInfo": {
                "configurationArn": (cluster.configuration_info or {}).get(
                    "arn", "string"
                ),
                "configurationRevision": (cluster.configuration_info or {}).get(
                    "revision", 1
                ),
                "kafkaVersion": cluster.kafka_version,
            },
            "currentVersion": cluster.current_version,
            "encryptionInfo": cluster.encryption_info or {},
            "enhancedMonitoring": cluster.enhanced_monitoring,
            "openMonitoring": cluster.open_monitoring or {},
            "loggingInfo": cluster.logging_info or {},
            "numberOfBrokerNodes": cluster.number_of_broker_nodes or 0,
            "state": cluster.state,
            "stateInfo": {
                "code": "string",
                "message": "Cluster state details.",
            },
            "tags": self.list_tags_for_resource(cluster.arn),
            "zookeeperConnectString": cluster.zookeeper_connect_string
            or "zookeeper.example.com:2181",
            "zookeeperConnectStringTls": cluster.zookeeper_connect_string_tls
            or "zookeeper.example.com:2181",
            "storageMode": cluster.storage_mode,
            "customerActionStatus": "NONE",
        }

    def list_clusters(
        self,
        cluster_name_filter: Optional[str],
        max_results: Optional[int],
        next_token: Optional[str],
    ) -> List[Dict[str, Any]]:
        cluster_info_list = [
            {
                "clusterArn": cluster.arn,
                "clusterName": cluster.cluster_name,
                "state": cluster.state,
                "creationTime": cluster.creation_time,
                "clusterType": cluster.cluster_type,
            }
            for cluster_arn, cluster in self.clusters.items()
        ]

        return cluster_info_list

    def delete_cluster(self, cluster_arn: str, current_version: str) -> Tuple[str, str]:
        cluster = self.clusters.pop(cluster_arn)
        return cluster_arn, cluster.state

    def list_tags_for_resource(self, resource_arn: str) -> Dict[str, str]:
        return self.tagger.get_tag_dict_for_resource(resource_arn)

    def tag_resource(self, resource_arn: str, tags: Dict[str, str]) -> None:
        tags_list = [{"Key": k, "Value": v} for k, v in tags.items()]
        self.tagger.tag_resource(resource_arn, tags_list)

    def untag_resource(self, resource_arn: str, tag_keys: List[str]) -> None:
        self.tagger.untag_resource_using_names(resource_arn, tag_keys)


kafka_backends = BackendDict(KafkaBackend, "kafka")
