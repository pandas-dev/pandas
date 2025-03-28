"""Handles incoming kafka requests, invokes methods, returns responses."""

import json
from urllib.parse import unquote

from moto.core.responses import BaseResponse

from .models import KafkaBackend, kafka_backends


class KafkaResponse(BaseResponse):
    """Handler for Kafka requests and responses."""

    def __init__(self) -> None:
        super().__init__(service_name="kafka")

    @property
    def kafka_backend(self) -> KafkaBackend:
        """Return backend instance specific for this region."""
        return kafka_backends[self.current_account][self.region]

    def create_cluster_v2(self) -> str:
        cluster_name = self._get_param("clusterName")
        tags = self._get_param("tags")
        provisioned = self._get_param("provisioned")
        serverless = self._get_param("serverless")
        cluster_arn, cluster_name, state, cluster_type = (
            self.kafka_backend.create_cluster_v2(
                cluster_name=cluster_name,
                tags=tags,
                provisioned=provisioned,
                serverless=serverless,
            )
        )
        return json.dumps(
            dict(
                clusterArn=cluster_arn,
                clusterName=cluster_name,
                state=state,
                clusterType=cluster_type,
            )
        )

    def describe_cluster_v2(self) -> str:
        cluster_arn = unquote(self.parsed_url.path.split("/clusters/")[-1])
        cluster_info = self.kafka_backend.describe_cluster_v2(
            cluster_arn=cluster_arn,
        )
        return json.dumps(dict(clusterInfo=cluster_info))

    def list_clusters_v2(self) -> str:
        cluster_name_filter = self._get_param("clusterNameFilter")
        cluster_type_filter = self._get_param("clusterTypeFilter")
        max_results = self._get_param("maxResults")
        next_token = self._get_param("nextToken")
        cluster_info_list, next_token = self.kafka_backend.list_clusters_v2(
            cluster_name_filter=cluster_name_filter,
            cluster_type_filter=cluster_type_filter,
            max_results=max_results,
            next_token=next_token,
        )
        return json.dumps(dict(clusterInfoList=cluster_info_list, nextToken=next_token))

    def list_tags_for_resource(self) -> str:
        resource_arn = unquote(self.parsed_url.path.split("/tags/")[-1])
        tags = self.kafka_backend.list_tags_for_resource(
            resource_arn=resource_arn,
        )
        return json.dumps(dict(tags=tags))

    def tag_resource(self) -> str:
        resource_arn = unquote(self._get_param("resourceArn"))
        tags = self._get_param("tags")
        self.kafka_backend.tag_resource(
            resource_arn=resource_arn,
            tags=tags,
        )
        return json.dumps(dict())

    def untag_resource(self) -> str:
        resource_arn = unquote(self._get_param("resourceArn"))
        tag_keys = self.__dict__["data"]["tagKeys"]
        self.kafka_backend.untag_resource(
            resource_arn=resource_arn,
            tag_keys=tag_keys,
        )
        return json.dumps(dict())

    def create_cluster(self) -> str:
        broker_node_group_info = self._get_param("brokerNodeGroupInfo")
        client_authentication = self._get_param("clientAuthentication")
        cluster_name = self._get_param("clusterName")
        configuration_info = self._get_param("configurationInfo")
        encryption_info = self._get_param("encryptionInfo")
        enhanced_monitoring = self._get_param("enhancedMonitoring")
        open_monitoring = self._get_param("openMonitoring")
        kafka_version = self._get_param("kafkaVersion")
        logging_info = self._get_param("loggingInfo")
        number_of_broker_nodes = self._get_param("numberOfBrokerNodes")
        tags = self._get_param("tags")
        storage_mode = self._get_param("storageMode")
        cluster_arn, cluster_name, state = self.kafka_backend.create_cluster(
            broker_node_group_info=broker_node_group_info,
            client_authentication=client_authentication,
            cluster_name=cluster_name,
            configuration_info=configuration_info,
            encryption_info=encryption_info,
            enhanced_monitoring=enhanced_monitoring,
            open_monitoring=open_monitoring,
            kafka_version=kafka_version,
            logging_info=logging_info,
            number_of_broker_nodes=number_of_broker_nodes,
            tags=tags,
            storage_mode=storage_mode,
        )
        return json.dumps(
            dict(clusterArn=cluster_arn, clusterName=cluster_name, state=state)
        )

    def describe_cluster(self) -> str:
        cluster_arn = unquote(self.parsed_url.path.split("/clusters/")[-1])
        cluster_info = self.kafka_backend.describe_cluster(
            cluster_arn=cluster_arn,
        )
        return json.dumps(dict(clusterInfo=cluster_info))

    def delete_cluster(self) -> str:
        cluster_arn = unquote(self.parsed_url.path.split("/clusters/")[-1])
        current_version = self._get_param("currentVersion")
        cluster_arn, state = self.kafka_backend.delete_cluster(
            cluster_arn=cluster_arn,
            current_version=current_version,
        )
        return json.dumps(dict(clusterArn=cluster_arn, state=state))

    def list_clusters(self) -> str:
        cluster_name_filter = self._get_param("clusterNameFilter")
        max_results = self._get_param("maxResults")
        next_token = self._get_param("nextToken")

        cluster_info_list = self.kafka_backend.list_clusters(
            cluster_name_filter=cluster_name_filter,
            max_results=max_results,
            next_token=next_token,
        )

        return json.dumps(dict(clusterInfoList=cluster_info_list, nextToken=next_token))
