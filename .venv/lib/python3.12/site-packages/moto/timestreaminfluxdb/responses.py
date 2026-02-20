import json

from moto.core.responses import BaseResponse

from .models import (
    DBStorageType,
    DeploymentType,
    NetworkType,
    TimestreamInfluxDBBackend,
    timestreaminfluxdb_backends,
)


class TimestreamInfluxDBResponse(BaseResponse):
    """Handler for TimestreamInfluxDB requests and responses."""

    def __init__(self) -> None:
        super().__init__(service_name="timestream-influxdb")

    @property
    def timestreaminfluxdb_backend(self) -> TimestreamInfluxDBBackend:
        return timestreaminfluxdb_backends[self.current_account][self.region]

    def create_db_instance(self) -> str:
        params = json.loads(self.body)
        name = params.get("name")
        username = params.get("username")
        password = params.get("password")
        organization = params.get("organization")
        bucket = params.get("bucket")
        db_instance_type = params.get("dbInstanceType")
        vpc_subnet_ids = params.get("vpcSubnetIds")
        vpc_security_group_ids = params.get("vpcSecurityGroupIds")
        publicly_accessible = params.get("publiclyAccessible", False)
        db_storage_type = params.get("dbStorageType", DBStorageType.InfluxIOIncludedT1)
        allocated_storage = params.get("allocatedStorage")
        db_parameter_group_identifier = params.get("dbParameterGroupIdentifier")
        deployment_type = params.get("deploymentType", DeploymentType.SINGLE_AZ)
        log_delivery_configuration = params.get("logDeliveryConfiguration", {})
        tags = params.get("tags", {})
        port = int(params.get("port", 8086))
        network_type = params.get("networkType", NetworkType.IPV4)

        created_instance = self.timestreaminfluxdb_backend.create_db_instance(
            name=name,
            username=username,
            password=password,
            organization=organization,
            bucket=bucket,
            db_instance_type=db_instance_type,
            vpc_subnet_ids=vpc_subnet_ids,
            vpc_security_group_ids=vpc_security_group_ids,
            publicly_accessible=publicly_accessible,
            db_storage_type=db_storage_type,
            allocated_storage=allocated_storage,
            db_parameter_group_identifier=db_parameter_group_identifier,
            deployment_type=deployment_type,
            log_delivery_configuration=log_delivery_configuration,
            tags=tags,
            port=port,
            network_type=network_type,
        )
        return json.dumps(created_instance.to_dict())

    def delete_db_instance(self) -> str:
        params = json.loads(self.body)
        id = params.get("identifier")
        deleted_instance = self.timestreaminfluxdb_backend.delete_db_instance(id=id)
        return json.dumps(deleted_instance.to_dict())

    def get_db_instance(self) -> str:
        params = json.loads(self.body)
        id = params.get("identifier")
        instance = self.timestreaminfluxdb_backend.get_db_instance(id=id)
        return json.dumps(instance.to_dict())

    def list_db_instances(self) -> str:
        """
        Pagination is not yet implemented
        """
        instances = self.timestreaminfluxdb_backend.list_db_instances()

        return json.dumps({"items": instances})

    def tag_resource(self) -> str:
        params = json.loads(self.body)
        arn = params.get("resourceArn")
        tags = params.get("tags")
        self.timestreaminfluxdb_backend.tag_resource(resource_arn=arn, tags=tags)
        return "{}"

    def untag_resource(self) -> str:
        params = json.loads(self.body)
        arn = params.get("resourceArn")
        tag_keys = params.get("tagKeys")
        self.timestreaminfluxdb_backend.untag_resource(
            resource_arn=arn, tag_keys=tag_keys
        )
        return "{}"

    def list_tags_for_resource(self) -> str:
        params = json.loads(self.body)
        arn = params.get("resourceArn")
        tags = self.timestreaminfluxdb_backend.list_tags_for_resource(resource_arn=arn)
        return json.dumps({"tags": tags})

    def create_db_parameter_group(self) -> str:
        params = json.loads(self.body)
        name = params.get("name")
        description = params.get("description")
        parameters = params.get("parameters")
        tags = params.get("tags")

        param_group = self.timestreaminfluxdb_backend.create_db_parameter_group(
            name=name,
            description=description,
            parameters=parameters,
            tags=tags,
        )

        return json.dumps(param_group.to_dict())

    def get_db_parameter_group(self) -> str:
        params = json.loads(self.body)
        identifier = params.get("identifier")

        param_group = self.timestreaminfluxdb_backend.get_db_parameter_group(
            identifier=identifier,
        )

        return json.dumps(param_group.to_dict())

    def list_db_parameter_groups(self) -> str:
        params = json.loads(self.body)
        next_token = params.get("nextToken")
        max_results = params.get("maxResults")

        items, next_token = self.timestreaminfluxdb_backend.list_db_parameter_groups(
            next_token=next_token,
            max_results=max_results,
        )

        response = {"items": items, "nextToken": next_token}
        return json.dumps(response)

    def list_db_clusters(self) -> str:
        params = json.loads(self.body)
        next_token = params.get("nextToken")
        max_results = params.get("maxResults")

        items, next_token = self.timestreaminfluxdb_backend.list_db_clusters(
            next_token=next_token,
            max_results=max_results,
        )

        response = {"items": items, "nextToken": next_token}
        return json.dumps(response)

    def get_db_cluster(self) -> str:
        params = json.loads(self.body)
        db_cluster_id = params.get("dbClusterId")

        cluster = self.timestreaminfluxdb_backend.get_db_cluster(
            db_cluster_id=db_cluster_id,
        )

        return json.dumps(cluster.to_dict())

    def create_db_cluster(self) -> str:
        params = json.loads(self.body)
        name = params.get("name")
        username = params.get("username")
        password = params.get("password")
        organization = params.get("organization")
        bucket = params.get("bucket")
        port = params.get("port")
        db_parameter_group_identifier = params.get("dbParameterGroupIdentifier")
        db_instance_type = params.get("dbInstanceType")
        db_storage_type = params.get("dbStorageType")
        allocated_storage = params.get("allocatedStorage")
        network_type = params.get("networkType")
        publicly_accessible = params.get("publiclyAccessible")
        vpc_subnet_ids = params.get("vpcSubnetIds")
        vpc_security_group_ids = params.get("vpcSecurityGroupIds")
        deployment_type = params.get("deploymentType")
        failover_mode = params.get("failoverMode")
        log_delivery_configuration = params.get("logDeliveryConfiguration")
        tags = params.get("tags")

        db_cluster_id, db_cluster_status = (
            self.timestreaminfluxdb_backend.create_db_cluster(
                name=name,
                username=username,
                password=password,
                organization=organization,
                bucket=bucket,
                port=port,
                db_parameter_group_identifier=db_parameter_group_identifier,
                db_instance_type=db_instance_type,
                db_storage_type=db_storage_type,
                allocated_storage=allocated_storage,
                network_type=network_type,
                publicly_accessible=publicly_accessible,
                vpc_subnet_ids=vpc_subnet_ids,
                vpc_security_group_ids=vpc_security_group_ids,
                deployment_type=deployment_type,
                failover_mode=failover_mode,
                log_delivery_configuration=log_delivery_configuration,
                tags=tags,
            )
        )

        return json.dumps(
            {"dbClusterId": db_cluster_id, "dbClusterStatus": db_cluster_status}
        )
