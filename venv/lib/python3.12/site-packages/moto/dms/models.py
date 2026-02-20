from collections.abc import Iterable
from datetime import datetime
from typing import Any, Optional

from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.common_models import BaseModel
from moto.core.utils import utcnow
from moto.moto_api._internal.managed_state_model import ManagedState
from moto.utilities.tagging_service import TaggingService
from moto.utilities.utils import get_partition

from .exceptions import (
    InvalidResourceStateFault,
    ResourceAlreadyExistsFault,
    ResourceNotFoundFault,
    ValidationError,
)
from .utils import filter_tasks, random_id


class DatabaseMigrationServiceBackend(BaseBackend):
    def __init__(self, region_name: str, account_id: str):
        super().__init__(region_name, account_id)
        self.replication_tasks: dict[str, FakeReplicationTask] = {}
        self.replication_instances: dict[str, FakeReplicationInstance] = {}
        self.endpoints: dict[str, Endpoint] = {}
        self.replication_subnet_groups: dict[str, FakeReplicationSubnetGroup] = {}
        self.connections: list[FakeConnection] = []
        self.tagger = TaggingService()

    def create_replication_task(
        self,
        replication_task_identifier: str,
        source_endpoint_arn: str,
        target_endpoint_arn: str,
        replication_instance_arn: str,
        migration_type: str,
        table_mappings: str,
        replication_task_settings: str,
        tags: Optional[list[dict[str, str]]] = None,
    ) -> "FakeReplicationTask":
        """
        The following parameters are not yet implemented:
        CDCStartTime, CDCStartPosition, CDCStopPosition, Tags, TaskData, ResourceIdentifier
        """
        replication_task = FakeReplicationTask(
            replication_task_identifier=replication_task_identifier,
            source_endpoint_arn=source_endpoint_arn,
            target_endpoint_arn=target_endpoint_arn,
            replication_instance_arn=replication_instance_arn,
            migration_type=migration_type,
            table_mappings=table_mappings,
            replication_task_settings=replication_task_settings,
            account_id=self.account_id,
            region_name=self.region_name,
        )

        if self.replication_tasks.get(replication_task.arn):
            raise ResourceAlreadyExistsFault(
                "The resource you are attempting to create already exists."
            )

        self.replication_tasks[replication_task.arn] = replication_task

        if tags:
            self.tagger.tag_resource(replication_task.arn, tags)

        return replication_task

    def start_replication_task(
        self, replication_task_arn: str
    ) -> "FakeReplicationTask":
        """
        The following parameters have not yet been implemented:
        StartReplicationTaskType, CDCStartTime, CDCStartPosition, CDCStopPosition
        """
        if not self.replication_tasks.get(replication_task_arn):
            raise ResourceNotFoundFault("Replication task could not be found.")

        return self.replication_tasks[replication_task_arn].start()

    def stop_replication_task(self, replication_task_arn: str) -> "FakeReplicationTask":
        if not self.replication_tasks.get(replication_task_arn):
            raise ResourceNotFoundFault("Replication task could not be found.")

        return self.replication_tasks[replication_task_arn].stop()

    def delete_replication_task(
        self, replication_task_arn: str
    ) -> "FakeReplicationTask":
        if not self.replication_tasks.get(replication_task_arn):
            raise ResourceNotFoundFault("Replication task could not be found.")

        task = self.replication_tasks[replication_task_arn]
        task.delete()
        self.replication_tasks.pop(replication_task_arn)

        return task

    def describe_replication_tasks(
        self, filters: list[dict[str, Any]], max_records: int
    ) -> Iterable["FakeReplicationTask"]:
        """
        The parameter WithoutSettings has not yet been implemented
        """
        replication_tasks = filter_tasks(list(self.replication_tasks.values()), filters)

        if max_records and max_records > 0:
            replication_tasks = replication_tasks[:max_records]

        return replication_tasks

    def create_replication_instance(
        self,
        replication_instance_identifier: str,
        replication_instance_class: str,
        allocated_storage: Optional[int] = None,
        vpc_security_group_ids: Optional[list[str]] = None,
        availability_zone: Optional[str] = None,
        replication_subnet_group_identifier: Optional[str] = None,
        preferred_maintenance_window: Optional[str] = None,
        multi_az: Optional[bool] = False,
        engine_version: Optional[str] = None,
        auto_minor_version_upgrade: Optional[bool] = True,
        tags: Optional[list[dict[str, str]]] = None,
        kms_key_id: Optional[str] = None,
        publicly_accessible: Optional[bool] = True,
        dns_name_servers: Optional[str] = None,
        resource_identifier: Optional[str] = None,
        network_type: Optional[str] = None,
        kerberos_authentication_settings: Optional[dict[str, str]] = None,
    ) -> "FakeReplicationInstance":
        replication_instance = FakeReplicationInstance(
            replication_instance_identifier=replication_instance_identifier,
            replication_instance_class=replication_instance_class,
            allocated_storage=allocated_storage,
            vpc_security_group_ids=vpc_security_group_ids or [],
            availability_zone=availability_zone,
            replication_subnet_group_identifier=replication_subnet_group_identifier,
            preferred_maintenance_window=preferred_maintenance_window,
            multi_az=multi_az,
            engine_version=engine_version,
            auto_minor_version_upgrade=auto_minor_version_upgrade,
            tags=tags or [],
            kms_key_id=kms_key_id,
            publicly_accessible=publicly_accessible,
            dns_name_servers=dns_name_servers,
            resource_identifier=resource_identifier,
            network_type=network_type,
            kerberos_authentication_settings=kerberos_authentication_settings or {},
            account_id=self.account_id,
            region_name=self.region_name,
        )

        if self.replication_instances.get(replication_instance.arn):
            raise ResourceAlreadyExistsFault(
                "The resource you are attempting to create already exists."
            )

        self.replication_instances[replication_instance.arn] = replication_instance
        if tags:
            self.tagger.tag_resource(replication_instance.arn, tags)

        return replication_instance

    def describe_replication_instances(
        self,
        filters: Optional[list[dict[str, Any]]] = None,
        max_records: Optional[int] = None,
        marker: Optional[str] = None,
    ) -> list["FakeReplicationInstance"]:
        """Get information about replication instances with optional filtering"""
        ### TODO: Implement pagination

        replication_instances = list(self.replication_instances.values())

        if filters:
            for filter_obj in filters:
                filter_name = filter_obj.get("Name", "")
                filter_values = filter_obj.get("Values", [])

                if filter_name == "replication-instance-id":
                    replication_instances = [
                        instance
                        for instance in replication_instances
                        if instance.id in filter_values
                    ]
                elif filter_name == "replication-instance-arn":
                    replication_instances = [
                        instance
                        for instance in replication_instances
                        if instance.arn in filter_values
                    ]
                elif filter_name == "replication-instance-class":
                    replication_instances = [
                        instance
                        for instance in replication_instances
                        if instance.replication_instance_class in filter_values
                    ]
                elif filter_name == "engine-version":
                    replication_instances = [
                        instance
                        for instance in replication_instances
                        if instance.engine_version in filter_values
                    ]
        for replication_instance in replication_instances:
            replication_instance.advance()
        return replication_instances

    def delete_replication_instance(
        self, replication_instance_arn: str
    ) -> "FakeReplicationInstance":
        if not self.replication_instances.get(replication_instance_arn):
            raise ResourceNotFoundFault("Replication instance could not be found.")

        replication_instance = self.replication_instances[replication_instance_arn]
        replication_instance.delete()
        self.replication_instances.pop(replication_instance_arn)

        # remove any connections we might have tested linked to this instance
        connections = [
            connection
            for connection in self.connections
            if connection.replication_instance_arn == replication_instance_arn
        ]
        for connection in connections:
            self.connections.remove(connection)

        return replication_instance

    def create_endpoint(
        self,
        endpoint_identifier: str,
        endpoint_type: str,
        engine_name: str,
        username: str,
        password: str,
        server_name: str,
        port: int,
        database_name: str,
        extra_connection_attributes: str,
        kms_key_id: str,
        tags: Optional[list[dict[str, str]]],
        certificate_arn: str,
        ssl_mode: str,
        service_access_role_arn: str,
        external_table_definition: str,
        dynamo_db_settings: Optional[dict[str, Any]],
        s3_settings: Optional[dict[str, Any]],
        dms_transfer_settings: Optional[dict[str, Any]],
        mongo_db_settings: Optional[dict[str, Any]],
        kinesis_settings: Optional[dict[str, Any]],
        kafka_settings: Optional[dict[str, Any]],
        elasticsearch_settings: Optional[dict[str, Any]],
        neptune_settings: Optional[dict[str, Any]],
        redshift_settings: Optional[dict[str, Any]],
        postgre_sql_settings: Optional[dict[str, Any]],
        my_sql_settings: Optional[dict[str, Any]],
        oracle_settings: Optional[dict[str, Any]],
        sybase_settings: Optional[dict[str, Any]],
        microsoft_sql_server_settings: Optional[dict[str, Any]],
        ibm_db2_settings: Optional[dict[str, Any]],
        resource_identifier: Optional[str],
        doc_db_settings: Optional[dict[str, Any]],
        redis_settings: Optional[dict[str, Any]],
        gcp_my_sql_settings: Optional[dict[str, Any]],
        timestream_settings: Optional[dict[str, Any]],
    ) -> "Endpoint":
        if endpoint_type not in ["source", "target"]:
            raise ValidationError("Invalid endpoint type")

        endpoint = Endpoint(
            endpoint_identifier=endpoint_identifier,
            endpoint_type=endpoint_type,
            engine_name=engine_name,
            engine_display_name=engine_name,  # Assuming engine_display_name is same as engine_name
            username=username,
            password=password,
            server_name=server_name,
            port=port,
            database_name=database_name,
            extra_connection_attributes=extra_connection_attributes,
            status="creating",  # Default status
            kms_key_id=kms_key_id,
            certificate_arn=certificate_arn,
            ssl_mode=ssl_mode,
            service_access_role_arn=service_access_role_arn,
            external_table_definition=external_table_definition,
            external_id=None,  # Default value
            dynamo_db_settings=dynamo_db_settings,
            s3_settings=s3_settings,
            dms_transfer_settings=dms_transfer_settings,
            mongo_db_settings=mongo_db_settings,
            kinesis_settings=kinesis_settings,
            kafka_settings=kafka_settings,
            elasticsearch_settings=elasticsearch_settings,
            neptune_settings=neptune_settings,
            redshift_settings=redshift_settings,
            postgre_sql_settings=postgre_sql_settings,
            my_sql_settings=my_sql_settings,
            oracle_settings=oracle_settings,
            sybase_settings=sybase_settings,
            microsoft_sql_server_settings=microsoft_sql_server_settings,
            ibm_db2_settings=ibm_db2_settings,
            resource_identifier=resource_identifier,
            doc_db_settings=doc_db_settings,
            redis_settings=redis_settings,
            gcp_my_sql_settings=gcp_my_sql_settings,
            timestream_settings=timestream_settings,
            account_id=self.account_id,
            region_name=self.region_name,
        )

        if tags:
            self.tagger.tag_resource(endpoint.endpoint_arn, tags)

        self.endpoints[endpoint.endpoint_identifier] = endpoint
        return endpoint

    # TODO implement pagination
    def describe_endpoints(
        self,
        filters: list[dict[str, Any]],
        max_records: Optional[int],
        marker: Optional[str],
    ) -> list["Endpoint"]:
        endpoints = list(self.endpoints.values())
        filter_map = {
            "endpoint-arn": "endpoint_arn",
            "endpoint-type": "endpoint_type",
            "endpoint-id": "endpoint_identifier",
            "engine-name": "engine_name",
        }

        for filter in filters:
            filter_name = filter.get("Name")
            filter_values = filter.get("Values", [])
            if filter_name not in filter_map:
                raise ValueError(f"Invalid filter name: {filter_name}")

            attribute = filter_map[filter_name]
            endpoints = [
                endpoint
                for endpoint in endpoints
                if getattr(endpoint, attribute) in filter_values
            ]
        return endpoints

    def list_tags_for_resource(
        self, resource_arn_list: list[str]
    ) -> list[dict[str, Any]]:
        result = []
        for resource_arn in resource_arn_list:
            tags = self.tagger.get_tag_dict_for_resource(resource_arn)
            for key, value in tags.items():
                result.append({"Key": key, "ResourceArn": resource_arn, "Value": value})
        return result

    def delete_endpoint(self, endpoint_arn: str) -> "Endpoint":
        endpoints = [
            endpoint
            for endpoint in list(self.endpoints.values())
            if endpoint.endpoint_arn == endpoint_arn
        ]

        if len(endpoints) == 0:
            raise ResourceNotFoundFault("Endpoint could not be found.")
        endpoint = endpoints[0]
        endpoint.delete()
        self.endpoints.pop(endpoint.endpoint_identifier)

        # remove any connections we might have tested linked to this instance
        connections = [
            connection
            for connection in self.connections
            if connection.endpoint_arn == endpoint_arn
        ]
        for connection in connections:
            self.connections.remove(connection)
        return endpoint

    def create_replication_subnet_group(
        self,
        replication_subnet_group_identifier: str,
        replication_subnet_group_description: str,
        subnet_ids: list[str],
        tags: Optional[list[dict[str, str]]] = None,
    ) -> "FakeReplicationSubnetGroup":
        replication_subnet_group = FakeReplicationSubnetGroup(
            replication_subnet_group_identifier=replication_subnet_group_identifier,
            replication_subnet_group_description=replication_subnet_group_description,
        )
        if tags:
            self.tagger.tag_resource(
                replication_subnet_group.replication_subnet_group_identifier, tags
            )
        if self.replication_subnet_groups.get(replication_subnet_group_identifier):
            raise ResourceAlreadyExistsFault(
                "The resource you are attempting to create already exists."
            )

        self.replication_subnet_groups[
            replication_subnet_group.replication_subnet_group_identifier
        ] = replication_subnet_group
        return replication_subnet_group

    # TODO implement pagination
    def describe_replication_subnet_groups(
        self,
        filters: list[dict[str, Any]],
        max_records: Optional[int],
        marker: Optional[str],
    ) -> list["FakeReplicationSubnetGroup"]:
        subnet_groups = list(self.replication_subnet_groups.values())
        filter_map = {
            "replication-subnet-group-id": "replication_subnet_group_identifier"
        }

        for filter in filters:
            filter_name = filter.get("Name")
            filter_values = filter.get("Values", [])
            if filter_name not in filter_map:
                raise ValueError(f"Invalid filter name: {filter_name}")

            attribute = filter_map[filter_name]
            subnet_groups = [
                subnet_group
                for subnet_group in subnet_groups
                if getattr(subnet_group, attribute) in filter_values
            ]

        return subnet_groups

    def delete_replication_subnet_group(
        self, replication_subnet_group_identifier: str
    ) -> "FakeReplicationSubnetGroup":
        if not self.replication_subnet_groups.get(replication_subnet_group_identifier):
            raise ResourceNotFoundFault("Replication subnet group could not be found.")

        replication_subnet_group = self.replication_subnet_groups[
            replication_subnet_group_identifier
        ]
        self.replication_subnet_groups.pop(replication_subnet_group_identifier)
        return replication_subnet_group

    def test_connection(
        self, replication_instance_arn: str, endpoint_arn: str
    ) -> "FakeConnection":
        if not self.replication_instances.get(replication_instance_arn):
            raise ResourceNotFoundFault("Replication instance could not be found.")
        replication_instance = self.replication_instances[replication_instance_arn]
        endpoints = [
            endpoint
            for endpoint in list(self.endpoints.values())
            if endpoint.endpoint_arn == endpoint_arn
        ]

        if len(endpoints) == 0:
            raise ResourceNotFoundFault("Endpoint could not be found.")
        endpoint = endpoints[0]
        connections = self.connections
        connections = [
            connection
            for connection in connections
            if connection.replication_instance_arn == replication_instance_arn
        ]
        connections = [
            connection
            for connection in connections
            if connection.endpoint_arn == endpoint_arn
        ]
        if len(connections) == 0:
            connection = FakeConnection(
                replication_instance_arn=replication_instance.arn,
                endpoint_arn=endpoint.endpoint_arn,
                replication_instance_identifier=replication_instance.id,
                endpoint_identifier=endpoint.endpoint_identifier,
            )
            self.connections.append(connection)
        return connection

    # TODO implement pagination
    def describe_connections(
        self,
        filters: list[dict[str, Any]],
        max_records: Optional[int],
        marker: Optional[str],
    ) -> list["FakeConnection"]:
        filter_map = {
            "endpoint-arn": "endpoint_arn",
            "replication-instance-arn": "replication_instance_arn",
        }
        connections = self.connections
        for filter in filters:
            filter_name = filter.get("Name")
            filter_values = filter.get("Values", [])
            if filter_name not in filter_map:
                raise ValueError(f"Invalid filter name: {filter_name}")

            attribute = filter_map[filter_name]
            connections = [
                connection
                for connection in self.connections
                if getattr(connection, attribute) in filter_values
            ]
        for connection in connections:
            connection.advance()
        return connections


class Endpoint(BaseModel):
    def __init__(
        self,
        endpoint_identifier: str,
        endpoint_type: str,
        engine_name: str,
        engine_display_name: str,
        username: str,
        password: str,
        server_name: str,
        port: int,
        database_name: str,
        extra_connection_attributes: str,
        status: str,
        kms_key_id: str,
        certificate_arn: Optional[str],
        ssl_mode: Optional[str],
        service_access_role_arn: Optional[str],
        external_table_definition: Optional[str],
        external_id: Optional[str],
        dynamo_db_settings: Optional[dict[str, Any]],
        s3_settings: Optional[dict[str, Any]],
        dms_transfer_settings: Optional[dict[str, Any]],
        mongo_db_settings: Optional[dict[str, Any]],
        kinesis_settings: Optional[dict[str, Any]],
        kafka_settings: Optional[dict[str, Any]],
        elasticsearch_settings: Optional[dict[str, Any]],
        neptune_settings: Optional[dict[str, Any]],
        redshift_settings: Optional[dict[str, Any]],
        postgre_sql_settings: Optional[dict[str, Any]],
        my_sql_settings: Optional[dict[str, Any]],
        oracle_settings: Optional[dict[str, Any]],
        sybase_settings: Optional[dict[str, Any]],
        microsoft_sql_server_settings: Optional[dict[str, Any]],
        ibm_db2_settings: Optional[dict[str, Any]],
        resource_identifier: Optional[str],
        doc_db_settings: Optional[dict[str, Any]],
        redis_settings: Optional[dict[str, Any]],
        gcp_my_sql_settings: Optional[dict[str, Any]],
        timestream_settings: Optional[dict[str, Any]],
        account_id: str,
        region_name: str,
    ):
        self.endpoint_identifier = endpoint_identifier
        self.endpoint_type = endpoint_type
        self.engine_name = engine_name
        self.engine_display_name = engine_display_name
        self.username = username
        self.password = password
        self.server_name = server_name
        self.port = port
        self.database_name = database_name
        self.extra_connection_attributes = extra_connection_attributes
        self.status = status

        # In real AWS, KMS key is used to encrypt connection parameters
        # but just storing it here for now
        self.kms_key_id = kms_key_id

        self.id = resource_identifier if resource_identifier else random_id(True, 27)
        self.endpoint_arn = f"arn:{get_partition(region_name)}:dms:{region_name}:{account_id}:endpoint:{self.id}"

        self.certificate_arn = certificate_arn
        self.ssl_mode = ssl_mode
        self.service_access_role_arn = service_access_role_arn
        self.external_table_definition = external_table_definition
        self.external_id = external_id
        self.dynamo_db_settings = dynamo_db_settings
        self.s3_settings = s3_settings
        self.dms_transfer_settings = dms_transfer_settings
        self.mongo_db_settings = mongo_db_settings
        self.kinesis_settings = kinesis_settings
        self.kafka_settings = kafka_settings
        self.elasticsearch_settings = elasticsearch_settings
        self.neptune_settings = neptune_settings
        self.redshift_settings = redshift_settings
        self.postgre_sql_settings = postgre_sql_settings
        self.my_sql_settings = my_sql_settings
        self.oracle_settings = oracle_settings
        self.sybase_settings = sybase_settings
        self.microsoft_sql_server_settings = microsoft_sql_server_settings
        self.ibm_db2_settings = ibm_db2_settings
        self.doc_db_settings = doc_db_settings
        self.redis_settings = redis_settings
        self.gcp_my_sql_settings = gcp_my_sql_settings
        self.timestream_settings = timestream_settings

    def to_dict(self) -> dict[str, Any]:
        return {
            "EndpointIdentifier": self.endpoint_identifier,
            "EndpointType": self.endpoint_type,
            "EngineName": self.engine_name,
            "EngineDisplayName": self.engine_display_name,
            "Username": self.username,
            "ServerName": self.server_name,
            "Port": self.port,
            "DatabaseName": self.database_name,
            "ExtraConnectionAttributes": self.extra_connection_attributes,
            "Status": self.status,
            "KmsKeyId": self.kms_key_id,
            "EndpointArn": self.endpoint_arn,
            "CertificateArn": self.certificate_arn,
            "SslMode": self.ssl_mode,
            "ServiceAccessRoleArn": self.service_access_role_arn,
            "ExternalTableDefinition": self.external_table_definition,
            "ExternalId": self.external_id,
            "DynamoDbSettings": self.dynamo_db_settings,
            "S3Settings": self.s3_settings,
            "DmsTransferSettings": self.dms_transfer_settings,
            "MongoDbSettings": self.mongo_db_settings,
            "KinesisSettings": self.kinesis_settings,
            "KafkaSettings": self.kafka_settings,
            "ElasticsearchSettings": self.elasticsearch_settings,
            "NeptuneSettings": self.neptune_settings,
            "RedshiftSettings": self.redshift_settings,
            "PostgreSQLSettings": self.postgre_sql_settings,
            "MySQLSettings": self.my_sql_settings,
            "OracleSettings": self.oracle_settings,
            "SybaseSettings": self.sybase_settings,
            "MicrosoftSQLServerSettings": self.microsoft_sql_server_settings,
            "IBMDb2Settings": self.ibm_db2_settings,
            "DocDbSettings": self.doc_db_settings,
            "RedisSettings": self.redis_settings,
            "GcpMySQLSettings": self.gcp_my_sql_settings,
            "TimestreamSettings": self.timestream_settings,
        }

    def delete(self) -> "Endpoint":
        self.status = "deleting"
        return self


class FakeReplicationTask(BaseModel):
    def __init__(
        self,
        replication_task_identifier: str,
        migration_type: str,
        replication_instance_arn: str,
        source_endpoint_arn: str,
        target_endpoint_arn: str,
        table_mappings: str,
        replication_task_settings: str,
        account_id: str,
        region_name: str,
    ):
        self.id = replication_task_identifier
        self.region = region_name
        self.migration_type = migration_type
        self.replication_instance_arn = replication_instance_arn
        self.source_endpoint_arn = source_endpoint_arn
        self.target_endpoint_arn = target_endpoint_arn
        self.table_mappings = table_mappings
        self.replication_task_settings = replication_task_settings

        self.arn = f"arn:{get_partition(region_name)}:dms:{region_name}:{account_id}:task:{self.id}"
        self.status = "creating"

        self.creation_date = utcnow()
        self.start_date: Optional[datetime] = None
        self.stop_date: Optional[datetime] = None

    def to_dict(self) -> dict[str, Any]:
        start_date = self.start_date.isoformat() if self.start_date else None
        stop_date = self.stop_date.isoformat() if self.stop_date else None

        return {
            "ReplicationTaskIdentifier": self.id,
            "SourceEndpointArn": self.source_endpoint_arn,
            "TargetEndpointArn": self.target_endpoint_arn,
            "ReplicationInstanceArn": self.replication_instance_arn,
            "MigrationType": self.migration_type,
            "TableMappings": self.table_mappings,
            "ReplicationTaskSettings": self.replication_task_settings,
            "Status": self.status,
            "ReplicationTaskCreationDate": self.creation_date.isoformat(),
            "ReplicationTaskStartDate": start_date,
            "ReplicationTaskArn": self.arn,
            "ReplicationTaskStats": {
                "FullLoadProgressPercent": 100,
                "ElapsedTimeMillis": 100,
                "TablesLoaded": 1,
                "TablesLoading": 0,
                "TablesQueued": 0,
                "TablesErrored": 0,
                "FreshStartDate": start_date,
                "StartDate": start_date,
                "StopDate": stop_date,
                "FullLoadStartDate": start_date,
                "FullLoadFinishDate": stop_date,
            },
        }

    def ready(self) -> "FakeReplicationTask":
        self.status = "ready"
        return self

    def start(self) -> "FakeReplicationTask":
        self.status = "starting"
        self.start_date = utcnow()
        self.run()
        return self

    def stop(self) -> "FakeReplicationTask":
        if self.status != "running":
            raise InvalidResourceStateFault("Replication task is not running")

        self.status = "stopped"
        self.stop_date = utcnow()
        return self

    def delete(self) -> "FakeReplicationTask":
        self.status = "deleting"
        return self

    def run(self) -> "FakeReplicationTask":
        self.status = "running"
        return self


class FakeReplicationInstance(ManagedState):
    def __init__(
        self,
        replication_instance_identifier: str,
        replication_instance_class: str,
        account_id: str,
        region_name: str,
        allocated_storage: Optional[int] = None,
        vpc_security_group_ids: Optional[list[str]] = None,
        availability_zone: Optional[str] = None,
        replication_subnet_group_identifier: Optional[str] = None,
        preferred_maintenance_window: Optional[str] = None,
        multi_az: Optional[bool] = False,
        engine_version: Optional[str] = None,
        auto_minor_version_upgrade: Optional[bool] = True,
        tags: Optional[list[dict[str, str]]] = None,
        kms_key_id: Optional[str] = None,
        publicly_accessible: Optional[bool] = True,
        dns_name_servers: Optional[str] = None,
        resource_identifier: Optional[str] = None,
        network_type: Optional[str] = None,
        kerberos_authentication_settings: Optional[dict[str, str]] = None,
    ):
        ManagedState.__init__(
            self,
            model_name="dms::replicationinstance",
            transitions=[("creating", "available")],
        )

        self.id = replication_instance_identifier
        self.replication_instance_class = replication_instance_class
        self.region = region_name
        self.allocated_storage = allocated_storage or 50
        self.vpc_security_groups = [
            {"VpcSecurityGroupId": sg_id, "Status": "active"}
            for sg_id in (vpc_security_group_ids or [])
        ]
        self.availability_zone = availability_zone
        self.replication_subnet_group_identifier = replication_subnet_group_identifier
        self.preferred_maintenance_window = preferred_maintenance_window
        self.multi_az = multi_az
        self.engine_version = engine_version
        self.auto_minor_version_upgrade = auto_minor_version_upgrade
        self.tags = tags or []
        self.kms_key_id = kms_key_id
        self.publicly_accessible = publicly_accessible
        self.dns_name_servers = dns_name_servers
        self.resource_identifier = resource_identifier
        self.network_type = network_type
        self.kerberos_authentication_settings = kerberos_authentication_settings or {}
        self.arn = f"arn:{get_partition(region_name)}:dms:{region_name}:{account_id}:rep:{self.id}"
        self.creation_date = utcnow()
        self.private_ip_addresses = ["10.0.0.1"]
        self.public_ip_addresses = ["54.0.0.1"] if publicly_accessible else []
        self.ipv6_addresses: list[str] = []

    def to_dict(self) -> dict[str, Any]:
        kerberos_settings = None
        if self.kerberos_authentication_settings:
            kerberos_settings = {
                "KeyCacheSecretId": self.kerberos_authentication_settings.get(
                    "KeyCacheSecretId"
                ),
                "KeyCacheSecretIamArn": self.kerberos_authentication_settings.get(
                    "KeyCacheSecretIamArn"
                ),
                "Krb5FileContents": self.kerberos_authentication_settings.get(
                    "Krb5FileContents"
                ),
            }

        subnet_group = None
        if self.replication_subnet_group_identifier:
            subnet_group = {
                "ReplicationSubnetGroupIdentifier": self.replication_subnet_group_identifier,
                "ReplicationSubnetGroupDescription": f"Subnet group for {self.id}",
                "VpcId": "vpc-12345",
                "SubnetGroupStatus": "Complete",
                "Subnets": [
                    {
                        "SubnetIdentifier": "subnet-12345",
                        "SubnetAvailabilityZone": {
                            "Name": self.availability_zone or "us-east-1a"
                        },
                        "SubnetStatus": "Active",
                    }
                ],
                "SupportedNetworkTypes": [self.network_type]
                if self.network_type
                else ["IPV4"],
            }

        return {
            "ReplicationInstanceIdentifier": self.id,
            "ReplicationInstanceClass": self.replication_instance_class,
            "ReplicationInstanceStatus": self.status,
            "AllocatedStorage": self.allocated_storage,
            "InstanceCreateTime": self.creation_date.isoformat(),
            "VpcSecurityGroups": self.vpc_security_groups,
            "AvailabilityZone": self.availability_zone,
            "ReplicationSubnetGroup": subnet_group,
            "PreferredMaintenanceWindow": self.preferred_maintenance_window,
            "PendingModifiedValues": {},
            "MultiAZ": self.multi_az,
            "EngineVersion": self.engine_version,
            "AutoMinorVersionUpgrade": self.auto_minor_version_upgrade,
            "KmsKeyId": self.kms_key_id,
            "ReplicationInstanceArn": self.arn,
            "ReplicationInstancePublicIpAddress": self.public_ip_addresses[0]
            if self.public_ip_addresses
            else None,
            "ReplicationInstancePrivateIpAddress": self.private_ip_addresses[0]
            if self.private_ip_addresses
            else None,
            "ReplicationInstancePublicIpAddresses": self.public_ip_addresses,
            "ReplicationInstancePrivateIpAddresses": self.private_ip_addresses,
            "ReplicationInstanceIpv6Addresses": self.ipv6_addresses,
            "PubliclyAccessible": self.publicly_accessible,
            "SecondaryAvailabilityZone": f"{self.availability_zone}b"
            if self.multi_az
            else None,
            "FreeUntil": None,
            "DnsNameServers": self.dns_name_servers,
            "NetworkType": self.network_type,
            "KerberosAuthenticationSettings": kerberos_settings,
        }

    def delete(self) -> "FakeReplicationInstance":
        self.status = "deleting"
        return self


class FakeReplicationSubnetGroup(BaseModel):
    def __init__(
        self,
        replication_subnet_group_identifier: str,
        replication_subnet_group_description: str,
        subnetIds: Optional[list[str]] = None,
        tags: Optional[list[dict[str, str]]] = None,
    ):
        self.replication_subnet_group_identifier = replication_subnet_group_identifier
        self.description = replication_subnet_group_description
        self.subnetIds = subnetIds
        self.tags = tags or []

    def to_dict(self) -> dict[str, Any]:
        return {
            "ReplicationSubnetGroupIdentifier": self.replication_subnet_group_identifier,
            "ReplicationSubnetGroupDescription": self.description,
            "VpcId": "vpc-12345",
            "SubnetGroupStatus": "Complete",
            "Subnets": [
                {
                    "SubnetIdentifier": "subnet-12345",
                    "SubnetAvailabilityZone": {"Name": "us-east-1a"},
                    "SubnetStatus": "Active",
                }
            ],
            "SupportedNetworkTypes": ["IPV4"],
        }


class FakeConnection(ManagedState):
    def __init__(
        self,
        replication_instance_arn: str,
        endpoint_arn: str,
        replication_instance_identifier: str,
        endpoint_identifier: str,
    ):
        ManagedState.__init__(
            self,
            model_name="dms::connection",
            transitions=[("testing", "successful")],
        )

        self.replication_instance_arn = replication_instance_arn
        self.endpoint_arn = endpoint_arn
        self.replication_instance_identifier = replication_instance_identifier
        self.endpoint_identifier = endpoint_identifier

    def to_dict(self) -> dict[str, Any]:
        return {
            "ReplicationInstanceArn": self.replication_instance_arn,
            "EndpointArn": self.endpoint_arn,
            "Status": self.status,
            "LastFailureMessage": "",
            "EndpointIdentifier": self.endpoint_identifier,
            "ReplicationInstanceIdentifier": self.replication_instance_identifier,
        }


dms_backends = BackendDict(DatabaseMigrationServiceBackend, "dms")
