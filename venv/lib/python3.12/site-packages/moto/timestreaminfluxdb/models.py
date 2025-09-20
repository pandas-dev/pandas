"""TimestreamInfluxDBBackend class with methods for supported APIs."""

from enum import Enum
from typing import Any, Dict, List, Optional, Tuple

from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.common_models import BaseModel
from moto.utilities.paginator import paginate
from moto.utilities.tagging_service import TaggingService

from .exceptions import (
    ConflictException,
    ResourceNotFoundException,
    ValidationException,
)
from .utils import random_id, validate_name

PAGINATION_MODEL = {
    "list_db_parameter_groups": {
        "input_token": "next_token",
        "limit_key": "max_results",
        "limit_default": 100,
        "unique_attribute": "id",
    },
    "list_db_clusters": {
        "input_token": "next_token",
        "limit_key": "max_results",
        "limit_default": 100,
        "unique_attribute": "id",
    },
}


class InstanceStatus(str, Enum):
    CREATING = "CREATING"
    AVAILABLE = "AVAILABLE"
    DELETING = "DELETING"
    MODIFYING = "MODIFYING"
    UPDATING = "UPDATING"
    DELETED = "DELETED"
    FAILED = "FAILED"
    UPDATING_DEPLOYMENT_TYPE = "UPDATING_DEPLOYMENT_TYPE"
    UPDATING_INSTANCE_TYPE = "UPDATING_INSTANCE_TYPE"


class NetworkType(str, Enum):
    IPV4 = "IPV4"
    DUAL = "DUAL"


class InstanceType(str, Enum):
    DB_INFLUX_MEDIUM = "db.influx.medium"
    DB_INFLUX_LARGE = "db.influx.large"
    DB_INFLUX_XLARGE = "db.influx.xlarge"
    DB_INFLUX_2XLARGE = "db.influx.2xlarge"
    DB_INFLUX_4XLARGE = "db.influx.4xlarge"
    DB_INFLUX_8XLARGE = "db.influx.8xlarge"
    DB_INFLUX_12XLARGE = "db.influx.12xlarge"
    DB_INFLUX_16XLARGE = "db.influx.16xlarge"


class DBStorageType(str, Enum):
    InfluxIOIncludedT1 = "InfluxIOIncludedT1"
    InfluxIOIncludedT2 = "InfluxIOIncludedT2"
    InfluxIOIncludedT3 = "InfluxIOIncludedT3"


class DeploymentType(str, Enum):
    SINGLE_AZ = "SINGLE_AZ"
    WITH_MULTIAZ_STANDBY = "WITH_MULTIAZ_STANDBY"


class ParameterGroup(BaseModel):
    def __init__(
        self,
        name: str,
        description: Optional[str] = None,
        parameters: Optional[Dict[str, Any]] = None,
        region_name: str = "",
        account_id: str = "",
    ):
        self.id = random_id()
        self.name = name
        self.description = description or ""
        self.parameters = parameters or {}
        self.arn = f"arn:aws:timestream-influxdb:{region_name}:{account_id}:db-parameter-group/{self.id}"

    def to_dict(self) -> Dict[str, Any]:
        return {
            "id": self.id,
            "name": self.name,
            "arn": self.arn,
            "description": self.description,
            "parameters": self.parameters,
        }

    def to_summary_dict(self) -> Dict[str, str]:
        return {
            "id": self.id,
            "name": self.name,
            "arn": self.arn,
            "description": self.description,
        }


class Cluster(BaseModel):
    def __init__(
        self,
        name: str,
        password: str,
        username: Optional[str] = None,
        organization: Optional[str] = None,
        bucket: Optional[str] = None,
        port: Optional[int] = None,
        db_parameter_group_identifier: Optional[str] = None,
        db_instance_type: Optional[str] = None,
        db_storage_type: Optional[str] = None,
        allocated_storage: Optional[int] = None,
        network_type: Optional[str] = None,
        publicly_accessible: Optional[bool] = None,
        vpc_subnet_ids: Optional[List[str]] = None,
        vpc_security_group_ids: Optional[List[str]] = None,
        deployment_type: Optional[str] = None,
        failover_mode: Optional[str] = None,
        log_delivery_configuration: Optional[Dict[str, Any]] = None,
        region_name: str = "",
        account_id: str = "",
        endpoint_id: str = "",
    ):
        self.id = random_id()
        self.name = name
        self.password = password
        self.username = username
        self.organization = organization
        self.bucket = bucket
        self.port = port or 8086
        self.db_parameter_group_identifier = db_parameter_group_identifier
        self.db_instance_type = db_instance_type or "db.influx.medium"
        self.db_storage_type = db_storage_type or DBStorageType.InfluxIOIncludedT1
        self.allocated_storage = allocated_storage or 100
        self.network_type = network_type or NetworkType.IPV4
        self.publicly_accessible = publicly_accessible or False
        self.vpc_subnet_ids = vpc_subnet_ids or ["subnet-default"]
        self.vpc_security_group_ids = vpc_security_group_ids or ["sg-default"]
        self.deployment_type = deployment_type or "MULTI_NODE_READ_REPLICAS"
        self.failover_mode = failover_mode or "AUTOMATIC"
        self.log_delivery_configuration = log_delivery_configuration or {}
        self.arn = f"arn:aws:timestream-influxdb:{region_name}:{account_id}:db-cluster/{self.id}"
        self.endpoint = (
            f"{self.id}-{endpoint_id}.timestream-influxdb.{region_name}.on.aws"
        )
        self.reader_endpoint = (
            f"{self.id}-{endpoint_id}.reader.timestream-influxdb.{region_name}.on.aws"
        )
        self.influx_auth_parameters_secret_arn = f"arn:aws:secretsmanager:{region_name}:{account_id}:secret:timestream-influxdb/{self.id}/auth-params-{random_id(6)}"
        self.status = "CREATING"

    def to_dict(self) -> Dict[str, Any]:
        return {
            "id": self.id,
            "name": self.name,
            "arn": self.arn,
            "status": self.status,
            "endpoint": self.endpoint,
            "readerEndpoint": self.reader_endpoint,
            "port": self.port,
            "deploymentType": self.deployment_type,
            "dbInstanceType": self.db_instance_type,
            "networkType": self.network_type,
            "dbStorageType": self.db_storage_type,
            "allocatedStorage": self.allocated_storage,
            "publiclyAccessible": self.publicly_accessible,
            "dbParameterGroupIdentifier": self.db_parameter_group_identifier,
            "logDeliveryConfiguration": self.log_delivery_configuration,
            "influxAuthParametersSecretArn": self.influx_auth_parameters_secret_arn,
            "vpcSubnetIds": self.vpc_subnet_ids,
            "vpcSecurityGroupIds": self.vpc_security_group_ids,
            "failoverMode": self.failover_mode,
        }

    def to_summary_dict(self) -> Dict[str, Any]:
        return {
            "id": self.id,
            "name": self.name,
            "arn": self.arn,
            "status": self.status,
            "endpoint": self.endpoint,
            "readerEndpoint": self.reader_endpoint,
            "port": self.port,
            "deploymentType": self.deployment_type,
            "dbInstanceType": self.db_instance_type,
            "networkType": self.network_type,
            "dbStorageType": self.db_storage_type,
            "allocatedStorage": self.allocated_storage,
        }


class DBInstance(BaseModel):
    def __init__(
        self,
        name: str,
        username: Optional[str],
        password: str,
        organization: str,
        bucket: str,
        dbInstanceType: str,
        vpcSubnetIds: List[str],
        vpcSecurityGroupIds: List[str],
        publiclyAccessible: bool,
        dbStorageType: str,
        allocatedStorage: int,
        dbParameterGroupIdentifier: Optional[str],
        deploymentType: str,
        logDeliveryConfiguration: Optional[Dict[str, Any]],
        tags: Optional[Dict[str, Any]],
        port: int,
        networkType: str,
        region_name: str,
        account_id: str,
        endpoint_id: str,
    ):
        # Generate a random id of size 10
        self.id = random_id()

        self.name = name
        self.username = username
        self.password = password
        self.organization = organization
        self.bucket = bucket
        self.db_instance_type = dbInstanceType
        self.vpc_subnet_ids = vpcSubnetIds
        self.vpc_security_group_ids = vpcSecurityGroupIds
        self.publicly_accessible = publiclyAccessible
        self.db_storage_type = dbStorageType
        self.allocated_storage = allocatedStorage
        self.db_parameter_group_id = dbParameterGroupIdentifier
        self.deployment_type = deploymentType
        self.log_delivery_configuration = logDeliveryConfiguration
        self.port = port
        self.network_type = networkType
        self.status = InstanceStatus.CREATING
        self.arn = f"arn:aws:timestream-influxdb:{region_name}:{account_id}:db-instance/{self.id}"
        self.endpoint = (
            f"{self.id}-{endpoint_id}.timestream-influxdb.{region_name}.on.aws"
        )
        # Before 12/09/2024, there was a different endpoint format.
        self.endpoint_old = (
            f"{self.name}-{endpoint_id}.timestream-influxdb.{region_name}.on.aws"
        )

        self.availability_zone = ""  # TODO implement this
        self.secondary_availability_zone = ""  # TODO implement this

    def to_dict(self) -> Dict[str, Any]:
        return {
            "id": self.id,
            "name": self.name,
            "arn": self.arn,
            "status": self.status,
            "endpoint": self.endpoint,
            "port": self.port,
            "networkType": self.network_type,
            "dbInstanceType": self.db_instance_type,
            "dbStorageType": self.db_storage_type,
            "allocatedStorage": self.allocated_storage,
            "deploymentType": self.deployment_type,
            "vpcSubnetIds": self.vpc_subnet_ids,
            "publiclyAccessible": self.publicly_accessible,
            "vpcSecurityGroupIds": self.vpc_security_group_ids,
            "dbParameterGroupIdentifier": self.db_parameter_group_id,  # TODO implement this
            "availabilityZone": self.availability_zone,  # TODO implement this
            "secondaryAvailabilityZone": self.secondary_availability_zone,  # TODO implement this
            "logDeliveryConfiguration": self.log_delivery_configuration,  # TODO implement this
            "influxAuthParametersSecretArn": "",  # TODO implement this
        }


class TimestreamInfluxDBBackend(BaseBackend):
    """Implementation of TimestreamInfluxDB APIs."""

    def __init__(self, region_name: str, account_id: str):
        super().__init__(region_name, account_id)

        # the endpoint identifier is unique per account and per region
        # https://docs.aws.amazon.com/timestream/latest/developerguide/timestream-for-influxdb.html
        self.endpoint_id: str = random_id(10)
        self.db_instances: Dict[str, DBInstance] = {}
        self.db_parameter_groups: Dict[str, ParameterGroup] = {}
        self.db_clusters: Dict[str, Cluster] = {}
        self.tagger = TaggingService()

    def create_db_instance(
        self,
        name: str,
        username: Optional[str],  # required if using InfluxDB UI though
        password: str,
        organization: str,
        bucket: str,
        db_instance_type: str,
        vpc_subnet_ids: List[str],
        vpc_security_group_ids: List[str],
        db_storage_type: str,
        publicly_accessible: bool,
        allocated_storage: int,
        db_parameter_group_identifier: str,
        deployment_type: str,
        log_delivery_configuration: Optional[Dict[str, Any]],
        tags: Optional[Dict[str, str]],
        port: int,
        network_type: str,
    ) -> DBInstance:
        """
        dbParameterGroupIdentifier argument is not yet handled
        deploymentType currently is auto set to 'SINGLE_AZ' if not passed in.
        publicAccessible is not yet handled
        logDeliveryConfiguration is not yet handled
        AvailabilityZone and SecondaryAvailabilityZone are not yet handled
        influxAuthParametersSecretArn is not yet handled
        """

        # Checks:
        for db_instance in self.db_instances.values():
            if db_instance.name == name:
                raise ConflictException(
                    f"A DB Instance with the name {name} already exists"
                )

        validate_name(name)

        if db_storage_type not in [t.value for t in DBStorageType]:
            raise ValidationException(f"Unknown DB storage type {db_storage_type}")

        if db_instance_type not in [t.value for t in InstanceType]:
            raise ValidationException(f"Unknown DB instance type {db_instance_type}")

        new_instance = DBInstance(
            name,
            username,
            password,
            organization,
            bucket,
            db_instance_type,
            vpc_subnet_ids,
            vpc_security_group_ids,
            publicly_accessible,
            db_storage_type,
            allocated_storage,
            db_parameter_group_identifier,
            deployment_type,
            log_delivery_configuration,
            tags,
            port,
            network_type,
            self.region_name,
            self.account_id,
            self.endpoint_id,
        )

        # add to the list
        self.db_instances[new_instance.id] = new_instance

        # add tags
        if tags:
            self.tag_resource(new_instance.arn, tags)

        return new_instance

    def delete_db_instance(self, id: str) -> DBInstance:
        if id not in self.db_instances:
            raise ResourceNotFoundException(f"DB Instance with id {id} not found")

        # mark as deleting
        self.db_instances[id].status = InstanceStatus.DELETING
        return self.db_instances.pop(id)

    def get_db_instance(self, id: str) -> DBInstance:
        if id not in self.db_instances:
            raise ResourceNotFoundException(f"DB Instance with id {id} not found")

        return self.db_instances[id]

    def list_db_instances(self) -> List[Dict[str, Any]]:
        """
        Pagination is not yet implemented
        """
        return [
            {
                "allocatedStorage": instance.allocated_storage,
                "arn": instance.arn,
                "dbInstanceType": instance.db_instance_type,
                "dbStorageType": instance.db_storage_type,
                "deploymentType": instance.deployment_type,
                "endpoint": instance.endpoint,
                "id": instance.id,
                "name": instance.name,
                "networkType": instance.network_type,
                "port": instance.port,
                "status": instance.status,
            }
            for instance in self.db_instances.values()
        ]

    def tag_resource(self, resource_arn: str, tags: Dict[str, str]) -> None:
        tag_list = self.tagger.convert_dict_to_tags_input(tags)
        errmsg = self.tagger.validate_tags(tag_list)
        if errmsg:
            raise ValidationException(errmsg)
        self.tagger.tag_resource(resource_arn, tag_list)

    def untag_resource(self, resource_arn: str, tag_keys: List[str]) -> None:
        self.tagger.untag_resource_using_names(resource_arn, tag_keys)

    def list_tags_for_resource(self, resource_arn: str) -> Dict[str, str]:
        return self.tagger.get_tag_dict_for_resource(resource_arn)

    def create_db_parameter_group(
        self,
        name: str,
        description: Optional[str] = None,
        parameters: Optional[Dict[str, Any]] = None,
        tags: Optional[Dict[str, str]] = None,
    ) -> ParameterGroup:
        validate_name(name)

        for param_group in self.db_parameter_groups.values():
            if param_group.name == name:
                raise ConflictException(
                    f"A DB parameter group with the name {name} already exists"
                )

        param_group = ParameterGroup(
            name=name,
            description=description,
            parameters=parameters,
            region_name=self.region_name,
            account_id=self.account_id,
        )

        self.db_parameter_groups[param_group.id] = param_group

        if tags:
            self.tag_resource(param_group.arn, tags)

        return param_group

    def get_db_parameter_group(self, identifier: str) -> ParameterGroup:
        param_group = self.db_parameter_groups.get(identifier)
        if not param_group:
            raise ResourceNotFoundException(
                f"DB parameter group with identifier {identifier} not found"
            )

        return param_group

    @paginate(pagination_model=PAGINATION_MODEL)
    def list_db_parameter_groups(self) -> List[Dict[str, str]]:
        if not self.db_parameter_groups:
            return []

        return [
            param_group.to_summary_dict()
            for param_group in self.db_parameter_groups.values()
        ]

    @paginate(pagination_model=PAGINATION_MODEL)
    def list_db_clusters(self) -> List[Dict[str, object]]:
        if not self.db_clusters:
            return []

        return [
            cluster.to_summary_dict()
            for cluster in self.db_clusters.values()
            if cluster.status != "DELETED"
        ]

    def get_db_cluster(self, db_cluster_id: str) -> Cluster:
        cluster = self.db_clusters.get(db_cluster_id)
        if not cluster:
            raise ResourceNotFoundException(
                f"DB cluster with ID {db_cluster_id} not found"
            )

        return cluster

    def create_db_cluster(
        self,
        name: str,
        password: str,
        username: Optional[str] = None,
        organization: Optional[str] = None,
        bucket: Optional[str] = None,
        port: Optional[int] = None,
        db_parameter_group_identifier: Optional[str] = None,
        db_instance_type: Optional[str] = None,
        db_storage_type: Optional[str] = None,
        allocated_storage: Optional[int] = None,
        network_type: Optional[str] = None,
        publicly_accessible: Optional[bool] = None,
        vpc_subnet_ids: Optional[List[str]] = None,
        vpc_security_group_ids: Optional[List[str]] = None,
        deployment_type: Optional[str] = None,
        failover_mode: Optional[str] = None,
        log_delivery_configuration: Optional[Dict[str, Any]] = None,
        tags: Optional[Dict[str, str]] = None,
    ) -> Tuple[str, str]:
        validate_name(name)

        for cluster in self.db_clusters.values():
            if cluster.name == name:
                raise ConflictException(
                    f"A DB cluster with the name {name} already exists"
                )

        new_cluster = Cluster(
            name=name,
            password=password,
            username=username,
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
            region_name=self.region_name,
            account_id=self.account_id,
            endpoint_id=self.endpoint_id,
        )

        new_cluster.status = "AVAILABLE"

        self.db_clusters[new_cluster.id] = new_cluster

        if tags:
            self.tag_resource(new_cluster.arn, tags)

        return (
            new_cluster.id,
            "AVAILABLE",
        )


timestreaminfluxdb_backends = BackendDict(
    TimestreamInfluxDBBackend,
    "timestream-influxdb",
    additional_regions=[
        "us-east-1",
        "us-east-2",
        "us-west-2",
        "eu-central-1",
        "eu-west-1",
        "ap-southeast-2",
        "ap-northeast-1",
    ],
)
