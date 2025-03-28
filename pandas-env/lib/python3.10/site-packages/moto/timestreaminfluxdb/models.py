"""TimestreamInfluxDBBackend class with methods for supported APIs."""

from enum import Enum
from typing import Any, Dict, List, Optional

from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.common_models import BaseModel
from moto.utilities.tagging_service import TaggingService

from .exceptions import (
    ConflictException,
    ResourceNotFoundException,
    ValidationException,
)
from .utils import random_id, validate_name


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
