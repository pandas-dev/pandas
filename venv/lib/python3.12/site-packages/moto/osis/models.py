"""OpenSearchIngestionBackend class with methods for supported APIs."""

from datetime import datetime
from typing import TYPE_CHECKING, Any, ClassVar, Dict, List, Optional

import yaml

from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.common_models import BaseModel
from moto.moto_api._internal import mock_random as random
from moto.moto_api._internal.managed_state_model import ManagedState
from moto.utilities.paginator import paginate
from moto.utilities.tagging_service import TaggingService
from moto.utilities.utils import get_partition

from .exceptions import (
    InvalidVPCOptionsException,
    PipelineAlreadyExistsException,
    PipelineInvalidStateException,
    PipelineNotFoundException,
    SecurityGroupNotFoundException,
    SubnetNotFoundException,
)

if TYPE_CHECKING:
    from moto.ec2.models import EC2Backend


class Pipeline(ManagedState, BaseModel):
    CREATING_REASON = "The pipeline is being created. It is not able to ingest data."
    ACTIVE_REASON = "The pipeline is ready to ingest data."
    DELETING_REASON = "The pipeline is being deleted"
    STOPPING_REASON = "The pipeline is being stopped"
    STOPPED_REASON = "The pipeline is stopped"
    STARTING_REASON = "The pipeline is starting. It is not able to ingest data"
    UPDATING_REASON = "An update was triggered for the pipeline. It is still available to ingest data."

    STATUS_REASON_MAP: ClassVar[Dict[str, str]] = {
        "CREATING": CREATING_REASON,
        "ACTIVE": ACTIVE_REASON,
        "STOPPING": STOPPING_REASON,
        "STOPPED": STOPPED_REASON,
        "STARTING": STARTING_REASON,
        "UPDATING": UPDATING_REASON,
        "DELETING": DELETING_REASON,
    }

    def __init__(
        self,
        pipeline_name: str,
        account_id: str,
        region: str,
        min_units: int,
        max_units: int,
        pipeline_configuration_body: str,
        log_publishing_options: Optional[Dict[str, Any]],
        vpc_options: Optional[Dict[str, Any]],
        buffer_options: Optional[Dict[str, Any]],
        encryption_at_rest_options: Optional[Dict[str, Any]],
        ingest_endpoint_urls: List[str],
        serverless: bool,
        vpc_endpoint_service: Optional[str],
        vpc_endpoint: Optional[str],
        vpc_id: Optional[str],
        backend: "OpenSearchIngestionBackend",
    ):
        ManagedState.__init__(
            self,
            model_name="osis::pipeline",
            transitions=[
                ("CREATING", "ACTIVE"),
                ("UPDATING", "ACTIVE"),
                ("DELETING", "DELETED"),
                ("STOPPING", "STOPPED"),
                ("STARTING", "ACTIVE"),
            ],
        )

        self.pipeline_name = pipeline_name
        self.account_id = account_id
        self.region = region
        self.min_units = min_units
        self.max_units = max_units
        self.pipeline_configuration_body_str = pipeline_configuration_body
        self.pipeline_configuration_body = yaml.safe_load(pipeline_configuration_body)
        self.log_publishing_options = log_publishing_options
        self.vpc_options = vpc_options
        self.buffer_options = buffer_options
        self.encryption_at_rest_options = encryption_at_rest_options
        self.ingest_endpoint_urls = ingest_endpoint_urls
        self.serverless = serverless
        self.vpc_endpoint_service = vpc_endpoint_service
        self.vpc_endpoint = vpc_endpoint
        self.vpc_id = vpc_id
        self.backend = backend

        self.status = "CREATING"
        self.arn = self._get_arn(self.pipeline_name)
        self.destinations = self._update_destinations()

        if (
            self.vpc_options is None
            or self.vpc_options.get("VpcEndpointManagement", "SERVICE") == "SERVICE"
        ):
            # Not returned in this case
            self.vpc_endpoint_service = None

        self.service_vpc_endpoints = self._get_service_vpc_endpoints()
        self.created_at: datetime = datetime.now()
        self.last_updated_at: datetime = datetime.now()

    def _get_arn(self, name: str) -> str:
        return f"arn:{get_partition(self.region)}:osis:{self.region}:{self.account_id}:pipeline/{name}"

    def _get_service_vpc_endpoints(self) -> Optional[List[Dict[str, str]]]:
        # ServiceVpcEndpoint.VpcEndpointId not implemented
        if self.serverless:
            return [{"ServiceName": "OPENSEARCH_SERVERLESS"}]
        else:
            return None

    def _update_destinations(self) -> List[Dict[str, str]]:
        destinations = []
        for sub_pipeline in self.pipeline_configuration_body:
            if sub_pipeline != "version":
                for sink in self.pipeline_configuration_body[sub_pipeline]["sink"]:
                    for sink_type, sink_config in sink.items():
                        if sink_type == "opensearch":
                            if sink_config["aws"].get("serverless") is True:
                                service_name = "OpenSearch_Serverless"
                            else:
                                service_name = "OpenSearch"
                            endpoint = sink_config["hosts"][0]
                        elif sink_type == "s3":
                            service_name = "S3"
                            endpoint = sink_config["bucket"]
                        else:
                            continue
                        destinations.append(
                            {"ServiceName": service_name, "Endpoint": endpoint}
                        )
        return destinations

    @staticmethod
    def is_serverless(pipeline_body: Dict[str, Any]) -> bool:  # type: ignore[misc]
        serverless = False
        for sub_pipeline in pipeline_body:
            if sub_pipeline != "version":
                for sink in pipeline_body[sub_pipeline]["sink"]:
                    for _, sink_config in sink.items():
                        serverless = (
                            sink_config.get("aws", {}).get("serverless", False)
                            or serverless
                        )
                source_type = list(pipeline_body[sub_pipeline]["source"].keys())[0]
                source_config = pipeline_body[sub_pipeline]["source"][source_type]
                serverless = (
                    source_config.get("aws", {}).get("serverless", False) or serverless
                )
        return serverless

    def delete(self) -> None:
        self.status = "DELETING"
        self.set_last_updated()

    def get_created_at(self) -> str:
        return self.created_at.astimezone().isoformat()

    def get_last_updated_at(self) -> str:
        return self.last_updated_at.astimezone().isoformat()

    def set_last_updated(self) -> None:
        self.last_updated_at = datetime.now()

    def start(self) -> None:
        self.status = "STARTING"
        self.set_last_updated()

    def stop(self) -> None:
        self.status = "STOPPING"
        self.set_last_updated()

    def update(
        self,
        min_units: Optional[int],
        max_units: Optional[int],
        pipeline_configuration_body: Optional[str],
        log_publishing_options: Optional[Dict[str, Any]],
        buffer_options: Optional[Dict[str, Any]],
        encryption_at_rest_options: Optional[Dict[str, Any]],
    ) -> None:
        if min_units is not None:
            self.min_units = min_units
        if max_units is not None:
            self.max_units = max_units
        if pipeline_configuration_body is not None:
            self.pipeline_configuration_body_str = pipeline_configuration_body
            self.pipeline_configuration_body = yaml.safe_load(
                pipeline_configuration_body
            )
        if log_publishing_options is not None:
            self.log_publishing_options = log_publishing_options
        if buffer_options is not None:
            self.buffer_options = buffer_options
        if encryption_at_rest_options is not None:
            self.encryption_at_rest_options = encryption_at_rest_options
        self.destinations = self._update_destinations()
        self.serverless = self.is_serverless(self.pipeline_configuration_body)
        self.service_vpc_endpoints = self._get_service_vpc_endpoints()
        self.status = "UPDATING"
        self.set_last_updated()

    def to_dict(self) -> Dict[str, Any]:
        return {
            "PipelineName": self.pipeline_name,
            "PipelineArn": self.arn,
            "MinUnits": self.min_units,
            "MaxUnits": self.max_units,
            "Status": self.status,
            "StatusReason": {
                "Description": self.STATUS_REASON_MAP.get(self.status or "", ""),
            },
            "PipelineConfigurationBody": self.pipeline_configuration_body_str,
            "CreatedAt": self.get_created_at(),
            "LastUpdatedAt": self.get_last_updated_at(),
            "IngestEndpointUrls": self.ingest_endpoint_urls,
            "LogPublishingOptions": self.log_publishing_options,
            "VpcEndpoints": None
            if self.vpc_options is None
            else [
                {
                    "VpcEndpointId": self.vpc_endpoint,
                    "VpcId": self.vpc_id,
                    "VpcOptions": self.vpc_options,
                }
            ],
            "BufferOptions": self.buffer_options,
            "EncryptionAtRestOptions": self.encryption_at_rest_options,
            "VpcEndpointService": self.vpc_endpoint_service,
            "ServiceVpcEndpoints": self.service_vpc_endpoints,
            "Destinations": self.destinations,
            "Tags": self.backend.list_tags_for_resource(self.arn)["Tags"],
        }

    def to_short_dict(self) -> Dict[str, Any]:
        return {
            "Status": self.status,
            "StatusReason": {
                "Description": self.STATUS_REASON_MAP.get(self.status or "", ""),
            },
            "PipelineName": self.pipeline_name,
            "PipelineArn": self.arn,
            "MinUnits": self.min_units,
            "MaxUnits": self.max_units,
            "CreatedAt": self.get_created_at(),
            "LastUpdatedAt": self.get_last_updated_at(),
            "Destinations": self.destinations,
            "Tags": self.backend.list_tags_for_resource(self.arn)["Tags"],
        }


class OpenSearchIngestionBackend(BaseBackend):
    """Implementation of OpenSearchIngestion APIs."""

    PAGINATION_MODEL = {
        "list_pipelines": {
            "input_token": "next_token",
            "limit_key": "max_results",
            "limit_default": 100,
            "unique_attribute": "PipelineName",
        },
    }

    PIPELINE_DELETE_VALID_STATES = [
        "UPDATE_FAILED",
        "ACTIVE",
        "START_FAILED",
        "STOPPED",
        "CREATE_FAILED",
    ]
    PIPELINE_STOP_VALID_STATES = ["UPDATE_FAILED", "ACTIVE"]
    PIPELINE_START_VALID_STATES = ["START_FAILED", "STOPPED"]
    PIPELINE_UPDATE_VALID_STATES = [
        "UPDATE_FAILED",
        "ACTIVE",
        "START_FAILED",
        "STOPPED",
    ]

    def __init__(self, region_name: str, account_id: str):
        super().__init__(region_name, account_id)
        self._pipelines: Dict[str, Pipeline] = dict()
        self.tagger = TaggingService()

    @property
    def ec2_backend(self) -> "EC2Backend":  # type: ignore[misc]
        from moto.ec2 import ec2_backends

        return ec2_backends[self.account_id][self.region_name]

    @property
    def pipelines(self) -> Dict[str, Pipeline]:
        self._pipelines = {
            name: pipeline
            for name, pipeline in self._pipelines.items()
            if pipeline.status != "DELETED"
        }
        return self._pipelines

    def _get_ingest_endpoint_urls(
        self, pipeline_name: str, endpoint_random_string: str
    ) -> List[str]:
        return [
            f"{pipeline_name}-{endpoint_random_string}.{self.region_name}.osis.amazonaws.com"
        ]

    def _get_random_endpoint_string(self) -> str:
        return random.get_random_string(length=26, lower_case=True)

    def _get_vpc_endpoint(
        self, vpc_id: str, vpc_options: Dict[str, Any], service_name: str
    ) -> Optional[str]:
        if vpc_options.get("VpcEndpointManagement", "SERVICE") == "SERVICE":
            service_managed_endpoint = self.ec2_backend.create_vpc_endpoint(
                vpc_id=vpc_id,
                service_name=service_name,
                endpoint_type="Interface",
                security_group_ids=vpc_options.get("SecurityGroupIds"),
                subnet_ids=vpc_options["SubnetIds"],
                private_dns_enabled=False,
                policy_document="OSIS Test Doc",
                route_table_ids=[],
                tags={"OSISManaged": "true"},
            )
            return service_managed_endpoint.id
        else:
            return None

    def _get_vpc_endpoint_service(
        self, pipeline_name: str, endpoint_random_string: str
    ) -> str:
        return f"com.amazonaws.osis.{self.region_name}.{pipeline_name}-{endpoint_random_string}"

    def _validate_and_get_vpc(self, vpc_options: Dict[str, Any]) -> str:
        from moto.ec2.exceptions import InvalidSubnetIdError

        vpc_id = ""
        for subnet_id in vpc_options["SubnetIds"]:
            try:
                subnet = self.ec2_backend.get_subnet(subnet_id)
            except InvalidSubnetIdError:
                # re-raising for more accurate error message
                raise SubnetNotFoundException(subnet_id)
            if vpc_id == "":
                vpc_id = subnet.vpc_id
            else:
                if subnet.vpc_id != vpc_id:
                    raise InvalidVPCOptionsException(
                        "All specified subnets must belong to the same VPC."
                    )

        for sg_id in vpc_options["SecurityGroupIds"]:
            sg = self.ec2_backend.get_security_group_from_id(sg_id)
            if sg is None:
                raise SecurityGroupNotFoundException(sg_id)

        return vpc_id

    def create_pipeline(
        self,
        pipeline_name: str,
        min_units: int,
        max_units: int,
        pipeline_configuration_body: str,
        log_publishing_options: Optional[Dict[str, Any]],
        vpc_options: Optional[Dict[str, Any]],
        buffer_options: Optional[Dict[str, bool]],
        encryption_at_rest_options: Optional[Dict[str, Any]],
        tags: List[Dict[str, str]],
    ) -> Pipeline:
        if pipeline_name in self.pipelines:
            raise PipelineAlreadyExistsException(pipeline_name)

        serverless = Pipeline.is_serverless(yaml.safe_load(pipeline_configuration_body))

        endpoint_random_string = self._get_random_endpoint_string()
        endpoint_service = self._get_vpc_endpoint_service(
            pipeline_name, endpoint_random_string
        )

        ingestion_endpoint_urls = self._get_ingest_endpoint_urls(
            pipeline_name, endpoint_random_string
        )

        vpc_endpoint = None
        vpc_id = None
        if vpc_options is not None:
            vpc_id = self._validate_and_get_vpc(vpc_options)
            vpc_endpoint = self._get_vpc_endpoint(vpc_id, vpc_options, endpoint_service)

        pipeline = Pipeline(
            pipeline_name,
            self.account_id,
            self.region_name,
            min_units,
            max_units,
            pipeline_configuration_body,
            log_publishing_options,
            vpc_options,
            buffer_options,
            encryption_at_rest_options,
            ingestion_endpoint_urls,
            serverless,
            endpoint_service,
            vpc_endpoint,
            vpc_id,
            backend=self,
        )
        self.pipelines[pipeline_name] = pipeline
        self.tag_resource(pipeline.arn, tags)
        return pipeline

    def delete_pipeline(self, pipeline_name: str) -> None:
        if pipeline_name not in self.pipelines:
            raise PipelineNotFoundException(pipeline_name)
        pipeline = self.pipelines[pipeline_name]
        if pipeline.status not in self.PIPELINE_DELETE_VALID_STATES:
            raise PipelineInvalidStateException(
                "deletion", self.PIPELINE_DELETE_VALID_STATES, pipeline.status
            )
        pipeline.delete()

    def start_pipeline(self, pipeline_name: str) -> Pipeline:
        if pipeline_name not in self.pipelines:
            raise PipelineNotFoundException(pipeline_name)
        pipeline = self.pipelines[pipeline_name]
        if pipeline.status not in self.PIPELINE_START_VALID_STATES:
            raise PipelineInvalidStateException(
                "starting", self.PIPELINE_START_VALID_STATES, pipeline.status
            )
        pipeline.start()
        return pipeline

    def stop_pipeline(self, pipeline_name: str) -> Pipeline:
        if pipeline_name not in self.pipelines:
            raise PipelineNotFoundException(pipeline_name)
        pipeline = self.pipelines[pipeline_name]
        if pipeline.status not in self.PIPELINE_STOP_VALID_STATES:
            raise PipelineInvalidStateException(
                "stopping", self.PIPELINE_STOP_VALID_STATES, pipeline.status
            )
        pipeline.stop()
        return pipeline

    def get_pipeline(self, pipeline_name: str) -> Pipeline:
        if pipeline_name not in self.pipelines:
            raise PipelineNotFoundException(pipeline_name)
        pipeline = self.pipelines[pipeline_name]
        pipeline.advance()
        return pipeline

    @paginate(pagination_model=PAGINATION_MODEL)  # type: ignore
    def list_pipelines(self) -> List[Pipeline]:
        for pipeline in self.pipelines.values():
            pipeline.advance()
        return [p for p in self.pipelines.values()]

    def list_tags_for_resource(self, arn: str) -> Dict[str, List[Dict[str, str]]]:
        return self.tagger.list_tags_for_resource(arn)

    def update_pipeline(
        self,
        pipeline_name: str,
        min_units: Optional[int],
        max_units: Optional[int],
        pipeline_configuration_body: Optional[str],
        log_publishing_options: Optional[Dict[str, Any]],
        buffer_options: Optional[Dict[str, Any]],
        encryption_at_rest_options: Optional[Dict[str, Any]],
    ) -> Pipeline:
        if pipeline_name not in self.pipelines:
            raise PipelineNotFoundException(pipeline_name)
        pipeline = self.pipelines[pipeline_name]
        if pipeline.status not in self.PIPELINE_UPDATE_VALID_STATES:
            raise PipelineInvalidStateException(
                "updates", self.PIPELINE_UPDATE_VALID_STATES, pipeline.status
            )
        pipeline.update(
            min_units,
            max_units,
            pipeline_configuration_body,
            log_publishing_options,
            buffer_options,
            encryption_at_rest_options,
        )
        return pipeline

    def tag_resource(self, arn: str, tags: List[Dict[str, str]]) -> None:
        self.tagger.tag_resource(arn, tags)

    def untag_resource(self, arn: str, tag_keys: List[str]) -> None:
        self.tagger.untag_resource_using_names(arn, tag_keys)


osis_backends = BackendDict(OpenSearchIngestionBackend, "osis")
