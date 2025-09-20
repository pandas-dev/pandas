import itertools
from typing import Any, Dict, List, Optional, Tuple

from moto.core.common_models import CloudFormationModel

from ..exceptions import (
    FlowLogAlreadyExists,
    InvalidAggregationIntervalParameterError,
    InvalidDependantParameterError,
    InvalidDependantParameterTypeError,
    InvalidFlowLogIdError,
)
from ..utils import (
    generic_filter,
    random_flow_log_id,
    utc_date_and_time,
)
from .core import TaggedEC2Resource


class FlowLogs(TaggedEC2Resource, CloudFormationModel):
    def __init__(
        self,
        ec2_backend: Any,
        flow_log_id: str,
        resource_id: str,
        traffic_type: str,
        log_destination: str,
        log_group_name: str,
        deliver_logs_permission_arn: str,
        max_aggregation_interval: str,
        log_destination_type: str,
        log_format: str,
        deliver_logs_status: str = "SUCCESS",
        deliver_logs_error_message: Optional[str] = None,
    ):
        self.ec2_backend = ec2_backend
        self.id = flow_log_id
        self.resource_id = resource_id
        self.traffic_type = traffic_type
        self.log_destination = log_destination
        self.log_group_name = log_group_name
        self.deliver_logs_permission_arn = deliver_logs_permission_arn
        self.deliver_logs_status = deliver_logs_status
        self.deliver_logs_error_message = deliver_logs_error_message
        self.max_aggregation_interval = max_aggregation_interval
        self.log_destination_type = log_destination_type
        self.log_format = log_format

        self.created_at = utc_date_and_time()

    @staticmethod
    def cloudformation_name_type() -> str:
        return ""

    @staticmethod
    def cloudformation_type() -> str:
        # https://docs.aws.amazon.com/AWSCloudFormation/latest/UserGuide/aws-resource-ec2-flowlog.html
        return "AWS::EC2::FlowLog"

    @classmethod
    def create_from_cloudformation_json(  # type: ignore[misc]
        cls,
        resource_name: str,
        cloudformation_json: Any,
        account_id: str,
        region_name: str,
        **kwargs: Any,
    ) -> "FlowLogs":
        from ..models import ec2_backends

        properties = cloudformation_json["Properties"]

        resource_type = properties.get("ResourceType")
        resource_id = [properties.get("ResourceId")]
        traffic_type = properties.get("TrafficType")
        deliver_logs_permission_arn = properties.get("DeliverLogsPermissionArn")
        log_destination_type = properties.get("LogDestinationType")
        log_destination = properties.get("LogDestination")
        log_group_name = properties.get("LogGroupName")
        log_format = properties.get("LogFormat")
        max_aggregation_interval = properties.get("MaxAggregationInterval")

        ec2_backend = ec2_backends[account_id][region_name]
        flow_log, _ = ec2_backend.create_flow_logs(
            resource_type,
            resource_id,
            traffic_type,
            deliver_logs_permission_arn,
            log_destination_type,
            log_destination,
            log_group_name,
            log_format,
            max_aggregation_interval,
        )
        for tag in properties.get("Tags", []):
            tag_key = tag["Key"]
            tag_value = tag["Value"]
            flow_log[0].add_tag(tag_key, tag_value)

        return flow_log[0]

    @property
    def physical_resource_id(self) -> str:
        return self.id

    def get_filter_value(
        self, filter_name: str, method_name: Optional[str] = None
    ) -> Any:
        """
        API Version 2016-11-15 defines the following filters for DescribeFlowLogs:

        * deliver-log-status
        * log-destination-type
        * flow-log-id
        * log-group-name
        * resource-id
        * traffic-type
        * tag:key=value
        * tag-key

        Taken from: https://docs.aws.amazon.com/AWSEC2/latest/APIReference/API_DescribeFlowLogs.html
        """
        if filter_name == "resource-id":
            return self.resource_id
        elif filter_name == "traffic-type":
            return self.traffic_type
        elif filter_name == "log-destination-type":
            return self.log_destination_type
        elif filter_name == "flow-log-id":
            return self.id
        elif filter_name == "log-group-name":
            return self.log_group_name
        elif filter_name == "deliver-log-status":
            return "SUCCESS"
        else:
            return super().get_filter_value(filter_name, "DescribeFlowLogs")


class FlowLogsBackend:
    def __init__(self) -> None:
        self.flow_logs: Dict[str, FlowLogs] = {}

    def _validate_request(
        self,
        log_group_name: str,
        log_destination: str,
        log_destination_type: str,
        max_aggregation_interval: str,
        deliver_logs_permission_arn: str,
    ) -> None:
        if log_group_name is None and log_destination is None:
            raise InvalidDependantParameterError(
                "LogDestination", "LogGroupName", "not provided"
            )

        if log_destination_type == "s3":
            if log_group_name is not None:
                raise InvalidDependantParameterTypeError(
                    "LogDestination", "cloud-watch-logs", "LogGroupName"
                )
        elif log_destination_type == "cloud-watch-logs":
            if deliver_logs_permission_arn is None:
                raise InvalidDependantParameterError(
                    "DeliverLogsPermissionArn", "LogDestinationType", "cloud-watch-logs"
                )

        if max_aggregation_interval not in ["60", "600"]:
            raise InvalidAggregationIntervalParameterError(
                "Flow Log Max Aggregation Interval"
            )

    def create_flow_logs(
        self,
        resource_type: str,
        resource_ids: List[str],
        traffic_type: str,
        deliver_logs_permission_arn: str,
        log_destination_type: str,
        log_destination: str,
        log_group_name: str,
        log_format: str,
        max_aggregation_interval: str,
    ) -> Tuple[List[FlowLogs], List[Any]]:
        # Guess it's best to put it here due to possible
        # lack of them in the CloudFormation template
        max_aggregation_interval = (
            "600" if max_aggregation_interval is None else max_aggregation_interval
        )
        log_destination_type = (
            "cloud-watch-logs" if log_destination_type is None else log_destination_type
        )
        log_format = (
            "${version} ${account-id} ${interface-id} ${srcaddr} ${dstaddr} ${srcport} ${dstport} ${protocol} ${packets} ${bytes} ${start} ${end} ${action} ${log-status}"
            if log_format is None
            else log_format
        )

        # Validate the requests paremeters
        self._validate_request(
            log_group_name,
            log_destination,
            log_destination_type,
            max_aggregation_interval,
            deliver_logs_permission_arn,
        )

        flow_logs_set = []
        unsuccessful = []

        for resource_id in resource_ids:
            deliver_logs_status = "SUCCESS"
            deliver_logs_error_message = None
            flow_log_id = random_flow_log_id()
            if resource_type == "VPC":
                # Validate VPCs exist
                self.get_vpc(resource_id)  # type: ignore[attr-defined]
            elif resource_type == "Subnet":
                # Validate Subnets exist
                self.get_subnet(resource_id)  # type: ignore[attr-defined]
            elif resource_type == "NetworkInterface":
                # Validate NetworkInterfaces exist
                self.get_network_interface(resource_id)  # type: ignore[attr-defined]

            if log_destination_type == "s3":
                from moto.s3.exceptions import MissingBucket
                from moto.s3.models import s3_backends

                arn = log_destination.split(":", 5)[5]
                try:
                    s3_backends[self.account_id][self.partition].get_bucket(arn)  # type: ignore[attr-defined]
                except MissingBucket:
                    # Instead of creating FlowLog report
                    # the unsuccessful status for the
                    # given resource_id
                    unsuccessful.append(
                        (resource_id, "400", f"LogDestination: {arn} does not exist.")
                    )
                    continue
            elif log_destination_type == "cloud-watch-logs":
                from moto.logs.exceptions import ResourceNotFoundException
                from moto.logs.models import logs_backends

                # API allows to create a FlowLog with a
                # non-existing LogGroup. It however later
                # on reports the FAILED delivery status.
                try:
                    # Need something easy to check the group exists.
                    # The list_tags_log_group seems to do the trick.
                    logs = logs_backends[self.account_id][self.region_name]  # type: ignore[attr-defined]
                    logs.list_tags_log_group(log_group_name)
                except ResourceNotFoundException:
                    deliver_logs_status = "FAILED"
                    deliver_logs_error_message = "Access error"

            all_flow_logs = self.describe_flow_logs()
            if any(
                (fl.resource_id, fl.log_group_name, fl.log_destination)
                == (resource_id, log_group_name, log_destination)
                for fl in all_flow_logs
            ):
                raise FlowLogAlreadyExists()
            flow_logs = FlowLogs(
                self,
                flow_log_id,
                resource_id,
                traffic_type,
                log_destination,
                log_group_name,
                deliver_logs_permission_arn,
                max_aggregation_interval,
                log_destination_type,
                log_format,
                deliver_logs_status,
                deliver_logs_error_message,
            )
            self.flow_logs[flow_log_id] = flow_logs
            flow_logs_set.append(flow_logs)

        return flow_logs_set, unsuccessful

    def describe_flow_logs(
        self, flow_log_ids: Optional[List[str]] = None, filters: Any = None
    ) -> List[FlowLogs]:
        matches = list(itertools.chain([i for i in self.flow_logs.values()]))
        if flow_log_ids:
            matches = [flow_log for flow_log in matches if flow_log.id in flow_log_ids]
        if filters:
            matches = generic_filter(filters, matches)
        return matches

    def delete_flow_logs(self, flow_log_ids: List[str]) -> None:
        non_existing = []
        for flow_log in flow_log_ids:
            if flow_log in self.flow_logs:
                self.flow_logs.pop(flow_log, None)
            else:
                non_existing.append(flow_log)

        if non_existing:
            raise InvalidFlowLogIdError(
                len(flow_log_ids), " ".join(x for x in flow_log_ids)
            )
