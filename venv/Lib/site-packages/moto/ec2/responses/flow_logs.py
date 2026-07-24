from moto.core.responses import ActionResult, EmptyResult
from moto.ec2.models import validate_resource_ids

from ._base_response import EC2BaseResponse


class FlowLogs(EC2BaseResponse):
    def create_flow_logs(self) -> ActionResult:
        resource_type = self._get_param("ResourceType")
        resource_ids = self._get_param("ResourceIds", [])
        traffic_type = self._get_param("TrafficType")
        deliver_logs_permission_arn = self._get_param("DeliverLogsPermissionArn")
        log_destination_type = self._get_param("LogDestinationType")
        log_destination = self._get_param("LogDestination")
        log_group_name = self._get_param("LogGroupName")
        log_format = self._get_param("LogFormat")
        max_aggregation_interval = self._get_int_param("MaxAggregationInterval")
        validate_resource_ids(resource_ids)

        tags = self._parse_tag_specification().get("vpc-flow-log", {})

        self.error_on_dryrun()

        flow_logs, errors = self.ec2_backend.create_flow_logs(
            resource_type=resource_type,
            resource_ids=resource_ids,
            traffic_type=traffic_type,
            deliver_logs_permission_arn=deliver_logs_permission_arn,
            log_destination_type=log_destination_type,
            log_destination=log_destination,
            log_group_name=log_group_name,
            log_format=log_format,
            max_aggregation_interval=max_aggregation_interval,
        )
        for fl in flow_logs:
            fl.add_tags(tags)

        unsuccessful = [
            {
                "ResourceId": error[0],
                "Error": {
                    "Code": error[1],
                    "Message": error[2],
                },
            }
            for error in errors
        ]

        result = {
            "FlowLogIds": [fl.id for fl in flow_logs],
            "Unsuccessful": unsuccessful,
        }
        return ActionResult(result)

    def describe_flow_logs(self) -> ActionResult:
        flow_log_ids = self._get_param("FlowLogIds", [])
        filters = self._filters_from_querystring()
        flow_logs = self.ec2_backend.describe_flow_logs(flow_log_ids, filters)

        self.error_on_dryrun()

        result = {"FlowLogs": flow_logs}
        return ActionResult(result)

    def delete_flow_logs(self) -> ActionResult:
        flow_log_ids = self._get_param("FlowLogIds", [])

        self.error_on_dryrun()

        self.ec2_backend.delete_flow_logs(flow_log_ids)
        return EmptyResult()
