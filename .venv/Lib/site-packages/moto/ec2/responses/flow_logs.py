from moto.ec2.models import validate_resource_ids

from ._base_response import EC2BaseResponse


class FlowLogs(EC2BaseResponse):
    def create_flow_logs(self) -> str:
        resource_type = self._get_param("ResourceType")
        resource_ids = self._get_multi_param("ResourceId")
        traffic_type = self._get_param("TrafficType")
        deliver_logs_permission_arn = self._get_param("DeliverLogsPermissionArn")
        log_destination_type = self._get_param("LogDestinationType")
        log_destination = self._get_param("LogDestination")
        log_group_name = self._get_param("LogGroupName")
        log_format = self._get_param("LogFormat")
        max_aggregation_interval = self._get_param("MaxAggregationInterval")
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
        template = self.response_template(CREATE_FLOW_LOGS_RESPONSE)
        return template.render(flow_logs=flow_logs, errors=errors)

    def describe_flow_logs(self) -> str:
        flow_log_ids = self._get_multi_param("FlowLogId")
        filters = self._filters_from_querystring()
        flow_logs = self.ec2_backend.describe_flow_logs(flow_log_ids, filters)

        self.error_on_dryrun()

        template = self.response_template(DESCRIBE_FLOW_LOGS_RESPONSE)
        return template.render(flow_logs=flow_logs)

    def delete_flow_logs(self) -> str:
        flow_log_ids = self._get_multi_param("FlowLogId")

        self.error_on_dryrun()

        self.ec2_backend.delete_flow_logs(flow_log_ids)
        return self.response_template(DELETE_FLOW_LOGS_RESPONSE).render()


CREATE_FLOW_LOGS_RESPONSE = """
<CreateFlowLogsResponse xmlns="http://ec2.amazonaws.com/doc/2016-11-15/">
  <requestId>2d96dae3-504b-4fc4-bf50-266EXAMPLE</requestId>
  <unsuccessful>
    {% for error in errors %}
      <item>
        <error>
          <code>{{ error.1 }}</code>
          <message>{{ error.2 }}</message>
        </error>
        <resourceId>{{ error.0 }}</resourceId>
      </item>
    {% endfor %}
  </unsuccessful>
  <flowLogIdSet>
    {% for flow_log in flow_logs %}
      <item>{{ flow_log.id }}</item>
    {% endfor %}
  </flowLogIdSet>
</CreateFlowLogsResponse>"""

DELETE_FLOW_LOGS_RESPONSE = """
<DeleteFlowLogsResponse xmlns="http://ec2.amazonaws.com/doc/2016-11-15/">
  <requestId>c5c4f51f-f4e9-42bc-8700-EXAMPLE</requestId>
  <unsuccessful/>
</DeleteFlowLogsResponse>"""

DESCRIBE_FLOW_LOGS_RESPONSE = """
<DescribeFlowLogsResponse xmlns="http://ec2.amazonaws.com/doc/2016-11-15/">
  <requestId>3cb46f23-099e-4bf0-891c-EXAMPLE</requestId>
  <flowLogSet>
  {% for flow_log in flow_logs %}
    <item>
      {% if flow_log.log_destination is not none %}
        <logDestination>{{ flow_log.log_destination }}</logDestination>
      {% endif %}
      <resourceId>{{ flow_log.resource_id }}</resourceId>
      <logDestinationType>{{ flow_log.log_destination_type }}</logDestinationType>
      <creationTime>{{ flow_log.created_at }}</creationTime>
      <trafficType>{{ flow_log.traffic_type }}</trafficType>
      <deliverLogsStatus>{{ flow_log.deliver_logs_status }}</deliverLogsStatus>
      {% if flow_log.deliver_logs_error_message is not none %}
        <deliverLogsErrorMessage>{{ flow_log.deliver_logs_error_message }}</deliverLogsErrorMessage>
      {% endif %}
      <logFormat>{{ flow_log.log_format }}</logFormat>
      <flowLogStatus>ACTIVE</flowLogStatus>
      <flowLogId>{{ flow_log.id }}</flowLogId>
      <maxAggregationInterval>{{ flow_log.max_aggregation_interval }}</maxAggregationInterval>
      {% if flow_log.deliver_logs_permission_arn is not none %}
        <deliverLogsPermissionArn>{{ flow_log.deliver_logs_permission_arn }}</deliverLogsPermissionArn>
      {% endif %}
      {% if flow_log.log_group_name is not none %}
        <logGroupName>{{ flow_log.log_group_name }}</logGroupName>
      {% endif %}
      {% if flow_log.get_tags() %}
        <tagSet>
          {% for tag in flow_log.get_tags() %}
            <item>
              <key>{{ tag.key }}</key>
              <value>{{ tag.value }}</value>
            </item>
          {% endfor %}
        </tagSet>
      {% endif %}
    </item>
  {% endfor %}
  </flowLogSet>
</DescribeFlowLogsResponse>"""
