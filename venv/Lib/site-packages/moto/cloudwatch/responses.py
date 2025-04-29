import json
from typing import Dict, Iterable, List, Tuple, Union

from dateutil.parser import parse as dtparse

from moto.core.responses import BaseResponse

from .exceptions import InvalidParameterCombination, ValidationError
from .models import (
    CloudWatchBackend,
    Dimension,
    FakeAlarm,
    Metric,
    MetricDataQuery,
    MetricStat,
    cloudwatch_backends,
)

ERROR_RESPONSE = Tuple[str, Dict[str, int]]


class CloudWatchResponse(BaseResponse):
    def __init__(self) -> None:
        super().__init__(service_name="cloudwatch")

    @property
    def cloudwatch_backend(self) -> CloudWatchBackend:
        return cloudwatch_backends[self.current_account][self.region]

    def _error(self, code: str, message: str, status: int = 400) -> ERROR_RESPONSE:
        template = self.response_template(ERROR_RESPONSE_TEMPLATE)
        return template.render(code=code, message=message), dict(status=status)

    def put_metric_alarm(self) -> str:
        name = self._get_param("AlarmName")
        namespace = self._get_param("Namespace")
        metric_name = self._get_param("MetricName")
        metrics = self._get_multi_param("Metrics.member", skip_result_conversion=True)
        metric_data_queries = None
        if metrics:
            metric_data_queries = []
            for metric in metrics:
                metric_dimensions = []
                dims = (
                    metric.get("MetricStat", {})
                    .get("Metric", {})
                    .get("Dimensions.member", [])
                )
                for dim in dims:
                    metric_dimensions.append(
                        Dimension(name=dim.get("Name"), value=dim.get("Value"))
                    )
                metric_stat = None
                stat_metric_name = (
                    metric.get("MetricStat", {}).get("Metric", {}).get("MetricName")
                )
                if stat_metric_name:
                    stat_details = metric.get("MetricStat", {})
                    stat_metric_ns = stat_details.get("Metric", {}).get("Namespace")
                    metric_stat = MetricStat(
                        metric=Metric(
                            metric_name=stat_metric_name,
                            namespace=stat_metric_ns,
                            dimensions=metric_dimensions,
                        ),
                        period=stat_details.get("Period"),
                        stat=stat_details.get("Stat"),
                        unit=stat_details.get("Unit"),
                    )
                metric_data_queries.append(
                    MetricDataQuery(
                        query_id=metric.get("Id"),
                        label=metric.get("Label"),
                        period=metric.get("Period"),
                        return_data=metric.get("ReturnData"),
                        expression=metric.get("Expression"),
                        metric_stat=metric_stat,
                    )
                )

        comparison_operator = self._get_param("ComparisonOperator")
        evaluation_periods = self._get_param("EvaluationPeriods")
        datapoints_to_alarm = self._get_param("DatapointsToAlarm")
        period = self._get_param("Period")
        threshold = self._get_param("Threshold")
        statistic = self._get_param("Statistic")
        extended_statistic = self._get_param("ExtendedStatistic")
        description = self._get_param("AlarmDescription")
        dimensions = self._get_list_prefix("Dimensions.member")
        alarm_actions = self._get_multi_param("AlarmActions.member")
        ok_actions = self._get_multi_param("OKActions.member")
        actions_enabled = self._get_bool_param("ActionsEnabled")
        insufficient_data_actions = self._get_multi_param(
            "InsufficientDataActions.member"
        )
        unit = self._get_param("Unit")
        treat_missing_data = self._get_param("TreatMissingData")
        evaluate_low_sample_count_percentile = self._get_param(
            "EvaluateLowSampleCountPercentile"
        )
        threshold_metric_id = self._get_param("ThresholdMetricId")
        # fetch AlarmRule to re-use this method for composite alarms as well
        rule = self._get_param("AlarmRule")
        tags = self._get_multi_param("Tags.member")
        alarm = self.cloudwatch_backend.put_metric_alarm(
            name=name,
            namespace=namespace,
            metric_name=metric_name,
            metric_data_queries=metric_data_queries,
            comparison_operator=comparison_operator,
            evaluation_periods=evaluation_periods,
            datapoints_to_alarm=datapoints_to_alarm,
            period=period,
            threshold=threshold,
            statistic=statistic,
            extended_statistic=extended_statistic,
            description=description,
            dimensions=dimensions,
            alarm_actions=alarm_actions,
            ok_actions=ok_actions,
            insufficient_data_actions=insufficient_data_actions,
            unit=unit,
            actions_enabled=actions_enabled,
            treat_missing_data=treat_missing_data,
            evaluate_low_sample_count_percentile=evaluate_low_sample_count_percentile,
            threshold_metric_id=threshold_metric_id,
            rule=rule,
            tags=tags,
        )
        template = self.response_template(PUT_METRIC_ALARM_TEMPLATE)
        return template.render(alarm=alarm)

    def describe_alarms(self) -> str:
        action_prefix = self._get_param("ActionPrefix")
        alarm_name_prefix = self._get_param("AlarmNamePrefix")
        alarm_names = self._get_multi_param("AlarmNames.member")
        state_value = self._get_param("StateValue")

        if action_prefix:
            alarms = self.cloudwatch_backend.get_alarms_by_action_prefix(action_prefix)
        elif alarm_name_prefix:
            alarms = self.cloudwatch_backend.get_alarms_by_alarm_name_prefix(
                alarm_name_prefix
            )
        elif alarm_names:
            alarms = self.cloudwatch_backend.get_alarms_by_alarm_names(alarm_names)
        elif state_value:
            alarms = self.cloudwatch_backend.get_alarms_by_state_value(state_value)
        else:
            alarms = self.cloudwatch_backend.describe_alarms()

        metric_alarms = [a for a in alarms if a.rule is None]
        composite_alarms = [a for a in alarms if a.rule is not None]

        template = self.response_template(DESCRIBE_ALARMS_TEMPLATE)
        return template.render(
            metric_alarms=metric_alarms, composite_alarms=composite_alarms
        )

    def delete_alarms(self) -> str:
        alarm_names = self._get_multi_param("AlarmNames.member")
        self.cloudwatch_backend.delete_alarms(alarm_names)
        template = self.response_template(DELETE_METRIC_ALARMS_TEMPLATE)
        return template.render()

    def put_metric_data(self) -> str:
        namespace = self._get_param("Namespace")
        metric_data = self._get_multi_param("MetricData.member")
        self.cloudwatch_backend.put_metric_data(namespace, metric_data)
        template = self.response_template(PUT_METRIC_DATA_TEMPLATE)
        return template.render()

    def get_metric_data(self) -> str:
        params = self._get_params()
        start = dtparse(params["StartTime"])
        end = dtparse(params["EndTime"])
        scan_by = params.get("ScanBy") or "TimestampDescending"

        queries = params.get("MetricDataQueries", [])
        for query in queries:
            if "MetricStat" not in query and "Expression" not in query:
                # AWS also returns the empty line
                raise ValidationError(
                    "The parameter MetricDataQueries.member.1.MetricStat is required.\n"
                )
        results = self.cloudwatch_backend.get_metric_data(
            start_time=start, end_time=end, queries=queries, scan_by=scan_by
        )

        template = self.response_template(GET_METRIC_DATA_TEMPLATE)
        return template.render(results=results)

    def get_metric_statistics(self) -> str:
        namespace = self._get_param("Namespace")
        metric_name = self._get_param("MetricName")
        start_time = dtparse(self._get_param("StartTime"))
        end_time = dtparse(self._get_param("EndTime"))
        period = int(self._get_param("Period"))
        statistics = self._get_multi_param("Statistics.member")
        dimensions = self._get_multi_param("Dimensions.member")

        # Unsupported Parameters (To Be Implemented)
        unit = self._get_param("Unit")
        extended_statistics = self._get_param("ExtendedStatistics")

        if not statistics and not extended_statistics:
            raise InvalidParameterCombination(
                "Must specify either Statistics or ExtendedStatistics"
            )

        datapoints = self.cloudwatch_backend.get_metric_statistics(
            namespace,
            metric_name,
            start_time,
            end_time,
            period,
            statistics,
            unit=unit,
            dimensions=dimensions,
        )
        template = self.response_template(GET_METRIC_STATISTICS_TEMPLATE)
        return template.render(label=metric_name, datapoints=datapoints)

    def list_metrics(self) -> str:
        namespace = self._get_param("Namespace")
        metric_name = self._get_param("MetricName")
        dimensions = self._get_params().get("Dimensions", [])
        next_token = self._get_param("NextToken")
        next_token, metrics = self.cloudwatch_backend.list_metrics(
            next_token, namespace, metric_name, dimensions
        )
        template = self.response_template(LIST_METRICS_TEMPLATE)
        return template.render(metrics=metrics, next_token=next_token)

    def delete_dashboards(self) -> Union[str, ERROR_RESPONSE]:
        dashboards = self._get_multi_param("DashboardNames.member")
        if not dashboards:
            return self._error("InvalidParameterValue", "Need at least 1 dashboard")

        error = self.cloudwatch_backend.delete_dashboards(dashboards)
        if error is not None:
            return self._error("ResourceNotFound", error)

        template = self.response_template(DELETE_DASHBOARD_TEMPLATE)
        return template.render()

    @staticmethod
    def filter_alarms(
        alarms: Iterable[FakeAlarm], metric_name: str, namespace: str
    ) -> List[FakeAlarm]:
        metric_filtered_alarms = []

        for alarm in alarms:
            if alarm.metric_name == metric_name and alarm.namespace == namespace:
                metric_filtered_alarms.append(alarm)
        return metric_filtered_alarms

    def describe_alarms_for_metric(self) -> str:
        alarms = self.cloudwatch_backend.describe_alarms()
        namespace = self._get_param("Namespace")
        metric_name = self._get_param("MetricName")
        filtered_alarms = self.filter_alarms(alarms, metric_name, namespace)
        template = self.response_template(DESCRIBE_METRIC_ALARMS_TEMPLATE)
        return template.render(alarms=filtered_alarms)

    def disable_alarm_actions(self) -> str:
        raise NotImplementedError()

    def enable_alarm_actions(self) -> str:
        raise NotImplementedError()

    def get_dashboard(self) -> Union[str, ERROR_RESPONSE]:
        dashboard_name = self._get_param("DashboardName")

        dashboard = self.cloudwatch_backend.get_dashboard(dashboard_name)
        if dashboard is None:
            return self._error("ResourceNotFound", "Dashboard does not exist")

        template = self.response_template(GET_DASHBOARD_TEMPLATE)
        return template.render(dashboard=dashboard)

    def list_dashboards(self) -> str:
        prefix = self._get_param("DashboardNamePrefix", "")

        dashboards = self.cloudwatch_backend.list_dashboards(prefix)

        template = self.response_template(LIST_DASHBOARD_RESPONSE)
        return template.render(dashboards=dashboards)

    def put_dashboard(self) -> Union[str, ERROR_RESPONSE]:
        name = self._get_param("DashboardName")
        body = self._get_param("DashboardBody")

        try:
            json.loads(body)
        except ValueError:
            return self._error("InvalidParameterInput", "Body is invalid JSON")

        self.cloudwatch_backend.put_dashboard(name, body)

        template = self.response_template(PUT_DASHBOARD_RESPONSE)
        return template.render()

    def set_alarm_state(self) -> str:
        alarm_name = self._get_param("AlarmName")
        reason = self._get_param("StateReason")
        reason_data = self._get_param("StateReasonData")
        state_value = self._get_param("StateValue")

        self.cloudwatch_backend.set_alarm_state(
            alarm_name, reason, reason_data, state_value
        )

        template = self.response_template(SET_ALARM_STATE_TEMPLATE)
        return template.render()

    def list_tags_for_resource(self) -> str:
        resource_arn = self._get_param("ResourceARN")

        tags = self.cloudwatch_backend.list_tags_for_resource(resource_arn)

        template = self.response_template(LIST_TAGS_FOR_RESOURCE_TEMPLATE)
        return template.render(tags=tags)

    def tag_resource(self) -> str:
        resource_arn = self._get_param("ResourceARN")
        tags = self._get_multi_param("Tags.member")

        self.cloudwatch_backend.tag_resource(resource_arn, tags)

        template = self.response_template(TAG_RESOURCE_TEMPLATE)
        return template.render()

    def untag_resource(self) -> str:
        resource_arn = self._get_param("ResourceARN")
        tag_keys = self._get_multi_param("TagKeys.member")

        self.cloudwatch_backend.untag_resource(resource_arn, tag_keys)

        template = self.response_template(UNTAG_RESOURCE_TEMPLATE)
        return template.render()


PUT_METRIC_ALARM_TEMPLATE = """<PutMetricAlarmResponse xmlns="http://monitoring.amazonaws.com/doc/2010-08-01/">
   <ResponseMetadata>
      <RequestId>
         {{ request_id }}
      </RequestId>
   </ResponseMetadata>
</PutMetricAlarmResponse>"""

DESCRIBE_ALARMS_TEMPLATE = """<DescribeAlarmsResponse xmlns="http://monitoring.amazonaws.com/doc/2010-08-01/">
    <DescribeAlarmsResult>
        {% for tag_name, alarms in (('MetricAlarms', metric_alarms), ('CompositeAlarms', composite_alarms)) %}
        <{{tag_name}}>
            {% for alarm in alarms %}
            <member>
                <ActionsEnabled>{{ "true" if alarm.actions_enabled else "false" }}</ActionsEnabled>
                <AlarmActions>
                    {% for action in alarm.alarm_actions %}
                    <member>{{ action }}</member>
                    {% endfor %}
                </AlarmActions>
                <AlarmArn>{{ alarm.alarm_arn }}</AlarmArn>
                <AlarmConfigurationUpdatedTimestamp>{{ alarm.configuration_updated_timestamp }}</AlarmConfigurationUpdatedTimestamp>
                <AlarmDescription>{{ alarm.description or '' }}</AlarmDescription>
                <AlarmName>{{ alarm.name }}</AlarmName>
                <ComparisonOperator>{{ alarm.comparison_operator }}</ComparisonOperator>
                {% if alarm.dimensions is not none %}
                    <Dimensions>
                        {% for dimension in alarm.dimensions %}
                        <member>
                            <Name>{{ dimension.name }}</Name>
                            <Value>{{ dimension.value }}</Value>
                        </member>
                        {% endfor %}
                    </Dimensions>
                {% endif %}
                <EvaluationPeriods>{{ alarm.evaluation_periods }}</EvaluationPeriods>
                {% if alarm.datapoints_to_alarm is not none %}
                <DatapointsToAlarm>{{ alarm.datapoints_to_alarm }}</DatapointsToAlarm>
                {% endif %}
                <InsufficientDataActions>
                    {% for action in alarm.insufficient_data_actions %}
                    <member>{{ action }}</member>
                    {% endfor %}
                </InsufficientDataActions>
                {% if alarm.metric_name is not none %}
                <MetricName>{{ alarm.metric_name }}</MetricName>
                {% endif %}
                {% if alarm.metric_data_queries is not none %}
                <Metrics>
                    {% for metric in alarm.metric_data_queries %}
                     <member>
                        <Id>{{ metric.id }}</Id>
                        {% if metric.label is not none %}
                        <Label>{{ metric.label }}</Label>
                        {% endif %}
                        {% if metric.expression is not none %}
                        <Expression>{{ metric.expression }}</Expression>
                        {% endif %}
                        {% if metric.metric_stat is not none %}
                        <MetricStat>
                            <Metric>
                                <Namespace>{{ metric.metric_stat.metric.namespace }}</Namespace>
                                <MetricName>{{ metric.metric_stat.metric.metric_name }}</MetricName>
                                <Dimensions>
                                {% for dim in metric.metric_stat.metric.dimensions %}
                                    <member>
                                        <Name>{{ dim.name }}</Name>
                                        <Value>{{ dim.value }}</Value>
                                    </member>
                                {% endfor %}
                                </Dimensions>
                            </Metric>
                            {% if metric.metric_stat.period is not none %}
                            <Period>{{ metric.metric_stat.period }}</Period>
                            {% endif %}
                            <Stat>{{ metric.metric_stat.stat }}</Stat>
                            {% if metric.metric_stat.unit is not none %}
                            <Unit>{{ metric.metric_stat.unit }}</Unit>
                            {% endif %}
                        </MetricStat>
                        {% endif %}
                        {% if metric.period is not none %}
                        <Period>{{ metric.period }}</Period>
                        {% endif %}
                        <ReturnData>{{ metric.return_data }}</ReturnData>
                    </member>
                    {% endfor %}
                </Metrics>
                {% endif %}
                {% if alarm.namespace is not none %}
                <Namespace>{{ alarm.namespace }}</Namespace>
                {% endif %}
                <OKActions>
                    {% for action in alarm.ok_actions %}
                    <member>{{ action }}</member>
                    {% endfor %}
                </OKActions>
                {% if alarm.period is not none %}
                <Period>{{ alarm.period }}</Period>
                {% endif %}
                <StateReason>{{ alarm.state_reason }}</StateReason>
                <StateReasonData>{{ alarm.state_reason_data }}</StateReasonData>
                <StateUpdatedTimestamp>{{ alarm.state_updated_timestamp }}</StateUpdatedTimestamp>
                <StateValue>{{ alarm.state_value }}</StateValue>
                {% if alarm.statistic is not none %}
                <Statistic>{{ alarm.statistic }}</Statistic>
                {% endif %}
                {% if alarm.extended_statistic is not none %}
                <ExtendedStatistic>{{ alarm.extended_statistic }}</ExtendedStatistic>
                {% endif %}
                {% if alarm.threshold is not none %}
                <Threshold>{{ alarm.threshold }}</Threshold>
                {% endif %}
                {% if alarm.unit is not none %}
                <Unit>{{ alarm.unit }}</Unit>
                {% endif %}
                {% if alarm.treat_missing_data is not none %}
                <TreatMissingData>{{ alarm.treat_missing_data }}</TreatMissingData>
                {% endif %}
                {% if alarm.evaluate_low_sample_count_percentile is not none %}
                <EvaluateLowSampleCountPercentile>{{ alarm.evaluate_low_sample_count_percentile }}</EvaluateLowSampleCountPercentile>
                {% endif %}
                {% if alarm.threshold_metric_id is not none %}
                <ThresholdMetricId>{{ alarm.threshold_metric_id }}</ThresholdMetricId>
                {% endif %}
                {% if alarm.rule is not none %}
                <AlarmRule>{{ alarm.rule }}</AlarmRule>
                {% endif %}
            </member>
            {% endfor %}
        </{{tag_name}}>
        {% endfor %}
    </DescribeAlarmsResult>
</DescribeAlarmsResponse>"""

DESCRIBE_METRIC_ALARMS_TEMPLATE = """<DescribeAlarmsForMetricResponse xmlns="http://monitoring.amazonaws.com/doc/2010-08-01/">
    <DescribeAlarmsForMetricResult>
        <MetricAlarms>
            {% for alarm in alarms %}
            <member>
                <ActionsEnabled>{{ "true" if alarm.actions_enabled else "false" }}</ActionsEnabled>
                <AlarmActions>
                    {% for action in alarm.alarm_actions %}
                    <member>{{ action }}</member>
                    {% endfor %}
                </AlarmActions>
                <AlarmArn>{{ alarm.alarm_arn }}</AlarmArn>
                <AlarmConfigurationUpdatedTimestamp>{{ alarm.configuration_updated_timestamp }}</AlarmConfigurationUpdatedTimestamp>
                <AlarmDescription>{{ alarm.description }}</AlarmDescription>
                <AlarmName>{{ alarm.name }}</AlarmName>
                <ComparisonOperator>{{ alarm.comparison_operator }}</ComparisonOperator>
                <Dimensions>
                    {% for dimension in alarm.dimensions %}
                    <member>
                        <Name>{{ dimension.name }}</Name>
                        <Value>{{ dimension.value }}</Value>
                    </member>
                    {% endfor %}
                </Dimensions>
                <EvaluationPeriods>{{ alarm.evaluation_periods }}</EvaluationPeriods>
                <InsufficientDataActions>
                    {% for action in alarm.insufficient_data_actions %}
                    <member>{{ action }}</member>
                    {% endfor %}
                </InsufficientDataActions>
                <MetricName>{{ alarm.metric_name }}</MetricName>
                <Namespace>{{ alarm.namespace }}</Namespace>
                <OKActions>
                    {% for action in alarm.ok_actions %}
                    <member>{{ action }}</member>
                    {% endfor %}
                </OKActions>
                <Period>{{ alarm.period }}</Period>
                <StateReason>{{ alarm.state_reason }}</StateReason>
                <StateReasonData>{{ alarm.state_reason_data }}</StateReasonData>
                <StateUpdatedTimestamp>{{ alarm.state_updated_timestamp }}</StateUpdatedTimestamp>
                <StateValue>{{ alarm.state_value }}</StateValue>
                <Statistic>{{ alarm.statistic }}</Statistic>
                {% if alarm.threshold is not none %}
                <Threshold>{{ alarm.threshold }}</Threshold>
                {% endif %}
                <Unit>{{ alarm.unit }}</Unit>
            </member>
            {% endfor %}
        </MetricAlarms>
    </DescribeAlarmsForMetricResult>
</DescribeAlarmsForMetricResponse>"""

DELETE_METRIC_ALARMS_TEMPLATE = """<DeleteMetricAlarmResponse xmlns="http://monitoring.amazonaws.com/doc/2010-08-01/">
   <ResponseMetadata>
      <RequestId>
         {{ request_id }}
      </RequestId>
   </ResponseMetadata>
</DeleteMetricAlarmResponse>"""

PUT_METRIC_DATA_TEMPLATE = """<PutMetricDataResponse xmlns="http://monitoring.amazonaws.com/doc/2010-08-01/">
   <ResponseMetadata>
      <RequestId>
         {{ request_id }}
      </RequestId>
   </ResponseMetadata>
</PutMetricDataResponse>"""

GET_METRIC_DATA_TEMPLATE = """<GetMetricDataResponse xmlns="http://monitoring.amazonaws.com/doc/2010-08-01/">
   <GetMetricDataResult>
       <MetricDataResults>
           {% for result in results %}
            <member>
                <Id>{{ result.id }}</Id>
                <Label>{{ result.label }}</Label>
                <StatusCode>Complete</StatusCode>
                <Timestamps>
                    {% for val in result.timestamps %}
                    <member>{{ val }}</member>
                    {% endfor %}
                </Timestamps>
                <Values>
                    {% for val in result.vals %}
                    <member>{{ val }}</member>
                    {% endfor %}
                </Values>
            </member>
            {% endfor %}
       </MetricDataResults>
   </GetMetricDataResult>
   <ResponseMetadata>
       <RequestId>
            {{ request_id }}
       </RequestId>
   </ResponseMetadata>
</GetMetricDataResponse>"""

GET_METRIC_STATISTICS_TEMPLATE = """<GetMetricStatisticsResponse xmlns="http://monitoring.amazonaws.com/doc/2010-08-01/">
  <GetMetricStatisticsResult>
      <Label>{{ label }}</Label>
      <Datapoints>
        {% for datapoint in datapoints %}
            <member>
              {% if datapoint.sum is not none %}
              <Sum>{{ datapoint.sum }}</Sum>
              {% endif %}

              {% if datapoint.average is not none %}
              <Average>{{ datapoint.average }}</Average>
              {% endif %}

              {% if datapoint.maximum is not none %}
              <Maximum>{{ datapoint.maximum }}</Maximum>
              {% endif %}

              {% if datapoint.minimum is not none %}
              <Minimum>{{ datapoint.minimum }}</Minimum>
              {% endif %}

              {% if datapoint.sample_count is not none %}
              <SampleCount>{{ datapoint.sample_count }}</SampleCount>
              {% endif %}

              {% if datapoint.extended_statistics is not none %}
              <ExtendedStatistics>{{ datapoint.extended_statistics }}</ExtendedStatistics>
              {% endif %}

              <Timestamp>{{ datapoint.timestamp }}</Timestamp>
              {% if datapoint.unit is not none %}
              <Unit>{{ datapoint.unit }}</Unit>
              {% endif %}
            </member>
        {% endfor %}
      </Datapoints>
    </GetMetricStatisticsResult>
    <ResponseMetadata>
      <RequestId>
        {{ request_id }}
      </RequestId>
    </ResponseMetadata>
</GetMetricStatisticsResponse>"""

LIST_METRICS_TEMPLATE = """<ListMetricsResponse xmlns="http://monitoring.amazonaws.com/doc/2010-08-01/">
    <ListMetricsResult>
        <Metrics>
            {% for metric in metrics %}
            <member>
                <Dimensions>
                    {% for dimension in metric.dimensions %}
                    <member>
                        <Name>{{ dimension.name }}</Name>
                        <Value>{{ dimension.value }}</Value>
                    </member>
                    {% endfor %}
                </Dimensions>
                <MetricName>{{ metric.name }}</MetricName>
                <Namespace>{{ metric.namespace }}</Namespace>
            </member>
            {% endfor %}
        </Metrics>
        {% if next_token is not none %}
        <NextToken>
            {{ next_token }}
        </NextToken>
        {% endif %}
    </ListMetricsResult>
</ListMetricsResponse>"""

PUT_DASHBOARD_RESPONSE = """<PutDashboardResponse xmlns="http://monitoring.amazonaws.com/doc/2010-08-01/">
  <PutDashboardResult>
    <DashboardValidationMessages/>
  </PutDashboardResult>
  <ResponseMetadata>
    <RequestId>{{ request_id }}</RequestId>
  </ResponseMetadata>
</PutDashboardResponse>"""

LIST_DASHBOARD_RESPONSE = """<ListDashboardsResponse xmlns="http://monitoring.amazonaws.com/doc/2010-08-01/">
  <ListDashboardsResult>
    <DashboardEntries>
      {% for dashboard in dashboards %}
      <member>
        <DashboardArn>{{ dashboard.arn }}</DashboardArn>
        <LastModified>{{ dashboard.last_modified_iso }}</LastModified>
        <Size>{{ dashboard.size }}</Size>
        <DashboardName>{{ dashboard.name }}</DashboardName>
      </member>
      {% endfor %}
    </DashboardEntries>
  </ListDashboardsResult>
  <ResponseMetadata>
    <RequestId>{{ request_id }}</RequestId>
  </ResponseMetadata>
</ListDashboardsResponse>"""

DELETE_DASHBOARD_TEMPLATE = """<DeleteDashboardsResponse xmlns="http://monitoring.amazonaws.com/doc/2010-08-01/">
  <DeleteDashboardsResult/>
  <ResponseMetadata>
    <RequestId>{{ request_id }}</RequestId>
  </ResponseMetadata>
</DeleteDashboardsResponse>"""

GET_DASHBOARD_TEMPLATE = """<GetDashboardResponse xmlns="http://monitoring.amazonaws.com/doc/2010-08-01/">
  <GetDashboardResult>
    <DashboardArn>{{ dashboard.arn }}</DashboardArn>
    <DashboardBody>{{ dashboard.body }}</DashboardBody>
    <DashboardName>{{ dashboard.name }}</DashboardName>
  </GetDashboardResult>
  <ResponseMetadata>
    <RequestId>{{ request_id }}</RequestId>
  </ResponseMetadata>
</GetDashboardResponse>
"""

SET_ALARM_STATE_TEMPLATE = """<SetAlarmStateResponse xmlns="http://monitoring.amazonaws.com/doc/2010-08-01/">
  <ResponseMetadata>
    <RequestId>{{ request_id }}</RequestId>
  </ResponseMetadata>
</SetAlarmStateResponse>"""

ERROR_RESPONSE_TEMPLATE = """<ErrorResponse xmlns="http://monitoring.amazonaws.com/doc/2010-08-01/">
  <Error>
    <Type>Sender</Type>
    <Code>{{ code }}</Code>
    <Message>{{ message }}</Message>
  </Error>
  <RequestId>{{ request_id }}</RequestId>
</ErrorResponse>"""

LIST_TAGS_FOR_RESOURCE_TEMPLATE = """<ListTagsForResourceResponse xmlns="http://monitoring.amazonaws.com/doc/2010-08-01/">
  <ListTagsForResourceResult>
    <Tags>
      {% for key, value in tags.items() %}
      <member>
        <Key>{{ key }}</Key>
        <Value>{{ value }}</Value>
      </member>
      {% endfor %}
    </Tags>
  </ListTagsForResourceResult>
  <ResponseMetadata>
    <RequestId>{{ request_id }}</RequestId>
  </ResponseMetadata>
</ListTagsForResourceResponse>
"""

TAG_RESOURCE_TEMPLATE = """<TagResourceResponse xmlns="http://monitoring.amazonaws.com/doc/2010-08-01/">
  <TagResourceResult/>
  <ResponseMetadata>
    <RequestId>{{ request_id }}</RequestId>
  </ResponseMetadata>
</TagResourceResponse>"""

UNTAG_RESOURCE_TEMPLATE = """<UntagResourceResponse xmlns="http://monitoring.amazonaws.com/doc/2010-08-01/">
  <UntagResourceResult/>
  <ResponseMetadata>
    <RequestId>{{ request_id }}</RequestId>
  </ResponseMetadata>
</UntagResourceResponse>"""
