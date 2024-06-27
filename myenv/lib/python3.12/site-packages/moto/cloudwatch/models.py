import json
from datetime import datetime, timedelta
from typing import Any, Dict, Iterable, List, Optional, SupportsFloat, Tuple

from dateutil import parser
from dateutil.tz import tzutc

from moto.core.base_backend import BaseBackend
from moto.core.common_models import BackendDict, BaseModel, CloudWatchMetricProvider
from moto.core.utils import (
    iso_8601_datetime_with_nanoseconds,
    iso_8601_datetime_without_milliseconds,
    utcnow,
)
from moto.moto_api._internal import mock_random

from ..utilities.tagging_service import TaggingService
from .exceptions import (
    InvalidFormat,
    InvalidParameterCombination,
    InvalidParameterValue,
    ResourceNotFound,
    ResourceNotFoundException,
    ValidationError,
)
from .metric_data_expression_parser import parse_expression
from .utils import make_arn_for_alarm, make_arn_for_dashboard

_EMPTY_LIST: Any = tuple()


class Dimension(object):
    def __init__(self, name: Optional[str], value: Optional[str]):
        self.name = name
        self.value = value

    def __eq__(self, item: Any) -> bool:
        if isinstance(item, Dimension):
            return self.name == item.name and (
                self.value is None or item.value is None or self.value == item.value
            )
        return False

    def __lt__(self, other: "Dimension") -> bool:
        return self.name < other.name and self.value < other.name  # type: ignore[operator]


class Metric(object):
    def __init__(self, metric_name: str, namespace: str, dimensions: List[Dimension]):
        self.metric_name = metric_name
        self.namespace = namespace
        self.dimensions = dimensions


class MetricStat(object):
    def __init__(self, metric: Metric, period: str, stat: str, unit: str):
        self.metric = metric
        self.period = period
        self.stat = stat
        self.unit = unit


class MetricDataQuery(object):
    def __init__(
        self,
        query_id: str,
        label: str,
        period: str,
        return_data: str,
        expression: Optional[str] = None,
        metric_stat: Optional[MetricStat] = None,
    ):
        self.id = query_id
        self.label = label
        self.period = period
        self.return_data = return_data
        self.expression = expression
        self.metric_stat = metric_stat


def daterange(
    start: datetime,
    stop: datetime,
    step: timedelta = timedelta(days=1),
    inclusive: bool = False,
) -> Iterable[datetime]:
    """
    This method will iterate from `start` to `stop` datetimes with a timedelta step of `step`
    (supports iteration forwards or backwards in time)

    :param start: start datetime
    :param stop: end datetime
    :param step: step size as a timedelta
    :param inclusive: if True, last item returned will be as step closest to `end` (or `end` if no remainder).
    """

    # inclusive=False to behave like range by default
    total_step_secs = step.total_seconds()
    assert total_step_secs != 0

    if total_step_secs > 0:
        while start < stop:
            yield start
            start = start + step
    else:
        while stop < start:
            yield start
            start = start + step

    if inclusive and start == stop:
        yield start


class FakeAlarm(BaseModel):
    def __init__(
        self,
        account_id: str,
        region_name: str,
        name: str,
        namespace: str,
        metric_name: str,
        metric_data_queries: Optional[List[MetricDataQuery]],
        comparison_operator: str,
        evaluation_periods: int,
        datapoints_to_alarm: Optional[int],
        period: int,
        threshold: float,
        statistic: str,
        extended_statistic: Optional[str],
        description: str,
        dimensions: List[Dict[str, str]],
        alarm_actions: List[str],
        ok_actions: Optional[List[str]],
        insufficient_data_actions: Optional[List[str]],
        unit: Optional[str],
        actions_enabled: bool,
        treat_missing_data: Optional[str],
        evaluate_low_sample_count_percentile: Optional[str],
        threshold_metric_id: Optional[str],
        rule: Optional[str],
    ):
        self.region_name = region_name
        self.name = name
        self.alarm_arn = make_arn_for_alarm(region_name, account_id, name)
        self.namespace = namespace
        self.metric_name = metric_name
        self.metric_data_queries = metric_data_queries or []
        self.comparison_operator = comparison_operator
        self.evaluation_periods = evaluation_periods
        self.datapoints_to_alarm = datapoints_to_alarm
        self.period = period
        self.threshold = threshold
        self.statistic = statistic
        self.extended_statistic = extended_statistic
        self.description = description
        self.dimensions = [
            Dimension(dimension["name"], dimension["value"]) for dimension in dimensions
        ]
        self.actions_enabled = True if actions_enabled is None else actions_enabled
        self.alarm_actions = alarm_actions
        self.ok_actions = ok_actions or []
        self.insufficient_data_actions = insufficient_data_actions or []
        self.unit = unit
        self.configuration_updated_timestamp = iso_8601_datetime_with_nanoseconds()
        self.treat_missing_data = treat_missing_data
        self.evaluate_low_sample_count_percentile = evaluate_low_sample_count_percentile
        self.threshold_metric_id = threshold_metric_id

        self.history: List[Any] = []

        self.state_reason = "Unchecked: Initial alarm creation"
        self.state_reason_data = "{}"
        self.state_value = "OK"
        self.state_updated_timestamp = iso_8601_datetime_with_nanoseconds()

        # only used for composite alarms
        self.rule = rule

    def update_state(self, reason: str, reason_data: str, state_value: str) -> None:
        # History type, that then decides what the rest of the items are, can be one of ConfigurationUpdate | StateUpdate | Action
        self.history.append(
            (
                "StateUpdate",
                self.state_reason,
                self.state_reason_data,
                self.state_value,
                self.state_updated_timestamp,
            )
        )

        self.state_reason = reason
        self.state_reason_data = reason_data
        self.state_value = state_value
        self.state_updated_timestamp = iso_8601_datetime_with_nanoseconds()


def are_dimensions_same(
    metric_dimensions: List[Dimension], dimensions: List[Dimension]
) -> bool:
    if len(metric_dimensions) != len(dimensions):
        return False
    for dimension in metric_dimensions:
        for new_dimension in dimensions:
            if (
                dimension.name != new_dimension.name
                or dimension.value != new_dimension.value
            ):
                return False
    return True


class MetricDatumBase(BaseModel):
    """
    Base class for Metrics Datum (represents value or statistics set by put-metric-data)
    """

    def __init__(
        self,
        namespace: str,
        name: str,
        dimensions: List[Dict[str, str]],
        timestamp: datetime,
        unit: Any = None,
    ):
        self.namespace = namespace
        self.name = name
        self.timestamp = timestamp or utcnow().replace(tzinfo=tzutc())
        self.dimensions = [
            Dimension(dimension["Name"], dimension["Value"]) for dimension in dimensions
        ]
        self.unit = unit

    def filter(
        self,
        namespace: Optional[str],
        name: Optional[str],
        dimensions: List[Dict[str, str]],
        already_present_metrics: Optional[List["MetricDatumBase"]] = None,
    ) -> bool:
        if namespace and namespace != self.namespace:
            return False
        if name and name != self.name:
            return False

        for metric in already_present_metrics or []:
            if (
                (
                    self.dimensions
                    and are_dimensions_same(metric.dimensions, self.dimensions)
                )
                and self.name == metric.name
                and self.namespace == metric.namespace
            ):  # should be considered as already present only when name, namespace and dimensions all three are same
                return False

        if dimensions and any(
            Dimension(d["Name"], d.get("Value")) not in self.dimensions
            for d in dimensions
        ):
            return False
        return True


class MetricDatum(MetricDatumBase):
    """
    Single Metric value, represents the "value" (or a single value from the list "values") used in put-metric-data
    """

    def __init__(
        self,
        namespace: str,
        name: str,
        value: float,
        dimensions: List[Dict[str, str]],
        timestamp: datetime,
        unit: Any = None,
    ):
        super().__init__(namespace, name, dimensions, timestamp, unit)
        self.value = value


class MetricAggregatedDatum(MetricDatumBase):
    """
    Metric Statistics, represents "statistics-values" used in put-metric-data
    """

    def __init__(
        self,
        namespace: str,
        name: str,
        min_stat: float,
        max_stat: float,
        sample_count: float,
        sum_stat: float,
        dimensions: List[Dict[str, str]],
        timestamp: datetime,
        unit: Any = None,
    ):
        super().__init__(namespace, name, dimensions, timestamp, unit)
        self.min = min_stat
        self.max = max_stat
        self.sample_count = sample_count
        self.sum = sum_stat


class Dashboard(BaseModel):
    def __init__(self, account_id: str, region_name: str, name: str, body: str):
        # Guaranteed to be unique for now as the name is also the key of a dictionary where they are stored
        self.arn = make_arn_for_dashboard(account_id, region_name, name)
        self.name = name
        self.body = body
        self.last_modified = datetime.now()

    @property
    def last_modified_iso(self) -> str:
        return self.last_modified.isoformat()

    @property
    def size(self) -> int:
        return len(self)

    def __len__(self) -> int:
        return len(self.body)

    def __repr__(self) -> str:
        return f"<CloudWatchDashboard {self.name}>"


class Statistics:
    """
    Helper class to calculate statics for a list of metrics (MetricDatum, or MetricAggregatedDatum)
    """

    def __init__(self, stats: List[str], dt: datetime, unit: Optional[str] = None):
        self.timestamp: str = iso_8601_datetime_without_milliseconds(dt or utcnow())
        self.metric_data: List[MetricDatumBase] = []
        self.stats = stats
        self.unit = unit

    def get_statistics_for_type(self, stat: str) -> Optional[SupportsFloat]:
        """Calculates the statistic for the metric_data provided

        :param stat: the statistic that should be returned, case-sensitive (Sum, Average, Minium, Maximum, SampleCount)
        :return: the statistic of the current 'metric_data' in this class, or 0
        """
        if stat == "Sum":
            return self.sum
        if stat == "Average":
            return self.average
        if stat == "Minimum":
            return self.minimum
        if stat == "Maximum":
            return self.maximum
        if stat == "SampleCount":
            return self.sample_count
        return None

    @property
    def metric_single_values_list(self) -> List[float]:
        """
        :return: list of all values for the MetricDatum instances of the metric_data list
        """
        return [m.value for m in self.metric_data or [] if isinstance(m, MetricDatum)]

    @property
    def metric_aggregated_list(self) -> List[MetricAggregatedDatum]:
        """
        :return: list of all MetricAggregatedDatum instances from the metric_data list
        """
        return [
            s for s in self.metric_data or [] if isinstance(s, MetricAggregatedDatum)
        ]

    @property
    def sample_count(self) -> Optional[SupportsFloat]:
        if "SampleCount" not in self.stats:
            return None

        return self.calc_sample_count()

    @property
    def sum(self) -> Optional[SupportsFloat]:
        if "Sum" not in self.stats:
            return None

        return self.calc_sum()

    @property
    def minimum(self) -> Optional[SupportsFloat]:
        if "Minimum" not in self.stats:
            return None
        if not self.metric_single_values_list and not self.metric_aggregated_list:
            return None

        metrics = self.metric_single_values_list + [
            s.min for s in self.metric_aggregated_list
        ]
        return min(metrics)

    @property
    def maximum(self) -> Optional[SupportsFloat]:
        if "Maximum" not in self.stats:
            return None

        if not self.metric_single_values_list and not self.metric_aggregated_list:
            return None

        metrics = self.metric_single_values_list + [
            s.max for s in self.metric_aggregated_list
        ]
        return max(metrics)

    @property
    def average(self) -> Optional[SupportsFloat]:
        if "Average" not in self.stats:
            return None

        sample_count = self.calc_sample_count()

        if not sample_count:
            return None

        return self.calc_sum() / sample_count

    def calc_sample_count(self) -> float:
        return len(self.metric_single_values_list) + sum(
            [s.sample_count for s in self.metric_aggregated_list]
        )

    def calc_sum(self) -> float:
        return sum(self.metric_single_values_list) + sum(
            [s.sum for s in self.metric_aggregated_list]
        )


class CloudWatchBackend(BaseBackend):
    def __init__(self, region_name: str, account_id: str):
        super().__init__(region_name, account_id)
        self.alarms: Dict[str, FakeAlarm] = {}
        self.dashboards: Dict[str, Dashboard] = {}
        self.metric_data: List[MetricDatumBase] = []
        self.paged_metric_data: Dict[str, List[MetricDatumBase]] = {}
        self.tagger = TaggingService()

    @property
    # Retrieve a list of all OOTB metrics that are provided by metrics providers
    # Computed on the fly
    def aws_metric_data(self) -> List[MetricDatumBase]:
        providers = CloudWatchMetricProvider.__subclasses__()
        md = []
        for provider in providers:
            md.extend(
                provider.get_cloudwatch_metrics(
                    self.account_id, region=self.region_name
                )
            )
        return md

    def put_metric_alarm(
        self,
        name: str,
        namespace: str,
        metric_name: str,
        comparison_operator: str,
        evaluation_periods: int,
        period: int,
        threshold: float,
        statistic: str,
        description: str,
        dimensions: List[Dict[str, str]],
        alarm_actions: List[str],
        metric_data_queries: Optional[List[MetricDataQuery]] = None,
        datapoints_to_alarm: Optional[int] = None,
        extended_statistic: Optional[str] = None,
        ok_actions: Optional[List[str]] = None,
        insufficient_data_actions: Optional[List[str]] = None,
        unit: Optional[str] = None,
        actions_enabled: bool = True,
        treat_missing_data: Optional[str] = None,
        evaluate_low_sample_count_percentile: Optional[str] = None,
        threshold_metric_id: Optional[str] = None,
        rule: Optional[str] = None,
        tags: Optional[List[Dict[str, str]]] = None,
    ) -> FakeAlarm:
        if extended_statistic and not extended_statistic.startswith("p"):
            raise InvalidParameterValue(
                f"The value {extended_statistic} for parameter ExtendedStatistic is not supported."
            )
        if (
            evaluate_low_sample_count_percentile
            and evaluate_low_sample_count_percentile not in ("evaluate", "ignore")
        ):
            raise ValidationError(
                f"Option {evaluate_low_sample_count_percentile} is not supported. "
                "Supported options for parameter EvaluateLowSampleCountPercentile are evaluate and ignore."
            )

        alarm = FakeAlarm(
            account_id=self.account_id,
            region_name=self.region_name,
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
        )

        self.alarms[name] = alarm
        if tags:
            self.tagger.tag_resource(alarm.alarm_arn, tags)

        return alarm

    def describe_alarms(self) -> Iterable[FakeAlarm]:
        return self.alarms.values()

    @staticmethod
    def _list_element_starts_with(items: List[str], needle: str) -> bool:
        """True of any of the list elements starts with needle"""
        for item in items:
            if item.startswith(needle):
                return True
        return False

    def get_alarms_by_action_prefix(self, action_prefix: str) -> Iterable[FakeAlarm]:
        return [
            alarm
            for alarm in self.alarms.values()
            if CloudWatchBackend._list_element_starts_with(
                alarm.alarm_actions, action_prefix
            )
        ]

    def get_alarms_by_alarm_name_prefix(self, name_prefix: str) -> Iterable[FakeAlarm]:
        return [
            alarm
            for alarm in self.alarms.values()
            if alarm.name.startswith(name_prefix)
        ]

    def get_alarms_by_alarm_names(self, alarm_names: List[str]) -> Iterable[FakeAlarm]:
        return [alarm for alarm in self.alarms.values() if alarm.name in alarm_names]

    def get_alarms_by_state_value(self, target_state: str) -> Iterable[FakeAlarm]:
        return filter(
            lambda alarm: alarm.state_value == target_state, self.alarms.values()
        )

    def delete_alarms(self, alarm_names: List[str]) -> None:
        for alarm_name in alarm_names:
            self.alarms.pop(alarm_name, None)

    def put_metric_data(
        self, namespace: str, metric_data: List[Dict[str, Any]]
    ) -> None:
        for i, metric in enumerate(metric_data):
            self._validate_parameters_put_metric_data(metric, i + 1)

        for metric_member in metric_data:
            # Preserve "datetime" for get_metric_statistics comparisons
            timestamp = metric_member.get("Timestamp")
            if timestamp is not None and type(timestamp) != datetime:
                timestamp = parser.parse(timestamp)
            metric_name = metric_member["MetricName"]
            dimension = metric_member.get("Dimensions.member", _EMPTY_LIST)
            unit = metric_member.get("Unit")

            # put_metric_data can include "value" as single value or "values" as a list
            if metric_member.get("Values.member"):
                values = metric_member["Values.member"]
                # value[i] should be added count[i] times (with default count 1)
                counts = metric_member.get("Counts.member") or ["1"] * len(values)
                for i in range(0, len(values)):
                    value = values[i]
                    timestamp = metric_member.get("Timestamp")
                    if timestamp is not None and type(timestamp) != datetime:
                        timestamp = parser.parse(timestamp)

                    # add the value count[i] times
                    for _ in range(0, int(float(counts[i]))):
                        self.metric_data.append(
                            MetricDatum(
                                namespace=namespace,
                                name=metric_name,
                                value=float(value),
                                dimensions=dimension,
                                timestamp=timestamp,
                                unit=unit,
                            )
                        )
            elif metric_member.get("StatisticValues"):
                stats = metric_member["StatisticValues"]
                self.metric_data.append(
                    MetricAggregatedDatum(
                        namespace=namespace,
                        name=metric_name,
                        sum_stat=float(stats["Sum"]),
                        min_stat=float(stats["Minimum"]),
                        max_stat=float(stats["Maximum"]),
                        sample_count=float(stats["SampleCount"]),
                        dimensions=dimension,
                        timestamp=timestamp,
                        unit=unit,
                    )
                )
            else:
                # there is only a single value
                self.metric_data.append(
                    MetricDatum(
                        namespace,
                        metric_name,
                        float(metric_member.get("Value", 0)),
                        dimension,
                        timestamp,
                        unit,
                    )
                )

    def get_metric_data(
        self,
        queries: List[Dict[str, Any]],
        start_time: datetime,
        end_time: datetime,
        scan_by: str = "TimestampAscending",
    ) -> List[Dict[str, Any]]:
        start_time = start_time.replace(microsecond=0)
        end_time = end_time.replace(microsecond=0)

        if start_time > end_time:
            raise ValidationError(
                "The parameter EndTime must be greater than StartTime."
            )
        if start_time == end_time:
            raise ValidationError(
                "The parameter StartTime must not equal parameter EndTime."
            )

        period_data = [
            md for md in self.get_all_metrics() if start_time <= md.timestamp < end_time
        ]

        results = []
        results_to_return = []
        metric_stat_queries = [q for q in queries if "MetricStat" in q]
        expression_queries = [q for q in queries if "Expression" in q]
        for query in metric_stat_queries:
            period_start_time = start_time
            metric_stat = query["MetricStat"]
            query_ns = metric_stat["Metric"]["Namespace"]
            query_name = metric_stat["Metric"]["MetricName"]
            delta = timedelta(seconds=int(metric_stat["Period"]))
            dimensions = [
                Dimension(name=d["Name"], value=d["Value"])
                for d in metric_stat["Metric"].get("Dimensions", [])
            ]
            unit = metric_stat.get("Unit")
            result_vals: List[SupportsFloat] = []
            timestamps: List[str] = []
            stat = metric_stat["Stat"]
            while period_start_time <= end_time:
                period_end_time = period_start_time + delta
                period_md = [
                    period_md
                    for period_md in period_data
                    if period_start_time <= period_md.timestamp < period_end_time
                ]

                query_period_data = [
                    md
                    for md in period_md
                    if md.namespace == query_ns and md.name == query_name
                ]
                if dimensions:
                    query_period_data = [
                        md
                        for md in period_md
                        if sorted(md.dimensions) == sorted(dimensions)
                        and md.name == query_name
                    ]
                # Filter based on unit value
                if unit:
                    query_period_data = [
                        md for md in query_period_data if md.unit == unit
                    ]

                if len(query_period_data) > 0:
                    stats = Statistics([stat], period_start_time)
                    stats.metric_data = query_period_data
                    result_vals.append(stats.get_statistics_for_type(stat))  # type: ignore[arg-type]

                    timestamps.append(stats.timestamp)
                period_start_time += delta
            if scan_by == "TimestampDescending" and len(timestamps) > 0:
                timestamps.reverse()
                result_vals.reverse()

            label = query.get("Label") or f"{query_name} {stat}"

            results.append(
                {
                    "id": query["Id"],
                    "label": label,
                    "vals": result_vals,
                    "timestamps": timestamps,
                }
            )
            if query.get("ReturnData", "true") == "true":
                results_to_return.append(
                    {
                        "id": query["Id"],
                        "label": label,
                        "vals": result_vals,
                        "timestamps": timestamps,
                    }
                )
        for query in expression_queries:
            label = query.get("Label") or f"{query_name} {stat}"
            result_vals, timestamps = parse_expression(query["Expression"], results)
            results_to_return.append(
                {
                    "id": query["Id"],
                    "label": label,
                    "vals": result_vals,
                    "timestamps": timestamps,
                }
            )
        return results_to_return

    def get_metric_statistics(
        self,
        namespace: str,
        metric_name: str,
        start_time: datetime,
        end_time: datetime,
        period: int,
        stats: List[str],
        dimensions: List[Dict[str, str]],
        unit: Optional[str] = None,
    ) -> List[Statistics]:
        start_time = start_time.replace(microsecond=0)
        end_time = end_time.replace(microsecond=0)

        if start_time >= end_time:
            raise InvalidParameterValue(
                "The parameter StartTime must be less than the parameter EndTime."
            )

        period_delta = timedelta(seconds=period)
        filtered_data = [
            md
            for md in self.get_all_metrics()
            if md.namespace == namespace
            and md.name == metric_name
            and start_time <= md.timestamp < end_time
        ]

        if unit:
            filtered_data = [md for md in filtered_data if md.unit == unit]
        if dimensions:
            filtered_data = [
                md for md in filtered_data if md.filter(None, None, dimensions)
            ]

        # earliest to oldest
        filtered_data = sorted(filtered_data, key=lambda x: x.timestamp)
        if not filtered_data:
            return []

        idx = 0
        data: List[Statistics] = list()
        for dt in daterange(
            filtered_data[0].timestamp,
            filtered_data[-1].timestamp + period_delta,
            period_delta,
        ):
            s = Statistics(stats, dt)
            while idx < len(filtered_data) and filtered_data[idx].timestamp < (
                dt + period_delta
            ):
                s.metric_data.append(filtered_data[idx])
                s.unit = filtered_data[idx].unit
                idx += 1

            if not s.metric_data:
                continue

            data.append(s)

        return data

    def get_all_metrics(self) -> List[MetricDatumBase]:
        return self.metric_data + self.aws_metric_data

    def put_dashboard(self, name: str, body: str) -> None:
        self.dashboards[name] = Dashboard(self.account_id, self.region_name, name, body)

    def list_dashboards(self, prefix: str = "") -> Iterable[Dashboard]:
        for key, value in self.dashboards.items():
            if key.startswith(prefix):
                yield value

    def delete_dashboards(self, dashboards: List[str]) -> Optional[str]:
        to_delete = set(dashboards)
        all_dashboards = set(self.dashboards.keys())

        left_over = to_delete - all_dashboards
        if len(left_over) > 0:
            # Some dashboards are not found
            db_list = ", ".join(left_over)
            return f"The specified dashboard does not exist. [{db_list}]"

        for dashboard in to_delete:
            del self.dashboards[dashboard]

        return None

    def get_dashboard(self, dashboard: str) -> Optional[Dashboard]:
        return self.dashboards.get(dashboard)

    def set_alarm_state(
        self, alarm_name: str, reason: str, reason_data: str, state_value: str
    ) -> None:
        try:
            if reason_data is not None:
                json.loads(reason_data)
        except ValueError:
            raise InvalidFormat("Unknown")

        if alarm_name not in self.alarms:
            raise ResourceNotFound

        if state_value not in ("OK", "ALARM", "INSUFFICIENT_DATA"):
            raise ValidationError(
                "1 validation error detected: "
                f"Value '{state_value}' at 'stateValue' failed to satisfy constraint: "
                "Member must satisfy enum value set: [INSUFFICIENT_DATA, ALARM, OK]"
            )

        self.alarms[alarm_name].update_state(reason, reason_data, state_value)

    def list_metrics(
        self,
        next_token: Optional[str],
        namespace: str,
        metric_name: str,
        dimensions: List[Dict[str, str]],
    ) -> Tuple[Optional[str], List[MetricDatumBase]]:
        if next_token:
            if next_token not in self.paged_metric_data:
                raise InvalidParameterValue("Request parameter NextToken is invalid")
            else:
                metrics = self.paged_metric_data[next_token]
                del self.paged_metric_data[next_token]  # Cant reuse same token twice
                return self._get_paginated(metrics)
        else:
            metrics = self.get_filtered_metrics(metric_name, namespace, dimensions)
            return self._get_paginated(metrics)

    def get_filtered_metrics(
        self, metric_name: str, namespace: str, dimensions: List[Dict[str, str]]
    ) -> List[MetricDatumBase]:
        metrics = self.get_all_metrics()
        new_metrics: List[MetricDatumBase] = []
        for md in metrics:
            if md.filter(
                namespace=namespace,
                name=metric_name,
                dimensions=dimensions,
                already_present_metrics=new_metrics,
            ):
                new_metrics.append(md)
        return new_metrics

    def list_tags_for_resource(self, arn: str) -> Dict[str, str]:
        return self.tagger.get_tag_dict_for_resource(arn)

    def tag_resource(self, arn: str, tags: List[Dict[str, str]]) -> None:
        # From boto3:
        # Currently, the only CloudWatch resources that can be tagged are alarms and Contributor Insights rules.
        all_arns = [alarm.alarm_arn for alarm in self.describe_alarms()]
        if arn not in all_arns:
            raise ResourceNotFoundException

        self.tagger.tag_resource(arn, tags)

    def untag_resource(self, arn: str, tag_keys: List[str]) -> None:
        if arn not in self.tagger.tags.keys():
            raise ResourceNotFoundException

        self.tagger.untag_resource_using_names(arn, tag_keys)

    def _get_paginated(
        self, metrics: List[MetricDatumBase]
    ) -> Tuple[Optional[str], List[MetricDatumBase]]:
        if len(metrics) > 500:
            next_token = str(mock_random.uuid4())
            self.paged_metric_data[next_token] = metrics[500:]
            return next_token, metrics[0:500]
        else:
            return None, metrics

    def _validate_parameters_put_metric_data(
        self, metric: Dict[str, Any], query_num: int
    ) -> None:
        """Runs some basic validation of the Metric Query

        :param metric: represents one metric query
        :param query_num: the query number (starting from 1)
        :returns: nothing if the validation passes, else an exception is thrown
        :raises: InvalidParameterValue
        :raises: InvalidParameterCombination
        """
        # basic validation of input
        if metric.get("Value") == "NaN":
            # single value
            raise InvalidParameterValue(
                f"The value NaN for parameter MetricData.member.{query_num}.Value is invalid."
            )
        if metric.get("Values.member"):
            # list of values
            if "Value" in metric:
                raise InvalidParameterValue(
                    f"The parameters MetricData.member.{query_num}.Value and MetricData.member.{query_num}.Values are mutually exclusive and you have specified both."
                )
            if metric.get("Counts.member"):
                if len(metric["Counts.member"]) != len(metric["Values.member"]):
                    raise InvalidParameterValue(
                        f"The parameters MetricData.member.{query_num}.Values and MetricData.member.{query_num}.Counts must be of the same size."
                    )
            for value in metric["Values.member"]:
                if value.lower() == "nan":
                    raise InvalidParameterValue(
                        f"The value {value} for parameter MetricData.member.{query_num}.Values is invalid."
                    )
        if metric.get("StatisticValues"):
            if metric.get("Value"):
                raise InvalidParameterCombination(
                    f"The parameters MetricData.member.{query_num}.Value and MetricData.member.{query_num}.StatisticValues are mutually exclusive and you have specified both."
                )

            # aggregated (statistic) for values, must contain sum, maximum, minimum and sample count
            statistic_values = metric["StatisticValues"]
            expected = ["Sum", "Maximum", "Minimum", "SampleCount"]
            for stat in expected:
                if stat not in statistic_values:
                    raise InvalidParameterValue(
                        f'Missing required parameter in MetricData[{query_num}].StatisticValues: "{stat}"'
                    )


cloudwatch_backends = BackendDict(CloudWatchBackend, "cloudwatch")
