import re
import time
from collections import OrderedDict
from enum import Enum, unique
from typing import TYPE_CHECKING, Any, Dict, List, Optional, Tuple, Union

from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.common_models import BaseModel
from moto.ecs import ecs_backends
from moto.moto_api._internal import mock_random
from moto.utilities.utils import ARN_PARTITION_REGEX, get_partition

from .exceptions import AWSValidationException

if TYPE_CHECKING:
    from moto.cloudwatch.models import FakeAlarm


@unique
class ResourceTypeExceptionValueSet(Enum):
    RESOURCE_TYPE = "ResourceType"
    # MSK currently only has the "broker-storage" resource type which is not part of the resource_id
    KAFKA_BROKER_STORAGE = "broker-storage"


@unique
class ServiceNamespaceValueSet(Enum):
    APPSTREAM = "appstream"
    RDS = "rds"
    LAMBDA = "lambda"
    CASSANDRA = "cassandra"
    DYNAMODB = "dynamodb"
    CUSTOM_RESOURCE = "custom-resource"
    ELASTICMAPREDUCE = "elasticmapreduce"
    EC2 = "ec2"
    COMPREHEND = "comprehend"
    ECS = "ecs"
    SAGEMAKER = "sagemaker"
    KAFKA = "kafka"


@unique
class ScalableDimensionValueSet(Enum):
    CASSANDRA_TABLE_READ_CAPACITY_UNITS = "cassandra:table:ReadCapacityUnits"
    CASSANDRA_TABLE_WRITE_CAPACITY_UNITS = "cassandra:table:WriteCapacityUnits"
    DYNAMODB_INDEX_READ_CAPACITY_UNITS = "dynamodb:index:ReadCapacityUnits"
    DYNAMODB_INDEX_WRITE_CAPACITY_UNITS = "dynamodb:index:WriteCapacityUnits"
    DYNAMODB_TABLE_READ_CAPACITY_UNITS = "dynamodb:table:ReadCapacityUnits"
    DYNAMODB_TABLE_WRITE_CAPACITY_UNITS = "dynamodb:table:WriteCapacityUnits"
    RDS_CLUSTER_READ_REPLICA_COUNT = "rds:cluster:ReadReplicaCount"
    RDS_CLUSTER_CAPACITY = "rds:cluster:Capacity"
    COMPREHEND_DOCUMENT_CLASSIFIER_ENDPOINT_DESIRED_INFERENCE_UNITS = (
        "comprehend:document-classifier-endpoint:DesiredInferenceUnits"
    )
    ELASTICMAPREDUCE_INSTANCE_FLEET_ON_DEMAND_CAPACITY = (
        "elasticmapreduce:instancefleet:OnDemandCapacity"
    )
    ELASTICMAPREDUCE_INSTANCE_FLEET_SPOT_CAPACITY = (
        "elasticmapreduce:instancefleet:SpotCapacity"
    )
    ELASTICMAPREDUCE_INSTANCE_GROUP_INSTANCE_COUNT = (
        "elasticmapreduce:instancegroup:InstanceCount"
    )
    LAMBDA_FUNCTION_PROVISIONED_CONCURRENCY = "lambda:function:ProvisionedConcurrency"
    APPSTREAM_FLEET_DESIRED_CAPACITY = "appstream:fleet:DesiredCapacity"
    CUSTOM_RESOURCE_RESOURCE_TYPE_PROPERTY = "custom-resource:ResourceType:Property"
    SAGEMAKER_VARIANT_DESIRED_INSTANCE_COUNT = "sagemaker:variant:DesiredInstanceCount"
    EC2_SPOT_FLEET_REQUEST_TARGET_CAPACITY = "ec2:spot-fleet-request:TargetCapacity"
    ECS_SERVICE_DESIRED_COUNT = "ecs:service:DesiredCount"
    KAFKA_BROKER_STORAGE_VOLUME_SIZE = "kafka:broker-storage:VolumeSize"


class ApplicationAutoscalingBackend(BaseBackend):
    def __init__(self, region_name: str, account_id: str) -> None:
        super().__init__(region_name, account_id)
        self.ecs_backend = ecs_backends[account_id][region_name]
        self.targets: Dict[str, Dict[str, FakeScalableTarget]] = OrderedDict()
        self.policies: Dict[str, FakeApplicationAutoscalingPolicy] = {}
        self.scheduled_actions: List[FakeScheduledAction] = list()

    def describe_scalable_targets(
        self, namespace: str, r_ids: Union[None, List[str]], dimension: Union[None, str]
    ) -> List["FakeScalableTarget"]:
        if r_ids is None:
            r_ids = []
        targets = self._flatten_scalable_targets(namespace)
        if dimension is not None:
            targets = [t for t in targets if t.scalable_dimension == dimension]
        if len(r_ids) > 0:
            targets = [t for t in targets if t.resource_id in r_ids]
        return targets

    def _flatten_scalable_targets(self, namespace: str) -> List["FakeScalableTarget"]:
        """Flatten scalable targets for a given service namespace down to a list."""
        targets = []
        for dimension in self.targets.keys():
            for resource_id in self.targets[dimension].keys():
                targets.append(self.targets[dimension][resource_id])
        targets = [t for t in targets if t.service_namespace == namespace]
        return targets

    def register_scalable_target(
        self,
        namespace: str,
        r_id: str,
        dimension: str,
        min_capacity: Optional[int],
        max_capacity: Optional[int],
        role_arn: str,
        suspended_state: str,
    ) -> "FakeScalableTarget":
        _ = _target_params_are_valid(namespace, r_id, dimension)
        if namespace == ServiceNamespaceValueSet.ECS.value:
            _ = self._ecs_service_exists_for_target(r_id)
        if self._scalable_target_exists(r_id, dimension):
            target = self.targets[dimension][r_id]
            target.update(min_capacity, max_capacity, suspended_state)
        else:
            target = FakeScalableTarget(
                self,
                namespace,
                r_id,
                dimension,
                min_capacity,
                max_capacity,
                role_arn,
                suspended_state,
            )
            self._add_scalable_target(target)
        return target

    def _scalable_target_exists(self, r_id: str, dimension: str) -> bool:
        return r_id in self.targets.get(dimension, [])

    def _ecs_service_exists_for_target(self, r_id: str) -> bool:
        """Raises a ValidationException if an ECS service does not exist
        for the specified resource ID.
        """
        _, cluster, service = r_id.split("/")
        result, _ = self.ecs_backend.describe_services(cluster, [service])
        if len(result) != 1:
            raise AWSValidationException(f"ECS service doesn't exist: {r_id}")
        return True

    def _add_scalable_target(
        self, target: "FakeScalableTarget"
    ) -> "FakeScalableTarget":
        if target.scalable_dimension not in self.targets:
            self.targets[target.scalable_dimension] = OrderedDict()
        if target.resource_id not in self.targets[target.scalable_dimension]:
            self.targets[target.scalable_dimension][target.resource_id] = target
        return target

    def deregister_scalable_target(
        self, namespace: str, r_id: str, dimension: str
    ) -> None:
        if self._scalable_target_exists(r_id, dimension):
            del self.targets[dimension][r_id]
        else:
            raise AWSValidationException(
                f"No scalable target found for service namespace: {namespace}, resource ID: {r_id}, scalable dimension: {dimension}"
            )

    def put_scaling_policy(
        self,
        policy_name: str,
        service_namespace: str,
        resource_id: str,
        scalable_dimension: str,
        policy_body: Dict[str, Any],
        policy_type: Optional[None],
    ) -> "FakeApplicationAutoscalingPolicy":
        policy_key = FakeApplicationAutoscalingPolicy.formulate_key(
            service_namespace, resource_id, scalable_dimension, policy_name
        )
        if policy_key in self.policies:
            old_policy = self.policies[policy_key]
            policy = FakeApplicationAutoscalingPolicy(
                account_id=self.account_id,
                region_name=self.region_name,
                policy_name=policy_name,
                service_namespace=service_namespace,
                resource_id=resource_id,
                scalable_dimension=scalable_dimension,
                policy_type=policy_type if policy_type else old_policy.policy_type,
                policy_body=policy_body if policy_body else old_policy._policy_body,
            )
        else:
            policy = FakeApplicationAutoscalingPolicy(
                account_id=self.account_id,
                region_name=self.region_name,
                policy_name=policy_name,
                service_namespace=service_namespace,
                resource_id=resource_id,
                scalable_dimension=scalable_dimension,
                policy_type=policy_type,
                policy_body=policy_body,
            )
        self.policies[policy_key] = policy
        return policy

    def describe_scaling_policies(
        self,
        service_namespace: str,
        resource_id: str,
        scalable_dimension: str,
        max_results: Optional[int],
        next_token: str,
    ) -> Tuple[Optional[str], List["FakeApplicationAutoscalingPolicy"]]:
        max_results = max_results or 100
        policies = [
            policy
            for policy in self.policies.values()
            if policy.service_namespace == service_namespace
        ]
        if resource_id:
            policies = [
                policy for policy in policies if policy.resource_id in resource_id
            ]
        if scalable_dimension:
            policies = [
                policy
                for policy in policies
                if policy.scalable_dimension in scalable_dimension
            ]
        starting_point = int(next_token) if next_token else 0
        ending_point = starting_point + max_results
        policies_page = policies[starting_point:ending_point]
        new_next_token = str(ending_point) if ending_point < len(policies) else None
        return new_next_token, policies_page

    def delete_scaling_policy(
        self,
        policy_name: str,
        service_namespace: str,
        resource_id: str,
        scalable_dimension: str,
    ) -> None:
        policy_key = FakeApplicationAutoscalingPolicy.formulate_key(
            service_namespace, resource_id, scalable_dimension, policy_name
        )
        if policy_key in self.policies:
            policy = self.policies[policy_key]
            policy.delete_alarms(self.account_id, self.region_name)
            del self.policies[policy_key]
        else:
            raise AWSValidationException(
                f"No scaling policy found for service namespace: {service_namespace}, resource ID: {resource_id}, scalable dimension: {scalable_dimension}, policy name: {policy_name}"
            )

    def delete_scheduled_action(
        self,
        service_namespace: str,
        scheduled_action_name: str,
        resource_id: str,
        scalable_dimension: str,
    ) -> None:
        self.scheduled_actions = [
            a
            for a in self.scheduled_actions
            if not (
                a.service_namespace == service_namespace
                and a.scheduled_action_name == scheduled_action_name
                and a.resource_id == resource_id
                and a.scalable_dimension == scalable_dimension
            )
        ]

    def describe_scheduled_actions(
        self,
        scheduled_action_names: str,
        service_namespace: str,
        resource_id: str,
        scalable_dimension: str,
    ) -> List["FakeScheduledAction"]:
        """
        Pagination is not yet implemented
        """
        result = [
            a
            for a in self.scheduled_actions
            if a.service_namespace == service_namespace
        ]
        if scheduled_action_names:
            result = [
                a for a in result if a.scheduled_action_name in scheduled_action_names
            ]
        if resource_id:
            result = [a for a in result if a.resource_id == resource_id]
        if scalable_dimension:
            result = [a for a in result if a.scalable_dimension == scalable_dimension]
        return result

    def put_scheduled_action(
        self,
        service_namespace: str,
        schedule: str,
        timezone: str,
        scheduled_action_name: str,
        resource_id: str,
        scalable_dimension: str,
        start_time: str,
        end_time: str,
        scalable_target_action: str,
    ) -> None:
        existing_action = next(
            (
                a
                for a in self.scheduled_actions
                if a.service_namespace == service_namespace
                and a.scheduled_action_name == scheduled_action_name
                and a.resource_id == resource_id
                and a.scalable_dimension == scalable_dimension
            ),
            None,
        )
        if existing_action:
            existing_action.update(
                schedule,
                timezone,
                scheduled_action_name,
                start_time,
                end_time,
                scalable_target_action,
            )
        else:
            action = FakeScheduledAction(
                service_namespace,
                schedule,
                timezone,
                scheduled_action_name,
                resource_id,
                scalable_dimension,
                start_time,
                end_time,
                scalable_target_action,
                self.account_id,
                self.region_name,
            )
            self.scheduled_actions.append(action)


def _target_params_are_valid(namespace: str, r_id: str, dimension: str) -> bool:
    """Check whether namespace, resource_id and dimension are valid and consistent with each other."""
    is_valid = True
    valid_namespaces = [n.value for n in ServiceNamespaceValueSet]
    if namespace not in valid_namespaces:
        is_valid = False
    if dimension is not None:
        try:
            valid_dimensions = [d.value for d in ScalableDimensionValueSet]
            resource_type_exceptions = [r.value for r in ResourceTypeExceptionValueSet]
            d_namespace, d_resource_type, _ = dimension.split(":")
            if d_resource_type not in resource_type_exceptions:
                resource_type = _get_resource_type_from_resource_id(r_id)
            else:
                resource_type = d_resource_type
            if (
                dimension not in valid_dimensions
                or d_namespace != namespace
                or resource_type != d_resource_type
            ):
                is_valid = False
        except ValueError:
            is_valid = False
    if not is_valid:
        raise AWSValidationException(
            "Unsupported service namespace, resource type or scalable dimension"
        )
    return is_valid


def _get_resource_type_from_resource_id(resource_id: str) -> str:
    # AWS Application Autoscaling resource_ids are multi-component (path-like) identifiers that vary in format,
    # depending on the type of resource it identifies.  resource_type is one of its components.
    #  resource_id format variations are described in
    #   https://docs.aws.amazon.com/autoscaling/application/APIReference/API_RegisterScalableTarget.html
    #  In a nutshell:
    #  - Most use slash separators, but some use colon separators.
    #  - The resource type is usually the first component of the resource_id...
    #    - ...except for sagemaker endpoints, dynamodb GSIs and keyspaces tables, where it's the third.
    #  - Comprehend uses an arn, with the resource type being the last element.

    if re.match(ARN_PARTITION_REGEX + ":comprehend", resource_id):
        resource_id = resource_id.split(":")[-1]
    resource_split = (
        resource_id.split("/") if "/" in resource_id else resource_id.split(":")
    )
    if (
        resource_split[0] == "endpoint"
        or (resource_split[0] == "table" and len(resource_split) > 2)
        or (resource_split[0] == "keyspace")
    ):
        resource_type = resource_split[2]
    else:
        resource_type = resource_split[0]
    return resource_type


class FakeScalableTarget(BaseModel):
    def __init__(
        self,
        backend: ApplicationAutoscalingBackend,
        service_namespace: str,
        resource_id: str,
        scalable_dimension: str,
        min_capacity: Optional[int],
        max_capacity: Optional[int],
        role_arn: str,
        suspended_state: str,
    ) -> None:
        self.applicationautoscaling_backend = backend
        self.service_namespace = service_namespace
        self.resource_id = resource_id
        self.scalable_dimension = scalable_dimension
        self.min_capacity = min_capacity
        self.max_capacity = max_capacity
        self.role_arn = role_arn
        self.suspended_state = suspended_state
        self.creation_time = time.time()

    def update(
        self,
        min_capacity: Optional[int],
        max_capacity: Optional[int],
        suspended_state: str,
    ) -> None:
        if min_capacity is not None:
            self.min_capacity = min_capacity
        if max_capacity is not None:
            self.max_capacity = max_capacity
        if suspended_state is not None:
            self.suspended_state = suspended_state


class FakeApplicationAutoscalingPolicy(BaseModel):
    def __init__(
        self,
        account_id: str,
        region_name: str,
        policy_name: str,
        service_namespace: str,
        resource_id: str,
        scalable_dimension: str,
        policy_type: Optional[str],
        policy_body: Dict[str, Any],
    ) -> None:
        self.step_scaling_policy_configuration = None
        self.target_tracking_scaling_policy_configuration = None

        if policy_type == "StepScaling":
            self.step_scaling_policy_configuration = policy_body
            self.target_tracking_scaling_policy_configuration = None
        elif policy_type == "TargetTrackingScaling":
            self.step_scaling_policy_configuration = None
            self.target_tracking_scaling_policy_configuration = policy_body
        else:
            raise AWSValidationException(
                f"1 validation error detected: Value '{policy_type}' at 'policyType' failed to satisfy constraint: Member must satisfy enum value set: [PredictiveScaling, StepScaling, TargetTrackingScaling]"
            )

        self._policy_body = policy_body
        self.service_namespace = service_namespace
        self.resource_id = resource_id
        self.scalable_dimension = scalable_dimension
        self.policy_name = policy_name
        self.policy_type = policy_type
        self._guid = mock_random.uuid4()
        self.policy_arn = f"arn:{get_partition(region_name)}:autoscaling:{region_name}:{account_id}:scalingPolicy:{self._guid}:resource/{self.service_namespace}/{self.resource_id}:policyName/{self.policy_name}"
        self.creation_time = time.time()
        self.alarms: List["FakeAlarm"] = []

        self.account_id = account_id
        self.region_name = region_name

        self.create_alarms()

    def create_alarms(self) -> None:
        if self.policy_type == "TargetTrackingScaling":
            if self.service_namespace == "dynamodb":
                self.alarms.extend(self._generate_dynamodb_alarms())
            if self.service_namespace == "ecs":
                self.alarms.extend(self._generate_ecs_alarms())

    def _generate_dynamodb_alarms(self) -> List["FakeAlarm"]:
        from moto.cloudwatch.models import CloudWatchBackend, cloudwatch_backends

        cloudwatch: CloudWatchBackend = cloudwatch_backends[self.account_id][
            self.region_name
        ]
        alarms = []
        table_name = self.resource_id.split("/")[-1]
        alarm_action = f"{self.policy_arn}:createdBy/{mock_random.uuid4()}"
        alarm1 = cloudwatch.put_metric_alarm(
            name=f"TargetTracking-table/{table_name}-AlarmHigh-{mock_random.uuid4()}",
            namespace="AWS/DynamoDB",
            metric_name="ConsumedReadCapacityUnits",
            metric_data_queries=[],
            comparison_operator="GreaterThanThreshold",
            evaluation_periods=2,
            period=60,
            threshold=42.0,
            statistic="Sum",
            description=f"DO NOT EDIT OR DELETE. For TargetTrackingScaling policy {alarm_action}",
            dimensions=[{"name": "TableName", "value": table_name}],
            alarm_actions=[alarm_action],
        )
        alarms.append(alarm1)
        alarm2 = cloudwatch.put_metric_alarm(
            name=f"TargetTracking-table/{table_name}-AlarmLow-{mock_random.uuid4()}",
            namespace="AWS/DynamoDB",
            metric_name="ConsumedReadCapacityUnits",
            metric_data_queries=[],
            comparison_operator="LessThanThreshold",
            evaluation_periods=15,
            period=60,
            threshold=30.0,
            statistic="Sum",
            description=f"DO NOT EDIT OR DELETE. For TargetTrackingScaling policy {alarm_action}",
            dimensions=[{"name": "TableName", "value": table_name}],
            alarm_actions=[alarm_action],
        )
        alarms.append(alarm2)
        alarm3 = cloudwatch.put_metric_alarm(
            name=f"TargetTracking-table/{table_name}-ProvisionedCapacityHigh-{mock_random.uuid4()}",
            namespace="AWS/DynamoDB",
            metric_name="ProvisionedReadCapacityUnits",
            metric_data_queries=[],
            comparison_operator="GreaterThanThreshold",
            evaluation_periods=2,
            period=300,
            threshold=1.0,
            statistic="Average",
            description=f"DO NOT EDIT OR DELETE. For TargetTrackingScaling policy {alarm_action}",
            dimensions=[{"name": "TableName", "value": table_name}],
            alarm_actions=[alarm_action],
        )
        alarms.append(alarm3)
        alarm4 = cloudwatch.put_metric_alarm(
            name=f"TargetTracking-table/{table_name}-ProvisionedCapacityLow-{mock_random.uuid4()}",
            namespace="AWS/DynamoDB",
            metric_name="ProvisionedReadCapacityUnits",
            metric_data_queries=[],
            comparison_operator="LessThanThreshold",
            evaluation_periods=3,
            period=300,
            threshold=1.0,
            statistic="Average",
            description=f"DO NOT EDIT OR DELETE. For TargetTrackingScaling policy {alarm_action}",
            dimensions=[{"name": "TableName", "value": table_name}],
            alarm_actions=[alarm_action],
        )
        alarms.append(alarm4)
        return alarms

    def _generate_ecs_alarms(self) -> List["FakeAlarm"]:
        from moto.cloudwatch.models import CloudWatchBackend, cloudwatch_backends

        cloudwatch: CloudWatchBackend = cloudwatch_backends[self.account_id][
            self.region_name
        ]
        alarms: List["FakeAlarm"] = []
        alarm_action = f"{self.policy_arn}:createdBy/{mock_random.uuid4()}"
        config = self.target_tracking_scaling_policy_configuration or {}
        metric_spec = config.get("PredefinedMetricSpecification", {})
        if "Memory" in metric_spec.get("PredefinedMetricType", ""):
            metric_name = "MemoryUtilization"
        else:
            metric_name = "CPUUtilization"
        _, cluster_name, service_name = self.resource_id.split("/")
        alarm1 = cloudwatch.put_metric_alarm(
            name=f"TargetTracking-{self.resource_id}-AlarmHigh-{mock_random.uuid4()}",
            namespace="AWS/ECS",
            metric_name=metric_name,
            metric_data_queries=[],
            comparison_operator="GreaterThanThreshold",
            evaluation_periods=3,
            period=60,
            threshold=6,
            unit="Percent",
            statistic="Average",
            description=f"DO NOT EDIT OR DELETE. For TargetTrackingScaling policy {alarm_action}",
            dimensions=[
                {"name": "ClusterName", "value": cluster_name},
                {"name": "ServiceName", "value": service_name},
            ],
            alarm_actions=[alarm_action],
        )
        alarms.append(alarm1)
        alarm2 = cloudwatch.put_metric_alarm(
            name=f"TargetTracking-{self.resource_id}-AlarmLow-{mock_random.uuid4()}",
            namespace="AWS/ECS",
            metric_name=metric_name,
            metric_data_queries=[],
            comparison_operator="LessThanThreshold",
            evaluation_periods=15,
            period=60,
            threshold=6,
            unit="Percent",
            statistic="Average",
            description=f"DO NOT EDIT OR DELETE. For TargetTrackingScaling policy {alarm_action}",
            dimensions=[
                {"name": "ClusterName", "value": cluster_name},
                {"name": "ServiceName", "value": service_name},
            ],
            alarm_actions=[alarm_action],
        )
        alarms.append(alarm2)
        return alarms

    def delete_alarms(self, account_id: str, region_name: str) -> None:
        from moto.cloudwatch.models import CloudWatchBackend, cloudwatch_backends

        cloudwatch: CloudWatchBackend = cloudwatch_backends[account_id][region_name]
        cloudwatch.delete_alarms([a.name for a in self.alarms])

    @staticmethod
    def formulate_key(
        service_namespace: str,
        resource_id: str,
        scalable_dimension: str,
        policy_name: str,
    ) -> str:
        return (
            f"{service_namespace}\t{resource_id}\t{scalable_dimension}\t{policy_name}"
        )


class FakeScheduledAction(BaseModel):
    def __init__(
        self,
        service_namespace: str,
        schedule: str,
        timezone: str,
        scheduled_action_name: str,
        resource_id: str,
        scalable_dimension: str,
        start_time: str,
        end_time: str,
        scalable_target_action: str,
        account_id: str,
        region: str,
    ) -> None:
        self.arn = f"arn:{get_partition(region)}:autoscaling:{region}:{account_id}:scheduledAction:{service_namespace}/{resource_id}:scheduledActionName/{scheduled_action_name}"
        self.service_namespace = service_namespace
        self.schedule = schedule
        self.timezone = timezone
        self.scheduled_action_name = scheduled_action_name
        self.resource_id = resource_id
        self.scalable_dimension = scalable_dimension
        self.start_time = start_time
        self.end_time = end_time
        self.scalable_target_action = scalable_target_action
        self.creation_time = time.time()

    def update(
        self,
        schedule: str,
        timezone: str,
        scheduled_action_name: str,
        start_time: str,
        end_time: str,
        scalable_target_action: str,
    ) -> None:
        if scheduled_action_name:
            self.scheduled_action_name = scheduled_action_name
        if schedule:
            self.schedule = schedule
        if timezone:
            self.timezone = timezone
        if scalable_target_action:
            self.scalable_target_action = scalable_target_action
        self.start_time = start_time
        self.end_time = end_time


applicationautoscaling_backends = BackendDict(
    ApplicationAutoscalingBackend, "application-autoscaling"
)
