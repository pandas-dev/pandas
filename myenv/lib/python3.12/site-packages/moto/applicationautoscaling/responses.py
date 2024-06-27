import json
from typing import Any, Dict, List

from moto.core.responses import BaseResponse

from .exceptions import AWSValidationException
from .models import (
    ApplicationAutoscalingBackend,
    FakeApplicationAutoscalingPolicy,
    FakeScalableTarget,
    FakeScheduledAction,
    ScalableDimensionValueSet,
    ServiceNamespaceValueSet,
    applicationautoscaling_backends,
)


class ApplicationAutoScalingResponse(BaseResponse):
    def __init__(self) -> None:
        super().__init__(service_name="application-autoscaling")

    @property
    def applicationautoscaling_backend(self) -> ApplicationAutoscalingBackend:
        return applicationautoscaling_backends[self.current_account][self.region]

    def describe_scalable_targets(self) -> str:
        self._validate_params()
        service_namespace = self._get_param("ServiceNamespace")
        resource_ids = self._get_param("ResourceIds")
        scalable_dimension = self._get_param("ScalableDimension")
        max_results = self._get_int_param("MaxResults") or 50
        marker = self._get_param("NextToken")
        all_scalable_targets = (
            self.applicationautoscaling_backend.describe_scalable_targets(
                service_namespace, resource_ids, scalable_dimension
            )
        )
        start = int(marker) + 1 if marker else 0
        next_token = None
        scalable_targets_resp = all_scalable_targets[start : start + max_results]
        if len(all_scalable_targets) > start + max_results:
            next_token = str(len(scalable_targets_resp) - 1)
        targets = [_build_target(t) for t in scalable_targets_resp]
        return json.dumps({"ScalableTargets": targets, "NextToken": next_token})

    def register_scalable_target(self) -> str:
        """Registers or updates a scalable target."""
        self._validate_params()
        self.applicationautoscaling_backend.register_scalable_target(
            self._get_param("ServiceNamespace"),
            self._get_param("ResourceId"),
            self._get_param("ScalableDimension"),
            min_capacity=self._get_int_param("MinCapacity"),
            max_capacity=self._get_int_param("MaxCapacity"),
            role_arn=self._get_param("RoleARN"),
            suspended_state=self._get_param("SuspendedState"),
        )
        return json.dumps({})

    def deregister_scalable_target(self) -> str:
        """Deregisters a scalable target."""
        self._validate_params()
        self.applicationautoscaling_backend.deregister_scalable_target(
            self._get_param("ServiceNamespace"),
            self._get_param("ResourceId"),
            self._get_param("ScalableDimension"),
        )
        return json.dumps({})

    def put_scaling_policy(self) -> str:
        policy = self.applicationautoscaling_backend.put_scaling_policy(
            policy_name=self._get_param("PolicyName"),
            service_namespace=self._get_param("ServiceNamespace"),
            resource_id=self._get_param("ResourceId"),
            scalable_dimension=self._get_param("ScalableDimension"),
            policy_type=self._get_param("PolicyType"),
            policy_body=self._get_param(
                "StepScalingPolicyConfiguration",
                self._get_param("TargetTrackingScalingPolicyConfiguration"),
            ),
        )
        return json.dumps(
            {"PolicyARN": policy.policy_arn, "Alarms": _build_alarms(policy)}
        )

    def describe_scaling_policies(self) -> str:
        (
            next_token,
            policy_page,
        ) = self.applicationautoscaling_backend.describe_scaling_policies(
            service_namespace=self._get_param("ServiceNamespace"),
            resource_id=self._get_param("ResourceId"),
            scalable_dimension=self._get_param("ScalableDimension"),
            max_results=self._get_int_param("MaxResults"),
            next_token=self._get_param("NextToken"),
        )
        response_obj = {
            "ScalingPolicies": [_build_policy(p) for p in policy_page],
            "NextToken": next_token,
        }
        return json.dumps(response_obj)

    def delete_scaling_policy(self) -> str:
        self.applicationautoscaling_backend.delete_scaling_policy(
            policy_name=self._get_param("PolicyName"),
            service_namespace=self._get_param("ServiceNamespace"),
            resource_id=self._get_param("ResourceId"),
            scalable_dimension=self._get_param("ScalableDimension"),
        )
        return json.dumps({})

    def _validate_params(self) -> None:
        """Validate parameters.
        TODO Integrate this validation with the validation in models.py
        """
        namespace = self._get_param("ServiceNamespace")
        dimension = self._get_param("ScalableDimension")
        messages = []
        dimensions = [d.value for d in ScalableDimensionValueSet]
        message = None
        if dimension is not None and dimension not in dimensions:
            messages.append(
                f"Value '{dimension}' at 'scalableDimension' failed to satisfy constraint: Member must satisfy enum value set: {dimensions}"
            )
        namespaces = [n.value for n in ServiceNamespaceValueSet]
        if namespace is not None and namespace not in namespaces:
            messages.append(
                f"Value '{namespace}' at 'serviceNamespace' failed to satisfy constraint: Member must satisfy enum value set: {namespaces}"
            )
        if len(messages) == 1:
            message = f"1 validation error detected: {messages[0]}"
        elif len(messages) > 1:
            message = (
                f'{len(messages)} validation errors detected: {"; ".join(messages)}'
            )
        if message:
            raise AWSValidationException(message)

    def delete_scheduled_action(self) -> str:
        params = json.loads(self.body)
        service_namespace = params.get("ServiceNamespace")
        scheduled_action_name = params.get("ScheduledActionName")
        resource_id = params.get("ResourceId")
        scalable_dimension = params.get("ScalableDimension")
        self.applicationautoscaling_backend.delete_scheduled_action(
            service_namespace=service_namespace,
            scheduled_action_name=scheduled_action_name,
            resource_id=resource_id,
            scalable_dimension=scalable_dimension,
        )
        return json.dumps(dict())

    def put_scheduled_action(self) -> str:
        params = json.loads(self.body)
        service_namespace = params.get("ServiceNamespace")
        schedule = params.get("Schedule")
        timezone = params.get("Timezone")
        scheduled_action_name = params.get("ScheduledActionName")
        resource_id = params.get("ResourceId")
        scalable_dimension = params.get("ScalableDimension")
        start_time = params.get("StartTime")
        end_time = params.get("EndTime")
        scalable_target_action = params.get("ScalableTargetAction")
        self.applicationautoscaling_backend.put_scheduled_action(
            service_namespace=service_namespace,
            schedule=schedule,
            timezone=timezone,
            scheduled_action_name=scheduled_action_name,
            resource_id=resource_id,
            scalable_dimension=scalable_dimension,
            start_time=start_time,
            end_time=end_time,
            scalable_target_action=scalable_target_action,
        )
        return json.dumps(dict())

    def describe_scheduled_actions(self) -> str:
        params = json.loads(self.body)
        scheduled_action_names = params.get("ScheduledActionNames")
        service_namespace = params.get("ServiceNamespace")
        resource_id = params.get("ResourceId")
        scalable_dimension = params.get("ScalableDimension")
        scheduled_actions = (
            self.applicationautoscaling_backend.describe_scheduled_actions(
                scheduled_action_names=scheduled_action_names,
                service_namespace=service_namespace,
                resource_id=resource_id,
                scalable_dimension=scalable_dimension,
            )
        )
        response_obj = {
            "ScheduledActions": [_build_scheduled_action(a) for a in scheduled_actions]
        }
        return json.dumps(response_obj)


def _build_target(t: FakeScalableTarget) -> Dict[str, Any]:
    return {
        "CreationTime": t.creation_time,
        "ServiceNamespace": t.service_namespace,
        "ResourceId": t.resource_id,
        "RoleARN": t.role_arn,
        "ScalableDimension": t.scalable_dimension,
        "MaxCapacity": t.max_capacity,
        "MinCapacity": t.min_capacity,
        "SuspendedState": t.suspended_state,
    }


def _build_alarms(policy: FakeApplicationAutoscalingPolicy) -> List[Dict[str, str]]:
    return [{"AlarmARN": a.alarm_arn, "AlarmName": a.name} for a in policy.alarms]


def _build_policy(p: FakeApplicationAutoscalingPolicy) -> Dict[str, Any]:
    response = {
        "PolicyARN": p.policy_arn,
        "PolicyName": p.policy_name,
        "ServiceNamespace": p.service_namespace,
        "ResourceId": p.resource_id,
        "ScalableDimension": p.scalable_dimension,
        "PolicyType": p.policy_type,
        "CreationTime": p.creation_time,
        "Alarms": _build_alarms(p),
    }
    if p.policy_type == "StepScaling":
        response["StepScalingPolicyConfiguration"] = p.step_scaling_policy_configuration
    elif p.policy_type == "TargetTrackingScaling":
        response["TargetTrackingScalingPolicyConfiguration"] = (
            p.target_tracking_scaling_policy_configuration
        )
    return response


def _build_scheduled_action(a: FakeScheduledAction) -> Dict[str, Any]:
    response = {
        "ScheduledActionName": a.scheduled_action_name,
        "ScheduledActionARN": a.arn,
        "ServiceNamespace": a.service_namespace,
        "Schedule": a.schedule,
        "Timezone": a.timezone,
        "ResourceId": a.resource_id,
        "ScalableDimension": a.scalable_dimension,
        "StartTime": a.start_time,
        "EndTime": a.end_time,
        "CreationTime": a.creation_time,
        "ScalableTargetAction": a.scalable_target_action,
    }
    return response
