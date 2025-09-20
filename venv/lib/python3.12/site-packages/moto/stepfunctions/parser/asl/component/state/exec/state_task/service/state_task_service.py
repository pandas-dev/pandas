from __future__ import annotations

import abc
import copy
import logging
from typing import Any, Final, List, Optional

from botocore.model import ListShape, Shape, StringShape, StructureShape
from botocore.response import StreamingBody

from moto.stepfunctions.parser.api import (
    HistoryEventExecutionDataDetails,
    HistoryEventType,
    TaskCredentials,
    TaskFailedEventDetails,
    TaskScheduledEventDetails,
    TaskStartedEventDetails,
    TaskSucceededEventDetails,
    TaskTimedOutEventDetails,
)
from moto.stepfunctions.parser.asl.component.common.error_name.failure_event import (
    FailureEvent,
    FailureEventException,
)
from moto.stepfunctions.parser.asl.component.common.error_name.states_error_name import (
    StatesErrorName,
)
from moto.stepfunctions.parser.asl.component.common.error_name.states_error_name_type import (
    StatesErrorNameType,
)
from moto.stepfunctions.parser.asl.component.state.exec.state_task.credentials import (
    StateCredentials,
)
from moto.stepfunctions.parser.asl.component.state.exec.state_task.mock_eval_utils import (
    eval_mocked_response,
)
from moto.stepfunctions.parser.asl.component.state.exec.state_task.service.resource import (
    ResourceRuntimePart,
    ServiceResource,
)
from moto.stepfunctions.parser.asl.component.state.exec.state_task.state_task import (
    StateTask,
)
from moto.stepfunctions.parser.asl.component.state.state_props import StateProps
from moto.stepfunctions.parser.asl.eval.environment import Environment
from moto.stepfunctions.parser.asl.eval.event.event_detail import EventDetails
from moto.stepfunctions.parser.asl.utils.encoding import to_json_str
from moto.stepfunctions.parser.mocking.mock_config import MockedResponse
from moto.stepfunctions.parser.quotas import is_within_size_quota
from moto.stepfunctions.parser.utils import (
    camel_to_snake_case,
    snake_to_camel_case,
    to_bytes,
    to_str,
)

LOG = logging.getLogger(__name__)


class StateTaskService(StateTask, abc.ABC):
    resource: ServiceResource

    _SERVICE_NAME_SFN_TO_BOTO_OVERRIDES: Final[dict[str, str]] = {
        "sfn": "stepfunctions",
        "states": "stepfunctions",
    }

    def from_state_props(self, state_props: StateProps) -> None:
        super().from_state_props(state_props=state_props)
        # Validate the service integration is supported on program creation.
        self._validate_service_integration_is_supported()

    def _validate_service_integration_is_supported(self):
        # Validate the service integration is supported.
        supported_parameters = self._get_supported_parameters()
        if supported_parameters is None:
            raise ValueError(
                f"The resource provided {self.resource.resource_arn} not recognized. "
                "The value is not a valid resource ARN, or the resource is not available in this region."
            )

    def _get_sfn_resource(self) -> str:
        return self.resource.api_action

    def _get_sfn_resource_type(self) -> str:
        return self.resource.service_name

    def _get_timed_out_failure_event(self, env: Environment) -> FailureEvent:
        return FailureEvent(
            env=env,
            error_name=StatesErrorName(typ=StatesErrorNameType.StatesTimeout),
            event_type=HistoryEventType.TaskTimedOut,
            event_details=EventDetails(
                taskTimedOutEventDetails=TaskTimedOutEventDetails(
                    resourceType=self._get_sfn_resource_type(),
                    resource=self._get_sfn_resource(),
                    error=StatesErrorNameType.StatesTimeout.to_name(),
                )
            ),
        )

    def _to_boto_request_value(self, request_value: Any, value_shape: Shape) -> Any:
        boto_request_value = request_value
        if isinstance(value_shape, StructureShape):
            self._to_boto_request(request_value, value_shape)
        elif isinstance(value_shape, ListShape) and isinstance(request_value, list):
            for request_list_value in request_value:
                self._to_boto_request_value(request_list_value, value_shape.member)  # noqa
        elif isinstance(value_shape, StringShape) and not isinstance(
            request_value, str
        ):
            boto_request_value = to_json_str(request_value)
        elif value_shape.type_name == "blob" and not isinstance(
            boto_request_value, bytes
        ):
            boto_request_value = to_json_str(request_value, separators=(",", ":"))
            boto_request_value = to_bytes(boto_request_value)
        return boto_request_value

    def _to_boto_request(
        self, parameters: dict, structure_shape: StructureShape
    ) -> None:
        if not isinstance(structure_shape, StructureShape):
            LOG.warning(
                "Step Functions could not normalise the request for integration '%s' due to the unexpected request template value of type '%s'",
                self.resource.resource_arn,
                type(structure_shape),
            )
            return
        shape_members = structure_shape.members
        norm_member_binds: dict[str, tuple[str, StructureShape]] = {
            camel_to_snake_case(member_key): (member_key, member_value)
            for member_key, member_value in shape_members.items()
        }
        parameters_bind_keys: List[str] = list(parameters.keys())
        for parameter_key in parameters_bind_keys:
            norm_parameter_key = camel_to_snake_case(parameter_key)
            norm_member_bind: Optional[tuple[str, Optional[StructureShape]]] = (
                norm_member_binds.get(norm_parameter_key)
            )
            if norm_member_bind is not None:
                norm_member_bind_key, norm_member_bind_shape = norm_member_bind
                parameter_value = parameters.pop(parameter_key)
                parameter_value = self._to_boto_request_value(
                    parameter_value, norm_member_bind_shape
                )
                parameters[norm_member_bind_key] = parameter_value

    @staticmethod
    def _to_sfn_cased(member_key: str) -> str:
        # Normalise the string to snake case, e.g. "HelloWorld_hello__world" -> "hello_world_hello_world"
        norm_member_key = camel_to_snake_case(member_key)
        # Normalise the snake case to camel case, e.g. "hello_world_hello_world" -> "HelloWorldHelloWorld"
        norm_member_key = snake_to_camel_case(norm_member_key)
        return norm_member_key

    @staticmethod
    def _from_boto_response_value(response_value: Any) -> Any:
        if isinstance(response_value, StreamingBody):
            body_str = to_str(response_value.read())
            return body_str
        return response_value

    def _from_boto_response(
        self, response: Any, structure_shape: StructureShape
    ) -> None:
        if not isinstance(response, dict):
            return

        if not isinstance(structure_shape, StructureShape):
            LOG.warning(
                "Step Functions could not normalise the response of integration '%s' due to the unexpected request template value of type '%s'",
                self.resource.resource_arn,
                type(structure_shape),
            )
            return

        shape_members = structure_shape.members
        response_bind_keys: List[str] = list(response.keys())
        for response_key in response_bind_keys:
            norm_response_key = self._to_sfn_cased(response_key)
            if response_key in shape_members:
                shape_member = shape_members[response_key]

                response_value = response.pop(response_key)
                response_value = self._from_boto_response_value(response_value)

                if isinstance(shape_member, StructureShape):
                    self._from_boto_response(response_value, shape_member)
                elif isinstance(shape_member, ListShape) and isinstance(
                    response_value, list
                ):
                    for response_value_member in response_value:
                        self._from_boto_response(
                            response_value_member, shape_member.member
                        )  # noqa

                response[norm_response_key] = response_value

    def _get_boto_service_name(self, boto_service_name: Optional[str] = None) -> str:
        api_name = boto_service_name or self.resource.api_name
        return self._SERVICE_NAME_SFN_TO_BOTO_OVERRIDES.get(api_name, api_name)

    def _get_boto_service_action(
        self, service_action_name: Optional[str] = None
    ) -> str:
        api_action = service_action_name or self.resource.api_action
        return camel_to_snake_case(api_action)

    def _normalise_parameters(
        self,
        parameters: dict,
        boto_service_name: Optional[str] = None,
        service_action_name: Optional[str] = None,
    ) -> None:
        pass

    def _normalise_response(
        self,
        response: Any,
        boto_service_name: Optional[str] = None,
        service_action_name: Optional[str] = None,
    ) -> None:
        pass

    def _verify_size_quota(self, env: Environment, value: Any) -> None:
        is_within: bool = is_within_size_quota(value)
        if is_within:
            return
        resource_type = self._get_sfn_resource_type()
        resource = self._get_sfn_resource()
        cause = (
            f"The state/task '{resource_type}' returned a result with a size "
            "exceeding the maximum number of bytes service limit."
        )
        raise FailureEventException(
            failure_event=FailureEvent(
                env=env,
                error_name=StatesErrorName(
                    typ=StatesErrorNameType.StatesStatesDataLimitExceeded
                ),
                event_type=HistoryEventType.TaskFailed,
                event_details=EventDetails(
                    taskFailedEventDetails=TaskFailedEventDetails(
                        error=StatesErrorNameType.StatesStatesDataLimitExceeded.to_name(),
                        cause=cause,
                        resourceType=resource_type,
                        resource=resource,
                    )
                ),
            )
        )

    @abc.abstractmethod
    def _eval_service_task(
        self,
        env: Environment,
        resource_runtime_part: ResourceRuntimePart,
        normalised_parameters: dict,
        state_credentials: StateCredentials,
    ): ...

    def _before_eval_execution(
        self,
        env: Environment,
        resource_runtime_part: ResourceRuntimePart,
        raw_parameters: dict,
        state_credentials: StateCredentials,
    ) -> None:
        parameters_str = to_json_str(raw_parameters)

        scheduled_event_details = TaskScheduledEventDetails(
            resource=self._get_sfn_resource(),
            resourceType=self._get_sfn_resource_type(),
            region=resource_runtime_part.region,
            parameters=parameters_str,
        )
        if not self.timeout.is_default_value():
            self.timeout.eval(env=env)
            timeout_seconds = env.stack.pop()
            scheduled_event_details["timeoutInSeconds"] = timeout_seconds
        if self.heartbeat is not None:
            self.heartbeat.eval(env=env)
            heartbeat_seconds = env.stack.pop()
            scheduled_event_details["heartbeatInSeconds"] = heartbeat_seconds
        if self.credentials:
            scheduled_event_details["taskCredentials"] = TaskCredentials(
                roleArn=state_credentials.role_arn
            )
        env.event_manager.add_event(
            context=env.event_history_context,
            event_type=HistoryEventType.TaskScheduled,
            event_details=EventDetails(
                taskScheduledEventDetails=scheduled_event_details
            ),
        )

        env.event_manager.add_event(
            context=env.event_history_context,
            event_type=HistoryEventType.TaskStarted,
            event_details=EventDetails(
                taskStartedEventDetails=TaskStartedEventDetails(
                    resource=self._get_sfn_resource(),
                    resourceType=self._get_sfn_resource_type(),
                )
            ),
        )

    def _after_eval_execution(
        self,
        env: Environment,
        resource_runtime_part: ResourceRuntimePart,
        normalised_parameters: dict,
        state_credentials: StateCredentials,
    ) -> None:
        output = env.stack[-1]
        self._verify_size_quota(env=env, value=output)
        env.event_manager.add_event(
            context=env.event_history_context,
            event_type=HistoryEventType.TaskSucceeded,
            event_details=EventDetails(
                taskSucceededEventDetails=TaskSucceededEventDetails(
                    resource=self._get_sfn_resource(),
                    resourceType=self._get_sfn_resource_type(),
                    output=to_json_str(output),
                    outputDetails=HistoryEventExecutionDataDetails(truncated=False),
                )
            ),
        )

    def _eval_execution(self, env: Environment) -> None:
        self.resource.eval(env=env)
        resource_runtime_part: ResourceRuntimePart = env.stack.pop()

        raw_parameters = self._eval_parameters(env=env)
        state_credentials = self._eval_state_credentials(env=env)

        self._before_eval_execution(
            env=env,
            resource_runtime_part=resource_runtime_part,
            raw_parameters=raw_parameters,
            state_credentials=state_credentials,
        )

        normalised_parameters = copy.deepcopy(raw_parameters)
        self._normalise_parameters(normalised_parameters)

        if env.is_mocked_mode():
            mocked_response: MockedResponse = env.get_current_mocked_response()
            eval_mocked_response(env=env, mocked_response=mocked_response)
        else:
            self._eval_service_task(
                env=env,
                resource_runtime_part=resource_runtime_part,
                normalised_parameters=normalised_parameters,
                state_credentials=state_credentials,
            )

        output_value = env.stack[-1]
        self._normalise_response(output_value)

        self._after_eval_execution(
            env=env,
            resource_runtime_part=resource_runtime_part,
            normalised_parameters=normalised_parameters,
            state_credentials=state_credentials,
        )
