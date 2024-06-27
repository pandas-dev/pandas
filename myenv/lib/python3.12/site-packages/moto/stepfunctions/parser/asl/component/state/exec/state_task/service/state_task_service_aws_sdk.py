from botocore.exceptions import ClientError

from moto.stepfunctions.parser.api import HistoryEventType, TaskFailedEventDetails
from moto.stepfunctions.parser.asl.component.common.error_name.failure_event import (
    FailureEvent,
)
from moto.stepfunctions.parser.asl.component.common.error_name.states_error_name import (
    StatesErrorName,
)
from moto.stepfunctions.parser.asl.component.common.error_name.states_error_name_type import (
    StatesErrorNameType,
)
from moto.stepfunctions.parser.asl.component.state.exec.state_task.service.resource import (
    ResourceRuntimePart,
)
from moto.stepfunctions.parser.asl.component.state.exec.state_task.service.state_task_service_callback import (
    StateTaskServiceCallback,
)
from moto.stepfunctions.parser.asl.component.state.state_props import StateProps
from moto.stepfunctions.parser.asl.eval.environment import Environment
from moto.stepfunctions.parser.asl.eval.event.event_detail import EventDetails
from moto.stepfunctions.parser.asl.utils.boto_client import boto_client_for


class StateTaskServiceAwsSdk(StateTaskServiceCallback):
    _NORMALISED_SERVICE_NAMES = {"dynamodb": "DynamoDb"}

    def from_state_props(self, state_props: StateProps) -> None:
        super().from_state_props(state_props=state_props)

    def _get_sfn_resource_type(self) -> str:
        return f"{self.resource.service_name}:{self.resource.api_name}"

    @staticmethod
    def _normalise_service_name(service_name: str) -> str:
        service_name_lower = service_name.lower()
        if service_name_lower in StateTaskServiceAwsSdk._NORMALISED_SERVICE_NAMES:
            return StateTaskServiceAwsSdk._NORMALISED_SERVICE_NAMES[service_name_lower]
        return service_name_lower

    @staticmethod
    def _normalise_exception_name(norm_service_name: str, ex: Exception) -> str:
        ex_name = ex.__class__.__name__
        norm_ex_name = f"{norm_service_name}.{norm_service_name if ex_name == 'ClientError' else ex_name}"
        if not norm_ex_name.endswith("Exception"):
            norm_ex_name += "Exception"
        return norm_ex_name

    def _get_task_failure_event(self, error: str, cause: str) -> FailureEvent:
        return FailureEvent(
            error_name=StatesErrorName(typ=StatesErrorNameType.StatesTaskFailed),
            event_type=HistoryEventType.TaskFailed,
            event_details=EventDetails(
                taskFailedEventDetails=TaskFailedEventDetails(
                    resource=self._get_sfn_resource(),
                    resourceType=self._get_sfn_resource_type(),
                    error=error,
                    cause=cause,
                )
            ),
        )

    def _from_error(self, env: Environment, ex: Exception) -> FailureEvent:
        if isinstance(ex, ClientError):
            norm_service_name: str = self._normalise_service_name(
                self.resource.api_name
            )
            error: str = self._normalise_exception_name(norm_service_name, ex)

            error_message: str = ex.response["Error"]["Message"]
            cause_details = [
                f"Service: {norm_service_name}",
                f"Status Code: {ex.response['ResponseMetadata']['HTTPStatusCode']}",
                f"Request ID: {ex.response['ResponseMetadata']['RequestId']}",
            ]
            if "HostId" in ex.response["ResponseMetadata"]:
                cause_details.append(
                    f'Extended Request ID: {ex.response["ResponseMetadata"]["HostId"]}'
                )

            cause: str = f"{error_message} ({', '.join(cause_details)})"
            failure_event = self._get_task_failure_event(error=error, cause=cause)
            return failure_event
        return super()._from_error(env=env, ex=ex)

    def _eval_service_task(
        self,
        env: Environment,
        resource_runtime_part: ResourceRuntimePart,
        normalised_parameters: dict,
    ):
        service_name = self._get_boto_service_name()
        api_action = self._get_boto_service_action()
        api_client = boto_client_for(
            region=resource_runtime_part.region,
            account=resource_runtime_part.account,
            service=service_name,
        )
        response = getattr(api_client, api_action)(**normalised_parameters) or dict()

        env.stack.append(response)
