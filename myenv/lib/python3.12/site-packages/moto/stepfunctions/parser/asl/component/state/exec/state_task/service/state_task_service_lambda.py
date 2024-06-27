from typing import Dict, Final, Optional, Set, Tuple

from botocore.exceptions import ClientError

from moto.stepfunctions.parser.api import HistoryEventType, TaskFailedEventDetails
from moto.stepfunctions.parser.asl.component.common.error_name.custom_error_name import (
    CustomErrorName,
)
from moto.stepfunctions.parser.asl.component.common.error_name.failure_event import (
    FailureEvent,
)
from moto.stepfunctions.parser.asl.component.state.exec.state_task import (
    lambda_eval_utils,
)
from moto.stepfunctions.parser.asl.component.state.exec.state_task.service.resource import (
    ResourceRuntimePart,
)
from moto.stepfunctions.parser.asl.component.state.exec.state_task.service.state_task_service_callback import (
    StateTaskServiceCallback,
)
from moto.stepfunctions.parser.asl.eval.environment import Environment
from moto.stepfunctions.parser.asl.eval.event.event_detail import EventDetails


class StateTaskServiceLambda(StateTaskServiceCallback):
    _SUPPORTED_API_PARAM_BINDINGS: Final[Dict[str, Set[str]]] = {
        "invoke": {
            "ClientContext",
            "FunctionName",
            "InvocationType",
            "Qualifier",
            "Payload",
            # Outside the specification, but supported in practice:
            "LogType",
        }
    }

    def _get_supported_parameters(self) -> Optional[Set[str]]:
        return self._SUPPORTED_API_PARAM_BINDINGS.get(self.resource.api_action.lower())

    @staticmethod
    def _error_cause_from_client_error(client_error: ClientError) -> Tuple[str, str]:
        error_code: str = client_error.response["Error"]["Code"]
        error_msg: str = client_error.response["Error"]["Message"]
        response_details = "; ".join(
            [
                "Service: AWSLambda",
                f"Status Code: {client_error.response['ResponseMetadata']['HTTPStatusCode']}",
                f"Error Code: {error_code}",
                f"Request ID: {client_error.response['ResponseMetadata']['RequestId']}",
                "Proxy: null",
            ]
        )
        error = f"Lambda.{error_code}"
        cause = f"{error_msg} ({response_details})"
        return error, cause

    def _from_error(self, env: Environment, ex: Exception) -> FailureEvent:
        if isinstance(ex, lambda_eval_utils.LambdaFunctionErrorException):
            error = "Exception"
            error_name = CustomErrorName(error)
            cause = ex.payload
        elif isinstance(ex, ClientError):
            error, cause = self._error_cause_from_client_error(ex)
            error_name = CustomErrorName(error)
        else:
            return super()._from_error(env=env, ex=ex)
        return FailureEvent(
            error_name=error_name,
            event_type=HistoryEventType.TaskFailed,
            event_details=EventDetails(
                taskFailedEventDetails=TaskFailedEventDetails(
                    error=error,
                    cause=cause,
                    resource=self._get_sfn_resource(),
                    resourceType=self._get_sfn_resource_type(),
                )
            ),
        )

    def _eval_service_task(
        self,
        env: Environment,
        resource_runtime_part: ResourceRuntimePart,
        normalised_parameters: dict,
    ):
        if "Payload" in normalised_parameters:
            normalised_parameters["Payload"] = lambda_eval_utils.to_payload_type(
                normalised_parameters["Payload"]
            )
        lambda_eval_utils.exec_lambda_function(
            env=env,
            parameters=normalised_parameters,
            region=resource_runtime_part.region,
            account=resource_runtime_part.account,
        )
