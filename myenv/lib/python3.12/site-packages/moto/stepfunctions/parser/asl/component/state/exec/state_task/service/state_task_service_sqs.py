from typing import Any, Dict, Final, Optional, Set

from botocore.exceptions import ClientError

from moto.stepfunctions.parser.api import HistoryEventType, TaskFailedEventDetails
from moto.stepfunctions.parser.asl.component.common.error_name.custom_error_name import (
    CustomErrorName,
)
from moto.stepfunctions.parser.asl.component.common.error_name.failure_event import (
    FailureEvent,
)
from moto.stepfunctions.parser.asl.component.state.exec.state_task.service.resource import (
    ResourceRuntimePart,
)
from moto.stepfunctions.parser.asl.component.state.exec.state_task.service.state_task_service_callback import (
    StateTaskServiceCallback,
)
from moto.stepfunctions.parser.asl.eval.environment import Environment
from moto.stepfunctions.parser.asl.eval.event.event_detail import EventDetails
from moto.stepfunctions.parser.asl.utils.boto_client import boto_client_for
from moto.stepfunctions.parser.asl.utils.encoding import to_json_str


class StateTaskServiceSqs(StateTaskServiceCallback):
    _ERROR_NAME_CLIENT: Final[str] = "SQS.SdkClientException"
    _ERROR_NAME_AWS: Final[str] = "SQS.AmazonSQSException"

    _SUPPORTED_API_PARAM_BINDINGS: Final[Dict[str, Set[str]]] = {
        "sendmessage": {
            "DelaySeconds",
            "MessageAttribute",
            "MessageBody",
            "MessageDeduplicationId",
            "MessageGroupId",
            "QueueUrl",
        }
    }

    def _get_supported_parameters(self) -> Optional[Set[str]]:
        return self._SUPPORTED_API_PARAM_BINDINGS.get(self.resource.api_action.lower())

    def _from_error(self, env: Environment, ex: Exception) -> FailureEvent:
        if isinstance(ex, ClientError):
            return FailureEvent(
                error_name=CustomErrorName(self._ERROR_NAME_CLIENT),
                event_type=HistoryEventType.TaskFailed,
                event_details=EventDetails(
                    taskFailedEventDetails=TaskFailedEventDetails(
                        error=self._ERROR_NAME_CLIENT,
                        cause=ex.response["Error"][
                            "Message"
                        ],  # TODO: update to report expected cause.
                        resource=self._get_sfn_resource(),
                        resourceType=self._get_sfn_resource_type(),
                    )
                ),
            )
        return super()._from_error(env=env, ex=ex)

    def _normalise_response(
        self,
        response: Any,
        boto_service_name: Optional[str] = None,
        service_action_name: Optional[str] = None,
    ) -> None:
        super()._normalise_response(
            response=response,
            boto_service_name=boto_service_name,
            service_action_name=service_action_name,
        )
        # Normalise output value key to SFN standard for Md5OfMessageBody.
        if response and "Md5OfMessageBody" in response:
            md5_message_body = response.pop("Md5OfMessageBody")
            response["MD5OfMessageBody"] = md5_message_body

    def _eval_service_task(
        self,
        env: Environment,
        resource_runtime_part: ResourceRuntimePart,
        normalised_parameters: dict,
    ):
        # TODO: Stepfunctions automatically dumps to json MessageBody's definitions.
        #  Are these other similar scenarios?
        if "MessageBody" in normalised_parameters:
            message_body = normalised_parameters["MessageBody"]
            if message_body is not None and not isinstance(message_body, str):
                normalised_parameters["MessageBody"] = to_json_str(message_body)

        service_name = self._get_boto_service_name()
        api_action = self._get_boto_service_action()
        sqs_client = boto_client_for(
            region=resource_runtime_part.region,
            account=resource_runtime_part.account,
            service=service_name,
        )
        response = getattr(sqs_client, api_action)(**normalised_parameters)
        response.pop("ResponseMetadata", None)
        env.stack.append(response)
