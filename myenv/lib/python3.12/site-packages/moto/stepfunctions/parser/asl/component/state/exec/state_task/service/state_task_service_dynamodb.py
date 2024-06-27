from typing import Dict, Final, Optional, Set, Tuple

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
from moto.stepfunctions.parser.asl.component.state.exec.state_task.service.state_task_service import (
    StateTaskService,
)
from moto.stepfunctions.parser.asl.eval.environment import Environment
from moto.stepfunctions.parser.asl.eval.event.event_detail import EventDetails
from moto.stepfunctions.parser.asl.utils.boto_client import boto_client_for


class StateTaskServiceDynamoDB(StateTaskService):
    _ERROR_NAME_AWS: Final[str] = "DynamoDB.AmazonDynamoDBException"

    _SUPPORTED_API_PARAM_BINDINGS: Final[Dict[str, Set[str]]] = {
        "getitem": {
            "Key",
            "TableName",
            "AttributesToGet",
            "ConsistentRead",
            "ExpressionAttributeNames",
            "ProjectionExpression",
            "ReturnConsumedCapacity",
        },
        "putitem": {
            "Item",
            "TableName",
            "ConditionalOperator",
            "ConditionExpression",
            "Expected",
            "ExpressionAttributeNames",
            "ExpressionAttributeValues",
            "ReturnConsumedCapacity",
            "ReturnItemCollectionMetrics",
            "ReturnValues",
        },
        "deleteitem": {
            "Key",
            "TableName",
            "ConditionalOperator",
            "ConditionExpression",
            "Expected",
            "ExpressionAttributeNames",
            "ExpressionAttributeValues",
            "ReturnConsumedCapacity",
            "ReturnItemCollectionMetrics",
            "ReturnValues",
        },
        "updateitem": {
            "Key",
            "TableName",
            "AttributeUpdates",
            "ConditionalOperator",
            "ConditionExpression",
            "Expected",
            "ExpressionAttributeNames",
            "ExpressionAttributeValues",
            "ReturnConsumedCapacity",
            "ReturnItemCollectionMetrics",
            "ReturnValues",
            "UpdateExpression",
        },
    }

    def _get_supported_parameters(self) -> Optional[Set[str]]:
        return self._SUPPORTED_API_PARAM_BINDINGS.get(self.resource.api_action.lower())

    @staticmethod
    def _error_cause_from_client_error(client_error: ClientError) -> Tuple[str, str]:
        error_code: str = client_error.response["Error"]["Code"]
        error_msg: str = client_error.response["Error"]["Message"]
        response_details = "; ".join(
            [
                "Service: AmazonDynamoDBv2",
                f"Status Code: {client_error.response['ResponseMetadata']['HTTPStatusCode']}",
                f"Error Code: {error_code}",
                f"Request ID: {client_error.response['ResponseMetadata']['RequestId']}",
                "Proxy: null",
            ]
        )
        error = f"DynamoDB.{error_code}"
        cause = f"{error_msg} ({response_details})"
        return error, cause

    def _from_error(self, env: Environment, ex: Exception) -> FailureEvent:
        if isinstance(ex, ClientError):
            error, cause = self._error_cause_from_client_error(ex)
            error_name = CustomErrorName(error)
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
        else:
            return FailureEvent(
                error_name=CustomErrorName(self._ERROR_NAME_AWS),
                event_type=HistoryEventType.TaskFailed,
                event_details=EventDetails(
                    taskFailedEventDetails=TaskFailedEventDetails(
                        error=self._ERROR_NAME_AWS,
                        cause=str(ex),  # TODO: update to report expected cause.
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
        service_name = self._get_boto_service_name()
        api_action = self._get_boto_service_action()
        dynamodb_client = boto_client_for(
            region=resource_runtime_part.region,
            account=resource_runtime_part.account,
            service=service_name,
        )
        response = getattr(dynamodb_client, api_action)(**normalised_parameters)
        response.pop("ResponseMetadata", None)
        env.stack.append(response)
