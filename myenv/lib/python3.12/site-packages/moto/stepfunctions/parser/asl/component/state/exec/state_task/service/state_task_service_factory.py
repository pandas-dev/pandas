from __future__ import annotations

from moto.stepfunctions.parser.asl.component.state.exec.state_task.service.state_task_service import (
    StateTaskService,
)
from moto.stepfunctions.parser.asl.component.state.exec.state_task.service.state_task_service_api_gateway import (
    StateTaskServiceApiGateway,
)
from moto.stepfunctions.parser.asl.component.state.exec.state_task.service.state_task_service_aws_sdk import (
    StateTaskServiceAwsSdk,
)
from moto.stepfunctions.parser.asl.component.state.exec.state_task.service.state_task_service_dynamodb import (
    StateTaskServiceDynamoDB,
)
from moto.stepfunctions.parser.asl.component.state.exec.state_task.service.state_task_service_events import (
    StateTaskServiceEvents,
)
from moto.stepfunctions.parser.asl.component.state.exec.state_task.service.state_task_service_lambda import (
    StateTaskServiceLambda,
)
from moto.stepfunctions.parser.asl.component.state.exec.state_task.service.state_task_service_sfn import (
    StateTaskServiceSfn,
)
from moto.stepfunctions.parser.asl.component.state.exec.state_task.service.state_task_service_sns import (
    StateTaskServiceSns,
)
from moto.stepfunctions.parser.asl.component.state.exec.state_task.service.state_task_service_sqs import (
    StateTaskServiceSqs,
)


# TODO: improve on factory constructor (don't use SubtypeManager: cannot reuse state task instances).
def state_task_service_for(service_name: str) -> StateTaskService:
    if service_name == "aws-sdk":
        return StateTaskServiceAwsSdk()
    if service_name == "lambda":
        return StateTaskServiceLambda()
    if service_name == "sqs":
        return StateTaskServiceSqs()
    if service_name == "states":
        return StateTaskServiceSfn()
    if service_name == "dynamodb":
        return StateTaskServiceDynamoDB()
    if service_name == "apigateway":
        return StateTaskServiceApiGateway()
    if service_name == "sns":
        return StateTaskServiceSns()
    if service_name == "events":
        return StateTaskServiceEvents()
    else:
        raise NotImplementedError(f"Unsupported service: '{service_name}'.")  # noqa
