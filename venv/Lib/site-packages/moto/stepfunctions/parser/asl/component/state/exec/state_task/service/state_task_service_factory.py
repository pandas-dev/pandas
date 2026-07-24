from __future__ import annotations

from typing import Final

from antlr4 import RecognitionException

from moto.stepfunctions.parser.asl.component.state.exec.state_task.service.state_task_service import (
    StateTaskService,
)
from moto.stepfunctions.parser.asl.component.state.exec.state_task.service.state_task_service_api_gateway import (
    StateTaskServiceApiGateway,
)
from moto.stepfunctions.parser.asl.component.state.exec.state_task.service.state_task_service_aws_sdk import (
    StateTaskServiceAwsSdk,
)
from moto.stepfunctions.parser.asl.component.state.exec.state_task.service.state_task_service_batch import (
    StateTaskServiceBatch,
)
from moto.stepfunctions.parser.asl.component.state.exec.state_task.service.state_task_service_dynamodb import (
    StateTaskServiceDynamoDB,
)
from moto.stepfunctions.parser.asl.component.state.exec.state_task.service.state_task_service_ecs import (
    StateTaskServiceEcs,
)
from moto.stepfunctions.parser.asl.component.state.exec.state_task.service.state_task_service_events import (
    StateTaskServiceEvents,
)
from moto.stepfunctions.parser.asl.component.state.exec.state_task.service.state_task_service_glue import (
    StateTaskServiceGlue,
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
from moto.stepfunctions.parser.asl.component.state.exec.state_task.service.state_task_service_unsupported import (
    StateTaskServiceUnsupported,
)

_UNSUPPORTED_SERVICE_NAMES: Final[set[str]] = {
    "athena",
    "bedrock",
    "codebuild",
    "eks",
    "elasticmapreduce",
    "emr-containers",
    "emr-serverless",
    "databrew",
    "mediaconvert",
    "sagemaker",
}


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
    if service_name == "ecs":
        return StateTaskServiceEcs()
    if service_name == "glue":
        return StateTaskServiceGlue()
    if service_name == "batch":
        return StateTaskServiceBatch()
    if service_name in _UNSUPPORTED_SERVICE_NAMES:
        return StateTaskServiceUnsupported()
    raise RecognitionException(f"Unknown service '{service_name}'")
