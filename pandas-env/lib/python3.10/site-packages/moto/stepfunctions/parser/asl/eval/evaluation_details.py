from typing import Any, Optional

from moto.stepfunctions.parser.api import Arn, Definition, LongArn, StateMachineType


class AWSExecutionDetails:
    account: str
    region: str
    role_arn: str

    def __init__(self, account: str, region: str, role_arn: str):
        self.account = account
        self.region = region
        self.role_arn = role_arn


class ExecutionDetails:
    arn: LongArn
    name: str
    role_arn: Arn
    inpt: Optional[Any]
    start_time: str

    def __init__(
        self,
        arn: LongArn,
        name: str,
        role_arn: Arn,
        inpt: Optional[Any],
        start_time: str,
    ):
        self.arn = arn
        self.name = name
        self.role_arn = role_arn
        self.inpt = inpt
        self.start_time = start_time


class StateMachineDetails:
    arn: Arn
    name: str
    typ: StateMachineType
    definition: Definition

    def __init__(self, arn: Arn, name: str, typ: StateMachineType, definition: str):
        self.arn = arn
        self.name = name
        self.typ = typ
        self.definition = definition


class EvaluationDetails:
    aws_execution_details: AWSExecutionDetails
    execution_details: ExecutionDetails
    state_machine_details: StateMachineDetails

    def __init__(
        self,
        aws_execution_details: AWSExecutionDetails,
        execution_details: ExecutionDetails,
        state_machine_details: StateMachineDetails,
    ):
        self.aws_execution_details = aws_execution_details
        self.execution_details = execution_details
        self.state_machine_details = state_machine_details
