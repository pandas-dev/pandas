import json
from typing import List, Optional

from moto.core.exceptions import JsonRESTError


class OpensearchIngestionExceptions(JsonRESTError):
    pass


class PipelineAlreadyExistsException(OpensearchIngestionExceptions):
    def __init__(self, pipeline_name: str):
        super().__init__(
            "ResourceAlreadyExistsException",
            f"Pipeline with given name {pipeline_name} already exists",
        )
        self.description = json.dumps({"message": self.message})


class PipelineInvalidStateException(OpensearchIngestionExceptions):
    def __init__(
        self, action: str, valid_states: List[str], current_state: Optional[str]
    ):
        super().__init__(
            "ConflictException",
            f"Only pipelines with one of the following statuses are eligible for {action}: {valid_states}. The current status is {current_state}.",
        )
        self.description = json.dumps({"message": self.message})


class PipelineNotFoundException(OpensearchIngestionExceptions):
    def __init__(self, pipeline_name: str):
        super().__init__(
            "ResourceNotFoundException", f"Pipeline {pipeline_name} could not be found."
        )
        self.description = json.dumps({"message": self.message})


class PipelineValidationException(OpensearchIngestionExceptions):
    def __init__(self, message: str):
        super().__init__("ValidationException", message)
        self.description = json.dumps({"message": self.message})


class InvalidVPCOptionsException(OpensearchIngestionExceptions):
    def __init__(self, message: str) -> None:
        super().__init__("ValidationException", f"Invalid VpcOptions: {message}")
        self.description = json.dumps({"message": self.message})


class SubnetNotFoundException(InvalidVPCOptionsException):
    def __init__(self, subnet_id: str) -> None:
        super().__init__(f"The subnet ID {subnet_id} does not exist")


class SecurityGroupNotFoundException(InvalidVPCOptionsException):
    def __init__(self, sg_id: str) -> None:
        super().__init__(f"The security group {sg_id} does not exist")
