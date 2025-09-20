from typing import List

from .generic_type import GenericType


class WorkflowType(GenericType):
    @property
    def _configuration_keys(self) -> List[str]:
        return [
            "defaultChildPolicy",
            "defaultExecutionStartToCloseTimeout",
            "defaultTaskStartToCloseTimeout",
            "defaultTaskPriority",
            "defaultLambdaRole",
        ]

    @property
    def kind(self) -> str:
        return "workflow"
