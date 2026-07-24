from .generic_type import GenericType


class WorkflowType(GenericType):
    @property
    def _configuration_keys(self) -> list[str]:
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
