from .generic_type import GenericType


class ActivityType(GenericType):
    @property
    def _configuration_keys(self) -> list[str]:
        return [
            "defaultTaskHeartbeatTimeout",
            "defaultTaskScheduleToCloseTimeout",
            "defaultTaskScheduleToStartTimeout",
            "defaultTaskStartToCloseTimeout",
        ]

    @property
    def kind(self) -> str:
        return "activity"
