from logging import Logger
from typing import Final

log: Logger
SERVICE_NAME: Final = "ecs"
ORIGIN: Final = "AWS::ECS::Container"

def initialize() -> None: ...
