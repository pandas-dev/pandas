from typing import Final

from moto.stepfunctions.parser.asl.component.component import Component


class ErrorDecl(Component):
    def __init__(self, error: str):
        self.error: Final[str] = error
