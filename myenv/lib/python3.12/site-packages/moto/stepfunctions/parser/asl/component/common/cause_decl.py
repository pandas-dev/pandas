from typing import Final

from moto.stepfunctions.parser.asl.component.component import Component


class CauseDecl(Component):
    def __init__(self, cause: str):
        self.cause: Final[str] = cause
