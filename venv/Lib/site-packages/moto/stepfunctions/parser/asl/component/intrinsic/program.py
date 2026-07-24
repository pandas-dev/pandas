from typing import Final

from moto.stepfunctions.parser.asl.component.intrinsic.component import Component


class Program(Component):
    def __init__(self):
        self.statements: Final[list[Component]] = []
