import abc

from moto.stepfunctions.parser.asl.component.component import Component


class FunctionName(Component, abc.ABC):
    name: str

    def __init__(self, name: str):
        self.name = name
