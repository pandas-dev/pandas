from moto.stepfunctions.parser.asl.component.component import Component
from moto.stepfunctions.parser.asl.component.state.state import CommonStateField


class States(Component):
    def __init__(self):
        self.states: dict[str, CommonStateField] = dict()
