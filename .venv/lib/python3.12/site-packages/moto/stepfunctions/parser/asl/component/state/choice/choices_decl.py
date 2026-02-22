from typing import Final

from moto.stepfunctions.parser.asl.component.component import Component
from moto.stepfunctions.parser.asl.component.state.choice.choice_rule import (
    ChoiceRule,
)


class ChoicesDecl(Component):
    def __init__(self, rules: list[ChoiceRule]):
        self.rules: Final[list[ChoiceRule]] = rules
