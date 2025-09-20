from typing import Final, List

from moto.stepfunctions.parser.asl.component.common.assign.assign_template_binding import (
    AssignTemplateBinding,
)
from moto.stepfunctions.parser.asl.component.common.assign.assign_template_value import (
    AssignTemplateValue,
)
from moto.stepfunctions.parser.asl.eval.environment import Environment


class AssignTemplateValueObject(AssignTemplateValue):
    bindings: Final[List[AssignTemplateBinding]]

    def __init__(self, bindings: List[AssignTemplateBinding]):
        self.bindings = bindings

    def _eval_body(self, env: Environment) -> None:
        env.stack.append(dict())
        for binding in self.bindings:
            binding.eval(env)
