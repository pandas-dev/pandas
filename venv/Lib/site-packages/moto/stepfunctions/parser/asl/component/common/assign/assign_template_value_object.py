from typing import Final

from moto.stepfunctions.parser.asl.component.common.assign.assign_template_binding import (
    AssignTemplateBinding,
)
from moto.stepfunctions.parser.asl.component.common.assign.assign_template_value import (
    AssignTemplateValue,
)
from moto.stepfunctions.parser.asl.eval.environment import Environment


class AssignTemplateValueObject(AssignTemplateValue):
    bindings: Final[list[AssignTemplateBinding]]

    def __init__(self, bindings: list[AssignTemplateBinding]):
        self.bindings = bindings

    def _eval_body(self, env: Environment) -> None:
        env.stack.append({})
        for binding in self.bindings:
            binding.eval(env)
