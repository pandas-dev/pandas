from typing import Final

from moto.stepfunctions.parser.asl.component.common.assign.assign_template_binding import (
    AssignTemplateBinding,
)
from moto.stepfunctions.parser.asl.component.eval_component import EvalComponent
from moto.stepfunctions.parser.asl.eval.environment import Environment


class AssignDeclBinding(EvalComponent):
    binding: Final[AssignTemplateBinding]

    def __init__(self, binding: AssignTemplateBinding):
        super().__init__()
        self.binding = binding

    def _eval_body(self, env: Environment) -> None:
        env.stack.append(dict())
        self.binding.eval(env=env)
