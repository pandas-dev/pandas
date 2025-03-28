from typing import List

from moto.stepfunctions.parser.asl.component.common.assign.assign_template_value import (
    AssignTemplateValue,
)
from moto.stepfunctions.parser.asl.eval.environment import Environment


class AssignTemplateValueArray(AssignTemplateValue):
    values: List[AssignTemplateValue]

    def __init__(self, values: List[AssignTemplateValue]):
        self.values = values

    def _eval_body(self, env: Environment) -> None:
        arr = list()
        for value in self.values:
            value.eval(env)
            arr.append(env.stack.pop())
        env.stack.append(arr)
