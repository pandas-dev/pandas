from typing import List

from moto.stepfunctions.parser.asl.component.common.jsonata.jsonata_template_value import (
    JSONataTemplateValue,
)
from moto.stepfunctions.parser.asl.eval.environment import Environment


class JSONataTemplateValueArray(JSONataTemplateValue):
    values: List[JSONataTemplateValue]

    def __init__(self, values: List[JSONataTemplateValue]):
        self.values = values

    def _eval_body(self, env: Environment) -> None:
        arr = list()
        for value in self.values:
            value.eval(env)
            arr.append(env.stack.pop())
        env.stack.append(arr)
