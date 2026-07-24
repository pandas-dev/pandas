import json

from moto.stepfunctions.parser.asl.component.eval_component import EvalComponent
from moto.stepfunctions.parser.asl.eval.environment import Environment


class Result(EvalComponent):
    result_obj: json

    def __init__(self, result_obj: json):
        self.result_obj = result_obj

    def _eval_body(self, env: Environment) -> None:
        env.stack.append(self.result_obj)
