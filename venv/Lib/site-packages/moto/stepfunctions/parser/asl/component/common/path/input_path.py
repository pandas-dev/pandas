import abc
import copy
from typing import Optional

from moto.stepfunctions.parser.asl.component.common.variable_sample import (
    VariableSample,
)
from moto.stepfunctions.parser.asl.component.eval_component import EvalComponent
from moto.stepfunctions.parser.asl.eval.environment import Environment
from moto.stepfunctions.parser.asl.utils.json_path import extract_json


class InputPath(EvalComponent, abc.ABC): ...


class InputPathBase(InputPath):
    DEFAULT_PATH: str = "$"

    path: Optional[str]

    def __init__(self, path: Optional[str]):
        self.path = path

    def _eval_body(self, env: Environment) -> None:
        if self.path is None:
            value = dict()
        elif self.path == self.DEFAULT_PATH:
            value = env.states.get_input()
        else:
            value = extract_json(self.path, env.states.get_input())
        env.stack.append(copy.deepcopy(value))


class InputPathContextObject(InputPathBase):
    def __init__(self, path: str):
        path_tail = path[1:]
        super().__init__(path=path_tail)

    def _eval_body(self, env: Environment) -> None:
        value = extract_json(self.path, env.states.context_object.context_object_data)
        env.stack.append(copy.deepcopy(value))


class InputPathVar(InputPath):
    variable_sample: VariableSample

    def __init__(self, variable_sample: VariableSample):
        self.variable_sample = variable_sample

    def _eval_body(self, env: Environment) -> None:
        self.variable_sample.eval(env=env)
