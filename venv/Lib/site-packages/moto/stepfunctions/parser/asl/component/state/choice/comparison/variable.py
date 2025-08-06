import abc
from typing import Final

from moto.stepfunctions.parser.asl.component.common.variable_sample import (
    VariableSample,
)
from moto.stepfunctions.parser.asl.component.eval_component import EvalComponent
from moto.stepfunctions.parser.asl.eval.environment import Environment
from moto.stepfunctions.parser.asl.utils.json_path import extract_json


class NoSuchVariable:
    def __init__(self, path: str):
        self.path: Final[str] = path


class Variable(EvalComponent, abc.ABC): ...


class VariableBase(Variable):
    def __init__(self, value: str):
        self.value: Final[str] = value

    def _eval_body(self, env: Environment) -> None:
        try:
            inp = env.stack[-1]
            value = extract_json(self.value, inp)
        except Exception as ex:
            value = NoSuchVariable(f"{self.value}, {ex}")
        env.stack.append(value)


class VariableContextObject(VariableBase):
    def __init__(self, value: str):
        value_tail = value[1:]
        super().__init__(value=value_tail)

    def _eval_body(self, env: Environment) -> None:
        try:
            value = extract_json(
                self.value, env.states.context_object.context_object_data
            )
        except Exception as ex:
            value = NoSuchVariable(f"{self.value}, {ex}")
        env.stack.append(value)


class VariableVar(Variable):
    variable_sample: VariableSample

    def __init__(self, variable_sample: VariableSample):
        self.variable_sample = variable_sample

    def _eval_body(self, env: Environment) -> None:
        self.variable_sample.eval(env=env)
