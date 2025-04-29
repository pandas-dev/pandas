from typing import Final, List

from moto.stepfunctions.parser.asl.component.eval_component import EvalComponent
from moto.stepfunctions.parser.asl.component.intrinsic.argument.function_argument import (
    FunctionArgument,
)
from moto.stepfunctions.parser.asl.eval.environment import Environment


class FunctionArgumentList(EvalComponent):
    def __init__(self, arg_list: List[FunctionArgument]):
        self.arg_list: Final[List[FunctionArgument]] = arg_list
        self.size: Final[int] = len(arg_list)

    def _eval_body(self, env: Environment) -> None:
        values = list()
        for arg in self.arg_list:
            arg.eval(env=env)
            values.append(env.stack.pop())
        env.stack.append(values)
