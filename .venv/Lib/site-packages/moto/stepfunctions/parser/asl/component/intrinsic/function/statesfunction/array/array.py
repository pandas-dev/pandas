from typing import Any, List

from moto.stepfunctions.parser.asl.component.intrinsic.argument.function_argument_list import (
    FunctionArgumentList,
)
from moto.stepfunctions.parser.asl.component.intrinsic.function.statesfunction.states_function import (
    StatesFunction,
)
from moto.stepfunctions.parser.asl.component.intrinsic.functionname.state_function_name_types import (
    StatesFunctionNameType,
)
from moto.stepfunctions.parser.asl.component.intrinsic.functionname.states_function_name import (
    StatesFunctionName,
)
from moto.stepfunctions.parser.asl.eval.environment import Environment


class Array(StatesFunction):
    def __init__(self, arg_list: FunctionArgumentList):
        super().__init__(
            states_name=StatesFunctionName(function_type=StatesFunctionNameType.Array),
            arg_list=arg_list,
        )

    def _eval_body(self, env: Environment) -> None:
        self.arg_list.eval(env=env)
        values: List[Any] = env.stack.pop()
        env.stack.append(values)
