from moto.moto_api._internal import mock_random
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


class StatesFunctionUUID(StatesFunction):
    def __init__(self, arg_list: FunctionArgumentList):
        super().__init__(
            states_name=StatesFunctionName(function_type=StatesFunctionNameType.UUID),
            arg_list=arg_list,
        )
        if len(arg_list.arg_list) != 0:
            raise ValueError(
                f"Expected no arguments for function type '{type(self)}', but got: '{arg_list}'."
            )

    def _eval_body(self, env: Environment) -> None:
        env.stack.append(str(mock_random.uuid4()))
