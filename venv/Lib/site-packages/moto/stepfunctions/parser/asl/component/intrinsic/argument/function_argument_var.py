from moto.stepfunctions.parser.asl.component.common.variable_sample import (
    VariableSample,
)
from moto.stepfunctions.parser.asl.component.intrinsic.argument.function_argument import (
    FunctionArgument,
)
from moto.stepfunctions.parser.asl.eval.environment import Environment


class FunctionArgumentVar(FunctionArgument):
    variable_sample: VariableSample

    def __init__(self, variable_sample: VariableSample):
        super().__init__()
        self.variable_sample = variable_sample

    def _eval_body(self, env: Environment) -> None:
        self.variable_sample.eval(env=env)
        self._value = env.stack.pop()
        super()._eval_body(env=env)
