from typing import Final

from moto.stepfunctions.parser.asl.component.state.choice.comparison.comparison import (
    Comparison,
)
from moto.stepfunctions.parser.asl.component.state.choice.comparison.comparison_func import (
    ComparisonFunc,
)
from moto.stepfunctions.parser.asl.component.state.choice.comparison.variable import (
    Variable,
)
from moto.stepfunctions.parser.asl.eval.environment import Environment


class ComparisonVariable(Comparison):
    variable: Final[Variable]
    comparison_function: Final[ComparisonFunc]

    def __init__(self, variable: Variable, func: ComparisonFunc):
        self.variable = variable
        self.comparison_function = func

    def _eval_body(self, env: Environment) -> None:
        variable: Variable = self.variable
        variable.eval(env)
        comparison_function: ComparisonFunc = self.comparison_function
        comparison_function.eval(env)
