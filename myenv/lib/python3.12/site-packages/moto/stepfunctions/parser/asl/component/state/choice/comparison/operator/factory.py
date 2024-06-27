from moto.stepfunctions.parser.asl.component.state.choice.comparison.comparison_operator_type import (
    ComparisonOperatorType,
)
from moto.stepfunctions.parser.asl.component.state.choice.comparison.operator.implementations.boolean_equals import *  # noqa
from moto.stepfunctions.parser.asl.component.state.choice.comparison.operator.implementations.is_operator import *  # noqa
from moto.stepfunctions.parser.asl.component.state.choice.comparison.operator.implementations.numeric import *  # noqa
from moto.stepfunctions.parser.asl.component.state.choice.comparison.operator.implementations.string_operators import *  # noqa
from moto.stepfunctions.parser.asl.component.state.choice.comparison.operator.implementations.timestamp_operators import *  # noqa
from moto.stepfunctions.parser.asl.component.state.choice.comparison.operator.operator import (
    Operator,
)


class OperatorFactory:
    @staticmethod
    def get(typ: ComparisonOperatorType) -> Operator:
        op = Operator.get((str(typ)))
        if op is None:
            raise NotImplementedError(f"{typ} is not supported.")
        return op
