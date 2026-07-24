import abc
from typing import Final

from moto.stepfunctions.parser.asl.component.eval_component import EvalComponent
from moto.stepfunctions.parser.asl.component.intrinsic.argument.argument import (
    ArgumentList,
)
from moto.stepfunctions.parser.asl.component.intrinsic.functionname.function_name import (
    FunctionName,
)


class Function(EvalComponent, abc.ABC):
    name: FunctionName
    argument_list: Final[ArgumentList]

    def __init__(self, name: FunctionName, argument_list: ArgumentList):
        self.name = name
        self.argument_list = argument_list
