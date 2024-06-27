from typing import Any, Final

from moto.stepfunctions.parser.asl.component.common.payload.payloadvalue.payloadbinding.payload_binding import (
    PayloadBinding,
)
from moto.stepfunctions.parser.asl.component.intrinsic.function.function import Function
from moto.stepfunctions.parser.asl.eval.environment import Environment
from moto.stepfunctions.parser.asl.parse.intrinsic.intrinsic_parser import (
    IntrinsicParser,
)


class PayloadBindingIntrinsicFunc(PayloadBinding):
    def __init__(self, field: str, intrinsic_func: str):
        super().__init__(field=field)
        self.src: Final[str] = intrinsic_func
        self.function: Final[Function] = IntrinsicParser.parse(self.src)

    @classmethod
    def from_raw(cls, string_dollar: str, intrinsic_func: str):
        field: str = string_dollar[:-2]
        return cls(field=field, intrinsic_func=intrinsic_func)

    def _eval_val(self, env: Environment) -> Any:
        self.function.eval(env=env)
        val = env.stack.pop()
        return val
