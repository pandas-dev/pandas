from typing import Any, Final

from moto.stepfunctions.parser.asl.component.common.payload.payloadvalue.payload_value import (
    PayloadValue,
)
from moto.stepfunctions.parser.asl.component.common.payload.payloadvalue.payloadbinding.payload_binding import (
    PayloadBinding,
)
from moto.stepfunctions.parser.asl.eval.environment import Environment


class PayloadBindingValue(PayloadBinding):
    def __init__(self, field: str, value: PayloadValue):
        super().__init__(field=field)
        self.value: Final[PayloadValue] = value

    def _eval_val(self, env: Environment) -> Any:
        self.value.eval(env)
        val: Any = env.stack.pop()
        return val
