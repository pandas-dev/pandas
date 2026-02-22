from typing import Final

from moto.stepfunctions.parser.asl.component.common.payload.payloadvalue.payload_value import (
    PayloadValue,
)
from moto.stepfunctions.parser.asl.component.common.payload.payloadvalue.payloadbinding.payload_binding import (
    PayloadBinding,
)
from moto.stepfunctions.parser.asl.eval.environment import Environment


class PayloadTmpl(PayloadValue):
    def __init__(self, payload_bindings: list[PayloadBinding]):
        self.payload_bindings: Final[list[PayloadBinding]] = payload_bindings

    def _eval_body(self, env: Environment) -> None:
        env.stack.append({})
        for payload_binding in self.payload_bindings:
            payload_binding.eval(env)
