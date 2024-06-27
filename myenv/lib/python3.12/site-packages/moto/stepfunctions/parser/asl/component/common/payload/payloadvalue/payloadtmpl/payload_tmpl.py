from typing import Final, List

from moto.stepfunctions.parser.asl.component.common.payload.payloadvalue.payload_value import (
    PayloadValue,
)
from moto.stepfunctions.parser.asl.component.common.payload.payloadvalue.payloadbinding.payload_binding import (
    PayloadBinding,
)
from moto.stepfunctions.parser.asl.eval.environment import Environment


class PayloadTmpl(PayloadValue):
    def __init__(self, payload_bindings: List[PayloadBinding]):
        self.payload_bindings: Final[List[PayloadBinding]] = payload_bindings

    def _eval_body(self, env: Environment) -> None:
        env.stack.append(dict())
        for payload_binding in self.payload_bindings:
            payload_binding.eval(env)
