import abc
from typing import Any

from moto.stepfunctions.parser.asl.component.common.payload.payloadvalue.payload_value import (
    PayloadValue,
)
from moto.stepfunctions.parser.asl.eval.environment import Environment


class PayloadValueLit(PayloadValue, abc.ABC):
    val: Any

    def __init__(self, val: Any):
        self.val: Any = val

    def _eval_body(self, env: Environment) -> None:
        env.stack.append(self.val)
