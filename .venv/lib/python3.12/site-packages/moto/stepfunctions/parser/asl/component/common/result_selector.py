from typing import Final

from moto.stepfunctions.parser.asl.component.common.payload.payloadvalue.payloadtmpl.payload_tmpl import (
    PayloadTmpl,
)
from moto.stepfunctions.parser.asl.component.eval_component import EvalComponent
from moto.stepfunctions.parser.asl.eval.environment import Environment


class ResultSelector(EvalComponent):
    payload_tmpl: Final[PayloadTmpl]

    def __init__(self, payload_tmpl: PayloadTmpl):
        self.payload_tmpl = payload_tmpl

    def _eval_body(self, env: Environment) -> None:
        self.payload_tmpl.eval(env=env)
