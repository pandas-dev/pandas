from moto.stepfunctions.parser.asl.component.common.payload.payloadvalue.payload_value import (
    PayloadValue,
)
from moto.stepfunctions.parser.asl.component.common.variable_sample import (
    VariableSample,
)
from moto.stepfunctions.parser.asl.eval.environment import Environment


class PayloadValueVariableSample(PayloadValue):
    variable_sample: VariableSample

    def __init__(self, variable_sample: VariableSample):
        super().__init__()
        self.variable_sample = variable_sample

    def _eval_body(self, env: Environment) -> None:
        self.variable_sample.eval(env=env)
