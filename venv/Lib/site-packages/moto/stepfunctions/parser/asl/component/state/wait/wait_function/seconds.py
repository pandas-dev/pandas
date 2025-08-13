from typing import Final

from moto.stepfunctions.parser.asl.component.common.jsonata.jsonata_template_value_terminal import (
    JSONataTemplateValueTerminalExpression,
)
from moto.stepfunctions.parser.asl.component.state.wait.wait_function.wait_function import (
    WaitFunction,
)
from moto.stepfunctions.parser.asl.eval.environment import Environment


class Seconds(WaitFunction):
    # Seconds
    # A time, in seconds, to state_wait before beginning the state specified in the Next
    # field. You must specify time as a positive, integer value.

    def __init__(self, seconds: int):
        self.seconds: int = seconds

    def _get_wait_seconds(self, env: Environment) -> int:
        return self.seconds


class SecondsJSONata(WaitFunction):
    jsonata_template_value_terminal_expression: Final[
        JSONataTemplateValueTerminalExpression
    ]

    def __init__(
        self,
        jsonata_template_value_terminal_expression: JSONataTemplateValueTerminalExpression,
    ):
        super().__init__()
        self.jsonata_template_value_terminal_expression = (
            jsonata_template_value_terminal_expression
        )

    def _get_wait_seconds(self, env: Environment) -> int:
        # TODO: add snapshot tests to verify AWS's behaviour about non integer values.
        self.jsonata_template_value_terminal_expression.eval(env=env)
        max_items: int = int(env.stack.pop())
        return max_items
