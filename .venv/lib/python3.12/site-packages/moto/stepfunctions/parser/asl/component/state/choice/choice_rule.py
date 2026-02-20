from typing import Final, Optional

from moto.stepfunctions.parser.asl.component.common.assign.assign_decl import AssignDecl
from moto.stepfunctions.parser.asl.component.common.comment import Comment
from moto.stepfunctions.parser.asl.component.common.flow.next import Next
from moto.stepfunctions.parser.asl.component.common.outputdecl import Output
from moto.stepfunctions.parser.asl.component.eval_component import EvalComponent
from moto.stepfunctions.parser.asl.component.state.choice.comparison.comparison_type import (
    Comparison,
)
from moto.stepfunctions.parser.asl.eval.environment import Environment


class ChoiceRule(EvalComponent):
    comparison: Final[Optional[Comparison]]
    next_stmt: Final[Optional[Next]]
    comment: Final[Optional[Comment]]
    assign: Final[Optional[AssignDecl]]
    output: Final[Optional[Output]]

    def __init__(
        self,
        comparison: Optional[Comparison],
        next_stmt: Optional[Next],
        comment: Optional[Comment],
        assign: Optional[AssignDecl],
        output: Optional[Output],
    ):
        self.comparison = comparison
        self.next_stmt = next_stmt
        self.comment = comment
        self.assign = assign
        self.output = output

    def _eval_body(self, env: Environment) -> None:
        self.comparison.eval(env)
        is_condition_true: bool = env.stack[-1]
        if not is_condition_true:
            return
        if self.assign:
            self.assign.eval(env=env)
        if self.output:
            self.output.eval(env=env)
