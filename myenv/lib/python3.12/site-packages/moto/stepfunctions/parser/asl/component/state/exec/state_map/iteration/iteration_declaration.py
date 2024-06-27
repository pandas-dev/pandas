from typing import Final, Optional

from moto.stepfunctions.parser.asl.component.common.comment import Comment
from moto.stepfunctions.parser.asl.component.common.flow.start_at import StartAt
from moto.stepfunctions.parser.asl.component.component import Component
from moto.stepfunctions.parser.asl.component.states import States


class IterationDecl(Component):
    comment: Final[Optional[Comment]]
    start_at: Final[StartAt]
    states: Final[States]

    def __init__(
        self,
        comment: Optional[Comment],
        start_at: StartAt,
        states: States,
    ):
        self.start_at = start_at
        self.comment = comment
        self.states = states
