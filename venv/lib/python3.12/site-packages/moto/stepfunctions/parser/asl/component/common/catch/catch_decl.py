from typing import Final, List

from moto.stepfunctions.parser.asl.component.common.catch.catch_outcome import (
    CatchOutcome,
)
from moto.stepfunctions.parser.asl.component.common.catch.catcher_decl import (
    CatcherDecl,
)
from moto.stepfunctions.parser.asl.component.common.catch.catcher_outcome import (
    CatcherOutcome,
    CatcherOutcomeCaught,
)
from moto.stepfunctions.parser.asl.component.eval_component import EvalComponent
from moto.stepfunctions.parser.asl.eval.environment import Environment


class CatchDecl(EvalComponent):
    def __init__(self, catchers: List[CatcherDecl]):
        self.catchers: Final[List[CatcherDecl]] = catchers

    def _eval_body(self, env: Environment) -> None:
        for catcher in self.catchers:
            catcher.eval(env)
            catcher_outcome: CatcherOutcome = env.stack.pop()

            if isinstance(catcher_outcome, CatcherOutcomeCaught):
                env.stack.append(CatchOutcome.Caught)
                return

        env.stack.append(CatchOutcome.NotCaught)
