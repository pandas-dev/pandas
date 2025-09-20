from typing import Final, List

from moto.stepfunctions.parser.asl.component.common.error_name.error_name import (
    ErrorName,
)
from moto.stepfunctions.parser.asl.component.common.retry.retrier_decl import (
    RetrierDecl,
)
from moto.stepfunctions.parser.asl.component.common.retry.retrier_outcome import (
    RetrierOutcome,
)
from moto.stepfunctions.parser.asl.component.common.retry.retry_outcome import (
    RetryOutcome,
)
from moto.stepfunctions.parser.asl.component.eval_component import EvalComponent
from moto.stepfunctions.parser.asl.eval.environment import Environment


class RetryDecl(EvalComponent):
    def __init__(self, retriers: List[RetrierDecl]):
        self.retriers: Final[List[RetrierDecl]] = retriers

    def _eval_body(self, env: Environment) -> None:
        error_name: ErrorName = env.stack.pop()

        for retrier in self.retriers:
            env.stack.append(error_name)
            retrier.eval(env)
            outcome: RetrierOutcome = env.stack.pop()

            if outcome == RetrierOutcome.Skipped:
                continue
            elif outcome == RetrierOutcome.Executed:
                env.stack.append(RetryOutcome.CanRetry)
                return
            elif outcome == RetrierOutcome.Failed:
                env.stack.append(RetryOutcome.CannotRetry)
                return

        env.stack.append(RetryOutcome.NoRetrier)
