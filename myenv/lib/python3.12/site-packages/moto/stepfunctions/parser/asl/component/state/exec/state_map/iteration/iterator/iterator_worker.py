import copy
import logging
from typing import Final, Optional

from moto.stepfunctions.parser.asl.component.common.parameters import Parameters
from moto.stepfunctions.parser.asl.component.state.exec.state_map.iteration.iteration_worker import (
    IterationWorker,
)
from moto.stepfunctions.parser.asl.component.state.exec.state_map.iteration.job import (
    JobPool,
)
from moto.stepfunctions.parser.asl.eval.environment import Environment

LOG = logging.getLogger(__name__)


class IteratorWorker(IterationWorker):
    _parameters: Final[Optional[Parameters]]

    def __init__(
        self,
        work_name: str,
        job_pool: JobPool,
        env: Environment,
        parameters: Optional[Parameters],
    ):
        super().__init__(work_name=work_name, job_pool=job_pool, env=env)
        self._parameters = parameters

    def _eval_input(self, env_frame: Environment) -> None:
        if self._parameters:
            map_state_input = self._env.stack[-1]
            env_frame.inp = copy.deepcopy(map_state_input)
            env_frame.stack.append(env_frame.inp)
            self._parameters.eval(env_frame)
            env_frame.inp = env_frame.stack.pop()
            env_frame.stack.append(env_frame.inp)
