import logging
import threading
from typing import Final

from moto.stepfunctions.parser.api import (
    ExecutionFailedEventDetails,
)
from moto.stepfunctions.parser.asl.component.common.error_name.failure_event import (
    FailureEventException,
)
from moto.stepfunctions.parser.asl.component.common.error_name.states_error_name import (
    StatesErrorName,
)
from moto.stepfunctions.parser.asl.component.common.error_name.states_error_name_type import (
    StatesErrorNameType,
)
from moto.stepfunctions.parser.asl.component.eval_component import EvalComponent
from moto.stepfunctions.parser.asl.component.state.state import CommonStateField
from moto.stepfunctions.parser.asl.eval.test_state.environment import (
    TestStateEnvironment,
)
from moto.stepfunctions.parser.asl.utils.encoding import to_json_str
from moto.stepfunctions.parser.utils import TMP_THREADS

LOG = logging.getLogger(__name__)

TEST_CASE_EXECUTION_TIMEOUT_SECONDS: Final[int] = 300  # 5 minutes.


class TestStateProgram(EvalComponent):
    test_state: Final[CommonStateField]

    def __init__(
        self,
        test_state: CommonStateField,
    ):
        self.test_state = test_state

    def eval(self, env: TestStateEnvironment) -> None:
        env.next_state_name = self.test_state.name
        worker_thread = threading.Thread(target=super().eval, args=(env,), daemon=True)
        TMP_THREADS.append(worker_thread)
        worker_thread.start()
        worker_thread.join(timeout=TEST_CASE_EXECUTION_TIMEOUT_SECONDS)
        is_timeout = worker_thread.is_alive()
        if is_timeout:
            env.set_timed_out()

    def _eval_body(self, env: TestStateEnvironment) -> None:
        try:
            env.inspection_data["input"] = to_json_str(env.states.get_input())
            self.test_state.eval(env=env)
        except FailureEventException as ex:
            env.set_error(error=ex.get_execution_failed_event_details())
        except Exception as ex:
            cause = f"{type(ex).__name__}({str(ex)})"
            LOG.error("Stepfunctions computation ended with exception '%s'.", cause)
            env.set_error(
                ExecutionFailedEventDetails(
                    error=StatesErrorName(
                        typ=StatesErrorNameType.StatesRuntime
                    ).error_name,
                    cause=cause,
                )
            )
