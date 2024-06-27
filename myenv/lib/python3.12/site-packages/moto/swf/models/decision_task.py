from typing import TYPE_CHECKING, Any, Dict, Optional

from moto.core.common_models import BaseModel
from moto.core.utils import unix_time, utcnow
from moto.moto_api._internal import mock_random

from ..exceptions import SWFWorkflowExecutionClosedError
from .timeout import Timeout

if TYPE_CHECKING:
    from .workflow_execution import WorkflowExecution


class DecisionTask(BaseModel):
    def __init__(
        self, workflow_execution: "WorkflowExecution", scheduled_event_id: int
    ):
        self.workflow_execution = workflow_execution
        self.workflow_type = workflow_execution.workflow_type
        self.task_token = str(mock_random.uuid4())
        self.scheduled_event_id = scheduled_event_id
        self.previous_started_event_id: Optional[int] = None
        self.started_event_id: Optional[int] = None
        self.started_timestamp: Optional[float] = None
        self.start_to_close_timeout = (
            self.workflow_execution.task_start_to_close_timeout
        )
        self.state = "SCHEDULED"
        # this is *not* necessarily coherent with workflow execution history,
        # but that shouldn't be a problem for tests
        self.scheduled_at = utcnow()
        self.timeout_type: Optional[str] = None

    @property
    def started(self) -> bool:
        return self.state == "STARTED"

    def _check_workflow_execution_open(self) -> None:
        if not self.workflow_execution.open:
            raise SWFWorkflowExecutionClosedError()

    def to_full_dict(self, reverse_order: bool = False) -> Dict[str, Any]:
        events = self.workflow_execution.events(reverse_order=reverse_order)
        hsh: Dict[str, Any] = {
            "events": [evt.to_dict() for evt in events],
            "taskToken": self.task_token,
            "workflowExecution": self.workflow_execution.to_short_dict(),
            "workflowType": self.workflow_type.to_short_dict(),
        }
        if self.previous_started_event_id is not None:
            hsh["previousStartedEventId"] = self.previous_started_event_id
        if self.started_event_id:
            hsh["startedEventId"] = self.started_event_id
        return hsh

    def start(
        self, started_event_id: int, previous_started_event_id: Optional[int] = None
    ) -> None:
        self.state = "STARTED"
        self.started_timestamp = unix_time()
        self.started_event_id = started_event_id
        self.previous_started_event_id = previous_started_event_id

    def complete(self) -> None:
        self._check_workflow_execution_open()
        self.state = "COMPLETED"

    def first_timeout(self) -> Optional[Timeout]:
        if not self.started or not self.workflow_execution.open:
            return None
        # TODO: handle the "NONE" case
        start_to_close_at = self.started_timestamp + int(self.start_to_close_timeout)  # type: ignore
        _timeout = Timeout(self, start_to_close_at, "START_TO_CLOSE")
        if _timeout.reached:
            return _timeout
        return None

    def process_timeouts(self) -> None:
        _timeout = self.first_timeout()
        if _timeout:
            self.timeout(_timeout)

    def timeout(self, _timeout: Timeout) -> None:
        self._check_workflow_execution_open()
        self.state = "TIMED_OUT"
        self.timeout_type = _timeout.kind
