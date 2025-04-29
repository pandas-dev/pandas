from typing import TYPE_CHECKING, Any, Dict, Optional

from moto.core.common_models import BaseModel
from moto.core.utils import unix_time, utcnow
from moto.moto_api._internal import mock_random

from ..exceptions import SWFWorkflowExecutionClosedError
from .timeout import Timeout

if TYPE_CHECKING:
    from .activity_type import ActivityType
    from .workflow_execution import WorkflowExecution


class ActivityTask(BaseModel):
    def __init__(
        self,
        activity_id: str,
        activity_type: "ActivityType",
        scheduled_event_id: int,
        workflow_execution: "WorkflowExecution",
        timeouts: Dict[str, Any],
        workflow_input: Any = None,
    ):
        self.activity_id = activity_id
        self.activity_type = activity_type
        self.details = None
        self.input = workflow_input
        self.last_heartbeat_timestamp = unix_time()
        self.scheduled_event_id = scheduled_event_id
        self.started_event_id: Optional[int] = None
        self.state = "SCHEDULED"
        self.task_token = str(mock_random.uuid4())
        self.timeouts = timeouts
        self.timeout_type: Optional[str] = None
        self.workflow_execution = workflow_execution
        # this is *not* necessarily coherent with workflow execution history,
        # but that shouldn't be a problem for tests
        self.scheduled_at = utcnow()

    def _check_workflow_execution_open(self) -> None:
        if not self.workflow_execution.open:
            raise SWFWorkflowExecutionClosedError()

    @property
    def open(self) -> bool:
        return self.state in ["SCHEDULED", "STARTED"]

    def to_full_dict(self) -> Dict[str, Any]:
        hsh: Dict[str, Any] = {
            "activityId": self.activity_id,
            "activityType": self.activity_type.to_short_dict(),
            "taskToken": self.task_token,
            "startedEventId": self.started_event_id,
            "workflowExecution": self.workflow_execution.to_short_dict(),
        }
        if self.input:
            hsh["input"] = self.input
        return hsh

    def start(self, started_event_id: int) -> None:
        self.state = "STARTED"
        self.started_event_id = started_event_id

    def complete(self) -> None:
        self._check_workflow_execution_open()
        self.state = "COMPLETED"

    def fail(self) -> None:
        self._check_workflow_execution_open()
        self.state = "FAILED"

    def reset_heartbeat_clock(self) -> None:
        self.last_heartbeat_timestamp = unix_time()

    def first_timeout(self) -> Optional[Timeout]:
        if not self.open or not self.workflow_execution.open:
            return None

        if self.timeouts["heartbeatTimeout"] == "NONE":
            return None

        heartbeat_timeout_at = self.last_heartbeat_timestamp + int(
            self.timeouts["heartbeatTimeout"]
        )
        _timeout = Timeout(self, heartbeat_timeout_at, "HEARTBEAT")
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
