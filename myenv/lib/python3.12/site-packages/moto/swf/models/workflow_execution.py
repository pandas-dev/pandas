from threading import Lock
from threading import Timer as ThreadingTimer
from typing import Any, Dict, Iterable, List, Optional

from moto.core.common_models import BaseModel
from moto.core.utils import camelcase_to_underscores, unix_time
from moto.moto_api._internal import mock_random

from ..constants import DECISIONS_FIELDS
from ..exceptions import (
    SWFDecisionValidationException,
    SWFDefaultUndefinedFault,
    SWFValidationException,
)
from ..utils import decapitalize
from .activity_task import ActivityTask
from .activity_type import ActivityType
from .decision_task import DecisionTask
from .domain import Domain
from .history_event import HistoryEvent
from .timeout import Timeout
from .timer import Timer
from .workflow_type import WorkflowType


# TODO: extract decision related logic into a Decision class
class WorkflowExecution(BaseModel):
    # NB: the list is ordered exactly as in SWF validation exceptions so we can
    # mimic error messages closely ; don't reorder it without checking SWF.
    KNOWN_DECISION_TYPES = [
        "CompleteWorkflowExecution",
        "StartTimer",
        "RequestCancelExternalWorkflowExecution",
        "SignalExternalWorkflowExecution",
        "CancelTimer",
        "RecordMarker",
        "ScheduleActivityTask",
        "ContinueAsNewWorkflowExecution",
        "ScheduleLambdaFunction",
        "FailWorkflowExecution",
        "RequestCancelActivityTask",
        "StartChildWorkflowExecution",
        "CancelWorkflowExecution",
    ]

    def __init__(
        self,
        domain: Domain,
        workflow_type: "WorkflowType",
        workflow_id: str,
        **kwargs: Any,
    ):
        self.domain = domain
        self.workflow_id = workflow_id
        self.run_id = mock_random.uuid4().hex
        # WorkflowExecutionInfo
        self.cancel_requested = False
        # TODO: check valid values among:
        # COMPLETED | FAILED | CANCELED | TERMINATED | CONTINUED_AS_NEW | TIMED_OUT
        # TODO: implement them all
        self.close_cause: Optional[str] = None
        self.close_status: Optional[str] = None
        self.close_timestamp: Optional[float] = None
        self.execution_status = "OPEN"
        self.latest_activity_task_timestamp: Optional[float] = None
        self.latest_execution_context: Optional[str] = None
        self.parent = None
        self.start_timestamp: Optional[float] = None
        self.tag_list = kwargs.get("tag_list", None) or []
        self.timeout_type: Optional[str] = None
        self.workflow_type = workflow_type
        # args processing
        # NB: the order follows boto/SWF order of exceptions appearance (if no
        # param is set, # SWF will raise DefaultUndefinedFault errors in the
        # same order as the few lines that follow)
        self.execution_start_to_close_timeout = self._get_from_kwargs_or_workflow_type(
            kwargs, "execution_start_to_close_timeout"
        )
        self.task_list = self._get_from_kwargs_or_workflow_type(
            kwargs, "task_list", "task_list"
        )
        self.task_start_to_close_timeout = self._get_from_kwargs_or_workflow_type(
            kwargs, "task_start_to_close_timeout"
        )
        self.child_policy = self._get_from_kwargs_or_workflow_type(
            kwargs, "child_policy"
        )
        self.input = kwargs.get("workflow_input")
        # counters
        self.open_counts = {
            "openTimers": 0,
            "openDecisionTasks": 0,
            "openActivityTasks": 0,
            "openChildWorkflowExecutions": 0,
            "openLambdaFunctions": 0,
        }
        # events
        self._events: List[HistoryEvent] = []
        # child workflows
        self.child_workflow_executions: List[WorkflowExecution] = []
        self._previous_started_event_id: Optional[int] = None
        # timers/thread utils
        self.threading_lock = Lock()
        self._timers: Dict[str, Timer] = {}

    def __repr__(self) -> str:
        return f"WorkflowExecution(run_id: {self.run_id})"

    def _get_from_kwargs_or_workflow_type(
        self,
        kwargs: Dict[str, Any],
        local_key: str,
        workflow_type_key: Optional[str] = None,
    ) -> Any:
        if workflow_type_key is None:
            workflow_type_key = "default_" + local_key
        value = kwargs.get(local_key)
        if not value and hasattr(self.workflow_type, workflow_type_key):
            value = getattr(self.workflow_type, workflow_type_key)
        if not value:
            raise SWFDefaultUndefinedFault(local_key)
        return value

    @property
    def _configuration_keys(self) -> List[str]:
        return [
            "executionStartToCloseTimeout",
            "childPolicy",
            "taskPriority",
            "taskStartToCloseTimeout",
        ]

    def to_short_dict(self) -> Dict[str, str]:
        return {"workflowId": self.workflow_id, "runId": self.run_id}

    def to_medium_dict(self) -> Dict[str, Any]:
        hsh: Dict[str, Any] = {
            "execution": self.to_short_dict(),
            "workflowType": self.workflow_type.to_short_dict(),
            "startTimestamp": 1420066800.123,
            "executionStatus": self.execution_status,
            "cancelRequested": self.cancel_requested,
        }
        if hasattr(self, "tag_list") and self.tag_list:
            hsh["tagList"] = self.tag_list
        return hsh

    def to_full_dict(self) -> Dict[str, Any]:
        hsh: Dict[str, Any] = {
            "executionInfo": self.to_medium_dict(),
            "executionConfiguration": {"taskList": {"name": self.task_list}},
        }
        # info
        if self.execution_status == "CLOSED":
            hsh["executionInfo"]["closeStatus"] = self.close_status
            hsh["executionInfo"]["closeTimestamp"] = self.close_timestamp
        # configuration
        for key in self._configuration_keys:
            attr = camelcase_to_underscores(key)
            if not hasattr(self, attr):
                continue
            if not getattr(self, attr):
                continue
            hsh["executionConfiguration"][key] = getattr(self, attr)
        # counters
        hsh["openCounts"] = self.open_counts
        # latest things
        if self.latest_execution_context:
            hsh["latestExecutionContext"] = self.latest_execution_context
        if self.latest_activity_task_timestamp:
            hsh["latestActivityTaskTimestamp"] = self.latest_activity_task_timestamp
        return hsh

    def to_list_dict(self) -> Dict[str, Any]:
        hsh: Dict[str, Any] = {
            "execution": {"workflowId": self.workflow_id, "runId": self.run_id},
            "workflowType": self.workflow_type.to_short_dict(),
            "startTimestamp": self.start_timestamp,
            "executionStatus": self.execution_status,
            "cancelRequested": self.cancel_requested,
        }
        if self.tag_list:
            hsh["tagList"] = self.tag_list
        if self.parent:
            hsh["parent"] = self.parent
        if self.close_status:
            hsh["closeStatus"] = self.close_status
        if self.close_timestamp:
            hsh["closeTimestamp"] = self.close_timestamp
        return hsh

    def _process_timeouts(self) -> None:
        """
        SWF timeouts can happen on different objects (workflow executions,
        activity tasks, decision tasks) and should be processed in order.

        A specific timeout can change the workflow execution state and have an
        impact on other timeouts: for instance, if the workflow execution
        timeouts, subsequent timeouts on activity or decision tasks are
        irrelevant ; if an activity task timeouts, other timeouts on this task
        are irrelevant, and a new decision is fired, which could well timeout
        before the end of the workflow.

        So the idea here is to find the earliest timeout that would have been
        triggered, process it, then make the workflow state progress and repeat
        the whole process.
        """
        # workflow execution timeout
        timeout_candidates_or_none = [self.first_timeout()]

        # decision tasks timeouts
        for d_task in self.decision_tasks:
            timeout_candidates_or_none.append(d_task.first_timeout())

        # activity tasks timeouts
        for a_task in self.activity_tasks:
            timeout_candidates_or_none.append(a_task.first_timeout())

        # remove blank values (foo.first_timeout() is a Timeout or None)
        timeout_candidates = list(filter(None, timeout_candidates_or_none))

        # now find the first timeout to process
        first_timeout = None
        if timeout_candidates:
            first_timeout = min(timeout_candidates, key=lambda t: t.timestamp)

        if first_timeout:
            should_schedule_decision_next = False
            if isinstance(first_timeout.obj, WorkflowExecution):
                self.timeout(first_timeout)
            elif isinstance(first_timeout.obj, DecisionTask):
                self.timeout_decision_task(first_timeout)
                should_schedule_decision_next = True
            elif isinstance(first_timeout.obj, ActivityTask):
                self.timeout_activity_task(first_timeout)
                should_schedule_decision_next = True
            else:
                raise NotImplementedError("Unhandled timeout object")

            # schedule decision task if needed
            if should_schedule_decision_next:
                self.schedule_decision_task()

            # the workflow execution progressed, let's see if another
            # timeout should be processed
            self._process_timeouts()

    def events(self, reverse_order: bool = False) -> Iterable[HistoryEvent]:
        if reverse_order:
            return reversed(self._events)
        else:
            return self._events

    def next_event_id(self) -> int:
        event_ids = [evt.event_id for evt in self._events]
        return max(event_ids or [0]) + 1

    def _add_event(self, *args: Any, **kwargs: Any) -> HistoryEvent:
        # lock here because the fire_timer function is called
        # async, and want to ensure uniqueness in event ids
        with self.threading_lock:
            evt = HistoryEvent(self.next_event_id(), *args, **kwargs)
            self._events.append(evt)
            return evt

    def start(self) -> None:
        self.start_timestamp = unix_time()
        self._add_event(
            "WorkflowExecutionStarted",
            child_policy=self.child_policy,
            execution_start_to_close_timeout=self.execution_start_to_close_timeout,
            # TODO: fix this hardcoded value
            parent_initiated_event_id=0,
            task_list=self.task_list,
            task_start_to_close_timeout=self.task_start_to_close_timeout,
            workflow_type=self.workflow_type,
            input=self.input,
        )
        self.schedule_decision_task()

    def _schedule_decision_task(self) -> None:
        has_scheduled_task = False
        has_started_task = False
        for task in self.decision_tasks:
            if task.state == "STARTED":
                has_started_task = True
            elif task.state == "SCHEDULED":
                has_scheduled_task = True
        # If a decision task is already running, we cannot schedule more than one additional task
        # See https://docs.aws.amazon.com/amazonswf/latest/developerguide/swf-dg-dev-deciders.html#swf-dg-deciders-launch
        if has_started_task and has_scheduled_task:
            return

        evt = self._add_event(
            "DecisionTaskScheduled",
            start_to_close_timeout=self.task_start_to_close_timeout,
            task_list=self.task_list,
        )
        self.domain.add_to_decision_task_list(
            self.task_list, DecisionTask(self, evt.event_id)
        )
        self.open_counts["openDecisionTasks"] += 1

    def schedule_decision_task(self) -> None:
        self._schedule_decision_task()

    # Shortcut for tests: helps having auto-starting decision tasks when needed
    def schedule_and_start_decision_task(self, identity: Optional[str] = None) -> None:
        self._schedule_decision_task()
        decision_task = self.decision_tasks[-1]
        self.start_decision_task(decision_task.task_token, identity=identity)

    @property
    def decision_tasks(self) -> List[DecisionTask]:
        return [t for t in self.domain.decision_tasks if t.workflow_execution == self]

    @property
    def activity_tasks(self) -> List[ActivityTask]:
        return [t for t in self.domain.activity_tasks if t.workflow_execution == self]

    def _find_decision_task(self, task_token: str) -> DecisionTask:
        for dt in self.decision_tasks:
            if dt.task_token == task_token:
                return dt
        raise ValueError(f"No decision task with token: {task_token}")

    def start_decision_task(
        self, task_token: str, identity: Optional[str] = None
    ) -> None:
        dt = self._find_decision_task(task_token)
        evt = self._add_event(
            "DecisionTaskStarted",
            scheduled_event_id=dt.scheduled_event_id,
            identity=identity,
        )
        dt.start(evt.event_id, self._previous_started_event_id)
        self._previous_started_event_id = evt.event_id

    def complete_decision_task(
        self,
        task_token: str,
        decisions: Optional[List[Dict[str, Any]]] = None,
        execution_context: Optional[str] = None,
    ) -> None:
        # 'decisions' can be None per boto.swf defaults, so replace it with something iterable
        if not decisions:
            decisions = []
        # In case of a malformed or invalid decision task, SWF will raise an error and
        # it won't perform any of the decisions in the decision set.
        self.validate_decisions(decisions)
        dt = self._find_decision_task(task_token)
        evt = self._add_event(
            "DecisionTaskCompleted",
            scheduled_event_id=dt.scheduled_event_id,
            started_event_id=dt.started_event_id,
            execution_context=execution_context,
        )
        dt.complete()
        self.should_schedule_decision_next = False
        self.handle_decisions(evt.event_id, decisions)
        if self.should_schedule_decision_next:
            self.schedule_decision_task()
        self.latest_execution_context = execution_context

    def _check_decision_attributes(
        self, kind: str, value: Dict[str, Any], decision_id: int
    ) -> List[Dict[str, str]]:
        problems = []
        constraints = DECISIONS_FIELDS.get(kind, {})
        for key, constraint in constraints.items():
            if constraint["required"] and not value.get(key):
                problems.append(
                    {
                        "type": "null_value",
                        "where": f"decisions.{decision_id}.member.{kind}.{key}",
                    }
                )
        return problems

    def validate_decisions(self, decisions: List[Dict[str, Any]]) -> None:
        """
        Performs some basic validations on decisions. The real SWF service
        seems to break early and *not* process any decision if there's a
        validation problem, such as a malformed decision for instance. I didn't
        find an explicit documentation for that though, so criticisms welcome.
        """
        problems = []

        # check close decision is last
        # (the real SWF service also works that way if you provide 2 close decision tasks)
        for dcs in decisions[:-1]:
            close_decision_types = [
                "CompleteWorkflowExecution",
                "FailWorkflowExecution",
                "CancelWorkflowExecution",
            ]
            if dcs["decisionType"] in close_decision_types:
                raise SWFValidationException("Close must be last decision in list")

        decision_number = 0
        for dcs in decisions:
            decision_number += 1
            # check decision types mandatory attributes
            # NB: the real SWF service seems to check attributes even for attributes list
            # that are not in line with the decisionType, so we do the same
            attrs_to_check = [d for d in dcs.keys() if d.endswith("DecisionAttributes")]
            if dcs["decisionType"] in self.KNOWN_DECISION_TYPES:
                decision_type = dcs["decisionType"]
                decision_attr = f"{decapitalize(decision_type)}DecisionAttributes"
                attrs_to_check.append(decision_attr)
            for attr in attrs_to_check:
                problems += self._check_decision_attributes(
                    attr, dcs.get(attr, {}), decision_number
                )
            # check decision type is correct
            if dcs["decisionType"] not in self.KNOWN_DECISION_TYPES:
                problems.append(
                    {
                        "type": "bad_decision_type",
                        "value": dcs["decisionType"],
                        "where": f"decisions.{decision_number}.member.decisionType",
                        "possible_values": ", ".join(self.KNOWN_DECISION_TYPES),
                    }
                )

        # raise if any problem
        if any(problems):
            raise SWFDecisionValidationException(problems)

    def handle_decisions(self, event_id: int, decisions: List[Dict[str, Any]]) -> None:
        """
        Handles a Decision according to SWF docs.
        See: http://docs.aws.amazon.com/amazonswf/latest/apireference/API_Decision.html
        """
        # handle each decision separately, in order
        for decision in decisions:
            decision_type = decision["decisionType"]
            attributes_key = f"{decapitalize(decision_type)}DecisionAttributes"
            attributes = decision.get(attributes_key, {})
            if decision_type == "CompleteWorkflowExecution":
                self.complete(event_id, attributes.get("result"))
            elif decision_type == "FailWorkflowExecution":
                self.fail(event_id, attributes.get("details"), attributes.get("reason"))
            elif decision_type == "ScheduleActivityTask":
                self.schedule_activity_task(event_id, attributes)
            elif decision_type == "RecordMarker":
                self.record_marker(event_id, attributes)
            elif decision_type == "StartTimer":
                self.start_timer(event_id, attributes)
            elif decision_type == "CancelTimer":
                self.cancel_timer(event_id, attributes["timerId"])
            elif decision_type == "CancelWorkflowExecution":
                self.cancel(event_id, attributes.get("details"))
            else:
                # TODO: implement Decision type: ContinueAsNewWorkflowExecution
                # TODO: implement Decision type: RequestCancelActivityTask
                # TODO: implement Decision type: RequestCancelExternalWorkflowExecution
                # TODO: implement Decision type: ScheduleLambdaFunction
                # TODO: implement Decision type: SignalExternalWorkflowExecution
                # TODO: implement Decision type: StartChildWorkflowExecution
                raise NotImplementedError(f"Cannot handle decision: {decision_type}")

        # finally decrement counter if and only if everything went well
        self.open_counts["openDecisionTasks"] -= 1

    def complete(self, event_id: int, result: Any = None) -> None:
        self.execution_status = "CLOSED"
        self.close_status = "COMPLETED"
        self.close_timestamp = unix_time()
        self._add_event(
            "WorkflowExecutionCompleted",
            decision_task_completed_event_id=event_id,
            result=result,
        )

    def fail(
        self, event_id: int, details: Any = None, reason: Optional[str] = None
    ) -> None:
        # TODO: implement length constraints on details/reason
        self.execution_status = "CLOSED"
        self.close_status = "FAILED"
        self.close_timestamp = unix_time()
        self._add_event(
            "WorkflowExecutionFailed",
            decision_task_completed_event_id=event_id,
            details=details,
            reason=reason,
        )

    def cancel(self, event_id: int, details: Any = None) -> None:
        # TODO: implement length constraints on details
        self.cancel_requested = True
        # Can only cancel if there are no other pending desicion tasks
        if self.open_counts["openDecisionTasks"] != 1:
            # TODO OPERATION_NOT_PERMITTED is a valid failure state
            self._add_event(
                "CancelWorkflowExecutionFailed",
                decision_task_completed_event_id=event_id,
                cause="UNHANDLED_DECISION",
            )
            return
        self.execution_status = "CLOSED"
        self.close_status = "CANCELED"
        self.close_timestamp = unix_time()
        self._add_event(
            "WorkflowExecutionCanceled",
            decision_task_completed_event_id=event_id,
            details=details,
        )

    def schedule_activity_task(self, event_id: int, attributes: Dict[str, Any]) -> None:
        # Helper function to avoid repeating ourselves in the next sections
        def fail_schedule_activity_task(_type: "ActivityType", _cause: str) -> None:
            # TODO: implement other possible failure mode: OPEN_ACTIVITIES_LIMIT_EXCEEDED
            # NB: some failure modes are not implemented and probably won't be implemented in
            # the future, such as ACTIVITY_CREATION_RATE_EXCEEDED or
            # OPERATION_NOT_PERMITTED
            self._add_event(
                "ScheduleActivityTaskFailed",
                activity_id=attributes["activityId"],
                activity_type=_type,
                cause=_cause,
                decision_task_completed_event_id=event_id,
            )
            self.should_schedule_decision_next = True

        activity_type: ActivityType = self.domain.get_type(  # type: ignore[assignment]
            "activity",
            attributes["activityType"]["name"],
            attributes["activityType"]["version"],
            ignore_empty=True,
        )
        if not activity_type:
            fake_type = ActivityType(
                attributes["activityType"]["name"],
                attributes["activityType"]["version"],
            )
            fail_schedule_activity_task(fake_type, "ACTIVITY_TYPE_DOES_NOT_EXIST")
            return
        if activity_type.status == "DEPRECATED":
            fail_schedule_activity_task(activity_type, "ACTIVITY_TYPE_DEPRECATED")
            return
        if any(
            at
            for at in self.activity_tasks
            if at.activity_id == attributes["activityId"]
        ):
            fail_schedule_activity_task(activity_type, "ACTIVITY_ID_ALREADY_IN_USE")
            return

        # find task list or default task list, else fail
        task_list = attributes.get("taskList", {}).get("name")
        if not task_list and activity_type.task_list:
            task_list = activity_type.task_list
        if not task_list:
            fail_schedule_activity_task(activity_type, "DEFAULT_TASK_LIST_UNDEFINED")
            return

        # find timeouts or default timeout, else fail
        timeouts = {}
        for _type in [
            "scheduleToStartTimeout",
            "scheduleToCloseTimeout",
            "startToCloseTimeout",
            "heartbeatTimeout",
        ]:
            default_key = "default_task_" + camelcase_to_underscores(_type)
            default_value = getattr(activity_type, default_key)
            timeouts[_type] = attributes.get(_type, default_value)
            if not timeouts[_type]:
                error_key = default_key.replace("default_task_", "default_")
                fail_schedule_activity_task(
                    activity_type, f"{error_key.upper()}_UNDEFINED"
                )
                return

        # Only add event and increment counters now that nothing went wrong
        evt = self._add_event(
            "ActivityTaskScheduled",
            activity_id=attributes["activityId"],
            activity_type=activity_type,
            control=attributes.get("control"),
            decision_task_completed_event_id=event_id,
            heartbeat_timeout=attributes.get("heartbeatTimeout"),
            input=attributes.get("input"),
            schedule_to_close_timeout=attributes.get("scheduleToCloseTimeout"),
            schedule_to_start_timeout=attributes.get("scheduleToStartTimeout"),
            start_to_close_timeout=attributes.get("startToCloseTimeout"),
            task_list=task_list,
            task_priority=attributes.get("taskPriority"),
        )
        task = ActivityTask(
            activity_id=attributes["activityId"],
            activity_type=activity_type,
            workflow_input=attributes.get("input"),
            scheduled_event_id=evt.event_id,
            workflow_execution=self,
            timeouts=timeouts,
        )
        self.domain.add_to_activity_task_list(task_list, task)
        self.open_counts["openActivityTasks"] += 1
        self.latest_activity_task_timestamp = unix_time()

    def _find_activity_task(self, task_token: str) -> ActivityTask:
        for task in self.activity_tasks:
            if task.task_token == task_token:
                return task
        raise ValueError(f"No activity task with token: {task_token}")

    def start_activity_task(self, task_token: str, identity: Any = None) -> None:
        task = self._find_activity_task(task_token)
        evt = self._add_event(
            "ActivityTaskStarted",
            scheduled_event_id=task.scheduled_event_id,
            identity=identity,
        )
        task.start(evt.event_id)

    def complete_activity_task(self, task_token: str, result: Any = None) -> None:
        task = self._find_activity_task(task_token)
        self._add_event(
            "ActivityTaskCompleted",
            scheduled_event_id=task.scheduled_event_id,
            started_event_id=task.started_event_id,
            result=result,
        )
        task.complete()
        self.open_counts["openActivityTasks"] -= 1
        # TODO: ensure we don't schedule multiple decisions at the same time!
        self.schedule_decision_task()

    def fail_activity_task(
        self, task_token: str, reason: Optional[str] = None, details: Any = None
    ) -> None:
        task = self._find_activity_task(task_token)
        self._add_event(
            "ActivityTaskFailed",
            scheduled_event_id=task.scheduled_event_id,
            started_event_id=task.started_event_id,
            reason=reason,
            details=details,
        )
        task.fail()
        self.open_counts["openActivityTasks"] -= 1
        # TODO: ensure we don't schedule multiple decisions at the same time!
        self.schedule_decision_task()

    def terminate(
        self,
        child_policy: Optional[str] = None,
        details: Optional[Dict[str, Any]] = None,
        reason: Optional[str] = None,
    ) -> None:
        # TODO: handle child policy for child workflows here
        # TODO: handle cause="CHILD_POLICY_APPLIED"
        # Until this, we set cause manually to "OPERATOR_INITIATED"
        cause = "OPERATOR_INITIATED"
        if not child_policy:
            child_policy = self.child_policy
        self._add_event(
            "WorkflowExecutionTerminated",
            cause=cause,
            child_policy=child_policy,
            details=details,
            reason=reason,
        )
        self.execution_status = "CLOSED"
        self.close_status = "TERMINATED"
        self.close_cause = "OPERATOR_INITIATED"

    def signal(self, signal_name: str, workflow_input: Dict[str, Any]) -> None:
        self._add_event(
            "WorkflowExecutionSignaled", signal_name=signal_name, input=workflow_input
        )
        self.schedule_decision_task()

    def first_timeout(self) -> Optional[Timeout]:
        if not self.open or not self.start_timestamp:
            return None
        start_to_close_at = self.start_timestamp + int(
            self.execution_start_to_close_timeout
        )
        _timeout = Timeout(self, start_to_close_at, "START_TO_CLOSE")
        if _timeout.reached:
            return _timeout
        return None

    def timeout(self, timeout: Timeout) -> None:
        # TODO: process child policy on child workflows here or in the
        # triggering function
        self.execution_status = "CLOSED"
        self.close_status = "TIMED_OUT"
        self.timeout_type = timeout.kind
        self._add_event(
            "WorkflowExecutionTimedOut",
            child_policy=self.child_policy,
            event_timestamp=timeout.timestamp,
            timeout_type=self.timeout_type,
        )

    def timeout_decision_task(self, _timeout: Timeout) -> None:
        task = _timeout.obj
        task.timeout(_timeout)
        self._add_event(
            "DecisionTaskTimedOut",
            event_timestamp=_timeout.timestamp,
            scheduled_event_id=task.scheduled_event_id,
            started_event_id=task.started_event_id,
            timeout_type=task.timeout_type,
        )

    def timeout_activity_task(self, _timeout: Timeout) -> None:
        task = _timeout.obj
        task.timeout(_timeout)
        self._add_event(
            "ActivityTaskTimedOut",
            details=task.details,
            event_timestamp=_timeout.timestamp,
            scheduled_event_id=task.scheduled_event_id,
            started_event_id=task.started_event_id,
            timeout_type=task.timeout_type,
        )

    def record_marker(self, event_id: int, attributes: Dict[str, Any]) -> None:
        self._add_event(
            "MarkerRecorded",
            decision_task_completed_event_id=event_id,
            details=attributes.get("details"),
            marker_name=attributes["markerName"],
        )

    def start_timer(self, event_id: int, attributes: Dict[str, Any]) -> None:
        timer_id = attributes["timerId"]
        existing_timer = self._timers.get(timer_id)
        if existing_timer and existing_timer.is_alive():
            # TODO 4 fail states are possible
            # TIMER_ID_ALREADY_IN_USE | OPEN_TIMERS_LIMIT_EXCEEDED | TIMER_CREATION_RATE_EXCEEDED | OPERATION_NOT_PERMITTED
            self._add_event(
                "StartTimerFailed",
                cause="TIMER_ID_ALREADY_IN_USE",
                decision_task_completed_event_id=event_id,
                timer_id=timer_id,
            )
            return

        time_to_wait = attributes["startToFireTimeout"]
        started_event_id = self._add_event(
            "TimerStarted",
            control=attributes.get("control"),
            decision_task_completed_event_id=event_id,
            start_to_fire_timeout=time_to_wait,
            timer_id=timer_id,
        ).event_id
        background_timer = ThreadingTimer(
            float(time_to_wait), self._fire_timer, args=(started_event_id, timer_id)
        )
        workflow_timer = Timer(background_timer, started_event_id)
        self._timers[timer_id] = workflow_timer
        workflow_timer.start()

    def _fire_timer(self, started_event_id: int, timer_id: str) -> None:
        self._add_event(
            "TimerFired", started_event_id=started_event_id, timer_id=timer_id
        )
        self._timers.pop(timer_id)
        self._schedule_decision_task()

    def cancel_timer(self, event_id: int, timer_id: str) -> None:
        requested_timer = self._timers.get(timer_id)
        if not requested_timer or not requested_timer.is_alive():
            # TODO there are 2 failure states
            # TIMER_ID_UNKNOWN | OPERATION_NOT_PERMITTED
            self._add_event(
                "CancelTimerFailed",
                cause="TIMER_ID_UNKNOWN",
                decision_task_completed_event_id=event_id,
            )
            return

        requested_timer.cancel()
        self._timers.pop(timer_id)
        self._add_event(
            "TimerCancelled",
            decision_task_completed_event_id=event_id,
            started_event_id=requested_timer.started_event_id,
            timer_id=timer_id,
        )

    @property
    def open(self) -> bool:
        return self.execution_status == "OPEN"
