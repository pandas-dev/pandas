import abc
import json
import time

from moto.stepfunctions.parser.api import (
    HistoryEventExecutionDataDetails,
    HistoryEventType,
    TaskFailedEventDetails,
    TaskSubmittedEventDetails,
)
from moto.stepfunctions.parser.asl.component.common.error_name.custom_error_name import (
    CustomErrorName,
)
from moto.stepfunctions.parser.asl.component.common.error_name.failure_event import (
    FailureEvent,
)
from moto.stepfunctions.parser.asl.component.state.exec.state_task.service.resource import (
    ResourceCondition,
    ResourceRuntimePart,
)
from moto.stepfunctions.parser.asl.component.state.exec.state_task.service.state_task_service import (
    StateTaskService,
)
from moto.stepfunctions.parser.asl.eval.callback.callback import (
    CallbackOutcomeFailure,
    CallbackOutcomeFailureError,
    CallbackOutcomeSuccess,
    CallbackTimeoutError,
    HeartbeatEndpoint,
    HeartbeatTimeoutError,
)
from moto.stepfunctions.parser.asl.eval.environment import Environment
from moto.stepfunctions.parser.asl.eval.event.event_detail import EventDetails
from moto.stepfunctions.parser.asl.utils.encoding import to_json_str


class StateTaskServiceCallback(StateTaskService, abc.ABC):
    def _get_sfn_resource(self) -> str:
        resource = super()._get_sfn_resource()
        if self.resource.condition is not None:
            resource += f".{self.resource.condition}"
        return resource

    def _wait_for_task_token(
        self,
        env: Environment,
        resource_runtime_part: ResourceRuntimePart,
        normalised_parameters: dict,
    ) -> None:
        # Discard the state evaluation output.
        env.stack.pop()

        callback_id = env.context_object_manager.context_object.get("Task", {}).get(
            "Token"
        )
        if callback_id is None:
            # We never generated a TaskToken, so there's no point in waiting until we receive it
            # AWS seems to just wait forever  in this scenario- so let's do the same
            timeout_seconds = self.timeout.timeout_seconds
            mustend_at = time.time() + timeout_seconds
            while time.time() < mustend_at:
                if not env.is_running():
                    return
                time.sleep(1)
            return

        callback_endpoint = env.callback_pool_manager.get(callback_id)

        # With Timeouts-only definition:
        if not self.heartbeat:
            self.timeout.eval(env=env)
            timeout_seconds = env.stack.pop()
            # Although the timeout is handled already be the superclass (ExecutionState),
            # the timeout value is specified here too, to allow this child process to terminate earlier even if
            # discarded by the main process.
            # Note: although this is the same timeout value, this can only decay strictly after the first timeout
            #       started as it is invoked strictly later.
            outcome = callback_endpoint.wait(timeout=timeout_seconds)
        else:
            self.heartbeat.eval(env=env)
            heartbeat_seconds = env.stack.pop()
            heartbeat_endpoint: HeartbeatEndpoint = (
                callback_endpoint.setup_heartbeat_endpoint(
                    heartbeat_seconds=heartbeat_seconds
                )
            )

            outcome = None
            while (
                env.is_running() and outcome is None
            ):  # Until subprocess hasn't timed out or result wasn't received.
                received = heartbeat_endpoint.clear_and_wait()
                if not received and env.is_running():  # Heartbeat timed out.
                    raise HeartbeatTimeoutError()
                outcome = callback_endpoint.get_outcome()

        if outcome is None:
            raise CallbackTimeoutError()
        if isinstance(outcome, CallbackOutcomeSuccess):
            outcome_output = json.loads(outcome.output)
            env.stack.append(outcome_output)
        elif isinstance(outcome, CallbackOutcomeFailure):
            raise CallbackOutcomeFailureError(callback_outcome_failure=outcome)
        else:
            raise NotImplementedError(
                f"Unsupported CallbackOutcome type '{type(outcome)}'."
            )

    def _sync(
        self,
        env: Environment,
        resource_runtime_part: ResourceRuntimePart,
        normalised_parameters: dict,
    ) -> None:
        raise RuntimeError(
            f"Unsupported .sync callback procedure in resource {self.resource.resource_arn}"
        )

    def _sync2(
        self,
        env: Environment,
        resource_runtime_part: ResourceRuntimePart,
        normalised_parameters: dict,
    ) -> None:
        raise RuntimeError(
            f"Unsupported .sync:2 callback procedure in resource {self.resource.resource_arn}"
        )

    def _is_condition(self):
        return self.resource.condition is not None

    def _get_callback_outcome_failure_event(
        self, ex: CallbackOutcomeFailureError
    ) -> FailureEvent:
        callback_outcome_failure: CallbackOutcomeFailure = ex.callback_outcome_failure
        error: str = callback_outcome_failure.error
        return FailureEvent(
            error_name=CustomErrorName(error_name=callback_outcome_failure.error),
            event_type=HistoryEventType.TaskFailed,
            event_details=EventDetails(
                taskFailedEventDetails=TaskFailedEventDetails(
                    resourceType=self._get_sfn_resource_type(),
                    resource=self._get_sfn_resource(),
                    error=error,
                    cause=callback_outcome_failure.cause,
                )
            ),
        )

    def _from_error(self, env: Environment, ex: Exception) -> FailureEvent:
        if isinstance(ex, CallbackOutcomeFailureError):
            return self._get_callback_outcome_failure_event(ex=ex)
        return super()._from_error(env=env, ex=ex)

    def _after_eval_execution(
        self,
        env: Environment,
        resource_runtime_part: ResourceRuntimePart,
        normalised_parameters: dict,
    ) -> None:
        if self._is_condition():
            output = env.stack[-1]
            env.event_history.add_event(
                context=env.event_history_context,
                hist_type_event=HistoryEventType.TaskSubmitted,
                event_detail=EventDetails(
                    taskSubmittedEventDetails=TaskSubmittedEventDetails(
                        resource=self._get_sfn_resource(),
                        resourceType=self._get_sfn_resource_type(),
                        output=to_json_str(output),
                        outputDetails=HistoryEventExecutionDataDetails(truncated=False),
                    )
                ),
            )
            if self.resource.condition == ResourceCondition.WaitForTaskToken:
                self._wait_for_task_token(
                    env=env,
                    resource_runtime_part=resource_runtime_part,
                    normalised_parameters=normalised_parameters,
                )
            elif self.resource.condition == ResourceCondition.Sync:
                self._sync(
                    env=env,
                    resource_runtime_part=resource_runtime_part,
                    normalised_parameters=normalised_parameters,
                )
            elif self.resource.condition == ResourceCondition.Sync2:
                self._sync2(
                    env=env,
                    resource_runtime_part=resource_runtime_part,
                    normalised_parameters=normalised_parameters,
                )
            else:
                raise NotImplementedError(
                    f"Unsupported callback type '{self.resource.condition}'."
                )

        super()._after_eval_execution(
            env=env,
            resource_runtime_part=resource_runtime_part,
            normalised_parameters=normalised_parameters,
        )
