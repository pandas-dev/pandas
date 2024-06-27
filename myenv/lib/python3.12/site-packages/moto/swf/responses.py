import json
from typing import Any, List

from moto.core.responses import BaseResponse

from .exceptions import SWFSerializationException, SWFValidationException
from .models import GenericType, SWFBackend, swf_backends


class SWFResponse(BaseResponse):
    def __init__(self) -> None:
        super().__init__(service_name="swf")

    @property
    def swf_backend(self) -> SWFBackend:
        return swf_backends[self.current_account][self.region]

    # SWF parameters are passed through a JSON body, so let's ease retrieval
    @property
    def _params(self) -> Any:  # type: ignore[misc]
        return json.loads(self.body)

    def _check_int(self, parameter: Any) -> None:
        if not isinstance(parameter, int):
            raise SWFSerializationException(parameter)

    def _check_float_or_int(self, parameter: Any) -> None:
        if not isinstance(parameter, float):
            if not isinstance(parameter, int):
                raise SWFSerializationException(parameter)

    def _check_none_or_string(self, parameter: Any) -> None:
        if parameter is not None:
            self._check_string(parameter)

    def _check_string(self, parameter: Any) -> None:
        if not isinstance(parameter, str):
            raise SWFSerializationException(parameter)

    def _check_none_or_list_of_strings(self, parameter: Any) -> None:
        if parameter is not None:
            self._check_list_of_strings(parameter)

    def _check_list_of_strings(self, parameter: Any) -> None:
        if not isinstance(parameter, list):
            raise SWFSerializationException(parameter)
        for i in parameter:
            if not isinstance(i, str):
                raise SWFSerializationException(parameter)

    def _check_exclusivity(self, **kwargs: Any) -> None:
        if list(kwargs.values()).count(None) >= len(kwargs) - 1:
            return
        keys = kwargs.keys()
        if len(keys) == 2:
            message = f"Cannot specify both a {keys[0]} and a {keys[1]}"  # type: ignore
        else:
            message = f"Cannot specify more than one exclusive filters in the same query: {keys}"
        raise SWFValidationException(message)

    def _list_types(self, kind: str) -> str:
        domain_name = self._params["domain"]
        status = self._params["registrationStatus"]
        reverse_order = self._params.get("reverseOrder", None)
        self._check_string(domain_name)
        self._check_string(status)
        types: List[GenericType] = self.swf_backend.list_types(
            kind, domain_name, status, reverse_order=reverse_order
        )
        return json.dumps({"typeInfos": [_type.to_medium_dict() for _type in types]})

    def _describe_type(self, kind: str) -> str:
        domain = self._params["domain"]
        _type_args = self._params[f"{kind}Type"]
        name = _type_args["name"]
        version = _type_args["version"]
        self._check_string(domain)
        self._check_string(name)
        self._check_string(version)
        _type = self.swf_backend.describe_type(kind, domain, name, version)

        return json.dumps(_type.to_full_dict())

    def _deprecate_type(self, kind: str) -> str:
        domain = self._params["domain"]
        _type_args = self._params[f"{kind}Type"]
        name = _type_args["name"]
        version = _type_args["version"]
        self._check_string(domain)
        self._check_string(name)
        self._check_string(version)
        self.swf_backend.deprecate_type(kind, domain, name, version)
        return ""

    def _undeprecate_type(self, kind: str) -> str:
        domain = self._params["domain"]
        _type_args = self._params[f"{kind}Type"]
        name = _type_args["name"]
        version = _type_args["version"]
        self._check_string(domain)
        self._check_string(name)
        self._check_string(version)
        self.swf_backend.undeprecate_type(kind, domain, name, version)
        return ""

    # TODO: implement pagination
    def list_domains(self) -> str:
        status = self._params["registrationStatus"]
        self._check_string(status)
        reverse_order = self._params.get("reverseOrder", None)
        domains = self.swf_backend.list_domains(status, reverse_order=reverse_order)
        return json.dumps(
            {"domainInfos": [domain.to_short_dict() for domain in domains]}
        )

    def list_closed_workflow_executions(self) -> str:
        domain = self._params["domain"]
        start_time_filter = self._params.get("startTimeFilter", None)
        close_time_filter = self._params.get("closeTimeFilter", None)
        execution_filter = self._params.get("executionFilter", None)
        workflow_id = execution_filter["workflowId"] if execution_filter else None
        maximum_page_size = self._params.get("maximumPageSize", 1000)
        reverse_order = self._params.get("reverseOrder", None)
        tag_filter = self._params.get("tagFilter", None)
        type_filter = self._params.get("typeFilter", None)
        close_status_filter = self._params.get("closeStatusFilter", None)

        self._check_string(domain)
        self._check_none_or_string(workflow_id)
        self._check_exclusivity(
            executionFilter=execution_filter,
            typeFilter=type_filter,
            tagFilter=tag_filter,
            closeStatusFilter=close_status_filter,
        )
        self._check_exclusivity(
            startTimeFilter=start_time_filter, closeTimeFilter=close_time_filter
        )
        if start_time_filter is None and close_time_filter is None:
            raise SWFValidationException("Must specify time filter")
        if start_time_filter:
            self._check_float_or_int(start_time_filter["oldestDate"])
            if "latestDate" in start_time_filter:
                self._check_float_or_int(start_time_filter["latestDate"])
        if close_time_filter:
            self._check_float_or_int(close_time_filter["oldestDate"])
            if "latestDate" in close_time_filter:
                self._check_float_or_int(close_time_filter["latestDate"])
        if tag_filter:
            self._check_string(tag_filter["tag"])
        if type_filter:
            self._check_string(type_filter["name"])
            self._check_string(type_filter["version"])
        if close_status_filter:
            self._check_string(close_status_filter["status"])
        self._check_int(maximum_page_size)

        workflow_executions = self.swf_backend.list_closed_workflow_executions(
            domain_name=domain,
            tag_filter=tag_filter,
            close_status_filter=close_status_filter,
            maximum_page_size=maximum_page_size,
            reverse_order=reverse_order,
        )

        return json.dumps(
            {"executionInfos": [wfe.to_list_dict() for wfe in workflow_executions]}
        )

    def list_open_workflow_executions(self) -> str:
        domain = self._params["domain"]
        start_time_filter = self._params["startTimeFilter"]
        execution_filter = self._params.get("executionFilter", None)
        workflow_id = execution_filter["workflowId"] if execution_filter else None
        maximum_page_size = self._params.get("maximumPageSize", 1000)
        reverse_order = self._params.get("reverseOrder", None)
        tag_filter = self._params.get("tagFilter", None)
        type_filter = self._params.get("typeFilter", None)

        self._check_string(domain)
        self._check_none_or_string(workflow_id)
        self._check_exclusivity(
            executionFilter=execution_filter,
            typeFilter=type_filter,
            tagFilter=tag_filter,
        )
        self._check_float_or_int(start_time_filter["oldestDate"])
        if "latestDate" in start_time_filter:
            self._check_float_or_int(start_time_filter["latestDate"])
        if tag_filter:
            self._check_string(tag_filter["tag"])
        if type_filter:
            self._check_string(type_filter["name"])
            self._check_string(type_filter["version"])
        self._check_int(maximum_page_size)

        workflow_executions = self.swf_backend.list_open_workflow_executions(
            domain_name=domain,
            maximum_page_size=maximum_page_size,
            tag_filter=tag_filter,
            reverse_order=reverse_order,
        )

        return json.dumps(
            {"executionInfos": [wfe.to_list_dict() for wfe in workflow_executions]}
        )

    def register_domain(self) -> str:
        name = self._params["name"]
        retention = self._params["workflowExecutionRetentionPeriodInDays"]
        description = self._params.get("description")
        self._check_string(retention)
        self._check_string(name)
        self._check_none_or_string(description)
        self.swf_backend.register_domain(name, retention, description=description)
        return ""

    def deprecate_domain(self) -> str:
        name = self._params["name"]
        self._check_string(name)
        self.swf_backend.deprecate_domain(name)
        return ""

    def undeprecate_domain(self) -> str:
        name = self._params["name"]
        self._check_string(name)
        self.swf_backend.undeprecate_domain(name)
        return ""

    def describe_domain(self) -> str:
        name = self._params["name"]
        self._check_string(name)
        domain = self.swf_backend.describe_domain(name)
        return json.dumps(domain.to_full_dict())  # type: ignore[union-attr]

    # TODO: implement pagination
    def list_activity_types(self) -> str:
        return self._list_types("activity")

    def register_activity_type(self) -> str:
        domain = self._params["domain"]
        name = self._params["name"]
        version = self._params["version"]
        default_task_list = self._params.get("defaultTaskList")
        if default_task_list:
            task_list = default_task_list.get("name")
        else:
            task_list = None
        default_task_heartbeat_timeout = self._params.get("defaultTaskHeartbeatTimeout")
        default_task_schedule_to_close_timeout = self._params.get(
            "defaultTaskScheduleToCloseTimeout"
        )
        default_task_schedule_to_start_timeout = self._params.get(
            "defaultTaskScheduleToStartTimeout"
        )
        default_task_start_to_close_timeout = self._params.get(
            "defaultTaskStartToCloseTimeout"
        )
        description = self._params.get("description")

        self._check_string(domain)
        self._check_string(name)
        self._check_string(version)
        self._check_none_or_string(task_list)
        self._check_none_or_string(default_task_heartbeat_timeout)
        self._check_none_or_string(default_task_schedule_to_close_timeout)
        self._check_none_or_string(default_task_schedule_to_start_timeout)
        self._check_none_or_string(default_task_start_to_close_timeout)
        self._check_none_or_string(description)

        # TODO: add defaultTaskPriority when boto gets to support it
        self.swf_backend.register_type(
            "activity",
            domain,
            name,
            version,
            task_list=task_list,
            default_task_heartbeat_timeout=default_task_heartbeat_timeout,
            default_task_schedule_to_close_timeout=default_task_schedule_to_close_timeout,
            default_task_schedule_to_start_timeout=default_task_schedule_to_start_timeout,
            default_task_start_to_close_timeout=default_task_start_to_close_timeout,
            description=description,
        )
        return ""

    def deprecate_activity_type(self) -> str:
        return self._deprecate_type("activity")

    def undeprecate_activity_type(self) -> str:
        return self._undeprecate_type("activity")

    def describe_activity_type(self) -> str:
        return self._describe_type("activity")

    def list_workflow_types(self) -> str:
        return self._list_types("workflow")

    def register_workflow_type(self) -> str:
        domain = self._params["domain"]
        name = self._params["name"]
        version = self._params["version"]
        default_task_list = self._params.get("defaultTaskList")
        if default_task_list:
            task_list = default_task_list.get("name")
        else:
            task_list = None
        default_child_policy = self._params.get("defaultChildPolicy")
        default_task_start_to_close_timeout = self._params.get(
            "defaultTaskStartToCloseTimeout"
        )
        default_execution_start_to_close_timeout = self._params.get(
            "defaultExecutionStartToCloseTimeout"
        )
        default_task_priority = self._params.get("defaultTaskPriority")
        default_lambda_role = self._params.get("defaultLambdaRole")
        description = self._params.get("description")

        self._check_string(domain)
        self._check_string(name)
        self._check_string(version)
        self._check_none_or_string(task_list)
        self._check_none_or_string(default_child_policy)
        self._check_none_or_string(default_task_start_to_close_timeout)
        self._check_none_or_string(default_execution_start_to_close_timeout)
        self._check_none_or_string(default_task_priority)
        self._check_none_or_string(default_lambda_role)
        self._check_none_or_string(description)

        self.swf_backend.register_type(
            "workflow",
            domain,
            name,
            version,
            task_list=task_list,
            default_child_policy=default_child_policy,
            default_task_start_to_close_timeout=default_task_start_to_close_timeout,
            default_execution_start_to_close_timeout=default_execution_start_to_close_timeout,
            default_task_priority=default_task_priority,
            default_lambda_role=default_lambda_role,
            description=description,
        )
        return ""

    def deprecate_workflow_type(self) -> str:
        return self._deprecate_type("workflow")

    def undeprecate_workflow_type(self) -> str:
        return self._undeprecate_type("workflow")

    def describe_workflow_type(self) -> str:
        return self._describe_type("workflow")

    def start_workflow_execution(self) -> str:
        domain = self._params["domain"]
        workflow_id = self._params["workflowId"]
        _workflow_type = self._params["workflowType"]
        workflow_name = _workflow_type["name"]
        workflow_version = _workflow_type["version"]
        _task_list = self._params.get("taskList")
        if _task_list:
            task_list = _task_list.get("name")
        else:
            task_list = None
        child_policy = self._params.get("childPolicy")
        execution_start_to_close_timeout = self._params.get(
            "executionStartToCloseTimeout"
        )
        input_ = self._params.get("input")
        tag_list = self._params.get("tagList")
        task_start_to_close_timeout = self._params.get("taskStartToCloseTimeout")

        self._check_string(domain)
        self._check_string(workflow_id)
        self._check_string(workflow_name)
        self._check_string(workflow_version)
        self._check_none_or_string(task_list)
        self._check_none_or_string(child_policy)
        self._check_none_or_string(execution_start_to_close_timeout)
        self._check_none_or_string(input_)
        self._check_none_or_list_of_strings(tag_list)
        self._check_none_or_string(task_start_to_close_timeout)

        wfe = self.swf_backend.start_workflow_execution(
            domain,
            workflow_id,
            workflow_name,
            workflow_version,
            task_list=task_list,
            child_policy=child_policy,
            execution_start_to_close_timeout=execution_start_to_close_timeout,
            workflow_input=input_,
            tag_list=tag_list,
            task_start_to_close_timeout=task_start_to_close_timeout,
        )

        return json.dumps({"runId": wfe.run_id})

    def describe_workflow_execution(self) -> str:
        domain_name = self._params["domain"]
        _workflow_execution = self._params["execution"]
        run_id = _workflow_execution["runId"]
        workflow_id = _workflow_execution["workflowId"]

        self._check_string(domain_name)
        self._check_string(run_id)
        self._check_string(workflow_id)

        wfe = self.swf_backend.describe_workflow_execution(
            domain_name, run_id, workflow_id
        )
        return json.dumps(wfe.to_full_dict())  # type: ignore[union-attr]

    def get_workflow_execution_history(self) -> str:
        domain_name = self._params["domain"]
        _workflow_execution = self._params["execution"]
        run_id = _workflow_execution["runId"]
        workflow_id = _workflow_execution["workflowId"]
        reverse_order = self._params.get("reverseOrder", None)
        wfe = self.swf_backend.describe_workflow_execution(
            domain_name, run_id, workflow_id
        )
        events = wfe.events(reverse_order=reverse_order)  # type: ignore[union-attr]
        return json.dumps({"events": [evt.to_dict() for evt in events]})

    def poll_for_decision_task(self) -> str:
        domain_name = self._params["domain"]
        task_list = self._params["taskList"]["name"]
        identity = self._params.get("identity")
        reverse_order = self._params.get("reverseOrder", None)

        self._check_string(domain_name)
        self._check_string(task_list)

        decision = self.swf_backend.poll_for_decision_task(
            domain_name, task_list, identity=identity
        )
        if decision:
            return json.dumps(decision.to_full_dict(reverse_order=reverse_order))
        else:
            return json.dumps({"previousStartedEventId": 0, "startedEventId": 0})

    def count_pending_decision_tasks(self) -> str:
        domain_name = self._params["domain"]
        task_list = self._params["taskList"]["name"]
        self._check_string(domain_name)
        self._check_string(task_list)
        count = self.swf_backend.count_pending_decision_tasks(domain_name, task_list)
        return json.dumps({"count": count, "truncated": False})

    def respond_decision_task_completed(self) -> str:
        task_token = self._params["taskToken"]
        execution_context = self._params.get("executionContext")
        decisions = self._params.get("decisions")
        self._check_string(task_token)
        self._check_none_or_string(execution_context)
        self.swf_backend.respond_decision_task_completed(
            task_token, decisions=decisions, execution_context=execution_context
        )
        return ""

    def poll_for_activity_task(self) -> str:
        domain_name = self._params["domain"]
        task_list = self._params["taskList"]["name"]
        identity = self._params.get("identity")
        self._check_string(domain_name)
        self._check_string(task_list)
        self._check_none_or_string(identity)
        activity_task = self.swf_backend.poll_for_activity_task(
            domain_name, task_list, identity=identity
        )
        if activity_task:
            return json.dumps(activity_task.to_full_dict())
        else:
            return json.dumps({"startedEventId": 0})

    def count_pending_activity_tasks(self) -> str:
        domain_name = self._params["domain"]
        task_list = self._params["taskList"]["name"]
        self._check_string(domain_name)
        self._check_string(task_list)
        count = self.swf_backend.count_pending_activity_tasks(domain_name, task_list)
        return json.dumps({"count": count, "truncated": False})

    def respond_activity_task_completed(self) -> str:
        task_token = self._params["taskToken"]
        result = self._params.get("result")
        self._check_string(task_token)
        self._check_none_or_string(result)
        self.swf_backend.respond_activity_task_completed(task_token, result=result)
        return ""

    def respond_activity_task_failed(self) -> str:
        task_token = self._params["taskToken"]
        reason = self._params.get("reason")
        details = self._params.get("details")
        self._check_string(task_token)
        # TODO: implement length limits on reason and details (common pb with
        # client libs)
        self._check_none_or_string(reason)
        self._check_none_or_string(details)
        self.swf_backend.respond_activity_task_failed(
            task_token, reason=reason, details=details
        )
        return ""

    def terminate_workflow_execution(self) -> str:
        domain_name = self._params["domain"]
        workflow_id = self._params["workflowId"]
        child_policy = self._params.get("childPolicy")
        details = self._params.get("details")
        reason = self._params.get("reason")
        run_id = self._params.get("runId")
        self._check_string(domain_name)
        self._check_string(workflow_id)
        self._check_none_or_string(child_policy)
        self._check_none_or_string(details)
        self._check_none_or_string(reason)
        self._check_none_or_string(run_id)
        self.swf_backend.terminate_workflow_execution(
            domain_name,
            workflow_id,
            child_policy=child_policy,
            details=details,
            reason=reason,
            run_id=run_id,
        )
        return ""

    def record_activity_task_heartbeat(self) -> str:
        task_token = self._params["taskToken"]
        details = self._params.get("details")
        self._check_string(task_token)
        self._check_none_or_string(details)
        self.swf_backend.record_activity_task_heartbeat(task_token, details=details)
        # TODO: make it dynamic when we implement activity tasks cancellation
        return json.dumps({"cancelRequested": False})

    def signal_workflow_execution(self) -> str:
        domain_name = self._params["domain"]
        signal_name = self._params["signalName"]
        workflow_id = self._params["workflowId"]
        _input = self._params["input"]
        run_id = self._params.get("runId")

        self._check_string(domain_name)
        self._check_string(signal_name)
        self._check_string(workflow_id)
        self._check_none_or_string(_input)
        self._check_none_or_string(run_id)

        self.swf_backend.signal_workflow_execution(
            domain_name, signal_name, workflow_id, _input, run_id
        )
        return ""
