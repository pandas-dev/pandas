from collections import defaultdict
from typing import TYPE_CHECKING, Any, Dict, List, Optional

from moto.core.common_models import BaseModel
from moto.utilities.utils import get_partition

from ..exceptions import (
    SWFUnknownResourceFault,
    SWFWorkflowExecutionAlreadyStartedFault,
)

if TYPE_CHECKING:
    from .activity_task import ActivityTask
    from .decision_task import DecisionTask
    from .generic_type import GenericType, TGenericType
    from .workflow_execution import WorkflowExecution


class Domain(BaseModel):
    def __init__(
        self,
        name: str,
        retention: int,
        account_id: str,
        region_name: str,
        description: Optional[str] = None,
    ):
        self.name = name
        self.retention = retention
        self.account_id = account_id
        self.region_name = region_name
        self.description = description
        self.status = "REGISTERED"
        self.types: Dict[str, Dict[str, Dict[str, GenericType]]] = {
            "activity": defaultdict(dict),
            "workflow": defaultdict(dict),
        }
        # Workflow executions have an id, which unicity is guaranteed
        # at domain level (not super clear in the docs, but I checked
        # that against SWF API) ; hence the storage method as a dict
        # of "workflow_id (client determined)" => WorkflowExecution()
        # here.
        self.workflow_executions: List["WorkflowExecution"] = []
        self.activity_task_lists: Dict[List[str], List["ActivityTask"]] = {}
        self.decision_task_lists: Dict[str, List["DecisionTask"]] = {}

    def __repr__(self) -> str:
        return f"Domain(name: {self.name}, status: {self.status})"

    def to_short_dict(self) -> Dict[str, str]:
        hsh = {"name": self.name, "status": self.status}
        if self.description:
            hsh["description"] = self.description
        hsh["arn"] = (
            f"arn:{get_partition(self.region_name)}:swf:{self.region_name}:{self.account_id}:/domain/{self.name}"
        )
        return hsh

    def to_full_dict(self) -> Dict[str, Any]:
        return {
            "domainInfo": self.to_short_dict(),
            "configuration": {"workflowExecutionRetentionPeriodInDays": self.retention},
        }

    def get_type(  # type: ignore
        self, kind: str, name: str, version: str, ignore_empty: bool = False
    ) -> "GenericType":
        try:
            return self.types[kind][name][version]
        except KeyError:
            if not ignore_empty:
                raise SWFUnknownResourceFault(
                    "type",
                    f"{kind.capitalize()}Type=[name={name}, version={version}]",
                )

    def add_type(self, _type: "TGenericType") -> None:
        self.types[_type.kind][_type.name][_type.version] = _type

    def find_types(self, kind: str, status: str) -> List["GenericType"]:
        _all = []
        for family in self.types[kind].values():
            for _type in family.values():
                if _type.status == status:
                    _all.append(_type)
        return _all

    def add_workflow_execution(self, workflow_execution: "WorkflowExecution") -> None:
        _id = workflow_execution.workflow_id
        if self.get_workflow_execution(_id, raise_if_none=False):
            raise SWFWorkflowExecutionAlreadyStartedFault()
        self.workflow_executions.append(workflow_execution)

    def get_workflow_execution(
        self,
        workflow_id: str,
        run_id: Optional[str] = None,
        raise_if_none: bool = True,
        raise_if_closed: bool = False,
    ) -> Optional["WorkflowExecution"]:
        # query
        if run_id:
            _all = [
                w
                for w in self.workflow_executions
                if w.workflow_id == workflow_id and w.run_id == run_id
            ]
        else:
            _all = [
                w
                for w in self.workflow_executions
                if w.workflow_id == workflow_id and w.open
            ]
        # reduce
        wfe = _all[0] if _all else None
        # raise if closed / none
        if raise_if_closed and wfe and wfe.execution_status == "CLOSED":
            wfe = None
        if not wfe and raise_if_none:
            if run_id:
                args = [
                    "execution",
                    f"WorkflowExecution=[workflowId={workflow_id}, runId={run_id}]",
                ]
            else:
                args = [f"execution, workflowId = {workflow_id}"]
            raise SWFUnknownResourceFault(*args)
        # at last return workflow execution
        return wfe

    def add_to_activity_task_list(
        self, task_list: List[str], obj: "ActivityTask"
    ) -> None:
        if task_list not in self.activity_task_lists:
            self.activity_task_lists[task_list] = []
        self.activity_task_lists[task_list].append(obj)

    @property
    def activity_tasks(self) -> List["ActivityTask"]:
        _all: List["ActivityTask"] = []
        for tasks in self.activity_task_lists.values():
            _all += tasks
        return _all

    def add_to_decision_task_list(self, task_list: str, obj: "DecisionTask") -> None:
        if task_list not in self.decision_task_lists:
            self.decision_task_lists[task_list] = []
        self.decision_task_lists[task_list].append(obj)

    @property
    def decision_tasks(self) -> List["DecisionTask"]:
        _all: List["DecisionTask"] = []
        for tasks in self.decision_task_lists.values():
            _all += tasks
        return _all
