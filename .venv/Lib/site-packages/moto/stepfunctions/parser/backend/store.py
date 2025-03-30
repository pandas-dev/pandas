from collections import OrderedDict
from typing import Dict

from localstack.services.stores import AccountRegionBundle, BaseStore, LocalAttribute

from moto.stepfunctions.parser.api import Arn
from moto.stepfunctions.parser.backend.activity import Activity
from moto.stepfunctions.parser.backend.execution import Execution
from moto.stepfunctions.parser.backend.state_machine import StateMachineInstance


class SFNStore(BaseStore):
    # Maps ARNs to state machines.
    state_machines: Dict[Arn, StateMachineInstance] = LocalAttribute(default=dict)
    # Maps Execution-ARNs to state machines.
    executions: Dict[Arn, Execution] = LocalAttribute(
        default=OrderedDict
    )  # TODO: when snapshot to pods stop execution(?)
    activities: OrderedDict[Arn, Activity] = LocalAttribute(default=dict)


sfn_stores: AccountRegionBundle = AccountRegionBundle("stepfunctions", SFNStore)
