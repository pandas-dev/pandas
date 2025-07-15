from typing import Any, Set

from moto.stepfunctions.parser.asl.component.common.flow.end import End
from moto.stepfunctions.parser.asl.component.common.flow.next import Next
from moto.stepfunctions.parser.asl.component.common.parargs import Parargs
from moto.stepfunctions.parser.asl.component.common.path.input_path import InputPath
from moto.stepfunctions.parser.asl.component.common.path.items_path import ItemsPath
from moto.stepfunctions.parser.asl.component.common.path.output_path import OutputPath
from moto.stepfunctions.parser.asl.component.common.timeouts.heartbeat import Heartbeat
from moto.stepfunctions.parser.asl.component.common.timeouts.timeout import Timeout
from moto.stepfunctions.parser.asl.component.state.choice.comparison.comparison_type import (
    Comparison,
)
from moto.stepfunctions.parser.asl.component.state.choice.comparison.variable import (
    Variable,
)
from moto.stepfunctions.parser.asl.component.state.exec.state_map.item_reader.reader_config.max_items_decl import (
    MaxItemsDecl,
)
from moto.stepfunctions.parser.asl.component.state.exec.state_map.items.items import (
    Items,
)
from moto.stepfunctions.parser.asl.component.state.exec.state_map.max_concurrency import (
    MaxConcurrencyDecl,
)
from moto.stepfunctions.parser.asl.component.state.exec.state_map.tolerated_failure import (
    ToleratedFailureCountDecl,
    ToleratedFailurePercentageDecl,
)
from moto.stepfunctions.parser.asl.component.state.exec.state_task.service.resource import (
    Resource,
)
from moto.stepfunctions.parser.asl.component.state.fail.cause_decl import CauseDecl
from moto.stepfunctions.parser.asl.component.state.fail.error_decl import ErrorDecl
from moto.stepfunctions.parser.asl.component.state.wait.wait_function.wait_function import (
    WaitFunction,
)
from moto.stepfunctions.parser.asl.parse.typed_props import TypedProps

UNIQUE_SUBINSTANCES: Set[type] = {
    InputPath,
    Items,
    ItemsPath,
    OutputPath,
    Resource,
    WaitFunction,
    Timeout,
    Heartbeat,
    MaxItemsDecl,
    MaxConcurrencyDecl,
    ToleratedFailureCountDecl,
    ToleratedFailurePercentageDecl,
    ErrorDecl,
    CauseDecl,
    Variable,
    Parargs,
    Comparison,
}


class StateProps(TypedProps):
    name: str

    def add(self, instance: Any) -> None:
        inst_type = type(instance)

        # End-Next conflicts:
        if inst_type == End and Next in self._instance_by_type:
            raise ValueError(
                f"End redefines Next, from '{self.get(Next)}' to '{instance}'."
            )
        if inst_type == Next and End in self._instance_by_type:
            raise ValueError(
                f"Next redefines End, from '{self.get(End)}' to '{instance}'."
            )

        # Subclasses
        for typ in UNIQUE_SUBINSTANCES:
            if issubclass(inst_type, typ):
                super()._add(typ, instance)
                return

        # Base and delegate to preprocessor.
        super().add(instance)
