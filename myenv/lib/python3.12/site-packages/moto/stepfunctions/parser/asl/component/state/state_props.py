from typing import Any

from moto.stepfunctions.parser.asl.component.common.flow.end import End
from moto.stepfunctions.parser.asl.component.common.flow.next import Next
from moto.stepfunctions.parser.asl.component.common.timeouts.heartbeat import Heartbeat
from moto.stepfunctions.parser.asl.component.common.timeouts.timeout import Timeout
from moto.stepfunctions.parser.asl.component.state.exec.state_map.item_reader.reader_config.max_items_decl import (
    MaxItemsDecl,
)
from moto.stepfunctions.parser.asl.component.state.exec.state_task.service.resource import (
    Resource,
)
from moto.stepfunctions.parser.asl.component.state.wait.wait_function.wait_function import (
    WaitFunction,
)
from moto.stepfunctions.parser.asl.parse.typed_props import TypedProps


class StateProps(TypedProps):
    _UNIQUE_SUBINSTANCES = {
        Resource,
        WaitFunction,
        Timeout,
        Heartbeat,
        MaxItemsDecl,
    }
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
        for typ in self._UNIQUE_SUBINSTANCES:
            if issubclass(inst_type, typ):
                super()._add(typ, instance)
                return

        # Base and delegate to preprocessor.
        super().add(instance)
