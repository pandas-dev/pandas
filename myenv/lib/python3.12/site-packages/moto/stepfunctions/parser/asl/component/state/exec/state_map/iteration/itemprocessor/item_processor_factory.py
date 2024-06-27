from moto.stepfunctions.parser.asl.component.state.exec.state_map.iteration.itemprocessor.distributed_item_processor import (
    DistributedItemProcessor,
)
from moto.stepfunctions.parser.asl.component.state.exec.state_map.iteration.itemprocessor.inline_item_processor import (
    InlineItemProcessor,
)
from moto.stepfunctions.parser.asl.component.state.exec.state_map.iteration.itemprocessor.item_processor_decl import (
    ItemProcessorDecl,
)
from moto.stepfunctions.parser.asl.component.state.exec.state_map.iteration.iteration_component import (
    IterationComponent,
)
from moto.stepfunctions.parser.asl.component.state.exec.state_map.mode import (
    Mode,
)


def from_item_processor_decl(
    item_processor_decl: ItemProcessorDecl,
) -> IterationComponent:
    if item_processor_decl.processor_config.mode == Mode.Inline:
        return InlineItemProcessor(
            start_at=item_processor_decl.start_at,
            states=item_processor_decl.states,
            comment=item_processor_decl.comment,
            processor_config=item_processor_decl.processor_config,
        )
    elif item_processor_decl.processor_config.mode == Mode.Distributed:
        return DistributedItemProcessor(
            start_at=item_processor_decl.start_at,
            states=item_processor_decl.states,
            comment=item_processor_decl.comment,
            processor_config=item_processor_decl.processor_config,
        )
    else:
        raise ValueError(
            f"Unknown Mode value '{item_processor_decl.processor_config.mode}'."
        )
