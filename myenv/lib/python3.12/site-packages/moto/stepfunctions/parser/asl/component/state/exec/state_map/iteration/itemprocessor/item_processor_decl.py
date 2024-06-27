from typing import Final, Optional

from moto.stepfunctions.parser.asl.component.common.comment import Comment
from moto.stepfunctions.parser.asl.component.common.flow.start_at import StartAt
from moto.stepfunctions.parser.asl.component.state.exec.state_map.iteration.itemprocessor.processor_config import (
    ProcessorConfig,
)
from moto.stepfunctions.parser.asl.component.state.exec.state_map.iteration.iteration_declaration import (
    IterationDecl,
)
from moto.stepfunctions.parser.asl.component.states import States


class ItemProcessorDecl(IterationDecl):
    processor_config: Final[ProcessorConfig]

    def __init__(
        self,
        comment: Optional[Comment],
        start_at: StartAt,
        states: States,
        processor_config: Optional[ProcessorConfig],
    ):
        super().__init__(start_at=start_at, comment=comment, states=states)
        self.processor_config = processor_config or ProcessorConfig()
