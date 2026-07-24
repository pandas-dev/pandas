from moto.stepfunctions.parser.asl.component.state.exec.state_map.iteration.iteration_component import (
    IterationComponent,
)
from moto.stepfunctions.parser.asl.component.state.exec.state_map.iteration.iterator.distributed_iterator import (
    DistributedIterator,
)
from moto.stepfunctions.parser.asl.component.state.exec.state_map.iteration.iterator.inline_iterator import (
    InlineIterator,
)
from moto.stepfunctions.parser.asl.component.state.exec.state_map.iteration.iterator.iterator_decl import (
    IteratorDecl,
)
from moto.stepfunctions.parser.asl.component.state.exec.state_map.mode import (
    Mode,
)


def from_iterator_decl(iterator_decl: IteratorDecl) -> IterationComponent:
    if iterator_decl.processor_config.mode == Mode.Inline:
        return InlineIterator(
            query_language=iterator_decl.query_language,
            start_at=iterator_decl.start_at,
            states=iterator_decl.states,
            comment=iterator_decl.comment,
            processor_config=iterator_decl.processor_config,
        )
    elif iterator_decl.processor_config.mode == Mode.Distributed:
        return DistributedIterator(
            query_language=iterator_decl.query_language,
            start_at=iterator_decl.start_at,
            states=iterator_decl.states,
            comment=iterator_decl.comment,
            processor_config=iterator_decl.processor_config,
        )
    else:
        raise ValueError(
            f"Unknown Map state processing mode: '{iterator_decl.processor_config.mode}'."
        )
