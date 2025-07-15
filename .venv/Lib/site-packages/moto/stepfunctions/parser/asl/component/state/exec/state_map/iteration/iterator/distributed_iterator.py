from __future__ import annotations

from typing import Optional

from moto.stepfunctions.parser.asl.component.common.comment import Comment
from moto.stepfunctions.parser.asl.component.common.flow.start_at import StartAt
from moto.stepfunctions.parser.asl.component.common.query_language import QueryLanguage
from moto.stepfunctions.parser.asl.component.program.states import States
from moto.stepfunctions.parser.asl.component.state.exec.state_map.iteration.distributed_iteration_component import (
    DistributedIterationComponent,
    DistributedIterationComponentEvalInput,
)
from moto.stepfunctions.parser.asl.component.state.exec.state_map.iteration.itemprocessor.processor_config import (
    ProcessorConfig,
)
from moto.stepfunctions.parser.asl.component.state.exec.state_map.iteration.iterator.distributed_iterator_worker import (
    DistributedIteratorWorker,
)
from moto.stepfunctions.parser.asl.eval.environment import Environment
from moto.stepfunctions.parser.asl.parse.typed_props import TypedProps


class DistributedIteratorEvalInput(DistributedIterationComponentEvalInput):
    pass


class DistributedIterator(DistributedIterationComponent):
    _eval_input: Optional[DistributedIteratorEvalInput]

    @classmethod
    def from_props(cls, props: TypedProps) -> DistributedIterator:
        item_processor = cls(
            query_language=props.get(QueryLanguage) or QueryLanguage(),
            start_at=props.get(
                typ=StartAt,
                raise_on_missing=ValueError(
                    f"Missing StartAt declaration in props '{props}'."
                ),
            ),
            states=props.get(
                typ=States,
                raise_on_missing=ValueError(
                    f"Missing States declaration in props '{props}'."
                ),
            ),
            comment=props.get(Comment),
            processor_config=props.get(ProcessorConfig),
        )
        return item_processor

    def _create_worker(self, env: Environment) -> DistributedIteratorWorker:
        return DistributedIteratorWorker(
            work_name=self._eval_input.state_name,
            job_pool=self._job_pool,
            env=env,
            parameters=self._eval_input.parameters,
            map_run_record=self._map_run_record,
            item_selector=self._eval_input.item_selector,
        )
