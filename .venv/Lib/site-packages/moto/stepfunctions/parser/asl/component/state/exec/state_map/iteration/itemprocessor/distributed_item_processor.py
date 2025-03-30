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
from moto.stepfunctions.parser.asl.component.state.exec.state_map.iteration.itemprocessor.distributed_item_processor_worker import (
    DistributedItemProcessorWorker,
)
from moto.stepfunctions.parser.asl.component.state.exec.state_map.iteration.itemprocessor.processor_config import (
    ProcessorConfig,
)
from moto.stepfunctions.parser.asl.eval.environment import Environment
from moto.stepfunctions.parser.asl.parse.typed_props import TypedProps


class DistributedItemProcessorEvalInput(DistributedIterationComponentEvalInput):
    pass


class DistributedItemProcessor(DistributedIterationComponent):
    _eval_input: Optional[DistributedItemProcessorEvalInput]

    @classmethod
    def from_props(cls, props: TypedProps) -> DistributedItemProcessor:
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
            processor_config=props.get(ProcessorConfig) or ProcessorConfig(),
        )
        return item_processor

    def _create_worker(self, env: Environment) -> DistributedItemProcessorWorker:
        return DistributedItemProcessorWorker(
            work_name=self._eval_input.state_name,
            job_pool=self._job_pool,
            env=env,
            item_reader=self._eval_input.item_reader,
            parameters=self._eval_input.parameters,
            item_selector=self._eval_input.item_selector,
            map_run_record=self._map_run_record,
        )
