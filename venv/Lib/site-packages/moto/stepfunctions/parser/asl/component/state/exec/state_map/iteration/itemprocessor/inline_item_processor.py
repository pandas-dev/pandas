from __future__ import annotations

import logging

from moto.stepfunctions.parser.asl.component.common.comment import Comment
from moto.stepfunctions.parser.asl.component.common.flow.start_at import StartAt
from moto.stepfunctions.parser.asl.component.common.query_language import QueryLanguage
from moto.stepfunctions.parser.asl.component.program.states import States
from moto.stepfunctions.parser.asl.component.state.exec.state_map.iteration.inline_iteration_component import (
    InlineIterationComponent,
    InlineIterationComponentEvalInput,
)
from moto.stepfunctions.parser.asl.component.state.exec.state_map.iteration.itemprocessor.inline_item_processor_worker import (
    InlineItemProcessorWorker,
)
from moto.stepfunctions.parser.asl.component.state.exec.state_map.iteration.itemprocessor.processor_config import (
    ProcessorConfig,
)
from moto.stepfunctions.parser.asl.component.state.exec.state_map.iteration.job import (
    JobPool,
)
from moto.stepfunctions.parser.asl.eval.environment import Environment
from moto.stepfunctions.parser.asl.parse.typed_props import TypedProps

LOG = logging.getLogger(__name__)


class InlineItemProcessorEvalInput(InlineIterationComponentEvalInput):
    pass


class InlineItemProcessor(InlineIterationComponent):
    @classmethod
    def from_props(cls, props: TypedProps) -> InlineItemProcessor:
        if not props.get(States):
            raise ValueError(f"Missing States declaration in props '{props}'.")
        if not props.get(StartAt):
            raise ValueError(f"Missing StartAt declaration in props '{props}'.")
        item_processor = cls(
            query_language=props.get(QueryLanguage) or QueryLanguage(),
            start_at=props.get(StartAt),
            states=props.get(States),
            comment=props.get(Comment),
            processor_config=props.get(ProcessorConfig) or ProcessorConfig(),
        )
        return item_processor

    def _create_worker(
        self,
        env: Environment,
        eval_input: InlineItemProcessorEvalInput,
        job_pool: JobPool,
    ) -> InlineItemProcessorWorker:
        return InlineItemProcessorWorker(
            work_name=eval_input.state_name,
            job_pool=job_pool,
            env=env,
            item_selector=eval_input.item_selector,
            parameters=eval_input.parameters,
        )
