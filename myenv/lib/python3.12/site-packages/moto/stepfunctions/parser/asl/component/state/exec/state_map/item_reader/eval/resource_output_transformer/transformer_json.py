import json

from moto.stepfunctions.parser.api import HistoryEventType, MapRunFailedEventDetails
from moto.stepfunctions.parser.asl.component.common.error_name.failure_event import (
    FailureEvent,
    FailureEventException,
)
from moto.stepfunctions.parser.asl.component.common.error_name.states_error_name import (
    StatesErrorName,
)
from moto.stepfunctions.parser.asl.component.common.error_name.states_error_name_type import (
    StatesErrorNameType,
)
from moto.stepfunctions.parser.asl.component.state.exec.state_map.item_reader.eval.resource_output_transformer.transformer import (
    ResourceOutputTransformer,
)
from moto.stepfunctions.parser.asl.component.state.exec.state_map.item_reader.reader_config.reader_config_decl import (
    ReaderConfigOutput,
)
from moto.stepfunctions.parser.asl.eval.environment import Environment
from moto.stepfunctions.parser.asl.eval.event.event_detail import EventDetails


class ResourceOutputTransformerJson(ResourceOutputTransformer):
    def _eval_body(self, env: Environment) -> None:
        _: ReaderConfigOutput = (
            env.stack.pop()
        )  # Not used, but expected by workflow (hence should consume the stack).
        resource_value: str = env.stack.pop()

        json_list = json.loads(resource_value)

        if not isinstance(json_list, list):
            error_name = StatesErrorName(typ=StatesErrorNameType.StatesItemReaderFailed)
            failure_event = FailureEvent(
                error_name=error_name,
                event_type=HistoryEventType.TaskFailed,
                event_details=EventDetails(
                    mapRunFailedEventDetails=MapRunFailedEventDetails(
                        error=error_name.error_name,
                        cause="Attempting to map over non-iterable node.",
                    )
                ),
            )
            raise FailureEventException(failure_event=failure_event)

        env.stack.append(json_list)
