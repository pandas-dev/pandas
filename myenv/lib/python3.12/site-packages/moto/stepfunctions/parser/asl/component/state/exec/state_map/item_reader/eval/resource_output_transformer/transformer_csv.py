import csv
import io
from collections import OrderedDict

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
    CSVHeaderLocationOutput,
    ReaderConfigOutput,
)
from moto.stepfunctions.parser.asl.eval.environment import Environment
from moto.stepfunctions.parser.asl.eval.event.event_detail import EventDetails


class ResourceOutputTransformerCSV(ResourceOutputTransformer):
    def _eval_body(self, env: Environment) -> None:
        reader_config: ReaderConfigOutput = env.stack.pop()
        resource_value: str = env.stack.pop()

        csv_file = io.StringIO(resource_value)
        csv_reader = csv.reader(csv_file)

        location = reader_config["CSVHeaderLocation"]
        if location == CSVHeaderLocationOutput.FIRST_ROW:
            headers = next(csv_reader)
        elif location == CSVHeaderLocationOutput.GIVEN:
            headers = reader_config["CSVHeaders"]
        else:
            raise ValueError(f"Unknown CSVHeaderLocation value '{location}'.")

        if len(set(headers)) < len(headers):
            error_name = StatesErrorName(typ=StatesErrorNameType.StatesItemReaderFailed)
            failure_event = FailureEvent(
                error_name=error_name,
                event_type=HistoryEventType.TaskFailed,
                event_details=EventDetails(
                    mapRunFailedEventDetails=MapRunFailedEventDetails(
                        error=error_name.error_name,
                        cause="CSV headers cannot contain duplicates.",
                    )
                ),
            )
            raise FailureEventException(failure_event=failure_event)

        transformed_outputs = list()
        for row in csv_reader:
            transformed_output = dict()
            for i, header in enumerate(headers):
                transformed_output[header] = row[i] if i < len(row) else ""
            transformed_outputs.append(
                OrderedDict(
                    sorted(
                        transformed_output.items(),
                        key=lambda item: (item[0].isalpha(), item[0]),
                    )
                )
            )

        env.stack.append(transformed_outputs)
