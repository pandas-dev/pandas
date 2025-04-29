from moto.stepfunctions.parser.asl.component.state.exec.state_map.item_reader.eval.resource_output_transformer.transformer import (
    ResourceOutputTransformer,
)
from moto.stepfunctions.parser.asl.component.state.exec.state_map.item_reader.eval.resource_output_transformer.transformer_csv import (
    ResourceOutputTransformerCSV,
)
from moto.stepfunctions.parser.asl.component.state.exec.state_map.item_reader.eval.resource_output_transformer.transformer_json import (
    ResourceOutputTransformerJson,
)
from moto.stepfunctions.parser.asl.component.state.exec.state_map.item_reader.reader_config.input_type import (
    InputType,
    InputTypeValue,
)


def resource_output_transformer_for(input_type: InputType) -> ResourceOutputTransformer:
    if input_type.input_type_value == InputTypeValue.CSV:
        return ResourceOutputTransformerCSV()
    elif input_type.input_type_value == InputTypeValue.JSON:
        return ResourceOutputTransformerJson()
    else:
        raise ValueError(f"Unknown InputType value: '{input_type.input_type_value}'.")
