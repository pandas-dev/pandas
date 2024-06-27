from moto.stepfunctions.parser.asl.component.intrinsic.argument.function_argument_list import (
    FunctionArgumentList,
)
from moto.stepfunctions.parser.asl.component.intrinsic.function.statesfunction.array import (
    array,
    array_contains,
    array_get_item,
    array_length,
    array_partition,
    array_range,
    array_unique,
)
from moto.stepfunctions.parser.asl.component.intrinsic.function.statesfunction.encoding_decoding import (
    base_64_decode,
    base_64_encode,
)
from moto.stepfunctions.parser.asl.component.intrinsic.function.statesfunction.generic import (
    string_format,
)
from moto.stepfunctions.parser.asl.component.intrinsic.function.statesfunction.hash_calculations import (
    hash_func,
)
from moto.stepfunctions.parser.asl.component.intrinsic.function.statesfunction.json_manipulation import (
    json_merge,
    json_to_string,
    string_to_json,
)
from moto.stepfunctions.parser.asl.component.intrinsic.function.statesfunction.math_operations import (
    math_add,
    math_random,
)
from moto.stepfunctions.parser.asl.component.intrinsic.function.statesfunction.states_function import (
    StatesFunction,
)
from moto.stepfunctions.parser.asl.component.intrinsic.function.statesfunction.string_operations import (
    string_split,
)
from moto.stepfunctions.parser.asl.component.intrinsic.function.statesfunction.unique_id_generation import (
    uuid,
)
from moto.stepfunctions.parser.asl.component.intrinsic.functionname.state_function_name_types import (
    StatesFunctionNameType,
)
from moto.stepfunctions.parser.asl.component.intrinsic.functionname.states_function_name import (
    StatesFunctionName,
)


# TODO: could use reflection on StatesFunctionNameType values.
class StatesFunctionFactory:
    @staticmethod
    def from_name(
        func_name: StatesFunctionName, arg_list: FunctionArgumentList
    ) -> StatesFunction:
        if func_name.function_type == StatesFunctionNameType.Array:
            return array.Array(arg_list=arg_list)
        if func_name.function_type == StatesFunctionNameType.ArrayPartition:
            return array_partition.ArrayPartition(arg_list=arg_list)
        if func_name.function_type == StatesFunctionNameType.ArrayContains:
            return array_contains.ArrayContains(arg_list=arg_list)
        if func_name.function_type == StatesFunctionNameType.ArrayRange:
            return array_range.ArrayRange(arg_list=arg_list)
        if func_name.function_type == StatesFunctionNameType.ArrayGetItem:
            return array_get_item.ArrayGetItem(arg_list=arg_list)
        if func_name.function_type == StatesFunctionNameType.ArrayLength:
            return array_length.ArrayLength(arg_list=arg_list)
        if func_name.function_type == StatesFunctionNameType.ArrayUnique:
            return array_unique.ArrayUnique(arg_list=arg_list)

        # JSON Manipulation
        if func_name.function_type == StatesFunctionNameType.JsonToString:
            return json_to_string.JsonToString(arg_list=arg_list)
        if func_name.function_type == StatesFunctionNameType.StringToJson:
            return string_to_json.StringToJson(arg_list=arg_list)
        if func_name.function_type == StatesFunctionNameType.JsonMerge:
            return json_merge.JsonMerge(arg_list=arg_list)

        # Unique Id Generation.
        if func_name.function_type == StatesFunctionNameType.UUID:
            return uuid.UUID(arg_list=arg_list)

        # String Operations.
        if func_name.function_type == StatesFunctionNameType.StringSplit:
            return string_split.StringSplit(arg_list=arg_list)

        # Hash Calculations.
        if func_name.function_type == StatesFunctionNameType.Hash:
            return hash_func.HashFunc(arg_list=arg_list)

        # Encoding and Decoding.
        if func_name.function_type == StatesFunctionNameType.Base64Encode:
            return base_64_encode.Base64Encode(arg_list=arg_list)
        if func_name.function_type == StatesFunctionNameType.Base64Decode:
            return base_64_decode.Base64Decode(arg_list=arg_list)

        # Math Operations.
        if func_name.function_type == StatesFunctionNameType.MathRandom:
            return math_random.MathRandom(arg_list=arg_list)
        if func_name.function_type == StatesFunctionNameType.MathAdd:
            return math_add.MathAdd(arg_list=arg_list)

        # Generic.
        if func_name.function_type == StatesFunctionNameType.Format:
            return string_format.StringFormat(arg_list=arg_list)

        # Unsupported.
        raise NotImplementedError(func_name.function_type)  # noqa
