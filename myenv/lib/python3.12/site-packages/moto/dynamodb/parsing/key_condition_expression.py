from enum import Enum
from typing import Any, Dict, List, Optional, Tuple, Union

from moto.dynamodb.exceptions import KeyIsEmptyStringException, MockValidationException
from moto.utilities.tokenizer import GenericTokenizer


class EXPRESSION_STAGES(Enum):
    INITIAL_STAGE = "INITIAL_STAGE"  # Can be a hash key, range key, or function
    KEY_NAME = "KEY_NAME"
    KEY_VALUE = "KEY_VALUE"
    COMPARISON = "COMPARISON"
    EOF = "EOF"


def get_key(schema: List[Dict[str, str]], key_type: str) -> Optional[str]:
    keys = [key for key in schema if key["KeyType"] == key_type]
    return keys[0]["AttributeName"] if keys else None


def parse_expression(
    key_condition_expression: str,
    expression_attribute_values: Dict[str, Dict[str, str]],
    expression_attribute_names: Dict[str, str],
    schema: List[Dict[str, str]],
) -> Tuple[Dict[str, Any], Optional[str], List[Dict[str, Any]]]:
    """
    Parse a KeyConditionExpression using the provided expression attribute names/values

    key_condition_expression:    hashkey = :id AND :sk = val
    expression_attribute_names:  {":sk": "sortkey"}
    expression_attribute_values: {":id": {"S": "some hash key"}}
    schema:                      [{'AttributeName': 'hashkey', 'KeyType': 'HASH'}, {"AttributeName": "sortkey", "KeyType": "RANGE"}]
    """

    current_stage: Optional[EXPRESSION_STAGES] = None
    current_phrase = ""
    key_name = comparison = ""
    key_values: List[Union[Dict[str, str], str]] = []
    results: List[Tuple[str, str, Any]] = []
    tokenizer = GenericTokenizer(key_condition_expression)
    for crnt_char in tokenizer:
        if crnt_char == " ":
            if current_stage == EXPRESSION_STAGES.INITIAL_STAGE:
                tokenizer.skip_white_space()
                if tokenizer.peek() == "(":
                    # begins_with(sk, :sk) and primary = :pk
                    #            ^
                    continue
                else:
                    # start_date < :sk and primary = :pk
                    #            ^
                    key_name = expression_attribute_names.get(
                        current_phrase, current_phrase
                    )
                    current_phrase = ""
                    current_stage = EXPRESSION_STAGES.COMPARISON
                    tokenizer.skip_white_space()
            elif current_stage == EXPRESSION_STAGES.KEY_VALUE:
                # job_id =          :id
                # job_id =          :id and  ...
                # pk=p and          x=y
                # pk=p and fn(x, y1, y1 )
                #                      ^ --> ^
                key_values.append(
                    expression_attribute_values.get(
                        current_phrase, {"S": current_phrase}
                    )
                )
                current_phrase = ""
                if comparison.upper() != "BETWEEN" or len(key_values) == 2:
                    results.append((key_name, comparison, key_values))
                    key_values = []
                tokenizer.skip_white_space()
                if tokenizer.peek() == ")":
                    tokenizer.skip_characters(")")
                    current_stage = EXPRESSION_STAGES.EOF
                    break
                elif tokenizer.is_eof():
                    break
                tokenizer.skip_characters("AND", case_sensitive=False)
                tokenizer.skip_white_space()
                if comparison.upper() == "BETWEEN":
                    # We can expect another key_value, i.e. BETWEEN x and y
                    # We should add some validation, to not allow BETWEEN x and y and z and ..
                    pass
                else:
                    current_stage = EXPRESSION_STAGES.INITIAL_STAGE
            elif current_stage == EXPRESSION_STAGES.COMPARISON:
                # hashkey = :id and sortkey       =      :sk
                # hashkey = :id and sortkey BETWEEN      x and y
                #                                  ^ --> ^
                comparison = current_phrase
                current_phrase = ""
                current_stage = EXPRESSION_STAGES.KEY_VALUE
            continue
        if crnt_char in ["=", "<", ">"] and current_stage in [
            EXPRESSION_STAGES.KEY_NAME,
            EXPRESSION_STAGES.INITIAL_STAGE,
            EXPRESSION_STAGES.COMPARISON,
        ]:
            if current_stage in [
                EXPRESSION_STAGES.KEY_NAME,
                EXPRESSION_STAGES.INITIAL_STAGE,
            ]:
                key_name = expression_attribute_names.get(
                    current_phrase, current_phrase
                )
            current_phrase = ""
            if crnt_char in ["<", ">"] and tokenizer.peek() == "=":
                comparison = crnt_char + tokenizer.__next__()
            else:
                comparison = crnt_char
            tokenizer.skip_white_space()
            current_stage = EXPRESSION_STAGES.KEY_VALUE
            continue
        if crnt_char in [","]:
            if current_stage == EXPRESSION_STAGES.KEY_NAME:
                # hashkey = :id and begins_with(sortkey,     :sk)
                #                                      ^ --> ^
                key_name = expression_attribute_names.get(
                    current_phrase, current_phrase
                )
                current_phrase = ""
                current_stage = EXPRESSION_STAGES.KEY_VALUE
                tokenizer.skip_white_space()
                continue
            else:
                raise MockValidationException(
                    f'Invalid KeyConditionExpression: Syntax error; token: "{current_phrase}"'
                )
        if crnt_char in [")"]:
            if current_stage == EXPRESSION_STAGES.KEY_VALUE:
                # hashkey = :id and begins_with(sortkey, :sk)
                #                                            ^
                value = expression_attribute_values.get(current_phrase, current_phrase)
                current_phrase = ""
                key_values.append(value)
                results.append((key_name, comparison, key_values))
                key_values = []
                tokenizer.skip_white_space()
                if tokenizer.is_eof() or tokenizer.peek() == ")":
                    break
                else:
                    tokenizer.skip_characters("AND", case_sensitive=False)
                    tokenizer.skip_white_space()
                    current_stage = EXPRESSION_STAGES.INITIAL_STAGE
                    continue
        if crnt_char in [""]:
            # hashkey =                   :id
            # hashkey = :id and sortkey = :sk
            #                                ^
            if current_stage == EXPRESSION_STAGES.KEY_VALUE:
                if current_phrase not in expression_attribute_values:
                    raise MockValidationException(
                        "Invalid condition in KeyConditionExpression: Multiple attribute names used in one condition"
                    )
                key_values.append(expression_attribute_values[current_phrase])
                results.append((key_name, comparison, key_values))
                break
        if crnt_char == "(":
            # hashkey = :id and begins_with(      sortkey,     :sk)
            #                              ^ --> ^
            # (hash_key = :id) and (sortkey = :sk)
            #                     ^
            if current_stage in [EXPRESSION_STAGES.INITIAL_STAGE]:
                if current_phrase not in ["begins_with", ""]:
                    raise MockValidationException(
                        f"Invalid KeyConditionExpression: Invalid function name; function: {current_phrase}"
                    )
                comparison = current_phrase
                current_phrase = ""
                tokenizer.skip_white_space()
                current_stage = EXPRESSION_STAGES.KEY_NAME
                continue
            if current_stage is None:
                # (hash_key = :id .. )
                # ^
                continue

        current_phrase += crnt_char
        if current_stage is None:
            current_stage = EXPRESSION_STAGES.INITIAL_STAGE

    hash_value, range_comparison, range_values = validate_schema(results, schema)

    return (
        hash_value,
        range_comparison.upper() if range_comparison else None,
        range_values,
    )


# Validate that the schema-keys are encountered in our query
def validate_schema(
    results: Any, schema: List[Dict[str, str]]
) -> Tuple[Dict[str, Any], Optional[str], List[Dict[str, Any]]]:
    index_hash_key = get_key(schema, "HASH")
    comparison, hash_value = next(
        (
            (comparison, value[0])
            for key, comparison, value in results
            if key == index_hash_key
        ),
        (None, None),
    )
    if hash_value is None:
        raise MockValidationException(
            f"Query condition missed key schema element: {index_hash_key}"
        )
    if comparison != "=":
        raise MockValidationException("Query key condition not supported")
    if "S" in hash_value and hash_value["S"] == "":
        raise KeyIsEmptyStringException(index_hash_key)  # type: ignore[arg-type]

    index_range_key = get_key(schema, "RANGE")
    range_key, range_comparison, range_values = next(
        (
            (key, comparison, values)
            for key, comparison, values in results
            if key == index_range_key
        ),
        (None, None, []),
    )
    if index_range_key:
        if len(results) > 1 and range_key != index_range_key:
            raise MockValidationException(
                f"Query condition missed key schema element: {index_range_key}"
            )
        if {"S": ""} in range_values:
            raise KeyIsEmptyStringException(index_range_key)

    return hash_value, range_comparison, range_values
