from enum import Enum
from typing import Any, Optional, Union

from moto.dynamodb.exceptions import KeyIsEmptyStringException, MockValidationException
from moto.utilities.tokenizer import GenericTokenizer


class EXPRESSION_STAGES(Enum):
    INITIAL_STAGE = "INITIAL_STAGE"  # Can be a hash key, range key, or function
    KEY_NAME = "KEY_NAME"
    KEY_VALUE = "KEY_VALUE"
    COMPARISON = "COMPARISON"
    EOF = "EOF"


def get_keys(schema: list[dict[str, str]], key_type: str) -> list[str]:
    """Get all keys of the specified type in schema order."""
    return [key["AttributeName"] for key in schema if key["KeyType"] == key_type]


def parse_expression(
    key_condition_expression: str,
    expression_attribute_values: dict[str, dict[str, str]],
    expression_attribute_names: dict[str, str],
    schema: list[dict[str, str]],
) -> tuple[
    list[tuple[str, dict[str, Any]]],
    list[tuple[str, str, list[dict[str, Any]]]],
    list[str],
]:
    """
    Parse a KeyConditionExpression using the provided expression attribute names/values

    key_condition_expression:    hashkey = :id AND :sk = val
    expression_attribute_names:  {":sk": "sortkey"}
    expression_attribute_values: {":id": {"S": "some hash key"}}
    schema:                      [{'AttributeName': 'hashkey', 'KeyType': 'HASH'}, {"AttributeName": "sortkey", "KeyType": "RANGE"}]

    Returns:
        - hash_key_conditions: List of (attr_name, value) for all hash keys
        - range_key_conditions: List of (attr_name, comparison, values) for range keys
        - expression_attribute_names_used: List of expression attribute names used
    """

    current_stage: Optional[EXPRESSION_STAGES] = None
    current_phrase = ""
    key_name = comparison = ""
    key_values: list[Union[dict[str, str], str]] = []
    expression_attribute_names_used: list[str] = []
    results: list[tuple[str, str, Any]] = []
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
                    if expression_attribute_names.get(current_phrase):
                        key_name = expression_attribute_names[current_phrase]
                        expression_attribute_names_used.append(current_phrase)
                    else:
                        key_name = current_phrase
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
                if expression_attribute_names.get(current_phrase):
                    key_name = expression_attribute_names[current_phrase]
                    expression_attribute_names_used.append(current_phrase)
                else:
                    key_name = current_phrase
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
                if expression_attribute_names.get(current_phrase):
                    key_name = expression_attribute_names[current_phrase]
                    expression_attribute_names_used.append(current_phrase)
                else:
                    key_name = current_phrase
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
                if not current_phrase:
                    # hashkey = :id and (begins_with(sortkey, :sk))
                    #                   ^
                    continue
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

    hash_key_conditions, range_key_conditions = validate_schema(results, schema)

    return (
        hash_key_conditions,
        range_key_conditions,
        expression_attribute_names_used,
    )


# Validate that the schema-keys are encountered in our query
def validate_schema(
    results: Any, schema: list[dict[str, str]]
) -> tuple[
    list[tuple[str, dict[str, Any]]],
    list[tuple[str, str, list[dict[str, Any]]]],
]:
    """
    Validate query conditions against the schema and extract key values.

    For multi-attribute keys:
    - ALL hash keys must be present with equality (=)
    - Range keys must be specified left-to-right (cannot skip)
    - Only the LAST range key in the query can use non-equality operators
    - Earlier range keys in the query must use equality (=)

    Returns:
        - hash_key_conditions: List of (attr_name, value) for all hash keys
        - range_key_conditions: List of (attr_name, comparison, values) for range keys
    """
    # Build a lookup from results: {key_name: (comparison, values)}
    results_by_key: dict[str, tuple[str, list[dict[str, Any]]]] = {
        key: (comparison, values) for key, comparison, values in results
    }

    # Get all hash and range keys from schema in order
    hash_keys = get_keys(schema, "HASH")
    range_keys = get_keys(schema, "RANGE")
    provided_keys = set(results_by_key.keys())
    schema_keys = set(hash_keys + range_keys)

    # === Validate HASH keys ===
    # All hash keys must be present with equality
    # Check missing hash keys FIRST (before checking for invalid keys)
    hash_key_conditions: list[tuple[str, dict[str, Any]]] = []

    for hash_key in hash_keys:
        if hash_key not in results_by_key:
            raise MockValidationException(
                f"Query condition missed key schema element: {hash_key}"
            )
        comparison, values = results_by_key[hash_key]
        if comparison != "=":
            raise MockValidationException("Query key condition not supported")
        value = values[0]
        if "S" in value and value["S"] == "":
            raise KeyIsEmptyStringException(hash_key)  # type: ignore[arg-type]
        hash_key_conditions.append((hash_key, value))

    # === Validate RANGE keys ===
    # Range keys must be specified left-to-right, only last can use non-equality
    range_key_conditions: list[tuple[str, str, list[dict[str, Any]]]] = []

    # Find which range keys are in the query (in schema order)
    range_keys_in_query: list[str] = [k for k in range_keys if k in results_by_key]

    # Check for unknown keys (keys in query but not in schema)
    unknown_keys = provided_keys - schema_keys
    if unknown_keys:
        # User provided keys not in schema - report first missing range key if any,
        # otherwise report generic error
        missing_range_keys = [k for k in range_keys if k not in results_by_key]
        if missing_range_keys:
            raise MockValidationException(
                f"Query condition missed key schema element: {missing_range_keys[0]}"
            )
        raise MockValidationException("Query key condition not supported")

    for i, range_key in enumerate(range_keys):
        if range_key not in results_by_key:
            # This range key is not in the query
            # Check if any later range keys ARE in the query (would be skipping)
            later_keys_in_query = [
                k for k in range_keys[i + 1 :] if k in results_by_key
            ]
            if later_keys_in_query:
                # User skipped this key but specified a later one - error
                raise MockValidationException(
                    f"RANGE key attributes {range_key} must have equality conditions "
                    f"specified in the query because a condition is present on key "
                    f"attribute {later_keys_in_query[0]}"
                )
            # No more range keys in query, we're done
            break

        comparison, values = results_by_key[range_key]
        is_last_in_query = range_key == range_keys_in_query[-1]

        if is_last_in_query:
            # Last range key in query can use any comparison
            if {"S": ""} in values:
                raise KeyIsEmptyStringException(range_key)
            range_key_conditions.append((range_key, comparison.upper(), values))
        else:
            # Not last in query - must be equality
            if comparison != "=":
                raise MockValidationException("Query key condition not supported")
            value = values[0]
            if "S" in value and value["S"] == "":
                raise KeyIsEmptyStringException(range_key)  # type: ignore[arg-type]
            range_key_conditions.append((range_key, "=", [value]))

    return hash_key_conditions, range_key_conditions
