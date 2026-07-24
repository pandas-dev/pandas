import base64
import json
import re
from datetime import datetime
from typing import Any


class DynamoJsonEncoder(json.JSONEncoder):
    def default(self, o: Any) -> Any:
        if isinstance(o, datetime):
            from moto.core.utils import unix_time

            return unix_time(o)
        elif hasattr(o, "to_json"):
            return o.to_json()
        elif isinstance(o, bytes):
            return base64.b64encode(o).decode("utf-8")


def dynamo_json_dump(dynamo_object: Any) -> str:
    return json.dumps(dynamo_object, cls=DynamoJsonEncoder)


def dynamo_to_dict(obj: Any) -> Any:
    """Recursively convert DynamoDB model objects to plain Python types.

    DynamoType objects are expanded to {type: value} dicts directly (without
    calling to_json()) so that binary values stay as raw bytes that the
    ActionResult serializer can base64-encode correctly.

    - DynamoType B/BS: base64-decoded to raw bytes for the serializer
    - DynamoType (other): expanded to {type: value} with recursive processing
    - Other objects with to_json(): recursively resolved via to_json()
    - datetime objects: preserved for ActionResult timestamp serialization
    - dicts and lists: recursively processed
    """
    # Lazy import to avoid circular dependency
    from moto.dynamodb.models.dynamo_type import DynamoType

    # Check datetime before to_json to prevent datetime subclasses
    # from being incorrectly routed through the to_json() path
    if isinstance(obj, datetime):
        return obj
    if isinstance(obj, DynamoType):
        if obj.type == "B":
            # Value may be raw bytes (after projection round-trip) or base64 str (from wire)
            if isinstance(obj.value, bytes):
                return {"B": obj.value}
            return {"B": base64.b64decode(obj.value)}
        if obj.type == "BS":
            return {
                "BS": [
                    v if isinstance(v, bytes) else base64.b64decode(v)
                    for v in obj.value
                ]
            }
        return {obj.type: dynamo_to_dict(obj.value)}
    if hasattr(obj, "to_json"):
        return dynamo_to_dict(obj.to_json())
    if isinstance(obj, dict):
        return {k: dynamo_to_dict(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [dynamo_to_dict(item) for item in obj]
    return obj


def bytesize(val: str) -> int:
    return len(val if isinstance(val, bytes) else val.encode("utf-8"))


def find_nested_key(
    keys: list[str],
    dct: dict[str, Any],
    processed_keys: list[str] | None = None,
    result: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """
    keys   : A list of keys that may be present in the provided dictionary
             ["level1", "level2"]
    dct    : A dictionary that we want to inspect
             {"level1": {"level2": "val", "irrelevant": ..}

    processed_keys:
        Should not be set by the caller, only by recursive invocations.
        Example value: ["level1"]
    result:
        Should not be set by the caller, only by recursive invocations
        Example value: {"level1": {}}

    returns: {"level1": {"level2": "val"}}
    """
    if result is None:
        result = {}
    if processed_keys is None:
        processed_keys = []

    # A key can refer to a list-item: 'level1[1].level2'
    is_list_expression = re.match(pattern=r"(.+)\[(\d+)\]$", string=keys[0])

    if len(keys) == 1:
        # Set 'current_key' and 'value'
        #   or return an empty dictionary if the key does not exist in our dictionary
        if is_list_expression:
            current_key = is_list_expression.group(1)
            idx = int(is_list_expression.group(2))
            if (
                current_key in dct
                and isinstance(dct[current_key], list)
                and len(dct[current_key]) > idx
            ):
                value = [dct[current_key][idx]]
            else:
                return {}
        elif keys[0] in dct:
            current_key = keys[0]
            value = dct[current_key]
        else:
            return {}

        # We may have already processed some keys
        # Dig into the result to find the appropriate key to append the value to
        #
        # result: {'level1': {'level2': {}}}
        # processed_keys: ['level1', 'level2']
        #     -->
        # result: {'level1': {'level2': value}}
        temp_result = result
        for key in processed_keys:
            if isinstance(temp_result, list):
                temp_result = temp_result[0][key]
            else:
                temp_result = temp_result[key]
        if isinstance(temp_result, list):
            temp_result.append({current_key: value})
        else:
            temp_result[current_key] = value
        return result
    else:
        # Set 'current_key'
        #   or return an empty dictionary if the key does not exist in our dictionary
        if is_list_expression:
            current_key = is_list_expression.group(1)
            idx = int(is_list_expression.group(2))
            if (
                current_key in dct
                and isinstance(dct[current_key], list)
                and len(dct[current_key]) > idx
            ):
                pass
            else:
                return {}
        elif keys[0] in dct:
            current_key = keys[0]
        else:
            return {}

        # Append the 'current_key' to the dictionary that is our result (so far)
        # {'level1': {}} --> {'level1': {current_key: {}}
        temp_result = result
        for key in processed_keys:
            temp_result = temp_result[key]
        if isinstance(temp_result, list):
            temp_result.append({current_key: [] if is_list_expression else {}})
        else:
            temp_result[current_key] = [] if is_list_expression else {}
        remaining_dct = (
            dct[current_key][idx] if is_list_expression else dct[current_key]
        )

        return find_nested_key(
            keys[1:],
            remaining_dct,
            processed_keys=processed_keys + [current_key],
            result=result,
        )
