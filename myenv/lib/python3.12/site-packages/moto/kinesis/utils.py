import base64
from datetime import datetime
from typing import Any, List, Optional

from .exceptions import InvalidArgumentError

encode_method = base64.encodebytes
decode_method = base64.decodebytes


PAGINATION_MODEL = {
    "list_shards": {
        "input_token": "next_token",
        "limit_key": "limit",
        "limit_default": 10000,
        "unique_attribute": "ShardId",
        "fail_on_invalid_token": False,
    },
}


PAGINATION_MODEL = {
    "list_shards": {
        "input_token": "next_token",
        "limit_key": "limit",
        "limit_default": 10000,
        "unique_attribute": "ShardId",
        "fail_on_invalid_token": False,
    },
}


def compose_new_shard_iterator(
    stream_name: Optional[str],
    shard: Any,
    shard_iterator_type: str,
    starting_sequence_number: int,
    at_timestamp: datetime,
) -> str:
    if shard_iterator_type == "AT_SEQUENCE_NUMBER":
        last_sequence_id = int(starting_sequence_number) - 1
    elif shard_iterator_type == "AFTER_SEQUENCE_NUMBER":
        last_sequence_id = int(starting_sequence_number)
    elif shard_iterator_type == "TRIM_HORIZON":
        last_sequence_id = 0
    elif shard_iterator_type == "LATEST":
        last_sequence_id = shard.get_max_sequence_number()
    elif shard_iterator_type == "AT_TIMESTAMP":
        last_sequence_id = shard.get_sequence_number_at(at_timestamp)
    else:
        raise InvalidArgumentError(f"Invalid ShardIteratorType: {shard_iterator_type}")
    return compose_shard_iterator(stream_name, shard, last_sequence_id)


def compose_shard_iterator(
    stream_name: Optional[str], shard: Any, last_sequence_id: int
) -> str:
    return encode_method(
        f"{stream_name}:{shard.shard_id}:{last_sequence_id}".encode("utf-8")
    ).decode("utf-8")


def decompose_shard_iterator(shard_iterator: str) -> List[str]:
    return decode_method(shard_iterator.encode("utf-8")).decode("utf-8").split(":")
