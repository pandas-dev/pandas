import base64
import json
import os
from typing import Any, Dict, Optional

from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.common_models import BaseModel
from moto.dynamodb.models import DynamoDBBackend, dynamodb_backends
from moto.dynamodb.models.table import StreamShard, Table
from moto.dynamodb.models.utilities import DynamoJsonEncoder


class ShardIterator(BaseModel):
    def __init__(
        self,
        streams_backend: "DynamoDBStreamsBackend",
        stream_shard: StreamShard,
        shard_iterator_type: str,
        sequence_number: Optional[int] = None,
    ):
        self.id = base64.b64encode(os.urandom(472)).decode("utf-8")
        self.streams_backend = streams_backend
        self.stream_shard = stream_shard
        self.shard_iterator_type = shard_iterator_type
        if shard_iterator_type == "TRIM_HORIZON":
            self.sequence_number = stream_shard.starting_sequence_number
        elif shard_iterator_type == "LATEST":
            self.sequence_number = stream_shard.starting_sequence_number + len(
                stream_shard.items
            )
        elif shard_iterator_type == "AT_SEQUENCE_NUMBER":
            self.sequence_number = sequence_number  # type: ignore[assignment]
        elif shard_iterator_type == "AFTER_SEQUENCE_NUMBER":
            self.sequence_number = sequence_number + 1  # type: ignore[operator]

    @property
    def arn(self) -> str:
        return f"{self.stream_shard.table.table_arn}/stream/{self.stream_shard.table.latest_stream_label}|1|{self.id}"

    def to_json(self) -> Dict[str, str]:
        return {"ShardIterator": self.arn}

    def get(self, limit: int = 1000) -> Dict[str, Any]:
        items = self.stream_shard.get(self.sequence_number, limit)
        try:
            last_sequence_number = max(
                int(i["dynamodb"]["SequenceNumber"]) for i in items
            )
            new_shard_iterator = ShardIterator(
                self.streams_backend,
                self.stream_shard,
                "AFTER_SEQUENCE_NUMBER",
                last_sequence_number,
            )
        except ValueError:
            new_shard_iterator = ShardIterator(
                self.streams_backend,
                self.stream_shard,
                "AT_SEQUENCE_NUMBER",
                self.sequence_number,
            )

        self.streams_backend.shard_iterators[new_shard_iterator.arn] = (
            new_shard_iterator
        )
        return {"NextShardIterator": new_shard_iterator.arn, "Records": items}


class DynamoDBStreamsBackend(BaseBackend):
    def __init__(self, region_name: str, account_id: str):
        super().__init__(region_name, account_id)
        self.shard_iterators: Dict[str, ShardIterator] = {}

    @property
    def dynamodb(self) -> DynamoDBBackend:
        return dynamodb_backends[self.account_id][self.region_name]

    def _get_table_from_arn(self, arn: str) -> Table:
        table_name = arn.split(":", 6)[5].split("/")[1]
        return self.dynamodb.get_table(table_name)

    def describe_stream(self, arn: str) -> str:
        table = self._get_table_from_arn(arn)
        resp = {
            "StreamDescription": {
                "StreamArn": arn,
                "StreamLabel": table.latest_stream_label,
                "StreamStatus": (
                    "ENABLED" if table.latest_stream_label else "DISABLED"
                ),
                "StreamViewType": table.stream_specification["StreamViewType"],  # type: ignore[index]
                "CreationRequestDateTime": table.stream_shard.created_on.isoformat(),  # type: ignore[union-attr]
                "TableName": table.name,
                "KeySchema": table.schema,
                "Shards": (
                    [table.stream_shard.to_json()] if table.stream_shard else []
                ),
            }
        }

        return json.dumps(resp)

    def list_streams(self, table_name: Optional[str] = None) -> str:
        streams = []
        for table in self.dynamodb.tables.values():
            if table_name is not None and table.name != table_name:
                continue
            if table.latest_stream_label:
                d = table.describe(base_key="Table")
                streams.append(
                    {
                        "StreamArn": d["Table"]["LatestStreamArn"],
                        "TableName": d["Table"]["TableName"],
                        "StreamLabel": d["Table"]["LatestStreamLabel"],
                    }
                )

        return json.dumps({"Streams": streams})

    def get_shard_iterator(
        self,
        arn: str,
        shard_id: str,
        shard_iterator_type: str,
        sequence_number: Optional[str] = None,
    ) -> str:
        table = self._get_table_from_arn(arn)
        assert table.stream_shard.id == shard_id  # type: ignore[union-attr]

        shard_iterator = ShardIterator(
            self,
            table.stream_shard,  # type: ignore[arg-type]
            shard_iterator_type,
            sequence_number,  # type: ignore[arg-type]
        )
        self.shard_iterators[shard_iterator.arn] = shard_iterator

        return json.dumps(shard_iterator.to_json())

    def get_records(self, iterator_arn: str, limit: int) -> str:
        shard_iterator = self.shard_iterators[iterator_arn]
        return json.dumps(shard_iterator.get(limit), cls=DynamoJsonEncoder)


dynamodbstreams_backends = BackendDict(DynamoDBStreamsBackend, "dynamodbstreams")
