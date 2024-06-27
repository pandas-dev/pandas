from moto.core.responses import BaseResponse

from .models import DynamoDBStreamsBackend, dynamodbstreams_backends


class DynamoDBStreamsHandler(BaseResponse):
    def __init__(self) -> None:
        super().__init__(service_name="dynamodb-streams")

    @property
    def backend(self) -> DynamoDBStreamsBackend:
        return dynamodbstreams_backends[self.current_account][self.region]

    def describe_stream(self) -> str:
        arn = self._get_param("StreamArn")
        return self.backend.describe_stream(arn)

    def list_streams(self) -> str:
        table_name = self._get_param("TableName")
        return self.backend.list_streams(table_name)

    def get_shard_iterator(self) -> str:
        arn = self._get_param("StreamArn")
        shard_id = self._get_param("ShardId")
        shard_iterator_type = self._get_param("ShardIteratorType")
        sequence_number = self._get_param("SequenceNumber")
        # according to documentation sequence_number param should be string
        if isinstance(sequence_number, str):
            sequence_number = int(sequence_number)

        return self.backend.get_shard_iterator(
            arn, shard_id, shard_iterator_type, sequence_number
        )

    def get_records(self) -> str:
        arn = self._get_param("ShardIterator")
        limit = self._get_param("Limit")
        if limit is None:
            limit = 1000
        return self.backend.get_records(arn, limit)
