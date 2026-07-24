from moto.core.responses import ActionResult, BaseResponse

from .models import DynamoDBStreamsBackend, dynamodbstreams_backends


class DynamoDBStreamsHandler(BaseResponse):
    def __init__(self) -> None:
        super().__init__(service_name="dynamodbstreams")

    @property
    def backend(self) -> DynamoDBStreamsBackend:
        return dynamodbstreams_backends[self.current_account][self.region]

    def describe_stream(self) -> ActionResult:
        arn = self._get_param("StreamArn")
        stream = self.backend.describe_stream(arn)
        return ActionResult({"StreamDescription": stream})

    def list_streams(self) -> ActionResult:
        table_name = self._get_param("TableName")
        streams = self.backend.list_streams(table_name)
        return ActionResult({"Streams": streams})

    def get_shard_iterator(self) -> ActionResult:
        arn = self._get_param("StreamArn")
        shard_id = self._get_param("ShardId")
        shard_iterator_type = self._get_param("ShardIteratorType")
        sequence_number = self._get_param("SequenceNumber")
        # according to documentation sequence_number param should be string
        if isinstance(sequence_number, str):
            sequence_number = int(sequence_number)

        iterator = self.backend.get_shard_iterator(
            arn, shard_id, shard_iterator_type, sequence_number
        )
        return ActionResult({"ShardIterator": iterator.arn})

    def get_records(self) -> ActionResult:
        arn = self._get_param("ShardIterator")
        limit = self._get_param("Limit")
        if limit is None:
            limit = 1000
        result = self.backend.get_records(arn, limit)
        return ActionResult(result)
