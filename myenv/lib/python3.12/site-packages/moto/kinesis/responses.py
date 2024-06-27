import json

from moto.core.responses import BaseResponse

from .models import KinesisBackend, kinesis_backends


class KinesisResponse(BaseResponse):
    def __init__(self) -> None:
        super().__init__(service_name="kinesis")

    @property
    def kinesis_backend(self) -> KinesisBackend:
        return kinesis_backends[self.current_account][self.region]

    def create_stream(self) -> str:
        stream_name = self._get_param("StreamName")
        shard_count = self._get_param("ShardCount")
        stream_mode = self._get_param("StreamModeDetails")
        self.kinesis_backend.create_stream(
            stream_name, shard_count, stream_mode=stream_mode
        )
        return ""

    def describe_stream(self) -> str:
        stream_name = self._get_param("StreamName")
        stream_arn = self._get_param("StreamARN")
        limit = self._get_param("Limit")
        stream = self.kinesis_backend.describe_stream(stream_arn, stream_name)
        return json.dumps(stream.to_json(shard_limit=limit))

    def describe_stream_summary(self) -> str:
        stream_arn = self._get_param("StreamARN")
        stream_name = self._get_param("StreamName")
        stream = self.kinesis_backend.describe_stream_summary(stream_arn, stream_name)
        return json.dumps(stream.to_json_summary())

    def list_streams(self) -> str:
        streams = self.kinesis_backend.list_streams()
        stream_names = [stream.stream_name for stream in streams]
        max_streams = self._get_param("Limit", 10)
        try:
            token = self._get_param("ExclusiveStartStreamName")
        except ValueError:
            token = self._get_param("ExclusiveStartStreamName")
        if token:
            start = stream_names.index(token) + 1
        else:
            start = 0
        streams_resp = stream_names[start : start + max_streams]
        has_more_streams = False
        if start + max_streams < len(stream_names):
            has_more_streams = True

        return json.dumps(
            {
                "HasMoreStreams": has_more_streams,
                "StreamNames": streams_resp,
                "StreamSummaries": [
                    stream.to_json_summary()["StreamDescriptionSummary"]
                    for stream in streams
                ],
            }
        )

    def delete_stream(self) -> str:
        stream_arn = self._get_param("StreamARN")
        stream_name = self._get_param("StreamName")
        self.kinesis_backend.delete_stream(stream_arn, stream_name)
        return ""

    def get_shard_iterator(self) -> str:
        stream_arn = self._get_param("StreamARN")
        stream_name = self._get_param("StreamName")
        shard_id = self._get_param("ShardId")
        shard_iterator_type = self._get_param("ShardIteratorType")
        starting_sequence_number = self._get_param("StartingSequenceNumber")
        at_timestamp = self._get_param("Timestamp")

        shard_iterator = self.kinesis_backend.get_shard_iterator(
            stream_arn,
            stream_name,
            shard_id,
            shard_iterator_type,
            starting_sequence_number,
            at_timestamp,
        )

        return json.dumps({"ShardIterator": shard_iterator})

    def get_records(self) -> str:
        stream_arn = self._get_param("StreamARN")
        shard_iterator = self._get_param("ShardIterator")
        limit = self._get_param("Limit")

        (
            next_shard_iterator,
            records,
            millis_behind_latest,
        ) = self.kinesis_backend.get_records(stream_arn, shard_iterator, limit)

        return json.dumps(
            {
                "NextShardIterator": next_shard_iterator,
                "Records": [record.to_json() for record in records],
                "MillisBehindLatest": millis_behind_latest,
            }
        )

    def put_record(self) -> str:
        stream_arn = self._get_param("StreamARN")
        stream_name = self._get_param("StreamName")
        partition_key = self._get_param("PartitionKey")
        explicit_hash_key = self._get_param("ExplicitHashKey")
        data = self._get_param("Data")

        sequence_number, shard_id = self.kinesis_backend.put_record(
            stream_arn,
            stream_name,
            partition_key,
            explicit_hash_key,
            data,
        )

        return json.dumps({"SequenceNumber": sequence_number, "ShardId": shard_id})

    def put_records(self) -> str:
        stream_arn = self._get_param("StreamARN")
        stream_name = self._get_param("StreamName")
        records = self._get_param("Records")

        response = self.kinesis_backend.put_records(stream_arn, stream_name, records)

        return json.dumps(response)

    def split_shard(self) -> str:
        stream_arn = self._get_param("StreamARN")
        stream_name = self._get_param("StreamName")
        shard_to_split = self._get_param("ShardToSplit")
        new_starting_hash_key = self._get_param("NewStartingHashKey")
        self.kinesis_backend.split_shard(
            stream_arn, stream_name, shard_to_split, new_starting_hash_key
        )
        return ""

    def merge_shards(self) -> str:
        stream_arn = self._get_param("StreamARN")
        stream_name = self._get_param("StreamName")
        shard_to_merge = self._get_param("ShardToMerge")
        adjacent_shard_to_merge = self._get_param("AdjacentShardToMerge")
        self.kinesis_backend.merge_shards(
            stream_arn, stream_name, shard_to_merge, adjacent_shard_to_merge
        )
        return ""

    def list_shards(self) -> str:
        stream_arn = self._get_param("StreamARN")
        stream_name = self._get_param("StreamName")
        next_token = self._get_param("NextToken")
        max_results = self._get_param("MaxResults", 10000)
        shards, token = self.kinesis_backend.list_shards(
            stream_arn=stream_arn,
            stream_name=stream_name,
            limit=max_results,
            next_token=next_token,
        )
        res = {"Shards": shards, "NextToken": token}
        return json.dumps(res)

    def update_shard_count(self) -> str:
        stream_arn = self._get_param("StreamARN")
        stream_name = self._get_param("StreamName")
        target_shard_count = self._get_param("TargetShardCount")
        current_shard_count = self.kinesis_backend.update_shard_count(
            stream_arn=stream_arn,
            stream_name=stream_name,
            target_shard_count=target_shard_count,
        )
        return json.dumps(
            dict(
                StreamName=stream_name,
                CurrentShardCount=current_shard_count,
                TargetShardCount=target_shard_count,
            )
        )

    def increase_stream_retention_period(self) -> str:
        stream_arn = self._get_param("StreamARN")
        stream_name = self._get_param("StreamName")
        retention_period_hours = self._get_param("RetentionPeriodHours")
        self.kinesis_backend.increase_stream_retention_period(
            stream_arn, stream_name, retention_period_hours
        )
        return ""

    def decrease_stream_retention_period(self) -> str:
        stream_arn = self._get_param("StreamARN")
        stream_name = self._get_param("StreamName")
        retention_period_hours = self._get_param("RetentionPeriodHours")
        self.kinesis_backend.decrease_stream_retention_period(
            stream_arn, stream_name, retention_period_hours
        )
        return ""

    def add_tags_to_stream(self) -> str:
        stream_arn = self._get_param("StreamARN")
        stream_name = self._get_param("StreamName")
        tags = self._get_param("Tags")
        self.kinesis_backend.add_tags_to_stream(stream_arn, stream_name, tags)
        return json.dumps({})

    def list_tags_for_stream(self) -> str:
        stream_arn = self._get_param("StreamARN")
        stream_name = self._get_param("StreamName")
        exclusive_start_tag_key = self._get_param("ExclusiveStartTagKey")
        limit = self._get_param("Limit")
        response = self.kinesis_backend.list_tags_for_stream(
            stream_arn, stream_name, exclusive_start_tag_key, limit
        )
        return json.dumps(response)

    def remove_tags_from_stream(self) -> str:
        stream_arn = self._get_param("StreamARN")
        stream_name = self._get_param("StreamName")
        tag_keys = self._get_param("TagKeys")
        self.kinesis_backend.remove_tags_from_stream(stream_arn, stream_name, tag_keys)
        return json.dumps({})

    def enable_enhanced_monitoring(self) -> str:
        stream_arn = self._get_param("StreamARN")
        stream_name = self._get_param("StreamName")
        shard_level_metrics = self._get_param("ShardLevelMetrics")
        arn, name, current, desired = self.kinesis_backend.enable_enhanced_monitoring(
            stream_arn=stream_arn,
            stream_name=stream_name,
            shard_level_metrics=shard_level_metrics,
        )
        return json.dumps(
            dict(
                StreamName=name,
                CurrentShardLevelMetrics=current,
                DesiredShardLevelMetrics=desired,
                StreamARN=arn,
            )
        )

    def disable_enhanced_monitoring(self) -> str:
        stream_arn = self._get_param("StreamARN")
        stream_name = self._get_param("StreamName")
        shard_level_metrics = self._get_param("ShardLevelMetrics")
        arn, name, current, desired = self.kinesis_backend.disable_enhanced_monitoring(
            stream_arn=stream_arn,
            stream_name=stream_name,
            to_be_disabled=shard_level_metrics,
        )
        return json.dumps(
            dict(
                StreamName=name,
                CurrentShardLevelMetrics=current,
                DesiredShardLevelMetrics=desired,
                StreamARN=arn,
            )
        )

    def list_stream_consumers(self) -> str:
        stream_arn = self._get_param("StreamARN")
        consumers = self.kinesis_backend.list_stream_consumers(stream_arn=stream_arn)
        return json.dumps(dict(Consumers=[c.to_json() for c in consumers]))

    def register_stream_consumer(self) -> str:
        stream_arn = self._get_param("StreamARN")
        consumer_name = self._get_param("ConsumerName")
        consumer = self.kinesis_backend.register_stream_consumer(
            stream_arn=stream_arn, consumer_name=consumer_name
        )
        return json.dumps(dict(Consumer=consumer.to_json()))

    def describe_stream_consumer(self) -> str:
        stream_arn = self._get_param("StreamARN")
        consumer_name = self._get_param("ConsumerName")
        consumer_arn = self._get_param("ConsumerARN")
        consumer = self.kinesis_backend.describe_stream_consumer(
            stream_arn=stream_arn,
            consumer_name=consumer_name,
            consumer_arn=consumer_arn,
        )
        return json.dumps(
            dict(ConsumerDescription=consumer.to_json(include_stream_arn=True))
        )

    def deregister_stream_consumer(self) -> str:
        stream_arn = self._get_param("StreamARN")
        consumer_name = self._get_param("ConsumerName")
        consumer_arn = self._get_param("ConsumerARN")
        self.kinesis_backend.deregister_stream_consumer(
            stream_arn=stream_arn,
            consumer_name=consumer_name,
            consumer_arn=consumer_arn,
        )
        return json.dumps(dict())

    def start_stream_encryption(self) -> str:
        stream_arn = self._get_param("StreamARN")
        stream_name = self._get_param("StreamName")
        encryption_type = self._get_param("EncryptionType")
        key_id = self._get_param("KeyId")
        self.kinesis_backend.start_stream_encryption(
            stream_arn=stream_arn,
            stream_name=stream_name,
            encryption_type=encryption_type,
            key_id=key_id,
        )
        return json.dumps(dict())

    def stop_stream_encryption(self) -> str:
        stream_arn = self._get_param("StreamARN")
        stream_name = self._get_param("StreamName")
        self.kinesis_backend.stop_stream_encryption(
            stream_arn=stream_arn, stream_name=stream_name
        )
        return json.dumps(dict())

    def update_stream_mode(self) -> str:
        stream_arn = self._get_param("StreamARN")
        stream_mode = self._get_param("StreamModeDetails")
        self.kinesis_backend.update_stream_mode(stream_arn, stream_mode)
        return "{}"
