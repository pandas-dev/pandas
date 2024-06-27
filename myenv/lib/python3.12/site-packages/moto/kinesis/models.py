import datetime
import io
import itertools
import json
import re
from base64 import b64decode, b64encode
from collections import OrderedDict
from gzip import GzipFile
from operator import attrgetter
from typing import TYPE_CHECKING, Any, Dict, Iterable, List, Optional, Tuple

from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.common_models import BaseModel, CloudFormationModel
from moto.core.utils import unix_time, utcnow
from moto.moto_api._internal import mock_random as random
from moto.utilities.paginator import paginate
from moto.utilities.utils import get_partition, md5_hash

from .exceptions import (
    ConsumerNotFound,
    InvalidArgumentError,
    InvalidDecreaseRetention,
    InvalidIncreaseRetention,
    InvalidRetentionPeriod,
    RecordSizeExceedsLimit,
    ResourceInUseError,
    ResourceNotFoundError,
    ShardNotFoundError,
    StreamCannotBeUpdatedError,
    StreamNotFoundError,
    TooManyRecords,
    TotalRecordsSizeExceedsLimit,
    ValidationException,
)
from .utils import (
    PAGINATION_MODEL,
    compose_new_shard_iterator,
    compose_shard_iterator,
    decompose_shard_iterator,
)

if TYPE_CHECKING:
    from moto.awslambda.models import EventSourceMapping


class Consumer(BaseModel):
    def __init__(
        self, consumer_name: str, account_id: str, region_name: str, stream_arn: str
    ):
        self.consumer_name = consumer_name
        self.created = unix_time()
        self.stream_arn = stream_arn
        stream_name = stream_arn.split("/")[-1]
        self.consumer_arn = f"arn:{get_partition(region_name)}:kinesis:{region_name}:{account_id}:stream/{stream_name}/consumer/{consumer_name}"

    def to_json(self, include_stream_arn: bool = False) -> Dict[str, Any]:
        resp = {
            "ConsumerName": self.consumer_name,
            "ConsumerARN": self.consumer_arn,
            "ConsumerStatus": "ACTIVE",
            "ConsumerCreationTimestamp": self.created,
        }
        if include_stream_arn:
            resp["StreamARN"] = self.stream_arn
        return resp


class Record(BaseModel):
    def __init__(
        self,
        partition_key: str,
        data: str,
        sequence_number: int,
        explicit_hash_key: str,
    ):
        self.partition_key = partition_key
        self.data = data
        self.sequence_number = sequence_number
        self.explicit_hash_key = explicit_hash_key
        self.created_at_datetime = utcnow()
        self.created_at = unix_time(self.created_at_datetime)

    def to_json(self) -> Dict[str, Any]:
        return {
            "Data": self.data,
            "PartitionKey": self.partition_key,
            "SequenceNumber": str(self.sequence_number),
            "ApproximateArrivalTimestamp": self.created_at,
        }


class Shard(BaseModel):
    def __init__(
        self,
        shard_id: int,
        starting_hash: int,
        ending_hash: int,
        parent: Optional[str] = None,
        adjacent_parent: Optional[str] = None,
    ):
        self._shard_id = shard_id
        self.starting_hash = starting_hash
        self.ending_hash = ending_hash
        self.records: Dict[int, Record] = OrderedDict()
        self.is_open = True
        self.parent = parent
        self.adjacent_parent = adjacent_parent

    @property
    def shard_id(self) -> str:
        return f"shardId-{str(self._shard_id).zfill(12)}"

    def get_records(
        self, last_sequence_id: str, limit: Optional[int]
    ) -> Tuple[List[Record], int, int]:
        last_sequence_int = int(last_sequence_id)
        results = []
        secs_behind_latest = 0.0

        for sequence_number, record in self.records.items():
            if sequence_number > last_sequence_int:
                results.append(record)
                last_sequence_int = sequence_number

                very_last_record = self.records[next(reversed(self.records))]
                secs_behind_latest = very_last_record.created_at - record.created_at

            if len(results) == limit:
                break

        millis_behind_latest = int(secs_behind_latest * 1000)
        return results, last_sequence_int, millis_behind_latest

    def put_record(self, partition_key: str, data: str, explicit_hash_key: str) -> str:
        # Note: this function is not safe for concurrency
        if self.records:
            last_sequence_number = self.get_max_sequence_number()
        else:
            last_sequence_number = 0
        sequence_number = last_sequence_number + 1
        self.records[sequence_number] = Record(
            partition_key, data, sequence_number, explicit_hash_key
        )
        return str(sequence_number)

    def get_min_sequence_number(self) -> int:
        if self.records:
            return list(self.records.keys())[0]
        return 0

    def get_max_sequence_number(self) -> int:
        if self.records:
            return list(self.records.keys())[-1]
        return 0

    def get_sequence_number_at(self, at_timestamp: float) -> int:
        if not self.records or at_timestamp < list(self.records.values())[0].created_at:
            return 0
        else:
            # find the last item in the list that was created before
            # at_timestamp
            r = next(
                (
                    r
                    for r in reversed(self.records.values())
                    if r.created_at < at_timestamp
                ),
                None,
            )
            return r.sequence_number  # type: ignore

    def to_json(self) -> Dict[str, Any]:
        response: Dict[str, Any] = {
            "HashKeyRange": {
                "EndingHashKey": str(self.ending_hash),
                "StartingHashKey": str(self.starting_hash),
            },
            "SequenceNumberRange": {
                "StartingSequenceNumber": str(self.get_min_sequence_number()),
            },
            "ShardId": self.shard_id,
        }
        if self.parent:
            response["ParentShardId"] = self.parent
        if self.adjacent_parent:
            response["AdjacentParentShardId"] = self.adjacent_parent
        if not self.is_open:
            response["SequenceNumberRange"]["EndingSequenceNumber"] = str(
                self.get_max_sequence_number()
            )
        return response


class Stream(CloudFormationModel):
    def __init__(
        self,
        stream_name: str,
        shard_count: int,
        stream_mode: Optional[Dict[str, str]],
        retention_period_hours: Optional[int],
        account_id: str,
        region_name: str,
    ):
        self.stream_name = stream_name
        self.creation_datetime = datetime.datetime.now().strftime(
            "%Y-%m-%dT%H:%M:%S.%f000"
        )
        self.region = region_name
        self.account_id = account_id
        self.arn = f"arn:{get_partition(region_name)}:kinesis:{region_name}:{account_id}:stream/{stream_name}"
        self.shards: Dict[str, Shard] = {}
        self.tags: Dict[str, str] = {}
        self.status = "ACTIVE"
        self.shard_count: Optional[int] = None
        self.stream_mode = stream_mode or {"StreamMode": "PROVISIONED"}
        if self.stream_mode.get("StreamMode", "") == "ON_DEMAND":
            shard_count = 4
        self.init_shards(shard_count)
        self.retention_period_hours = retention_period_hours or 24
        self.shard_level_metrics: List[str] = []
        self.encryption_type = "NONE"
        self.key_id: Optional[str] = None
        self.consumers: List[Consumer] = []
        self.lambda_event_source_mappings: Dict[str, "EventSourceMapping"] = {}

    def delete_consumer(self, consumer_arn: str) -> None:
        self.consumers = [c for c in self.consumers if c.consumer_arn != consumer_arn]

    def get_consumer_by_arn(self, consumer_arn: str) -> Optional[Consumer]:
        return next((c for c in self.consumers if c.consumer_arn == consumer_arn), None)

    def init_shards(self, shard_count: int) -> None:
        self.shard_count = shard_count

        step = 2**128 // shard_count
        hash_ranges = itertools.chain(
            map(lambda i: (i, i * step, (i + 1) * step - 1), range(shard_count - 1)),
            [(shard_count - 1, (shard_count - 1) * step, 2**128)],
        )
        for index, start, end in hash_ranges:
            shard = Shard(index, start, end)
            self.shards[shard.shard_id] = shard

    def split_shard(self, shard_to_split: str, new_starting_hash_key: str) -> None:
        new_starting_hash_int = int(new_starting_hash_key)

        shard = self.shards[shard_to_split]

        if shard.starting_hash < new_starting_hash_int < shard.ending_hash:
            pass
        else:
            raise InvalidArgumentError(
                message=f"NewStartingHashKey {new_starting_hash_int} used in SplitShard() on shard {shard_to_split} in stream {self.stream_name} under account {self.account_id} is not both greater than one plus the shard's StartingHashKey {shard.starting_hash} and less than the shard's EndingHashKey {(shard.ending_hash - 1)}."
            )

        if not shard.is_open:
            raise InvalidArgumentError(
                message=f"Shard {shard.shard_id} in stream {self.stream_name} under account {self.account_id} has already been merged or split, and thus is not eligible for merging or splitting."
            )

        last_id = sorted(self.shards.values(), key=attrgetter("_shard_id"))[
            -1
        ]._shard_id

        # Create two new shards
        new_shard_1 = Shard(
            last_id + 1,
            starting_hash=shard.starting_hash,
            ending_hash=new_starting_hash_int - 1,
            parent=shard.shard_id,
        )
        new_shard_2 = Shard(
            last_id + 2,
            starting_hash=new_starting_hash_int,
            ending_hash=shard.ending_hash,
            parent=shard.shard_id,
        )
        self.shards[new_shard_1.shard_id] = new_shard_1
        self.shards[new_shard_2.shard_id] = new_shard_2
        shard.is_open = False

        records = shard.records
        shard.records = OrderedDict()

        for index in records:
            record = records[index]
            self.put_record(record.partition_key, record.explicit_hash_key, record.data)

    def merge_shards(self, shard_to_merge: str, adjacent_shard_to_merge: str) -> None:
        shard1 = self.shards[shard_to_merge]
        shard2 = self.shards[adjacent_shard_to_merge]

        # Validate the two shards are adjacent
        if shard1.ending_hash == (shard2.starting_hash - 1):
            pass
        elif shard2.ending_hash == (shard1.starting_hash + 1):
            pass
        else:
            raise InvalidArgumentError(adjacent_shard_to_merge)

        # Create a new shard
        last_id = sorted(self.shards.values(), key=attrgetter("_shard_id"))[
            -1
        ]._shard_id
        new_shard = Shard(
            last_id + 1,
            starting_hash=shard1.starting_hash,
            ending_hash=shard2.ending_hash,
            parent=shard1.shard_id,
            adjacent_parent=shard2.shard_id,
        )
        self.shards[new_shard.shard_id] = new_shard

        # Close the merged shards
        shard1.is_open = False
        shard2.is_open = False

        # Move all data across
        for record in shard1.records.values():
            new_shard.put_record(
                record.partition_key, record.data, record.explicit_hash_key
            )
        for record in shard2.records.values():
            new_shard.put_record(
                record.partition_key, record.data, record.explicit_hash_key
            )

    def update_shard_count(self, target_shard_count: int) -> None:
        if self.stream_mode.get("StreamMode", "") == "ON_DEMAND":
            raise StreamCannotBeUpdatedError(
                stream_name=self.stream_name, account_id=self.account_id
            )
        current_shard_count = len([s for s in self.shards.values() if s.is_open])
        if current_shard_count == target_shard_count:
            return

        # Split shards until we have enough shards
        # AWS seems to split until we have (current * 2) shards, and then merge until we reach the target
        # That's what observable at least - the actual algorithm is probably more advanced
        #
        if current_shard_count < target_shard_count:
            open_shards = [
                (shard_id, shard)
                for shard_id, shard in self.shards.items()
                if shard.is_open
            ]
            for shard_id, shard in open_shards:
                # Split the current shard
                new_starting_hash_key = str(
                    int((shard.ending_hash + shard.starting_hash) / 2)
                )
                self.split_shard(shard_id, new_starting_hash_key)

        current_shard_count = len([s for s in self.shards.values() if s.is_open])

        # If we need to reduce the shard count, merge shards until we get there
        while current_shard_count > target_shard_count:
            # Keep track of how often we need to merge to get to the target shard count
            required_shard_merges = current_shard_count - target_shard_count
            # Get a list of pairs of adjacent shards
            shard_list = sorted(
                [s for s in self.shards.values() if s.is_open],
                key=lambda x: x.starting_hash,
            )
            adjacent_shards = zip(
                [s for s in shard_list[0:-1:2]], [s for s in shard_list[1::2]]
            )

            for shard, adjacent in adjacent_shards:
                self.merge_shards(shard.shard_id, adjacent.shard_id)
                required_shard_merges -= 1
                if required_shard_merges == 0:
                    break

            current_shard_count = len([s for s in self.shards.values() if s.is_open])

        self.shard_count = target_shard_count

    def get_shard(self, shard_id: str) -> Shard:
        if shard_id in self.shards:
            return self.shards[shard_id]
        else:
            raise ShardNotFoundError(shard_id, stream="", account_id=self.account_id)

    def get_shard_for_key(
        self, partition_key: str, explicit_hash_key: str
    ) -> Optional[Shard]:
        if not isinstance(partition_key, str):
            raise InvalidArgumentError("partition_key")
        if len(partition_key) > 256:
            raise InvalidArgumentError("partition_key")

        if explicit_hash_key:
            if not isinstance(explicit_hash_key, str):
                raise InvalidArgumentError("explicit_hash_key")

            int_key = int(explicit_hash_key)

            if int_key >= 2**128:
                raise InvalidArgumentError("explicit_hash_key")

        else:
            int_key = int(md5_hash(partition_key.encode("utf-8")).hexdigest(), 16)

        for shard in self.shards.values():
            if shard.starting_hash <= int_key < shard.ending_hash:
                return shard
        return None

    def put_record(
        self, partition_key: str, explicit_hash_key: str, data: str
    ) -> Tuple[str, str]:
        shard: Shard = self.get_shard_for_key(partition_key, explicit_hash_key)  # type: ignore

        sequence_number = shard.put_record(partition_key, data, explicit_hash_key)

        from moto.awslambda.utils import get_backend

        for arn, esm in self.lambda_event_source_mappings.items():
            region = arn.split(":")[3]

            get_backend(self.account_id, region).send_kinesis_message(
                function_name=esm.function_arn,
                kinesis_stream=self.arn,
                kinesis_data=data,
                kinesis_shard_id=shard.shard_id,
                kinesis_partition_key=partition_key,
                kinesis_sequence_number=sequence_number,
            )

        return sequence_number, shard.shard_id

    def to_json(self, shard_limit: Optional[int] = None) -> Dict[str, Any]:
        all_shards = list(self.shards.values())
        requested_shards = all_shards[0 : shard_limit or len(all_shards)]
        return {
            "StreamDescription": {
                "StreamARN": self.arn,
                "StreamName": self.stream_name,
                "StreamCreationTimestamp": self.creation_datetime,
                "StreamStatus": self.status,
                "HasMoreShards": len(requested_shards) != len(all_shards),
                "RetentionPeriodHours": self.retention_period_hours,
                "EnhancedMonitoring": [{"ShardLevelMetrics": self.shard_level_metrics}],
                "EncryptionType": self.encryption_type,
                "KeyId": self.key_id,
                "Shards": [shard.to_json() for shard in requested_shards],
            }
        }

    def to_json_summary(self) -> Dict[str, Any]:
        return {
            "StreamDescriptionSummary": {
                "StreamARN": self.arn,
                "StreamName": self.stream_name,
                "StreamStatus": self.status,
                "StreamModeDetails": self.stream_mode,
                "RetentionPeriodHours": self.retention_period_hours,
                "StreamCreationTimestamp": self.creation_datetime,
                "EnhancedMonitoring": [{"ShardLevelMetrics": self.shard_level_metrics}],
                "OpenShardCount": self.shard_count,
                "EncryptionType": self.encryption_type,
                "KeyId": self.key_id,
            }
        }

    @staticmethod
    def cloudformation_name_type() -> str:
        return "Name"

    @staticmethod
    def cloudformation_type() -> str:
        # https://docs.aws.amazon.com/AWSCloudFormation/latest/UserGuide/aws-resource-kinesis-stream.html
        return "AWS::Kinesis::Stream"

    @classmethod
    def create_from_cloudformation_json(  # type: ignore[misc]
        cls,
        resource_name: str,
        cloudformation_json: Any,
        account_id: str,
        region_name: str,
        **kwargs: Any,
    ) -> "Stream":
        properties = cloudformation_json.get("Properties", {})
        shard_count = properties.get("ShardCount", 1)
        retention_period_hours = properties.get("RetentionPeriodHours", resource_name)
        tags = {
            tag_item["Key"]: tag_item["Value"]
            for tag_item in properties.get("Tags", [])
        }

        backend: KinesisBackend = kinesis_backends[account_id][region_name]
        stream = backend.create_stream(
            resource_name, shard_count, retention_period_hours=retention_period_hours
        )
        if any(tags):
            backend.add_tags_to_stream(
                stream_arn=None, stream_name=stream.stream_name, tags=tags
            )
        return stream

    @classmethod
    def update_from_cloudformation_json(  # type: ignore[misc]
        cls,
        original_resource: Any,
        new_resource_name: str,
        cloudformation_json: Any,
        account_id: str,
        region_name: str,
    ) -> "Stream":
        properties = cloudformation_json["Properties"]

        if Stream.is_replacement_update(properties):
            resource_name_property = cls.cloudformation_name_type()
            if resource_name_property not in properties:
                properties[resource_name_property] = new_resource_name
            new_resource = cls.create_from_cloudformation_json(
                resource_name=properties[resource_name_property],
                cloudformation_json=cloudformation_json,
                account_id=account_id,
                region_name=region_name,
            )
            properties[resource_name_property] = original_resource.name
            cls.delete_from_cloudformation_json(
                resource_name=original_resource.name,
                cloudformation_json=cloudformation_json,
                account_id=account_id,
                region_name=region_name,
            )
            return new_resource

        else:  # No Interruption
            if "ShardCount" in properties:
                original_resource.update_shard_count(properties["ShardCount"])
            if "RetentionPeriodHours" in properties:
                original_resource.retention_period_hours = properties[
                    "RetentionPeriodHours"
                ]
            if "Tags" in properties:
                original_resource.tags = {
                    tag_item["Key"]: tag_item["Value"]
                    for tag_item in properties.get("Tags", [])
                }
            return original_resource

    @classmethod
    def delete_from_cloudformation_json(  # type: ignore[misc]
        cls,
        resource_name: str,
        cloudformation_json: Any,
        account_id: str,
        region_name: str,
    ) -> None:
        backend: KinesisBackend = kinesis_backends[account_id][region_name]
        backend.delete_stream(stream_arn=None, stream_name=resource_name)

    @staticmethod
    def is_replacement_update(properties: List[str]) -> bool:
        properties_requiring_replacement_update = ["BucketName", "ObjectLockEnabled"]
        return any(
            [
                property_requiring_replacement in properties
                for property_requiring_replacement in properties_requiring_replacement_update
            ]
        )

    @classmethod
    def has_cfn_attr(cls, attr: str) -> bool:
        return attr in ["Arn"]

    def get_cfn_attribute(self, attribute_name: str) -> str:
        from moto.cloudformation.exceptions import UnformattedGetAttTemplateException

        if attribute_name == "Arn":
            return self.arn
        raise UnformattedGetAttTemplateException()

    @property
    def physical_resource_id(self) -> str:
        return self.stream_name


class KinesisBackend(BaseBackend):
    def __init__(self, region_name: str, account_id: str):
        super().__init__(region_name, account_id)
        self.streams: Dict[str, Stream] = OrderedDict()

    @staticmethod
    def default_vpc_endpoint_service(
        service_region: str, zones: List[str]
    ) -> List[Dict[str, str]]:
        """Default VPC endpoint service."""
        return BaseBackend.default_vpc_endpoint_service_factory(
            service_region, zones, "kinesis", special_service_name="kinesis-streams"
        )

    def create_stream(
        self,
        stream_name: str,
        shard_count: int,
        stream_mode: Optional[Dict[str, str]] = None,
        retention_period_hours: Optional[int] = None,
    ) -> Stream:
        if stream_name in self.streams:
            raise ResourceInUseError(stream_name)
        stream = Stream(
            stream_name,
            shard_count,
            stream_mode=stream_mode,
            retention_period_hours=retention_period_hours,
            account_id=self.account_id,
            region_name=self.region_name,
        )
        self.streams[stream_name] = stream
        return stream

    def describe_stream(
        self, stream_arn: Optional[str], stream_name: Optional[str]
    ) -> Stream:
        if stream_name and stream_name in self.streams:
            return self.streams[stream_name]
        if stream_arn:
            for stream in self.streams.values():
                if stream.arn == stream_arn:
                    return stream
        if stream_arn:
            stream_name = stream_arn.split("/")[1]
        raise StreamNotFoundError(stream_name, self.account_id)  # type: ignore

    def describe_stream_summary(
        self, stream_arn: Optional[str], stream_name: Optional[str]
    ) -> Stream:
        return self.describe_stream(stream_arn=stream_arn, stream_name=stream_name)

    def list_streams(self) -> Iterable[Stream]:
        return self.streams.values()

    def delete_stream(
        self, stream_arn: Optional[str], stream_name: Optional[str]
    ) -> Stream:
        stream = self.describe_stream(stream_arn=stream_arn, stream_name=stream_name)
        return self.streams.pop(stream.stream_name)

    def get_shard_iterator(
        self,
        stream_arn: Optional[str],
        stream_name: Optional[str],
        shard_id: str,
        shard_iterator_type: str,
        starting_sequence_number: int,
        at_timestamp: datetime.datetime,
    ) -> str:
        # Validate params
        stream = self.describe_stream(stream_arn=stream_arn, stream_name=stream_name)
        try:
            shard = stream.get_shard(shard_id)
        except ShardNotFoundError:
            raise ResourceNotFoundError(
                message=f"Shard {shard_id} in stream {stream.stream_name} under account {self.account_id} does not exist"
            )

        shard_iterator = compose_new_shard_iterator(
            stream_name,
            shard,
            shard_iterator_type,
            starting_sequence_number,
            at_timestamp,
        )
        return shard_iterator

    def get_records(
        self, stream_arn: Optional[str], shard_iterator: str, limit: Optional[int]
    ) -> Tuple[str, List[Record], int]:
        decomposed = decompose_shard_iterator(shard_iterator)
        stream_name, shard_id, last_sequence_id = decomposed

        stream = self.describe_stream(stream_arn=stream_arn, stream_name=stream_name)
        shard = stream.get_shard(shard_id)

        records, last_sequence_id, millis_behind_latest = shard.get_records(  # type: ignore
            last_sequence_id, limit
        )

        next_shard_iterator = compose_shard_iterator(
            stream_name,
            shard,
            last_sequence_id,  # type: ignore
        )

        return next_shard_iterator, records, millis_behind_latest

    def put_record(
        self,
        stream_arn: str,
        stream_name: str,
        partition_key: str,
        explicit_hash_key: str,
        data: str,
    ) -> Tuple[str, str]:
        stream = self.describe_stream(stream_arn=stream_arn, stream_name=stream_name)

        sequence_number, shard_id = stream.put_record(
            partition_key, explicit_hash_key, data
        )

        return sequence_number, shard_id

    def put_records(
        self, stream_arn: str, stream_name: str, records: List[Dict[str, Any]]
    ) -> Dict[str, Any]:
        stream = self.describe_stream(stream_arn=stream_arn, stream_name=stream_name)

        response: Dict[str, Any] = {"FailedRecordCount": 0, "Records": []}

        if len(records) > 500:
            raise TooManyRecords

        data_sizes = [
            len(b64decode(r.get("Data", ""))) + len(r.get("PartitionKey", ""))
            for r in records
        ]
        if sum(data_sizes) > 5242880:
            raise TotalRecordsSizeExceedsLimit
        idx_over_limit = next(
            (idx for idx, x in enumerate(data_sizes) if x > 1048576), None
        )
        if idx_over_limit is not None:
            raise RecordSizeExceedsLimit(position=idx_over_limit + 1)

        for record in records:
            partition_key = record.get("PartitionKey")
            explicit_hash_key = record.get("ExplicitHashKey")
            data = record.get("Data")

            sequence_number, shard_id = stream.put_record(
                partition_key,  # type: ignore[arg-type]
                explicit_hash_key,  # type: ignore[arg-type]
                data,  # type: ignore[arg-type]
            )
            response["Records"].append(
                {"SequenceNumber": sequence_number, "ShardId": shard_id}
            )

        return response

    def split_shard(
        self,
        stream_arn: str,
        stream_name: str,
        shard_to_split: str,
        new_starting_hash_key: str,
    ) -> None:
        stream = self.describe_stream(stream_arn=stream_arn, stream_name=stream_name)

        if not re.match("[a-zA-Z0-9_.-]+", shard_to_split):
            raise ValidationException(
                value=shard_to_split,
                position="shardToSplit",
                regex_to_match="[a-zA-Z0-9_.-]+",
            )

        if shard_to_split not in stream.shards:
            raise ShardNotFoundError(
                shard_id=shard_to_split, stream=stream_name, account_id=self.account_id
            )

        if not re.match(r"0|([1-9]\d{0,38})", new_starting_hash_key):
            raise ValidationException(
                value=new_starting_hash_key,
                position="newStartingHashKey",
                regex_to_match=r"0|([1-9]\d{0,38})",
            )

        stream.split_shard(shard_to_split, new_starting_hash_key)

    def merge_shards(
        self,
        stream_arn: str,
        stream_name: str,
        shard_to_merge: str,
        adjacent_shard_to_merge: str,
    ) -> None:
        stream = self.describe_stream(stream_arn=stream_arn, stream_name=stream_name)

        if shard_to_merge not in stream.shards:
            raise ShardNotFoundError(
                shard_to_merge, stream=stream.stream_name, account_id=self.account_id
            )

        if adjacent_shard_to_merge not in stream.shards:
            raise ShardNotFoundError(
                adjacent_shard_to_merge,
                stream=stream.stream_name,
                account_id=self.account_id,
            )

        stream.merge_shards(shard_to_merge, adjacent_shard_to_merge)

    def update_shard_count(
        self, stream_arn: str, stream_name: str, target_shard_count: int
    ) -> int:
        stream = self.describe_stream(stream_arn=stream_arn, stream_name=stream_name)
        current_shard_count = len([s for s in stream.shards.values() if s.is_open])

        stream.update_shard_count(target_shard_count)

        return current_shard_count

    @paginate(pagination_model=PAGINATION_MODEL)  # type: ignore[misc]
    def list_shards(
        self, stream_arn: Optional[str], stream_name: Optional[str]
    ) -> List[Dict[str, Any]]:
        stream = self.describe_stream(stream_arn=stream_arn, stream_name=stream_name)
        shards = sorted(stream.shards.values(), key=lambda x: x.shard_id)
        return [shard.to_json() for shard in shards]

    def increase_stream_retention_period(
        self,
        stream_arn: Optional[str],
        stream_name: Optional[str],
        retention_period_hours: int,
    ) -> None:
        stream = self.describe_stream(stream_arn=stream_arn, stream_name=stream_name)
        if retention_period_hours < 24:
            raise InvalidRetentionPeriod(retention_period_hours, too_short=True)
        if retention_period_hours > 8760:
            raise InvalidRetentionPeriod(retention_period_hours, too_short=False)
        if retention_period_hours < stream.retention_period_hours:
            raise InvalidIncreaseRetention(
                name=stream_name,
                requested=retention_period_hours,
                existing=stream.retention_period_hours,
            )
        stream.retention_period_hours = retention_period_hours

    def decrease_stream_retention_period(
        self,
        stream_arn: Optional[str],
        stream_name: Optional[str],
        retention_period_hours: int,
    ) -> None:
        stream = self.describe_stream(stream_arn=stream_arn, stream_name=stream_name)
        if retention_period_hours < 24:
            raise InvalidRetentionPeriod(retention_period_hours, too_short=True)
        if retention_period_hours > 8760:
            raise InvalidRetentionPeriod(retention_period_hours, too_short=False)
        if retention_period_hours > stream.retention_period_hours:
            raise InvalidDecreaseRetention(
                name=stream_name,
                requested=retention_period_hours,
                existing=stream.retention_period_hours,
            )
        stream.retention_period_hours = retention_period_hours

    def list_tags_for_stream(
        self,
        stream_arn: str,
        stream_name: str,
        exclusive_start_tag_key: Optional[str] = None,
        limit: Optional[int] = None,
    ) -> Dict[str, Any]:
        stream = self.describe_stream(stream_arn=stream_arn, stream_name=stream_name)

        tags: List[Dict[str, str]] = []
        result: Dict[str, Any] = {"HasMoreTags": False, "Tags": tags}
        for key, val in sorted(stream.tags.items(), key=lambda x: x[0]):
            if limit and len(tags) >= limit:
                result["HasMoreTags"] = True
                break
            if exclusive_start_tag_key and key < exclusive_start_tag_key:
                continue

            tags.append({"Key": key, "Value": val})

        return result

    def add_tags_to_stream(
        self,
        stream_arn: Optional[str],
        stream_name: Optional[str],
        tags: Dict[str, str],
    ) -> None:
        stream = self.describe_stream(stream_arn=stream_arn, stream_name=stream_name)
        stream.tags.update(tags)

    def remove_tags_from_stream(
        self, stream_arn: Optional[str], stream_name: Optional[str], tag_keys: List[str]
    ) -> None:
        stream = self.describe_stream(stream_arn=stream_arn, stream_name=stream_name)
        for key in tag_keys:
            if key in stream.tags:
                del stream.tags[key]

    def enable_enhanced_monitoring(
        self,
        stream_arn: Optional[str],
        stream_name: Optional[str],
        shard_level_metrics: List[str],
    ) -> Tuple[str, str, List[str], List[str]]:
        stream = self.describe_stream(stream_arn=stream_arn, stream_name=stream_name)
        current_shard_level_metrics = stream.shard_level_metrics
        desired_metrics = list(set(current_shard_level_metrics + shard_level_metrics))
        stream.shard_level_metrics = desired_metrics
        return (
            stream.arn,
            stream.stream_name,
            current_shard_level_metrics,
            desired_metrics,
        )

    def disable_enhanced_monitoring(
        self,
        stream_arn: Optional[str],
        stream_name: Optional[str],
        to_be_disabled: List[str],
    ) -> Tuple[str, str, List[str], List[str]]:
        stream = self.describe_stream(stream_arn=stream_arn, stream_name=stream_name)
        current_metrics = stream.shard_level_metrics
        if "ALL" in to_be_disabled:
            desired_metrics = []
        else:
            desired_metrics = [
                metric for metric in current_metrics if metric not in to_be_disabled
            ]
        stream.shard_level_metrics = desired_metrics
        return stream.arn, stream.stream_name, current_metrics, desired_metrics

    def _find_stream_by_arn(self, stream_arn: str) -> Stream:  # type: ignore[return]
        for stream in self.streams.values():
            if stream.arn == stream_arn:
                return stream

    def list_stream_consumers(self, stream_arn: str) -> List[Consumer]:
        """
        Pagination is not yet implemented
        """
        stream = self._find_stream_by_arn(stream_arn)
        return stream.consumers

    def register_stream_consumer(self, stream_arn: str, consumer_name: str) -> Consumer:
        consumer = Consumer(
            consumer_name, self.account_id, self.region_name, stream_arn
        )
        stream = self._find_stream_by_arn(stream_arn)
        stream.consumers.append(consumer)
        return consumer

    def describe_stream_consumer(
        self, stream_arn: str, consumer_name: str, consumer_arn: str
    ) -> Consumer:
        if stream_arn:
            stream = self._find_stream_by_arn(stream_arn)
            for consumer in stream.consumers:
                if consumer_name and consumer.consumer_name == consumer_name:
                    return consumer
        if consumer_arn:
            for stream in self.streams.values():
                _consumer = stream.get_consumer_by_arn(consumer_arn)
                if _consumer:
                    return _consumer
        raise ConsumerNotFound(
            consumer=consumer_name or consumer_arn, account_id=self.account_id
        )

    def deregister_stream_consumer(
        self, stream_arn: str, consumer_name: str, consumer_arn: str
    ) -> None:
        if stream_arn:
            stream = self._find_stream_by_arn(stream_arn)
            stream.consumers = [
                c for c in stream.consumers if c.consumer_name == consumer_name
            ]
        if consumer_arn:
            for stream in self.streams.values():
                # Only one stream will actually have this consumer
                # It will be a noop for other streams
                stream.delete_consumer(consumer_arn)

    def start_stream_encryption(
        self, stream_arn: str, stream_name: str, encryption_type: str, key_id: str
    ) -> None:
        stream = self.describe_stream(stream_arn=stream_arn, stream_name=stream_name)
        stream.encryption_type = encryption_type
        stream.key_id = key_id

    def stop_stream_encryption(self, stream_arn: str, stream_name: str) -> None:
        stream = self.describe_stream(stream_arn=stream_arn, stream_name=stream_name)
        stream.encryption_type = "NONE"
        stream.key_id = None

    def update_stream_mode(self, stream_arn: str, stream_mode: Dict[str, str]) -> None:
        stream = self._find_stream_by_arn(stream_arn)
        stream.stream_mode = stream_mode

    """Send log events to a Stream after encoding and gzipping it."""

    def send_log_event(
        self,
        delivery_stream_arn: str,
        filter_name: str,
        log_group_name: str,
        log_stream_name: str,
        log_events: List[Dict[str, Any]],
    ) -> None:
        data = {
            "logEvents": log_events,
            "logGroup": log_group_name,
            "logStream": log_stream_name,
            "messageType": "DATA_MESSAGE",
            "owner": self.account_id,
            "subscriptionFilters": [filter_name],
        }

        output = io.BytesIO()
        with GzipFile(fileobj=output, mode="w") as fhandle:
            fhandle.write(json.dumps(data, separators=(",", ":")).encode("utf-8"))
        gzipped_payload = b64encode(output.getvalue()).decode("UTF-8")

        stream = self.describe_stream(stream_arn=delivery_stream_arn, stream_name=None)
        random_partition_key = random.get_random_string(length=32, lower_case=True)
        stream.put_record(
            partition_key=random_partition_key,
            data=gzipped_payload,
            explicit_hash_key="",
        )


kinesis_backends = BackendDict(KinesisBackend, "kinesis")
