import copy
from collections import defaultdict
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union

from moto.core.common_models import BaseModel, CloudFormationModel
from moto.core.utils import unix_time, unix_time_millis, utcnow
from moto.dynamodb.comparisons import get_expected, get_filter_expression
from moto.dynamodb.exceptions import (
    ConditionalCheckFailed,
    HashKeyTooLong,
    InvalidAttributeTypeError,
    InvalidConversion,
    InvalidIndexNameError,
    MockValidationException,
    RangeKeyTooLong,
    SerializationException,
)
from moto.dynamodb.limits import HASH_KEY_MAX_LENGTH, RANGE_KEY_MAX_LENGTH
from moto.dynamodb.models.dynamo_type import DynamoType, Item
from moto.dynamodb.models.utilities import dynamo_json_dump
from moto.moto_api._internal import mock_random
from moto.utilities.utils import get_partition

RESULT_SIZE_LIMIT = 1000000  # DynamoDB has a 1MB size limit


class SecondaryIndex(BaseModel):
    def __init__(
        self,
        index_name: str,
        schema: List[Dict[str, str]],
        projection: Dict[str, Any],
        table_key_attrs: List[str],
    ):
        self.name = index_name
        self.schema = schema
        self.table_key_attrs = table_key_attrs
        self.projection = projection
        self.schema_key_attrs = [k["AttributeName"] for k in schema]

    def project(self, item: Item) -> Item:
        """
        Enforces the ProjectionType of this Index (LSI/GSI)
        Removes any non-wanted attributes from the item
        :param item:
        :return:
        """
        if self.projection:
            projection_type = self.projection.get("ProjectionType", None)
            key_attributes = self.table_key_attrs + [
                key["AttributeName"] for key in self.schema
            ]

            if projection_type == "KEYS_ONLY":
                # 'project' expects lists of lists of strings
                # project([["attr1"], ["nested", "attr2"]]
                #
                # In our case, we need to convert
                # ["key1", "key2"]
                # into
                # [["key1"], ["key2"]]
                item = item.project([[attr] for attr in key_attributes])
            elif projection_type == "INCLUDE":
                allowed_attributes = key_attributes
                allowed_attributes.extend(self.projection.get("NonKeyAttributes", []))
                item = item.project([[attr] for attr in allowed_attributes])
            # ALL is handled implicitly by not filtering
        return item


class LocalSecondaryIndex(SecondaryIndex):
    def describe(self) -> Dict[str, Any]:
        return {
            "IndexName": self.name,
            "KeySchema": self.schema,
            "Projection": self.projection,
        }

    @staticmethod
    def create(  # type: ignore[misc]
        dct: Dict[str, Any], table_key_attrs: List[str]
    ) -> "LocalSecondaryIndex":
        return LocalSecondaryIndex(
            index_name=dct["IndexName"],
            schema=dct["KeySchema"],
            projection=dct["Projection"],
            table_key_attrs=table_key_attrs,
        )


class GlobalSecondaryIndex(SecondaryIndex):
    def __init__(
        self,
        index_name: str,
        schema: List[Dict[str, str]],
        projection: Dict[str, Any],
        table_key_attrs: List[str],
        status: str = "ACTIVE",
        throughput: Optional[Dict[str, Any]] = None,
    ):
        super().__init__(index_name, schema, projection, table_key_attrs)
        self.status = status
        self.throughput = throughput or {
            "ReadCapacityUnits": 0,
            "WriteCapacityUnits": 0,
        }

    def describe(self) -> Dict[str, Any]:
        return {
            "IndexName": self.name,
            "KeySchema": self.schema,
            "Projection": self.projection,
            "IndexStatus": self.status,
            "ProvisionedThroughput": self.throughput,
        }

    @staticmethod
    def create(  # type: ignore[misc]
        dct: Dict[str, Any], table_key_attrs: List[str]
    ) -> "GlobalSecondaryIndex":
        return GlobalSecondaryIndex(
            index_name=dct["IndexName"],
            schema=dct["KeySchema"],
            projection=dct["Projection"],
            table_key_attrs=table_key_attrs,
            throughput=dct.get("ProvisionedThroughput", None),
        )

    def update(self, u: Dict[str, Any]) -> None:
        self.name = u.get("IndexName", self.name)
        self.schema = u.get("KeySchema", self.schema)
        self.projection = u.get("Projection", self.projection)
        self.throughput = u.get("ProvisionedThroughput", self.throughput)


class StreamRecord(BaseModel):
    def __init__(
        self,
        table: "Table",
        stream_type: str,
        event_name: str,
        old: Optional[Item],
        new: Optional[Item],
        seq: int,
    ):
        old_a = old.to_json()["Attributes"] if old is not None else {}
        new_a = new.to_json()["Attributes"] if new is not None else {}

        rec = old if old is not None else new
        keys = {table.hash_key_attr: rec.hash_key.to_json()}  # type: ignore[union-attr]
        if table.range_key_attr is not None and rec is not None:
            keys[table.range_key_attr] = rec.range_key.to_json()  # type: ignore

        self.record: Dict[str, Any] = {
            "eventID": mock_random.uuid4().hex,
            "eventName": event_name,
            "eventSource": "aws:dynamodb",
            "eventVersion": "1.0",
            "awsRegion": "us-east-1",
            "dynamodb": {
                "StreamViewType": stream_type,
                "ApproximateCreationDateTime": utcnow().isoformat(),
                "SequenceNumber": str(seq),
                "SizeBytes": 1,
                "Keys": keys,
            },
        }

        if stream_type in ("NEW_IMAGE", "NEW_AND_OLD_IMAGES"):
            self.record["dynamodb"]["NewImage"] = new_a
        if stream_type in ("OLD_IMAGE", "NEW_AND_OLD_IMAGES"):
            self.record["dynamodb"]["OldImage"] = old_a

        # This is a substantial overestimate but it's the easiest to do now
        self.record["dynamodb"]["SizeBytes"] = len(
            dynamo_json_dump(self.record["dynamodb"])
        )

    def to_json(self) -> Dict[str, Any]:
        return self.record


class StreamShard(BaseModel):
    def __init__(self, account_id: str, table: "Table"):
        self.account_id = account_id
        self.table = table
        self.id = "shardId-00000001541626099285-f35f62ef"
        self.starting_sequence_number = 1100000000017454423009
        self.items: List[StreamRecord] = []
        self.created_on = utcnow()

    def to_json(self) -> Dict[str, Any]:
        return {
            "ShardId": self.id,
            "SequenceNumberRange": {
                "StartingSequenceNumber": str(self.starting_sequence_number)
            },
        }

    def add(self, old: Optional[Item], new: Optional[Item]) -> None:
        t = self.table.stream_specification["StreamViewType"]  # type: ignore
        if old is None:
            event_name = "INSERT"
        elif new is None:
            event_name = "REMOVE"
        else:
            event_name = "MODIFY"
        seq = len(self.items) + self.starting_sequence_number
        self.items.append(StreamRecord(self.table, t, event_name, old, new, seq))
        result = None
        from moto.awslambda.utils import get_backend

        for arn, esm in self.table.lambda_event_source_mappings.items():
            region = arn.split(":")[3]

            result = get_backend(self.account_id, region).send_dynamodb_items(
                arn, self.items, esm.event_source_arn
            )

        if result:
            self.items = []

    def get(self, start: int, quantity: int) -> List[Dict[str, Any]]:
        start -= self.starting_sequence_number
        assert start >= 0
        end = start + quantity
        return [i.to_json() for i in self.items[start:end]]


class Table(CloudFormationModel):
    def __init__(
        self,
        table_name: str,
        account_id: str,
        region: str,
        schema: List[Dict[str, Any]],
        attr: List[Dict[str, str]],
        throughput: Optional[Dict[str, int]] = None,
        billing_mode: Optional[str] = None,
        indexes: Optional[List[Dict[str, Any]]] = None,
        global_indexes: Optional[List[Dict[str, Any]]] = None,
        streams: Optional[Dict[str, Any]] = None,
        sse_specification: Optional[Dict[str, Any]] = None,
        tags: Optional[List[Dict[str, str]]] = None,
        deletion_protection_enabled: Optional[bool] = False,
    ):
        self.name = table_name
        self.account_id = account_id
        self.region_name = region
        self.attr = attr
        self.schema = schema
        self.range_key_attr: Optional[str] = None
        self.hash_key_attr: str = ""
        self.range_key_type: Optional[str] = None
        self.hash_key_type: str = ""
        for elem in schema:
            attr_type = [
                a["AttributeType"]
                for a in attr
                if a["AttributeName"] == elem["AttributeName"]
            ][0]
            if elem["KeyType"] == "HASH":
                self.hash_key_attr = elem["AttributeName"]
                self.hash_key_type = attr_type
            elif elem["KeyType"] == "RANGE":
                self.range_key_attr = elem["AttributeName"]
                self.range_key_type = attr_type
        self.table_key_attrs = [
            key for key in (self.hash_key_attr, self.range_key_attr) if key is not None
        ]
        self.billing_mode = billing_mode
        if throughput is None:
            self.throughput = {"WriteCapacityUnits": 0, "ReadCapacityUnits": 0}
        else:
            self.throughput = throughput
        self.throughput["NumberOfDecreasesToday"] = 0
        self.indexes = [
            LocalSecondaryIndex.create(i, self.table_key_attrs)
            for i in (indexes if indexes else [])
        ]
        self.global_indexes = [
            GlobalSecondaryIndex.create(i, self.table_key_attrs)
            for i in (global_indexes if global_indexes else [])
        ]
        self.created_at = utcnow()
        self.items = defaultdict(dict)  # type: ignore  # [hash: DynamoType] or [hash: [range: DynamoType]]
        self.table_arn = self._generate_arn(table_name)
        self.tags = tags or []
        self.ttl = {
            "TimeToLiveStatus": "DISABLED"  # One of 'ENABLING'|'DISABLING'|'ENABLED'|'DISABLED',
            # 'AttributeName': 'string'  # Can contain this
        }
        self.stream_specification: Optional[Dict[str, Any]] = {"StreamEnabled": False}
        self.latest_stream_label: Optional[str] = None
        self.stream_shard: Optional[StreamShard] = None
        self.set_stream_specification(streams)
        self.lambda_event_source_mappings: Dict[str, Any] = {}
        self.continuous_backups: Dict[str, Any] = {
            "ContinuousBackupsStatus": "ENABLED",  # One of 'ENABLED'|'DISABLED', it's enabled by default
            "PointInTimeRecoveryDescription": {
                "PointInTimeRecoveryStatus": "DISABLED"  # One of 'ENABLED'|'DISABLED'
            },
        }
        self.sse_specification = sse_specification
        if self.sse_specification and "KMSMasterKeyId" not in self.sse_specification:
            self.sse_specification["KMSMasterKeyId"] = self._get_default_encryption_key(
                account_id, region
            )
        self.deletion_protection_enabled = deletion_protection_enabled

    def _get_default_encryption_key(self, account_id: str, region: str) -> str:
        from moto.kms import kms_backends

        # https://aws.amazon.com/kms/features/#AWS_Service_Integration
        # An AWS managed CMK is created automatically when you first create
        # an encrypted resource using an AWS service integrated with KMS.
        kms = kms_backends[account_id][region]
        ddb_alias = "alias/aws/dynamodb"
        if not kms.alias_exists(ddb_alias):
            key = kms.create_key(
                policy="",
                key_usage="ENCRYPT_DECRYPT",
                key_spec="SYMMETRIC_DEFAULT",
                description="Default master key that protects my DynamoDB table storage",
                tags=None,
            )
            kms.create_alias(key.id, ddb_alias)
        ebs_key = kms.describe_key(ddb_alias)
        return ebs_key.arn

    @classmethod
    def has_cfn_attr(cls, attr: str) -> bool:
        return attr in ["Arn", "StreamArn"]

    def get_cfn_attribute(self, attribute_name: str) -> Any:
        from moto.cloudformation.exceptions import UnformattedGetAttTemplateException

        if attribute_name == "Arn":
            return self.table_arn
        elif attribute_name == "StreamArn" and self.stream_specification:
            return self.describe()["TableDescription"]["LatestStreamArn"]

        raise UnformattedGetAttTemplateException()

    @property
    def physical_resource_id(self) -> str:
        return self.name

    @property
    def attribute_keys(self) -> List[str]:
        # A set of all the hash or range attributes for all indexes
        def keys_from_index(idx: SecondaryIndex) -> List[str]:
            schema = idx.schema
            return [attr["AttributeName"] for attr in schema]

        fieldnames = copy.copy(self.table_key_attrs)
        for idx in self.indexes + self.global_indexes:
            fieldnames += keys_from_index(idx)
        return fieldnames

    @staticmethod
    def cloudformation_name_type() -> str:
        return "TableName"

    @staticmethod
    def cloudformation_type() -> str:
        # https://docs.aws.amazon.com/AWSCloudFormation/latest/UserGuide/aws-resource-dynamodb-table.html
        return "AWS::DynamoDB::Table"

    @classmethod
    def create_from_cloudformation_json(  # type: ignore[misc]
        cls,
        resource_name: str,
        cloudformation_json: Dict[str, Any],
        account_id: str,
        region_name: str,
        **kwargs: Any,
    ) -> "Table":
        from moto.dynamodb.models import dynamodb_backends

        properties = cloudformation_json["Properties"]
        params = {}

        if "KeySchema" in properties:
            params["schema"] = properties["KeySchema"]
        if "AttributeDefinitions" in properties:
            params["attr"] = properties["AttributeDefinitions"]
        params["global_indexes"] = properties.get("GlobalSecondaryIndexes", [])
        params["throughput"] = properties.get("ProvisionedThroughput")
        params["indexes"] = properties.get("LocalSecondaryIndexes", [])
        params["streams"] = properties.get("StreamSpecification")
        params["tags"] = properties.get("Tags")
        params["deletion_protection_enabled"] = properties.get(
            "DeletionProtectionEnabled", False
        )
        params["sse_specification"] = properties.get("SSESpecification")

        billing_mode = (
            "PAY_PER_REQUEST"
            if properties.get("BillingMode") == "PAY_PER_REQUEST"
            else "PROVISIONED"
        )
        params["billing_mode"] = billing_mode

        table = dynamodb_backends[account_id][region_name].create_table(
            name=resource_name, **params
        )
        return table

    @classmethod
    def delete_from_cloudformation_json(  # type: ignore[misc]
        cls,
        resource_name: str,
        cloudformation_json: Dict[str, Any],
        account_id: str,
        region_name: str,
    ) -> None:
        from moto.dynamodb.models import dynamodb_backends

        dynamodb_backends[account_id][region_name].delete_table(name=resource_name)

    def _generate_arn(self, name: str) -> str:
        return f"arn:{get_partition(self.region_name)}:dynamodb:{self.region_name}:{self.account_id}:table/{name}"

    def set_stream_specification(self, streams: Optional[Dict[str, Any]]) -> None:
        self.stream_specification = streams
        if (
            self.stream_specification
            and streams
            and (streams.get("StreamEnabled") or streams.get("StreamViewType"))
        ):
            self.stream_specification["StreamEnabled"] = True
            self.latest_stream_label = utcnow().isoformat()
            self.stream_shard = StreamShard(self.account_id, self)
        else:
            self.stream_specification = {"StreamEnabled": False}

    def describe(self, base_key: str = "TableDescription") -> Dict[str, Any]:
        results: Dict[str, Any] = {
            base_key: {
                "AttributeDefinitions": self.attr,
                "ProvisionedThroughput": self.throughput,
                "BillingModeSummary": {"BillingMode": self.billing_mode},
                "TableSizeBytes": 0,
                "TableName": self.name,
                "TableStatus": "ACTIVE",
                "TableArn": self.table_arn,
                "KeySchema": self.schema,
                "ItemCount": len(self),
                "CreationDateTime": unix_time(self.created_at),
                "GlobalSecondaryIndexes": [
                    index.describe() for index in self.global_indexes
                ],
                "LocalSecondaryIndexes": [index.describe() for index in self.indexes],
                "DeletionProtectionEnabled": self.deletion_protection_enabled,
            }
        }
        if self.latest_stream_label:
            results[base_key]["LatestStreamLabel"] = self.latest_stream_label
            results[base_key]["LatestStreamArn"] = (
                f"{self.table_arn}/stream/{self.latest_stream_label}"
            )
        if self.stream_specification and self.stream_specification["StreamEnabled"]:
            results[base_key]["StreamSpecification"] = self.stream_specification
        if self.sse_specification and self.sse_specification.get("Enabled") is True:
            results[base_key]["SSEDescription"] = {
                "Status": "ENABLED",
                "SSEType": "KMS",
                "KMSMasterKeyArn": self.sse_specification.get("KMSMasterKeyId"),
            }
        return results

    def __len__(self) -> int:
        return sum(
            [(len(value) if self.has_range_key else 1) for value in self.items.values()]
        )

    @property
    def hash_key_names(self) -> List[str]:
        keys = [self.hash_key_attr]
        for index in self.global_indexes:
            for key in index.schema:
                if key["KeyType"] == "HASH":
                    keys.append(key["AttributeName"])
        return keys

    @property
    def range_key_names(self) -> List[str]:
        keys = [self.range_key_attr] if self.has_range_key else []
        for index in self.global_indexes:
            for key in index.schema:
                if key["KeyType"] == "RANGE":
                    keys.append(key["AttributeName"])
        return keys  # type: ignore[return-value]

    def _validate_key_sizes(self, item_attrs: Dict[str, Any]) -> None:
        for hash_name in self.hash_key_names:
            hash_value = item_attrs.get(hash_name)
            if hash_value:
                if DynamoType(hash_value).size() > HASH_KEY_MAX_LENGTH:
                    raise HashKeyTooLong
        for range_name in self.range_key_names:
            range_value = item_attrs.get(range_name)
            if range_value:
                if DynamoType(range_value).size() > RANGE_KEY_MAX_LENGTH:
                    raise RangeKeyTooLong

    def _validate_item_types(
        self, item_attrs: Dict[str, Any], attr: Optional[str] = None
    ) -> None:
        for key, value in item_attrs.items():
            if isinstance(value, dict):
                self._validate_item_types(value, attr=key if attr is None else key)
            elif isinstance(value, int) and key == "N":
                raise InvalidConversion
            if key == "S":
                # This scenario is usually caught by boto3, but the user can disable parameter validation
                # Which is why we need to catch it 'server-side' as well
                if isinstance(value, int):
                    raise SerializationException(
                        "NUMBER_VALUE cannot be converted to String"
                    )
                if attr and attr in self.table_key_attrs and isinstance(value, dict):
                    raise SerializationException(
                        "Start of structure or map found where not expected"
                    )

    def put_item(
        self,
        item_attrs: Dict[str, Any],
        expected: Optional[Dict[str, Any]] = None,
        condition_expression: Optional[str] = None,
        expression_attribute_names: Optional[Dict[str, str]] = None,
        expression_attribute_values: Optional[Dict[str, Any]] = None,
        overwrite: bool = False,
        return_values_on_condition_check_failure: Optional[str] = None,
    ) -> Item:
        if self.hash_key_attr not in item_attrs.keys():
            raise MockValidationException(
                "One or more parameter values were invalid: Missing the key "
                + self.hash_key_attr
                + " in the item"
            )
        hash_value = DynamoType(item_attrs[self.hash_key_attr])
        if self.range_key_attr is not None:
            if self.range_key_attr not in item_attrs.keys():
                raise MockValidationException(
                    f"One or more parameter values were invalid: Missing the key {self.range_key_attr} in the item"
                )
            range_value = DynamoType(item_attrs[self.range_key_attr])
        else:
            range_value = None

        if hash_value.type != self.hash_key_type:
            raise InvalidAttributeTypeError(
                self.hash_key_attr,
                expected_type=self.hash_key_type,
                actual_type=hash_value.type,
            )
        if range_value and range_value.type != self.range_key_type:
            raise InvalidAttributeTypeError(
                self.range_key_attr,
                expected_type=self.range_key_type,
                actual_type=range_value.type,
            )

        self._validate_item_types(item_attrs)
        self._validate_key_sizes(item_attrs)

        if expected is None:
            expected = {}
            lookup_range_value = range_value
        else:
            expected_range_value = expected.get(self.range_key_attr, {}).get("Value")  # type: ignore
            if expected_range_value is None:
                lookup_range_value = range_value
            else:
                lookup_range_value = DynamoType(expected_range_value)
        current = self.get_item(hash_value, lookup_range_value)
        item = Item(hash_value, range_value, item_attrs)

        if not overwrite:
            if not get_expected(expected).expr(current):
                raise ConditionalCheckFailed
            condition_op = get_filter_expression(
                condition_expression,
                expression_attribute_names,
                expression_attribute_values,
            )
            if not condition_op.expr(current):
                if (
                    return_values_on_condition_check_failure == "ALL_OLD"
                    and current is not None
                ):
                    raise ConditionalCheckFailed(item=current.to_json()["Attributes"])
                else:
                    raise ConditionalCheckFailed

        if range_value:
            self.items[hash_value][range_value] = item
        else:
            self.items[hash_value] = item  # type: ignore[assignment]

        if self.stream_shard is not None:
            self.stream_shard.add(current, item)

        return item

    def __nonzero__(self) -> bool:
        return True

    def __bool__(self) -> bool:
        return self.__nonzero__()

    @property
    def has_range_key(self) -> bool:
        return self.range_key_attr is not None

    def get_item(
        self,
        hash_key: DynamoType,
        range_key: Optional[DynamoType] = None,
        projection_expression: Optional[List[List[str]]] = None,
    ) -> Optional[Item]:
        if self.has_range_key and not range_key:
            raise MockValidationException(
                "Table has a range key, but no range key was passed into get_item"
            )
        try:
            result = None

            if range_key:
                result = self.items[hash_key][range_key]
            elif hash_key in self.items:
                result = self.items[hash_key]

            if projection_expression and result:
                result = result.project(projection_expression)

            return result
        except KeyError:
            return None

    def delete_item(
        self, hash_key: DynamoType, range_key: Optional[DynamoType]
    ) -> Optional[Item]:
        try:
            if range_key:
                item = self.items[hash_key].pop(range_key)
            else:
                item = self.items.pop(hash_key)

            if self.stream_shard is not None:
                self.stream_shard.add(item, None)

            return item
        except KeyError:
            return None

    def query(
        self,
        hash_key: DynamoType,
        range_comparison: Optional[str],
        range_objs: List[DynamoType],
        limit: int,
        exclusive_start_key: Dict[str, Any],
        scan_index_forward: bool,
        projection_expressions: Optional[List[List[str]]],
        index_name: Optional[str] = None,
        consistent_read: bool = False,
        filter_expression: Any = None,
        **filter_kwargs: Any,
    ) -> Tuple[List[Item], int, Optional[Dict[str, Any]]]:
        # FIND POSSIBLE RESULTS
        if index_name:
            all_indexes = self.all_indexes()
            indexes_by_name = dict((i.name, i) for i in all_indexes)
            if index_name not in indexes_by_name:
                all_names = ", ".join(indexes_by_name.keys())
                raise MockValidationException(
                    f"Invalid index: {index_name} for table: {self.name}. Available indexes are: {all_names}"
                )

            index = indexes_by_name[index_name]

            if consistent_read and index in self.global_indexes:
                raise MockValidationException(
                    "Consistent reads are not supported on global secondary indexes"
                )

            try:
                index_hash_key = [
                    key for key in index.schema if key["KeyType"] == "HASH"
                ][0]
            except IndexError:
                raise MockValidationException(
                    f"Missing Hash Key. KeySchema: {index.name}"
                )

            try:
                index_range_key = [
                    key for key in index.schema if key["KeyType"] == "RANGE"
                ][0]
            except IndexError:
                if isinstance(index, GlobalSecondaryIndex) and self.range_key_attr:
                    # If we're querying a GSI that does not have a range key, the main range key acts as a range key
                    index_range_key = {"AttributeName": self.range_key_attr}
                else:
                    # If we don't have a range key on the main table either, the hash key acts as a range key
                    index_range_key = {"AttributeName": self.hash_key_attr}
                if range_comparison:
                    raise ValueError(
                        f"Range Key comparison but no range key found for index: {index_name}"
                    )

            hash_attrs = [index_hash_key["AttributeName"], self.hash_key_attr]
            if index_range_key:
                range_attrs = [
                    index_range_key["AttributeName"],
                    self.range_key_attr,
                ]
            else:
                range_attrs = [self.range_key_attr]

            possible_results = []
            for item in self.all_items():
                if not isinstance(item, Item):
                    continue
                item_hash_key = item.attrs.get(hash_attrs[0])
                if len(range_attrs) == 1:
                    if item_hash_key and item_hash_key == hash_key:
                        possible_results.append(item)
                else:
                    item_range_key = item.attrs.get(range_attrs[0])  # type: ignore
                    if item_hash_key and item_hash_key == hash_key and item_range_key:
                        possible_results.append(item)
        else:
            hash_attrs = [self.hash_key_attr]
            range_attrs = [self.range_key_attr]

            possible_results = [
                item
                for item in self.all_items()
                if isinstance(item, Item) and item.hash_key == hash_key
            ]

        # SORT
        if index_name:
            possible_results = self.sorted_items(
                hash_attrs, range_attrs, possible_results
            )
        else:
            possible_results.sort(key=lambda item: item.range_key)  # type: ignore

        if scan_index_forward is False:
            possible_results.reverse()

        # FILTER
        results: List[Item] = []
        result_size = 0
        scanned_count = 0
        last_evaluated_key = None
        processing_previous_page = exclusive_start_key is not None
        for result in possible_results:
            # Cycle through the previous page of results
            # When we encounter our start key, we know we've reached the end of the previous page
            if processing_previous_page:
                if self._item_comes_before_dct(
                    result,
                    exclusive_start_key,
                    hash_attrs,
                    range_attrs,
                    scan_index_forward,
                ):
                    continue
                else:
                    processing_previous_page = False

            # Check wether we've reached the limit of our result set
            # That can be either in number, or in size
            reached_length_limit = len(results) == limit
            reached_size_limit = (result_size + result.size()) > RESULT_SIZE_LIMIT
            if reached_length_limit or reached_size_limit:
                last_evaluated_key = self._get_last_evaluated_key(
                    results[-1], index_name
                )
                break

            if not range_comparison and not filter_kwargs:
                # If we're not filtering on range key or on an index
                results.append(result)
                result_size += result.size()
                scanned_count += 1

            if range_comparison:
                if (
                    index_name
                    and index_range_key
                    and result.attrs.get(index_range_key["AttributeName"])
                ):
                    if result.attrs.get(index_range_key["AttributeName"]).compare(  # type: ignore
                        range_comparison, range_objs
                    ):
                        results.append(result)
                        result_size += result.size()
                        scanned_count += 1
                else:
                    if result.range_key.compare(range_comparison, range_objs):  # type: ignore[union-attr]
                        results.append(result)
                        result_size += result.size()
                        scanned_count += 1

            if filter_kwargs:
                for field, value in filter_kwargs.items():
                    dynamo_types = [
                        DynamoType(ele) for ele in value["AttributeValueList"]
                    ]
                    if result.attrs.get(field).compare(  # type: ignore[union-attr]
                        value["ComparisonOperator"], dynamo_types
                    ):
                        results.append(result)
                        result_size += result.size()
                scanned_count += 1

        results = copy.deepcopy(results)
        if index_name:
            index = self.get_index(index_name)
            results = [index.project(r) for r in results]

        if filter_expression is not None:
            results = [item for item in results if filter_expression.expr(item)]

        if projection_expressions:
            results = [r.project(projection_expressions) for r in results]

        return results, scanned_count, last_evaluated_key

    def all_items(self) -> List[Item]:
        items: List[Item] = []
        for hash_set in self.items.values():
            if self.range_key_attr:
                for item in hash_set.values():
                    items.append(item)
            else:
                items.append(hash_set)  # type: ignore
        return sorted(
            items,
            key=lambda x: (x.hash_key, x.range_key) if x.range_key else x.hash_key,
        )

    def all_indexes(self) -> Sequence[SecondaryIndex]:
        return (self.global_indexes or []) + (self.indexes or [])  # type: ignore

    def get_index(self, index_name: str, error_if_not: bool = False) -> SecondaryIndex:
        all_indexes = self.all_indexes()
        indexes_by_name = dict((i.name, i) for i in all_indexes)
        if error_if_not and index_name not in indexes_by_name:
            raise InvalidIndexNameError(
                f"The table does not have the specified index: {index_name}"
            )
        return indexes_by_name[index_name]

    def has_idx_items(self, index_name: str) -> List[Item]:
        idx = self.get_index(index_name)
        idx_col_set = set([i["AttributeName"] for i in idx.schema])

        items: List[Item] = []

        for hash_set in self.items.values():
            if self.range_key_attr:
                for item in hash_set.values():
                    if idx_col_set.issubset(set(item.attrs)):
                        items.append(item)
            else:
                if idx_col_set.issubset(set(hash_set.attrs)):  # type: ignore
                    items.append(hash_set)  # type: ignore
        return items

    def scan(
        self,
        filters: Dict[str, Any],
        limit: int,
        exclusive_start_key: Dict[str, Any],
        filter_expression: Any = None,
        index_name: Optional[str] = None,
        consistent_read: bool = False,
        projection_expression: Optional[List[List[str]]] = None,
        segments: Union[Tuple[None, None], Tuple[int, int]] = (None, None),
    ) -> Tuple[List[Item], int, Optional[Dict[str, Any]]]:
        results: List[Item] = []
        result_size = 0
        scanned_count = 0

        if index_name:
            index = self.get_index(index_name, error_if_not=True)

            if consistent_read and index in self.global_indexes:
                raise MockValidationException(
                    "Consistent reads are not supported on global secondary indexes"
                )

            try:
                index_hash_key = [
                    key for key in index.schema if key["KeyType"] == "HASH"
                ][0]
            except IndexError:
                raise MockValidationException(
                    f"Missing Hash Key. KeySchema: {index.name}"
                )

            try:
                index_range_key = [
                    key for key in index.schema if key["KeyType"] == "RANGE"
                ][0]
            except IndexError:
                index_range_key = None

            hash_attrs = [index_hash_key["AttributeName"], self.hash_key_attr]
            if index_range_key:
                range_attrs = [index_range_key["AttributeName"], self.range_key_attr]
            else:
                range_attrs = [self.range_key_attr]

            items = self.has_idx_items(index_name)
            items = self.sorted_items(hash_attrs, range_attrs, items)
        else:
            hash_attrs = [self.hash_key_attr]
            range_attrs = [self.range_key_attr]

            items = self.all_items()

        last_evaluated_key = None
        processing_previous_page = exclusive_start_key is not None
        for item in items:
            if not item.is_within_segment(segments):
                continue

            # Cycle through the previous page of results
            # When we encounter our start key, we know we've reached the end of the previous page
            if processing_previous_page:
                if self._item_comes_before_dct(
                    item,
                    exclusive_start_key,
                    hash_attrs,
                    range_attrs,
                    True,
                ):
                    continue
                else:
                    processing_previous_page = False

            # Check whether we've reached the limit of our result set
            # That can be either in number, or in size
            reached_length_limit = len(results) == limit
            reached_size_limit = (result_size + item.size()) > RESULT_SIZE_LIMIT
            if reached_length_limit or reached_size_limit:
                last_evaluated_key = self._get_last_evaluated_key(
                    results[-1], index_name
                )
                break

            passes_all_conditions = True
            for attribute_name in filters:
                attribute = item.attrs.get(attribute_name)
                (comparison_operator, comparison_objs) = filters[attribute_name]

                if attribute:
                    # Attribute found
                    if not attribute.compare(comparison_operator, comparison_objs):
                        passes_all_conditions = False
                        break
                elif comparison_operator == "NULL":
                    # Comparison is NULL and we don't have the attribute
                    continue
                else:
                    # No attribute found and comparison is no NULL. This item
                    # fails
                    passes_all_conditions = False
                    break

            if passes_all_conditions:
                if index_name:
                    index = self.get_index(index_name)
                    results.append(index.project(copy.deepcopy(item)))
                else:
                    results.append(copy.deepcopy(item))
                result_size += item.size()

            scanned_count += 1

        # https://docs.aws.amazon.com/amazondynamodb/latest/developerguide/Query.html#Query.FilterExpression
        # the filter expression should be evaluated after the query.
        if filter_expression is not None:
            results = [item for item in results if filter_expression.expr(item)]

        if projection_expression:
            results = [r.project(projection_expression) for r in results]

        return results, scanned_count, last_evaluated_key

    def _item_comes_before_dct(
        self,
        item: Item,
        dct: Dict[str, Any],
        hash_key_attrs: List[str],
        range_key_attrs: List[Optional[str]],
        scan_index_forward: bool,
    ) -> bool:
        """
        Does item appear before or at dct relative to sort options?

        hash_key_attrs: The list of hash keys.
        Includes the key of the GSI (first, if it exists) and the key of the main table.

        range_key_attrs: A list of range keys (RK).
        Includes the RK-name of the GSI (first, if it exists) and the RK-name of the main table (second - can be None).
        When sorting a GSI, we'll try to sort by the GSI RK first.
        However, because GSI RK's are not unique (by design), item and dct can have the same RK-value
        If that is the case, we compare by the RK of the main table instead
        Related: https://github.com/getmoto/moto/issues/7761
        """
        attrs_to_sort_by = self._generate_attr_to_sort_by(
            hash_key_attrs, range_key_attrs
        )
        for attr in attrs_to_sort_by:
            if attr in item.attrs and item.attrs[attr] != DynamoType(dct.get(attr)):  # type: ignore
                return (
                    (item.attrs[attr] < DynamoType(dct.get(attr)))  # type: ignore
                    == scan_index_forward
                )
        # Keys were equal, items are identical
        return True

    def sorted_items(
        self,
        hash_key_attrs: List[str],
        range_key_attrs: List[Optional[str]],
        items: List[Item],
    ) -> List[Item]:
        attrs_to_sort_by = self._generate_attr_to_sort_by(
            hash_key_attrs, range_key_attrs
        )
        items.sort(
            key=lambda x: tuple([x.attrs[key] for key in attrs_to_sort_by]),
        )
        return items

    def _generate_attr_to_sort_by(
        self, hash_key_attrs: List[str], range_key_attrs: List[Optional[str]]
    ) -> List[str]:
        gsi_hash_key = hash_key_attrs[0] if len(hash_key_attrs) == 2 else None
        table_hash_key = str(
            hash_key_attrs[0] if gsi_hash_key is None else hash_key_attrs[1]
        )
        gsi_range_key = range_key_attrs[0] if len(range_key_attrs) == 2 else None
        table_range_key = str(
            range_key_attrs[0] if gsi_range_key is None else range_key_attrs[1]
        )
        # Gets the GSI and table hash and range keys in the order to try sorting by
        attrs_to_sort_by = [
            gsi_hash_key,
            gsi_range_key,
            table_hash_key,
            table_range_key,
        ]
        return [
            attr for attr in attrs_to_sort_by if attr is not None and attr != "None"
        ]

    def _get_last_evaluated_key(
        self, last_result: Item, index_name: Optional[str]
    ) -> Dict[str, Any]:
        last_evaluated_key = {self.hash_key_attr: last_result.hash_key}
        if self.range_key_attr is not None and last_result.range_key is not None:
            last_evaluated_key[self.range_key_attr] = last_result.range_key
        if index_name:
            index = self.get_index(index_name)
            idx_col_list = [i["AttributeName"] for i in index.schema]
            for col in idx_col_list:
                last_evaluated_key[col] = last_result.attrs[col]
        return last_evaluated_key

    def delete(self, account_id: str, region_name: str) -> None:
        from moto.dynamodb.models import dynamodb_backends

        dynamodb_backends[account_id][region_name].delete_table(self.name)


class Backup:
    def __init__(
        self,
        account_id: str,
        region_name: str,
        name: str,
        table: Table,
        status: Optional[str] = None,
        type_: Optional[str] = None,
    ):
        self.region_name = region_name
        self.account_id = account_id
        self.name = name
        self.table = copy.deepcopy(table)
        self.status = status or "AVAILABLE"
        self.type = type_ or "USER"
        self.creation_date_time = utcnow()
        self.identifier = self._make_identifier()

    def _make_identifier(self) -> str:
        timestamp = int(unix_time_millis(self.creation_date_time))
        timestamp_padded = str("0" + str(timestamp))[-16:16]
        guid = str(mock_random.uuid4())
        guid_shortened = guid[:8]
        return f"{timestamp_padded}-{guid_shortened}"

    @property
    def arn(self) -> str:
        return f"arn:{get_partition(self.region_name)}:dynamodb:{self.region_name}:{self.account_id}:table/{self.table.name}/backup/{self.identifier}"

    @property
    def details(self) -> Dict[str, Any]:  # type: ignore[misc]
        return {
            "BackupArn": self.arn,
            "BackupName": self.name,
            "BackupSizeBytes": 123,
            "BackupStatus": self.status,
            "BackupType": self.type,
            "BackupCreationDateTime": unix_time(self.creation_date_time),
        }

    @property
    def summary(self) -> Dict[str, Any]:  # type: ignore[misc]
        return {
            "TableName": self.table.name,
            # 'TableId': 'string',
            "TableArn": self.table.table_arn,
            "BackupArn": self.arn,
            "BackupName": self.name,
            "BackupCreationDateTime": unix_time(self.creation_date_time),
            # 'BackupExpiryDateTime': datetime(2015, 1, 1),
            "BackupStatus": self.status,
            "BackupType": self.type,
            "BackupSizeBytes": 123,
        }

    @property
    def description(self) -> Dict[str, Any]:  # type: ignore[misc]
        source_table_details = self.table.describe()["TableDescription"]
        source_table_details["TableCreationDateTime"] = source_table_details[
            "CreationDateTime"
        ]
        description = {
            "BackupDetails": self.details,
            "SourceTableDetails": source_table_details,
        }
        return description


class RestoredTable(Table):
    def __init__(self, name: str, account_id: str, region: str, backup: "Backup"):
        params = self._parse_params_from_backup(backup)
        super().__init__(name, account_id=account_id, region=region, **params)
        self.indexes = copy.deepcopy(backup.table.indexes)
        self.global_indexes = copy.deepcopy(backup.table.global_indexes)
        self.items = copy.deepcopy(backup.table.items)
        # Restore Attrs
        self.source_backup_arn = backup.arn
        self.source_table_arn = backup.table.table_arn
        self.restore_date_time = self.created_at

    def _parse_params_from_backup(self, backup: "Backup") -> Dict[str, Any]:
        return {
            "schema": copy.deepcopy(backup.table.schema),
            "attr": copy.deepcopy(backup.table.attr),
            "throughput": copy.deepcopy(backup.table.throughput),
        }

    def describe(self, base_key: str = "TableDescription") -> Dict[str, Any]:
        result = super().describe(base_key=base_key)
        result[base_key]["RestoreSummary"] = {
            "SourceBackupArn": self.source_backup_arn,
            "SourceTableArn": self.source_table_arn,
            "RestoreDateTime": unix_time(self.restore_date_time),
            "RestoreInProgress": False,
        }
        return result


class RestoredPITTable(Table):
    def __init__(self, name: str, account_id: str, region: str, source: Table):
        params = self._parse_params_from_table(source)
        super().__init__(name, account_id=account_id, region=region, **params)
        self.indexes = copy.deepcopy(source.indexes)
        self.global_indexes = copy.deepcopy(source.global_indexes)
        self.items = copy.deepcopy(source.items)
        # Restore Attrs
        self.source_table_arn = source.table_arn
        self.restore_date_time = self.created_at

    def _parse_params_from_table(self, table: Table) -> Dict[str, Any]:
        return {
            "schema": copy.deepcopy(table.schema),
            "attr": copy.deepcopy(table.attr),
            "throughput": copy.deepcopy(table.throughput),
        }

    def describe(self, base_key: str = "TableDescription") -> Dict[str, Any]:
        result = super().describe(base_key=base_key)
        result[base_key]["RestoreSummary"] = {
            "SourceTableArn": self.source_table_arn,
            "RestoreDateTime": unix_time(self.restore_date_time),
            "RestoreInProgress": False,
        }
        return result


class ResourcePolicy:
    def __init__(self, resource_arn: str, policy_doc: str):
        self.resource_arn = resource_arn
        self.policy_doc = policy_doc
        self.revision_id = str(int(unix_time_millis()))
