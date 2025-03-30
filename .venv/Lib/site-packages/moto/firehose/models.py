"""FirehoseBackend class with methods for supported APIs.

Incomplete list of unfinished items:
  - Data record size and number of transactions are ignored.
  - Better validation of delivery destination parameters, e.g.,
    validation of the url for an http endpoint (boto3 does this).
  - Better handling of the put_record_batch() API.  Not only is
    the existing logic bare bones, but for the ElasticSearch and
    RedShift destinations, the data is just ignored.
  - put_record_batch() handling of errors is minimal and no errors
    are reported back to the user.  Instead an exception is raised.
  - put_record(), put_record_batch() always set "Encrypted" to False.
"""

import io
import json
import warnings
from base64 import b64decode, b64encode
from datetime import datetime, timezone
from gzip import GzipFile
from time import time
from typing import Any, Dict, List, Optional, Tuple

import requests

from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.common_models import BaseModel
from moto.core.utils import utcnow
from moto.firehose.exceptions import (
    ConcurrentModificationException,
    InvalidArgumentException,
    LimitExceededException,
    ResourceInUseException,
    ResourceNotFoundException,
    ValidationException,
)
from moto.moto_api._internal import mock_random
from moto.s3.models import s3_backends
from moto.utilities.tagging_service import TaggingService
from moto.utilities.utils import get_partition

MAX_TAGS_PER_DELIVERY_STREAM = 50

DESTINATION_TYPES_TO_NAMES = {
    "s3": "S3",
    "extended_s3": "ExtendedS3",
    "http_endpoint": "HttpEndpoint",
    "elasticsearch": "Elasticsearch",
    "redshift": "Redshift",
    "snowflake": "Snowflake",
    "splunk": "Splunk",  # Unimplemented
}


def find_destination_config_in_args(api_args: Dict[str, Any]) -> Tuple[str, Any]:
    """Return (config_name, config) tuple for destination config.

    Determines which destination config(s) have been specified.  The
    alternative is to use a bunch of 'if' statements to check each
    destination configuration.  If more than one destination config is
    specified, than an exception is raised.

    A logical name for the destination type is returned along with the
    destination config as it's useful way to compare current and replacement
    destinations.
    """
    destination_names = DESTINATION_TYPES_TO_NAMES.keys()
    configs = []
    for arg_name, arg_value in api_args.items():
        # Ignore arguments that are not destination configs.
        if "_destination" not in arg_name:
            continue

        # If the destination config value is non-null, save it.
        name = arg_name.split("_destination")[0]
        if name in destination_names and arg_value:
            configs.append((DESTINATION_TYPES_TO_NAMES[name], arg_value))

    # One and only one destination configuration is allowed.
    if len(configs) != 1:
        raise InvalidArgumentException(
            "Exactly one destination configuration is supported for a Firehose"
        )

    return configs[0]


def create_s3_destination_config(
    extended_s3_destination_config: Dict[str, Any],
) -> Dict[str, Any]:
    """Return dict with selected fields copied from ExtendedS3 config.

    When an ExtendedS3 config is chosen, AWS tacks on a S3 config as
    well.  When the same field names for S3 and ExtendedS3 exists,
    the ExtendedS3 fields are copied to the added S3 destination.
    """
    fields_not_needed = [
        "S3BackupMode",
        "S3Description",
        "DataFormatconversionConfiguration",
        "DynamicPartitionConfiguration",
    ]
    destination = {}
    for field, value in extended_s3_destination_config.items():
        if field in fields_not_needed:
            continue
        destination[field] = value
    return destination


class DeliveryStream(BaseModel):  # pylint: disable=too-few-public-methods,too-many-instance-attributes
    """Represents a delivery stream, its source and destination configs."""

    STATES = {"CREATING", "ACTIVE", "CREATING_FAILED"}

    MAX_STREAMS_PER_REGION = 50

    ALTERNATIVE_FIELD_NAMES = [
        ("S3Configuration", "S3DestinationDescription"),
        ("S3Update", "S3DestinationDescription"),
        ("S3BackupConfiguration", "S3BackupDescription"),
        ("S3BackupUpdate", "S3BackupDescription"),
    ]

    def __init__(
        self,
        account_id: str,
        region: str,
        delivery_stream_name: str,
        delivery_stream_type: str,
        encryption: Dict[str, Any],
        kinesis_stream_source_configuration: Dict[str, Any],
        destination_name: str,
        destination_config: Dict[str, Any],
    ):  # pylint: disable=too-many-arguments
        self.delivery_stream_status = "CREATING"
        self.delivery_stream_name = delivery_stream_name
        self.delivery_stream_type = (
            delivery_stream_type if delivery_stream_type else "DirectPut"
        )
        self.delivery_stream_encryption_configuration = encryption

        self.source = kinesis_stream_source_configuration
        self.destinations: List[Dict[str, Any]] = [
            {
                "destination_id": "destinationId-000000000001",
                destination_name: destination_config,
            }
        ]

        if destination_name == "ExtendedS3":
            # Add a S3 destination as well, minus a few ExtendedS3 fields.
            self.destinations[0]["S3"] = create_s3_destination_config(
                destination_config
            )

        # S3Configuration becomes S3DestinationDescription for the
        # other destinations. Same for S3Backup
        for old, new in DeliveryStream.ALTERNATIVE_FIELD_NAMES:
            if old in destination_config:
                self.destinations[0][destination_name][new] = destination_config[old]
                del self.destinations[0][destination_name][old]

        self.delivery_stream_status = "ACTIVE"
        self.delivery_stream_arn = f"arn:{get_partition(region)}:firehose:{region}:{account_id}:deliverystream/{delivery_stream_name}"

        self.create_timestamp = datetime.now(timezone.utc).isoformat()
        self.version_id = "1"  # Used to track updates of destination configs

        # I believe boto3 only adds this field after an update ...
        self.last_update_timestamp = datetime.now(timezone.utc).isoformat()


class FirehoseBackend(BaseBackend):
    """Implementation of Firehose APIs."""

    def __init__(self, region_name: str, account_id: str):
        super().__init__(region_name, account_id)
        self.delivery_streams: Dict[str, DeliveryStream] = {}
        self.tagger = TaggingService()

    @staticmethod
    def default_vpc_endpoint_service(
        service_region: str, zones: List[str]
    ) -> List[Dict[str, str]]:
        """Default VPC endpoint service."""
        return BaseBackend.default_vpc_endpoint_service_factory(
            service_region, zones, "firehose", special_service_name="kinesis-firehose"
        )

    def create_delivery_stream(  # pylint: disable=unused-argument
        self,
        region: str,
        delivery_stream_name: str,
        delivery_stream_type: str,
        kinesis_stream_source_configuration: Dict[str, Any],
        delivery_stream_encryption_configuration_input: Dict[str, Any],
        s3_destination_configuration: Dict[str, Any],
        extended_s3_destination_configuration: Dict[str, Any],
        redshift_destination_configuration: Dict[str, Any],
        elasticsearch_destination_configuration: Dict[str, Any],
        splunk_destination_configuration: Dict[str, Any],
        http_endpoint_destination_configuration: Dict[str, Any],
        snowflake_destination_configuration: Dict[str, Any],
        tags: List[Dict[str, str]],
    ) -> str:
        """Create a Kinesis Data Firehose delivery stream."""
        (destination_name, destination_config) = find_destination_config_in_args(
            locals()
        )

        if delivery_stream_name in self.delivery_streams:
            raise ResourceInUseException(
                f"Firehose {delivery_stream_name} under accountId {self.account_id} "
                f"already exists"
            )

        if len(self.delivery_streams) == DeliveryStream.MAX_STREAMS_PER_REGION:
            raise LimitExceededException(
                f"You have already consumed your firehose quota of "
                f"{DeliveryStream.MAX_STREAMS_PER_REGION} hoses. Firehose "
                f"names: {list(self.delivery_streams.keys())}"
            )

        # Rule out situations that are not yet implemented.
        if destination_name == "Splunk":
            warnings.warn("A Splunk destination delivery stream is not yet implemented")

        if (
            kinesis_stream_source_configuration
            and delivery_stream_type != "KinesisStreamAsSource"
        ):
            raise InvalidArgumentException(
                "KinesisSourceStreamConfig is only applicable for "
                "KinesisStreamAsSource stream type"
            )

        # Validate the tags before proceeding.
        errmsg = self.tagger.validate_tags(tags or [])
        if errmsg:
            raise ValidationException(errmsg)

        if tags and len(tags) > MAX_TAGS_PER_DELIVERY_STREAM:
            raise ValidationException(
                f"1 validation error detected: Value '{tags}' at 'tags' "
                f"failed to satisify contstraint: Member must have length "
                f"less than or equal to {MAX_TAGS_PER_DELIVERY_STREAM}"
            )

        # Create a DeliveryStream instance that will be stored and indexed
        # by delivery stream name.  This instance will update the state and
        # create the ARN.
        delivery_stream = DeliveryStream(
            account_id=self.account_id,
            region=region,
            delivery_stream_name=delivery_stream_name,
            delivery_stream_type=delivery_stream_type,
            encryption=delivery_stream_encryption_configuration_input,
            kinesis_stream_source_configuration=kinesis_stream_source_configuration,
            destination_name=destination_name,
            destination_config=destination_config,
        )
        self.tagger.tag_resource(delivery_stream.delivery_stream_arn, tags or [])

        self.delivery_streams[delivery_stream_name] = delivery_stream
        return self.delivery_streams[delivery_stream_name].delivery_stream_arn

    def delete_delivery_stream(self, delivery_stream_name: str) -> None:
        """Delete a delivery stream and its data.

        AllowForceDelete option is ignored as we only superficially
        apply state.
        """
        delivery_stream = self.delivery_streams.get(delivery_stream_name)
        if not delivery_stream:
            raise ResourceNotFoundException(
                f"Firehose {delivery_stream_name} under account {self.account_id} "
                f"not found."
            )

        self.tagger.delete_all_tags_for_resource(delivery_stream.delivery_stream_arn)

        delivery_stream.delivery_stream_status = "DELETING"
        self.delivery_streams.pop(delivery_stream_name)

    def describe_delivery_stream(self, delivery_stream_name: str) -> Dict[str, Any]:
        """Return description of specified delivery stream and its status.

        Note:  the 'limit' and 'exclusive_start_destination_id' parameters
        are not currently processed/implemented.
        """
        delivery_stream = self.delivery_streams.get(delivery_stream_name)
        if not delivery_stream:
            raise ResourceNotFoundException(
                f"Firehose {delivery_stream_name} under account {self.account_id} "
                f"not found."
            )

        result: Dict[str, Any] = {
            "DeliveryStreamDescription": {"HasMoreDestinations": False}
        }
        for attribute, attribute_value in vars(delivery_stream).items():
            if not attribute_value:
                continue

            # Convert from attribute's snake case to camel case for outgoing
            # JSON.
            name = "".join([x.capitalize() for x in attribute.split("_")])

            # Fooey ... always an exception to the rule:
            if name == "DeliveryStreamArn":
                name = "DeliveryStreamARN"

            if name != "Destinations":
                if name == "Source":
                    result["DeliveryStreamDescription"][name] = {
                        "KinesisStreamSourceDescription": attribute_value
                    }
                else:
                    result["DeliveryStreamDescription"][name] = attribute_value
                continue

            result["DeliveryStreamDescription"]["Destinations"] = []
            for destination in attribute_value:
                description = {}
                for key, value in destination.items():
                    if key == "destination_id":
                        description["DestinationId"] = value
                    else:
                        description[f"{key}DestinationDescription"] = value

                result["DeliveryStreamDescription"]["Destinations"].append(description)

        return result

    def list_delivery_streams(
        self,
        limit: Optional[int],
        delivery_stream_type: str,
        exclusive_start_delivery_stream_name: str,
    ) -> Dict[str, Any]:
        """Return list of delivery streams in alphabetic order of names."""
        result = {"DeliveryStreamNames": [], "HasMoreDeliveryStreams": False}
        if not self.delivery_streams:
            return result

        # If delivery_stream_type is specified, filter out any stream that's
        # not of that type.
        stream_list = list(self.delivery_streams.keys())
        if delivery_stream_type:
            stream_list = [
                x
                for x in stream_list
                if self.delivery_streams[x].delivery_stream_type == delivery_stream_type
            ]

        # The list is sorted alphabetically, not alphanumerically.
        sorted_list = sorted(stream_list)

        # Determine the limit or number of names to return in the list.
        limit = limit or DeliveryStream.MAX_STREAMS_PER_REGION

        # If a starting delivery stream name is given, find the index into
        # the sorted list, then add one to get the name following it.  If the
        # exclusive_start_delivery_stream_name doesn't exist, it's ignored.
        start = 0
        if exclusive_start_delivery_stream_name:
            if self.delivery_streams.get(exclusive_start_delivery_stream_name):
                start = sorted_list.index(exclusive_start_delivery_stream_name) + 1

        result["DeliveryStreamNames"] = sorted_list[start : start + limit]
        if len(sorted_list) > (start + limit):
            result["HasMoreDeliveryStreams"] = True
        return result

    def list_tags_for_delivery_stream(
        self,
        delivery_stream_name: str,
        exclusive_start_tag_key: str,
        limit: Optional[int],
    ) -> Dict[str, Any]:
        """Return list of tags."""
        result = {"Tags": [], "HasMoreTags": False}
        delivery_stream = self.delivery_streams.get(delivery_stream_name)
        if not delivery_stream:
            raise ResourceNotFoundException(
                f"Firehose {delivery_stream_name} under account {self.account_id} not found."
            )

        tags = self.tagger.list_tags_for_resource(delivery_stream.delivery_stream_arn)[
            "Tags"
        ]
        keys = self.tagger.extract_tag_names(tags)

        # If a starting tag is given and can be found, find the index into
        # tags, then add one to get the tag following it.
        start = 0
        if exclusive_start_tag_key:
            if exclusive_start_tag_key in keys:
                start = keys.index(exclusive_start_tag_key) + 1

        limit = limit or MAX_TAGS_PER_DELIVERY_STREAM
        result["Tags"] = tags[start : start + limit]
        if len(tags) > (start + limit):
            result["HasMoreTags"] = True
        return result

    def put_record(
        self, delivery_stream_name: str, record: Dict[str, bytes]
    ) -> Dict[str, Any]:
        """Write a single data record into a Kinesis Data firehose stream."""
        result = self.put_record_batch(delivery_stream_name, [record])
        return {
            "RecordId": result["RequestResponses"][0]["RecordId"],
            "Encrypted": False,
        }

    @staticmethod
    def put_http_records(  # type: ignore[misc]
        http_destination: Dict[str, Any], records: List[Dict[str, bytes]]
    ) -> List[Dict[str, str]]:
        """Put records to a HTTP destination."""
        # Mostly copied from localstack
        url = http_destination["EndpointConfiguration"]["Url"]
        headers = {"Content-Type": "application/json"}
        record_to_send = {
            "requestId": str(mock_random.uuid4()),
            "timestamp": int(time()),
            "records": [{"data": record["Data"]} for record in records],
        }
        try:
            requests.post(url, json=record_to_send, headers=headers)
        except Exception as exc:
            # This could be better ...
            raise RuntimeError(
                "Firehose PutRecord(Batch) to HTTP destination failed"
            ) from exc
        return [{"RecordId": str(mock_random.uuid4())} for _ in range(len(records))]

    @staticmethod
    def _format_s3_object_path(
        delivery_stream_name: str, version_id: str, prefix: str
    ) -> str:
        """Return a S3 object path in the expected format."""
        # Taken from LocalStack's firehose logic, with minor changes.
        # See https://docs.aws.amazon.com/firehose/latest/dev/basic-deliver.html#s3-object-name
        # Path prefix pattern: myApp/YYYY/MM/DD/HH/
        # Object name pattern:
        # DeliveryStreamName-DeliveryStreamVersion-YYYY-MM-DD-HH-MM-SS-RandomString
        prefix = f"{prefix}{'' if prefix.endswith('/') or prefix == '' else '/'}"
        now = utcnow()
        return (
            f"{prefix}{now.strftime('%Y/%m/%d/%H')}/"
            f"{delivery_stream_name}-{version_id}-"
            f"{now.strftime('%Y-%m-%d-%H-%M-%S')}-{str(mock_random.uuid4())}"
        )

    def put_s3_records(
        self,
        delivery_stream_name: str,
        version_id: str,
        s3_destination: Dict[str, Any],
        records: List[Dict[str, bytes]],
    ) -> List[Dict[str, str]]:
        """Put records to a ExtendedS3 or S3 destination."""
        # Taken from LocalStack's firehose logic, with minor changes.
        bucket_name = s3_destination["BucketARN"].split(":")[-1]
        prefix = s3_destination.get("Prefix", "")
        object_path = self._format_s3_object_path(
            delivery_stream_name, version_id, prefix
        )

        batched_data = b"".join([b64decode(r["Data"]) for r in records])
        try:
            s3_backends[self.account_id][self.partition].put_object(
                bucket_name, object_path, batched_data
            )
        except Exception as exc:
            # This could be better ...
            raise RuntimeError(
                "Firehose PutRecord(Batch to S3 destination failed"
            ) from exc
        return [{"RecordId": str(mock_random.uuid4())} for _ in range(len(records))]

    def put_record_batch(
        self, delivery_stream_name: str, records: List[Dict[str, bytes]]
    ) -> Dict[str, Any]:
        """Write multiple data records into a Kinesis Data firehose stream."""
        delivery_stream = self.delivery_streams.get(delivery_stream_name)
        if not delivery_stream:
            raise ResourceNotFoundException(
                f"Firehose {delivery_stream_name} under account {self.account_id} "
                f"not found."
            )

        request_responses = []
        for destination in delivery_stream.destinations:
            if "ExtendedS3" in destination:
                # ExtendedS3 will be handled like S3,but in the future
                # this will probably need to be revisited.  This destination
                # must be listed before S3 otherwise both destinations will
                # be processed instead of just ExtendedS3.
                request_responses = self.put_s3_records(
                    delivery_stream_name,
                    delivery_stream.version_id,
                    destination["ExtendedS3"],
                    records,
                )
            elif "S3" in destination:
                request_responses = self.put_s3_records(
                    delivery_stream_name,
                    delivery_stream.version_id,
                    destination["S3"],
                    records,
                )
            elif "HttpEndpoint" in destination:
                request_responses = self.put_http_records(
                    destination["HttpEndpoint"], records
                )
            elif {"Elasticsearch", "Redshift", "Snowflake"} & set(destination):
                # This isn't implemented as these services aren't implemented,
                # so ignore the data, but return a "proper" response.
                request_responses = [
                    {"RecordId": str(mock_random.uuid4())} for _ in range(len(records))
                ]

        return {
            "FailedPutCount": 0,
            "Encrypted": False,
            "RequestResponses": request_responses,
        }

    def tag_delivery_stream(
        self, delivery_stream_name: str, tags: List[Dict[str, str]]
    ) -> None:
        """Add/update tags for specified delivery stream."""
        delivery_stream = self.delivery_streams.get(delivery_stream_name)
        if not delivery_stream:
            raise ResourceNotFoundException(
                f"Firehose {delivery_stream_name} under account {self.account_id} "
                f"not found."
            )

        if len(tags) > MAX_TAGS_PER_DELIVERY_STREAM:
            raise ValidationException(
                f"1 validation error detected: Value '{tags}' at 'tags' "
                f"failed to satisify contstraint: Member must have length "
                f"less than or equal to {MAX_TAGS_PER_DELIVERY_STREAM}"
            )

        errmsg = self.tagger.validate_tags(tags)
        if errmsg:
            raise ValidationException(errmsg)

        self.tagger.tag_resource(delivery_stream.delivery_stream_arn, tags)

    def untag_delivery_stream(
        self, delivery_stream_name: str, tag_keys: List[str]
    ) -> None:
        """Removes tags from specified delivery stream."""
        delivery_stream = self.delivery_streams.get(delivery_stream_name)
        if not delivery_stream:
            raise ResourceNotFoundException(
                f"Firehose {delivery_stream_name} under account {self.account_id} not found."
            )

        # If a tag key doesn't exist for the stream, boto3 ignores it.
        self.tagger.untag_resource_using_names(
            delivery_stream.delivery_stream_arn, tag_keys
        )

    def start_delivery_stream_encryption(
        self, stream_name: str, encryption_config: Dict[str, Any]
    ) -> None:
        delivery_stream = self.delivery_streams.get(stream_name)
        if not delivery_stream:
            raise ResourceNotFoundException(
                f"Firehose {stream_name} under account {self.account_id} not found."
            )

        delivery_stream.delivery_stream_encryption_configuration = encryption_config
        delivery_stream.delivery_stream_encryption_configuration["Status"] = "ENABLED"

    def stop_delivery_stream_encryption(self, stream_name: str) -> None:
        delivery_stream = self.delivery_streams.get(stream_name)
        if not delivery_stream:
            raise ResourceNotFoundException(
                f"Firehose {stream_name} under account {self.account_id} not found."
            )

        delivery_stream.delivery_stream_encryption_configuration["Status"] = "DISABLED"

    def update_destination(  # pylint: disable=unused-argument
        self,
        delivery_stream_name: str,
        current_delivery_stream_version_id: str,
        destination_id: str,
        s3_destination_update: Dict[str, Any],
        extended_s3_destination_update: Dict[str, Any],
        s3_backup_mode: str,
        redshift_destination_update: Dict[str, Any],
        elasticsearch_destination_update: Dict[str, Any],
        splunk_destination_update: Dict[str, Any],
        http_endpoint_destination_update: Dict[str, Any],
        snowflake_destination_configuration: Dict[str, Any],
    ) -> None:
        (dest_name, dest_config) = find_destination_config_in_args(locals())

        delivery_stream = self.delivery_streams.get(delivery_stream_name)
        if not delivery_stream:
            raise ResourceNotFoundException(
                f"Firehose {delivery_stream_name} under accountId {self.account_id} not found."
            )

        if dest_name == "Splunk":
            warnings.warn("A Splunk destination delivery stream is not yet implemented")

        if delivery_stream.version_id != current_delivery_stream_version_id:
            raise ConcurrentModificationException(
                f"Cannot update firehose: {delivery_stream_name} since the "
                f"current version id: {delivery_stream.version_id} and "
                f"specified version id: {current_delivery_stream_version_id} "
                f"do not match"
            )

        destination = {}
        destination_idx = 0
        for destination in delivery_stream.destinations:
            if destination["destination_id"] == destination_id:
                break
            destination_idx += 1
        else:
            raise InvalidArgumentException("Destination Id {destination_id} not found")

        # Switching between Amazon ES and other services is not supported.
        # For an Amazon ES destination, you can only update to another Amazon
        # ES destination.  Same with HTTP.  Didn't test Splunk.
        if (dest_name == "Elasticsearch" and "Elasticsearch" not in destination) or (
            dest_name == "HttpEndpoint" and "HttpEndpoint" not in destination
        ):
            raise InvalidArgumentException(
                f"Changing the destination type to or from {dest_name} "
                f"is not supported at this time."
            )

        # If this is a different type of destination configuration,
        # the existing configuration is reset first.
        if dest_name in destination:
            delivery_stream.destinations[destination_idx][dest_name].update(dest_config)
        else:
            delivery_stream.destinations[destination_idx] = {
                "destination_id": destination_id,
                dest_name: dest_config,
            }

        # Some names in the Update-request differ from the names used by Describe
        for old, new in DeliveryStream.ALTERNATIVE_FIELD_NAMES:
            if old in dest_config:
                if new not in delivery_stream.destinations[destination_idx][dest_name]:
                    delivery_stream.destinations[destination_idx][dest_name][new] = {}
                delivery_stream.destinations[destination_idx][dest_name][new].update(
                    dest_config[old]
                )

        # Once S3 is updated to an ExtendedS3 destination, both remain in
        # the destination.  That means when one is updated, the other needs
        # to be updated as well.  The problem is that they don't have the
        # same fields.
        if dest_name == "ExtendedS3":
            delivery_stream.destinations[destination_idx]["S3"] = (
                create_s3_destination_config(dest_config)
            )
        elif dest_name == "S3" and "ExtendedS3" in destination:
            destination["ExtendedS3"] = {
                k: v
                for k, v in destination["S3"].items()
                if k in destination["ExtendedS3"]
            }

        # Increment version number and update the timestamp.
        delivery_stream.version_id = str(int(current_delivery_stream_version_id) + 1)
        delivery_stream.last_update_timestamp = datetime.now(timezone.utc).isoformat()

        # Unimplemented: processing of the "S3BackupMode" parameter.  Per the
        # documentation:  "You can update a delivery stream to enable Amazon
        # S3 backup if it is disabled.  If backup is enabled, you can't update
        # the delivery stream to disable it."

    def lookup_name_from_arn(self, arn: str) -> Optional[DeliveryStream]:
        """Given an ARN, return the associated delivery stream name."""
        return self.delivery_streams.get(arn.split("/")[-1])

    def send_log_event(
        self,
        delivery_stream_arn: str,
        filter_name: str,
        log_group_name: str,
        log_stream_name: str,
        log_events: List[Dict[str, Any]],
    ) -> None:
        """Send log events to a S3 bucket after encoding and gzipping it."""
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
        gzipped_payload = b64encode(output.getvalue())

        delivery_stream: DeliveryStream = self.lookup_name_from_arn(delivery_stream_arn)  # type: ignore[assignment]
        self.put_s3_records(
            delivery_stream.delivery_stream_name,
            delivery_stream.version_id,
            delivery_stream.destinations[0]["S3"],
            [{"Data": gzipped_payload}],
        )


firehose_backends = BackendDict(FirehoseBackend, "firehose")
