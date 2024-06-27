import base64
import hashlib
import json
import re
import string
import struct
from copy import deepcopy
from threading import Condition
from typing import TYPE_CHECKING, Any, Dict, List, Optional, Set, Tuple
from urllib.parse import ParseResult
from xml.sax.saxutils import escape

from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.common_models import BaseModel, CloudFormationModel
from moto.core.exceptions import RESTError
from moto.core.utils import (
    camelcase_to_underscores,
    tags_from_cloudformation_tags_list,
    unix_time,
    unix_time_millis,
)
from moto.moto_api._internal import mock_random as random
from moto.utilities.utils import get_partition, md5_hash

from .constants import MAXIMUM_VISIBILITY_TIMEOUT
from .exceptions import (
    BatchEntryIdsNotDistinct,
    BatchRequestTooLong,
    InvalidAttributeName,
    InvalidAttributeValue,
    InvalidBatchEntryId,
    InvalidParameterValue,
    MessageAttributesInvalid,
    MissingParameter,
    OverLimit,
    QueueAlreadyExists,
    QueueDoesNotExist,
    ReceiptHandleIsInvalid,
    TooManyEntriesInBatchRequest,
)
from .utils import generate_receipt_handle

if TYPE_CHECKING:
    from moto.awslambda.models import EventSourceMapping

DEFAULT_SENDER_ID = "AIDAIT2UOQQY3AUEKVGXU"

MAXIMUM_MESSAGE_LENGTH = 262144  # 256 KiB

MAXIMUM_MESSAGE_SIZE_ATTR_LOWER_BOUND = 1024
MAXIMUM_MESSAGE_SIZE_ATTR_UPPER_BOUND = MAXIMUM_MESSAGE_LENGTH

MAXIMUM_MESSAGE_DELAY = 900

TRANSPORT_TYPE_ENCODINGS = {
    "String": b"\x01",
    "Binary": b"\x02",
    "Number": b"\x01",
    "String.custom": b"\x01",
}

STRING_TYPE_FIELD_INDEX = 1
BINARY_TYPE_FIELD_INDEX = 2
STRING_LIST_TYPE_FIELD_INDEX = 3
BINARY_LIST_TYPE_FIELD_INDEX = 4

# Valid attribute name rules can found at
# https://docs.aws.amazon.com/AWSSimpleQueueService/latest/SQSDeveloperGuide/sqs-message-metadata.html
ATTRIBUTE_NAME_PATTERN = re.compile("^([a-z]|[A-Z]|[0-9]|[_.\\-])+$")

DEDUPLICATION_TIME_IN_SECONDS = 300


class Message(BaseModel):
    def __init__(
        self,
        message_id: str,
        body: str,
        system_attributes: Optional[Dict[str, Any]] = None,
    ):
        self.id = message_id
        self._body = body
        self.message_attributes: Dict[str, Any] = {}
        self.receipt_handle: Optional[str] = None
        self._old_receipt_handles: List[str] = []
        self.sender_id = DEFAULT_SENDER_ID
        self.sent_timestamp = None
        self.approximate_first_receive_timestamp: Optional[int] = None
        self.approximate_receive_count = 0
        self.deduplication_id: Optional[str] = None
        self.group_id: Optional[str] = None
        self.sequence_number: Optional[str] = None
        self.visible_at = 0.0
        self.delayed_until = 0.0
        self.system_attributes = system_attributes or {}

    @property
    def body_md5(self) -> str:
        md5 = md5_hash()
        md5.update(self._body.encode("utf-8"))
        return md5.hexdigest()

    @property
    def attribute_md5(self) -> str:
        md5 = md5_hash()

        for attrName in sorted(self.message_attributes.keys()):
            self.validate_attribute_name(attrName)
            attrValue = self.message_attributes[attrName]
            # Encode name
            self.update_binary_length_and_value(md5, self.utf8(attrName))
            # Encode type
            self.update_binary_length_and_value(md5, self.utf8(attrValue["data_type"]))

            if attrValue.get("string_value"):
                md5.update(bytearray([STRING_TYPE_FIELD_INDEX]))
                self.update_binary_length_and_value(
                    md5, self.utf8(attrValue.get("string_value"))
                )
            elif attrValue.get("binary_value"):
                md5.update(bytearray([BINARY_TYPE_FIELD_INDEX]))
                decoded_binary_value = base64.b64decode(attrValue.get("binary_value"))
                self.update_binary_length_and_value(md5, decoded_binary_value)
            # string_list_value type is not implemented, reserved for the future use.
            # See https://docs.aws.amazon.com/AWSSimpleQueueService/latest/APIReference/API_MessageAttributeValue.html
            elif len(attrValue["string_list_value"]) > 0:
                md5.update(bytearray([STRING_LIST_TYPE_FIELD_INDEX]))
                for strListMember in attrValue["string_list_value"]:
                    self.update_binary_length_and_value(md5, self.utf8(strListMember))
            # binary_list_value type is not implemented, reserved for the future use.
            # See https://docs.aws.amazon.com/AWSSimpleQueueService/latest/APIReference/API_MessageAttributeValue.html
            elif len(attrValue["binary_list_value"]) > 0:
                md5.update(bytearray([BINARY_LIST_TYPE_FIELD_INDEX]))
                for strListMember in attrValue["binary_list_value"]:
                    decoded_binary_value = base64.b64decode(strListMember)
                    self.update_binary_length_and_value(md5, decoded_binary_value)

        return md5.hexdigest()

    @staticmethod
    def update_binary_length_and_value(md5: Any, value: bytes) -> None:  # type: ignore[misc]
        length_bytes = struct.pack("!I".encode("ascii"), len(value))
        md5.update(length_bytes)
        md5.update(value)

    @staticmethod
    def validate_attribute_name(name: str) -> None:
        if not ATTRIBUTE_NAME_PATTERN.match(name):
            raise MessageAttributesInvalid(
                f"The message attribute name '{name}' is invalid. "
                "Attribute name can contain A-Z, a-z, 0-9, "
                "underscore (_), hyphen (-), and period (.) characters."
            )

    @staticmethod
    def utf8(value: Any) -> bytes:  # type: ignore[misc]
        if isinstance(value, str):
            return value.encode("utf-8")
        return value

    @property
    def body(self) -> str:
        return escape(self._body).replace('"', "&quot;").replace("\r", "&#xD;")

    @property
    def original_body(self) -> str:
        return self._body

    def mark_sent(self, delay_seconds: Optional[int] = None) -> None:
        self.sent_timestamp = int(unix_time_millis())  # type: ignore
        if delay_seconds:
            self.delay(delay_seconds=delay_seconds)

    def mark_received(self, visibility_timeout: Optional[int] = None) -> None:
        """
        When a message is received we will set the first receive timestamp,
        tap the ``approximate_receive_count`` and the ``visible_at`` time.
        """
        if visibility_timeout:
            visibility_timeout = int(visibility_timeout)
        else:
            visibility_timeout = 0

        if not self.approximate_first_receive_timestamp:
            self.approximate_first_receive_timestamp = int(unix_time_millis())

        self.approximate_receive_count += 1

        # Make message visible again in the future unless its
        # destroyed.
        if visibility_timeout:
            self.change_visibility(visibility_timeout)

        self._old_receipt_handles.append(self.receipt_handle)  # type: ignore
        self.receipt_handle = generate_receipt_handle()

    def change_visibility(self, visibility_timeout: int) -> None:
        # We're dealing with milliseconds internally
        visibility_timeout_msec = int(visibility_timeout) * 1000
        self.visible_at = unix_time_millis() + visibility_timeout_msec

    def delay(self, delay_seconds: int) -> None:
        delay_msec = int(delay_seconds) * 1000
        self.delayed_until = unix_time_millis() + delay_msec

    @property
    def visible(self) -> bool:
        current_time = unix_time_millis()
        if current_time > self.visible_at:
            return True
        return False

    @property
    def delayed(self) -> bool:
        current_time = unix_time_millis()
        if current_time < self.delayed_until:
            return True
        return False

    @property
    def all_receipt_handles(self) -> List[Optional[str]]:
        return [self.receipt_handle] + self._old_receipt_handles  # type: ignore

    def had_receipt_handle(self, receipt_handle: str) -> bool:
        """
        Check if this message ever had this receipt_handle in the past
        """
        return receipt_handle in self.all_receipt_handles


class Queue(CloudFormationModel):
    BASE_ATTRIBUTES = [
        "ApproximateNumberOfMessages",
        "ApproximateNumberOfMessagesDelayed",
        "ApproximateNumberOfMessagesNotVisible",
        "CreatedTimestamp",
        "DelaySeconds",
        "LastModifiedTimestamp",
        "MaximumMessageSize",
        "MessageRetentionPeriod",
        "QueueArn",
        "Policy",
        "RedrivePolicy",
        "ReceiveMessageWaitTimeSeconds",
        "VisibilityTimeout",
        "SqsManagedSseEnabled",
    ]
    FIFO_ATTRIBUTES = [
        "ContentBasedDeduplication",
        "DeduplicationScope",
        "FifoQueue",
        "FifoThroughputLimit",
    ]
    KMS_ATTRIBUTES = ["KmsDataKeyReusePeriodSeconds", "KmsMasterKeyId"]
    ALLOWED_PERMISSIONS = (
        "*",
        "ChangeMessageVisibility",
        "DeleteMessage",
        "GetQueueAttributes",
        "GetQueueUrl",
        "ListDeadLetterSourceQueues",
        "PurgeQueue",
        "ReceiveMessage",
        "SendMessage",
    )

    def __init__(self, name: str, region: str, account_id: str, **kwargs: Any):
        self.name = name
        self.region = region
        self.account_id = account_id
        self.tags: Dict[str, str] = {}
        self.permissions: Dict[str, Any] = {}

        self._messages: List[Message] = []
        self._pending_messages: Set[Message] = set()
        self.deleted_messages: Set[str] = set()
        self._messages_lock = Condition()

        now = unix_time()
        self.created_timestamp = now
        self.queue_arn = f"arn:{get_partition(region)}:sqs:{region}:{account_id}:{name}"
        self.dead_letter_queue: Optional["Queue"] = None
        self.fifo_queue = False

        self.lambda_event_source_mappings: Dict[str, "EventSourceMapping"] = {}

        # default settings for a non fifo queue
        defaults = {
            "ContentBasedDeduplication": "false",
            "DeduplicationScope": "queue",
            "DelaySeconds": 0,
            "FifoQueue": "false",
            "FifoThroughputLimit": "perQueue",
            "KmsDataKeyReusePeriodSeconds": 300,  # five minutes
            "KmsMasterKeyId": None,
            "MaximumMessageSize": MAXIMUM_MESSAGE_LENGTH,
            "MessageRetentionPeriod": 86400 * 4,  # four days
            "Policy": None,
            "ReceiveMessageWaitTimeSeconds": 0,
            "RedrivePolicy": None,
            "VisibilityTimeout": 30,
            "SqsManagedSseEnabled": True,
        }

        defaults.update(kwargs)
        self._set_attributes(defaults, now)

        # Check some conditions
        if self.fifo_queue and not self.name.endswith(".fifo"):
            raise InvalidParameterValue("Queue name must end in .fifo for FIFO queues")
        if (
            self.maximum_message_size < MAXIMUM_MESSAGE_SIZE_ATTR_LOWER_BOUND  # type: ignore
            or self.maximum_message_size > MAXIMUM_MESSAGE_SIZE_ATTR_UPPER_BOUND  # type: ignore
        ):
            raise InvalidAttributeValue("MaximumMessageSize")

    @property
    def name(self) -> str:
        return self._name

    @name.setter
    def name(self, value: str) -> None:
        # https://docs.aws.amazon.com/AWSSimpleQueueService/latest/APIReference/API_CreateQueue.html#SQS-CreateQueue-request-QueueName
        if re.match(r"^[a-zA-Z0-9\_-]{1,80}$", value) or re.match(
            r"^[a-zA-Z0-9\_-]{1,75}(\.fifo)?$", value
        ):
            self._name = value
        else:
            raise InvalidParameterValue(
                "Can only include alphanumeric characters, hyphens, or underscores. 1 to 80 in length"
            )

    @property
    def pending_messages(self) -> Set[Message]:
        return self._pending_messages

    @property
    def pending_message_groups(self) -> Set[str]:
        return set(
            message.group_id
            for message in self._pending_messages
            if message.group_id is not None
        )

    def _set_attributes(
        self, attributes: Dict[str, Any], now: Optional[float] = None
    ) -> None:
        if not now:
            now = unix_time()

        integer_fields = (
            "DelaySeconds",
            "KmsDataKeyreusePeriodSeconds",
            "MaximumMessageSize",
            "MessageRetentionPeriod",
            "ReceiveMessageWaitTime",
            "VisibilityTimeout",
        )
        bool_fields = ("ContentBasedDeduplication", "FifoQueue")

        for key, value in attributes.items():
            if key in integer_fields:
                value = int(value)
            if key in bool_fields:
                value = str(value).lower() == "true"

            if key in ["Policy", "RedrivePolicy"] and value is not None:
                continue

            setattr(self, camelcase_to_underscores(key), value)

        if attributes.get("RedrivePolicy", None) is not None:
            self._setup_dlq(attributes["RedrivePolicy"])

        self.policy = attributes.get("Policy")

        self.last_modified_timestamp = now

    @staticmethod
    def _is_empty_redrive_policy(policy: Any) -> bool:  # type: ignore[misc]
        if isinstance(policy, str):
            if policy == "" or len(json.loads(policy)) == 0:
                return True
        elif isinstance(policy, dict) and len(policy) == 0:
            return True

        return False

    def _setup_dlq(self, policy: Any) -> None:
        if Queue._is_empty_redrive_policy(policy):
            self.redrive_policy = None
            self.dead_letter_queue = None
            return

        if isinstance(policy, str):
            try:
                self.redrive_policy = json.loads(policy)
            except ValueError:
                raise RESTError(
                    "InvalidParameterValue",
                    "Redrive policy is not a dict or valid json",
                )
        elif isinstance(policy, dict):
            self.redrive_policy = policy
        else:
            raise RESTError(
                "InvalidParameterValue", "Redrive policy is not a dict or valid json"
            )

        if "deadLetterTargetArn" not in self.redrive_policy:
            raise RESTError(
                "InvalidParameterValue",
                "Redrive policy does not contain deadLetterTargetArn",
            )
        if "maxReceiveCount" not in self.redrive_policy:
            raise RESTError(
                "InvalidParameterValue",
                "Redrive policy does not contain maxReceiveCount",
            )

        # 'maxReceiveCount' is stored as int
        self.redrive_policy["maxReceiveCount"] = int(
            self.redrive_policy["maxReceiveCount"]
        )

        sqs_backend = sqs_backends[self.account_id][self.region]
        for queue in sqs_backend.queues.values():
            if queue.queue_arn == self.redrive_policy["deadLetterTargetArn"]:
                self.dead_letter_queue = queue

                if self.fifo_queue and not queue.fifo_queue:
                    raise RESTError(
                        "InvalidParameterCombination",
                        "Fifo queues cannot use non fifo dead letter queues",
                    )
                break
        else:
            raise RESTError(
                "AWS.SimpleQueueService.NonExistentQueue",
                f"Could not find DLQ for {self.redrive_policy['deadLetterTargetArn']}",
            )

    @staticmethod
    def cloudformation_name_type() -> str:
        return "QueueName"

    @staticmethod
    def cloudformation_type() -> str:
        # https://docs.aws.amazon.com/AWSCloudFormation/latest/UserGuide/aws-resource-sqs-queue.html
        return "AWS::SQS::Queue"

    @classmethod
    def create_from_cloudformation_json(  # type: ignore[misc]
        cls,
        resource_name: str,
        cloudformation_json: Any,
        account_id: str,
        region_name: str,
        **kwargs: Any,
    ) -> "Queue":
        properties = deepcopy(cloudformation_json["Properties"])
        # remove Tags from properties and convert tags list to dict
        tags = properties.pop("Tags", [])
        tags_dict = tags_from_cloudformation_tags_list(tags)

        # Could be passed as an integer - just treat it as a string
        # take first 80 characters of the q name, more is invalid
        # https://docs.aws.amazon.com/AWSSimpleQueueService/latest/APIReference/API_CreateQueue.html#SQS-CreateQueue-request-QueueName
        resource_name = str(resource_name)[0:80]

        sqs_backend = sqs_backends[account_id][region_name]
        return sqs_backend.create_queue(
            name=resource_name, tags=tags_dict, region=region_name, **properties
        )

    @classmethod
    def update_from_cloudformation_json(  # type: ignore[misc]
        cls,
        original_resource: Any,
        new_resource_name: str,
        cloudformation_json: Any,
        account_id: str,
        region_name: str,
    ) -> "Queue":
        properties = cloudformation_json["Properties"]
        queue_name = original_resource.name

        sqs_backend = sqs_backends[account_id][region_name]
        queue = sqs_backend.get_queue(queue_name)
        if "VisibilityTimeout" in properties:
            queue.visibility_timeout = int(properties["VisibilityTimeout"])  # type: ignore[attr-defined]

        if "ReceiveMessageWaitTimeSeconds" in properties:
            queue.receive_message_wait_time_seconds = int(  # type: ignore[attr-defined]
                properties["ReceiveMessageWaitTimeSeconds"]
            )
        return queue

    @classmethod
    def delete_from_cloudformation_json(  # type: ignore[misc]
        cls,
        resource_name: str,
        cloudformation_json: Any,
        account_id: str,
        region_name: str,
    ) -> None:
        # ResourceName will be the full queue URL - we only need the name
        # https://sqs.us-west-1.amazonaws.com/123456789012/queue_name
        queue_name = resource_name.split("/")[-1]
        sqs_backend = sqs_backends[account_id][region_name]
        sqs_backend.delete_queue(queue_name)

    @property
    def approximate_number_of_messages_delayed(self) -> int:
        return len([m for m in self._messages if m.delayed])

    @property
    def approximate_number_of_messages_not_visible(self) -> int:
        return len([m for m in self._messages if not m.visible])

    @property
    def approximate_number_of_messages(self) -> int:
        return len(self.messages)

    @property
    def physical_resource_id(self) -> str:
        return f"https://sqs.{self.region}.amazonaws.com/{self.account_id}/{self.name}"

    @property
    def attributes(self) -> Dict[str, Any]:  # type: ignore[misc]
        result: Dict[str, Any] = {}

        for attribute in self.BASE_ATTRIBUTES:
            attr = getattr(self, camelcase_to_underscores(attribute))
            result[attribute] = attr

        if self.fifo_queue:
            for attribute in self.FIFO_ATTRIBUTES:
                attr = getattr(self, camelcase_to_underscores(attribute))
                result[attribute] = attr

        if self.kms_master_key_id:  # type: ignore
            for attribute in self.KMS_ATTRIBUTES:
                attr = getattr(self, camelcase_to_underscores(attribute))
                result[attribute] = attr

        if self.policy:
            result["Policy"] = self.policy

        if self.redrive_policy:
            result["RedrivePolicy"] = json.dumps(self.redrive_policy)

        for key in result:
            if isinstance(result[key], bool):
                result[key] = str(result[key]).lower()

        return result

    def url(self, request_url: ParseResult) -> str:
        return (
            f"{request_url.scheme}://{request_url.netloc}/{self.account_id}/{self.name}"
        )

    @property
    def messages(self) -> List[Message]:
        # TODO: This can become very inefficient if a large number of messages are in-flight
        return [
            message
            for message in self._messages
            if message.visible and not message.delayed
        ]

    def add_message(self, message: Message) -> None:
        if self.fifo_queue:
            # the cases in which we dedupe fifo messages
            # from https://docs.aws.amazon.com/AWSSimpleQueueService/latest/SQSDeveloperGuide/using-messagededuplicationid-property.html
            # https://docs.aws.amazon.com/AWSSimpleQueueService/latest/APIReference/API_SendMessage.html
            if (
                self.attributes.get("ContentBasedDeduplication") == "true"
                or message.deduplication_id
            ):
                for m in self._messages:
                    if m.deduplication_id == message.deduplication_id:
                        diff = message.sent_timestamp - m.sent_timestamp  # type: ignore
                        # if a duplicate message is received within the deduplication time then it should
                        # not be added to the queue
                        if diff / 1000 < DEDUPLICATION_TIME_IN_SECONDS:
                            return

        with self._messages_lock:
            self._messages.append(message)
            self._messages_lock.notify_all()

        for arn, esm in self.lambda_event_source_mappings.items():
            backend = sqs_backends[self.account_id][self.region]

            """
            Lambda polls the queue and invokes your function synchronously with an event
            that contains queue messages. Lambda reads messages in batches and invokes
            your function once for each batch. When your function successfully processes
            a batch, Lambda deletes its messages from the queue.
            """
            messages = backend.receive_message(
                self.name,
                esm.batch_size,
                self.receive_message_wait_time_seconds,  # type: ignore
                self.visibility_timeout,  # type: ignore
            )

            from moto.awslambda.utils import get_backend

            result = get_backend(self.account_id, self.region).send_sqs_batch(
                arn, messages, self.queue_arn
            )

            if result:
                for m in messages:
                    backend.delete_message(self.name, m.receipt_handle)  # type: ignore[arg-type]
            else:
                # Make messages visible again
                for m in messages:
                    backend.change_message_visibility(
                        self.name,
                        m.receipt_handle,  # type: ignore[arg-type]
                        visibility_timeout=0,
                    )

    def delete_message(self, receipt_handle: str) -> None:
        if receipt_handle in self.deleted_messages:
            # Already deleted - gracefully handle deleting it again
            return

        if not any(
            message.had_receipt_handle(receipt_handle) for message in self._messages
        ):
            raise ReceiptHandleIsInvalid()

        # Delete message from queue regardless of pending state
        new_messages = []
        for message in self._messages:
            if message.had_receipt_handle(receipt_handle):
                self.pending_messages.discard(message)
                self.deleted_messages.update(message.all_receipt_handles)  # type: ignore
                continue
            new_messages.append(message)
        self._messages = new_messages

    def wait_for_messages(self, timeout: int) -> None:
        with self._messages_lock:
            self._messages_lock.wait_for(lambda: self.messages, timeout=timeout)

    @classmethod
    def has_cfn_attr(cls, attr: str) -> bool:
        return attr in ["Arn", "QueueName"]

    def get_cfn_attribute(self, attribute_name: str) -> str:
        from moto.cloudformation.exceptions import UnformattedGetAttTemplateException

        if attribute_name == "Arn":
            return self.queue_arn
        elif attribute_name == "QueueName":
            return self.name
        raise UnformattedGetAttTemplateException()

    @property
    def policy(self) -> Any:  # type: ignore[misc]
        if self._policy_json.get("Statement"):
            return json.dumps(self._policy_json)
        else:
            return None

    @policy.setter
    def policy(self, policy: Any) -> None:
        if policy:
            self._policy_json = json.loads(policy)
        else:
            self._policy_json = {
                "Version": "2012-10-17",
                "Id": f"{self.queue_arn}/SQSDefaultPolicy",
                "Statement": [],
            }


def _filter_message_attributes(
    message: Message, input_message_attributes: List[str]
) -> None:
    filtered_message_attributes = {}
    return_all = "All" in input_message_attributes
    for key, value in message.message_attributes.items():
        if return_all or key in input_message_attributes:
            filtered_message_attributes[key] = value
    message.message_attributes = filtered_message_attributes


class SQSBackend(BaseBackend):
    def __init__(self, region_name: str, account_id: str):
        super().__init__(region_name, account_id)
        self.queues: Dict[str, Queue] = {}

    def create_queue(
        self, name: str, tags: Optional[Dict[str, str]] = None, **kwargs: Any
    ) -> Queue:
        queue = self.queues.get(name)
        if queue:
            try:
                kwargs.pop("region")
            except KeyError:
                pass

            new_queue = Queue(
                name, region=self.region_name, account_id=self.account_id, **kwargs
            )

            queue_attributes = queue.attributes
            new_queue_attributes = new_queue.attributes

            # only the attributes which are being sent for the queue
            # creation have to be compared if the queue is existing.
            for key in kwargs:
                if queue_attributes.get(key) != new_queue_attributes.get(key):
                    raise QueueAlreadyExists("The specified queue already exists.")
        else:
            try:
                kwargs.pop("region")
            except KeyError:
                pass
            queue = Queue(
                name, region=self.region_name, account_id=self.account_id, **kwargs
            )
            self.queues[name] = queue

        if tags:
            queue.tags = tags

        return queue

    def get_queue_url(self, queue_name: str) -> Queue:
        return self.get_queue(queue_name)

    def list_queues(self, queue_name_prefix: str) -> List[Queue]:
        re_str = ".*"
        if queue_name_prefix:
            re_str = f"^{queue_name_prefix}.*"
        prefix_re = re.compile(re_str)
        qs = []
        for name, q in self.queues.items():
            if prefix_re.search(name):
                qs.append(q)
        return qs[:1000]

    def get_queue(self, queue_name: str) -> Queue:
        queue = self.queues.get(queue_name)
        if queue is None:
            raise QueueDoesNotExist()
        return queue

    def delete_queue(self, queue_name: str) -> None:
        self.get_queue(queue_name)

        del self.queues[queue_name]

    def get_queue_attributes(
        self, queue_name: str, attribute_names: List[str]
    ) -> Dict[str, Any]:
        queue = self.get_queue(queue_name)
        if not attribute_names:
            return {}

        valid_names = (
            ["All"]
            + queue.BASE_ATTRIBUTES
            + queue.FIFO_ATTRIBUTES
            + queue.KMS_ATTRIBUTES
        )
        invalid_name = next(
            (name for name in attribute_names if name not in valid_names), None
        )

        if invalid_name or invalid_name == "":
            raise InvalidAttributeName(invalid_name)

        attributes = {}

        if "All" in attribute_names:
            attributes = queue.attributes
        else:
            for name in (name for name in attribute_names if name in queue.attributes):
                if queue.attributes.get(name) is not None:
                    attributes[name] = queue.attributes.get(name)

        return attributes

    def set_queue_attributes(
        self, queue_name: str, attributes: Dict[str, Any]
    ) -> Queue:
        queue = self.get_queue(queue_name)
        queue._set_attributes(attributes)
        return queue

    def _validate_message(
        self,
        queue: Queue,
        message_body: str,
        delay_seconds: int,
        deduplication_id: Optional[str] = None,
        group_id: Optional[str] = None,
    ) -> None:
        if queue.fifo_queue:
            if (
                queue.attributes.get("ContentBasedDeduplication") == "false"
                and not group_id
            ):
                msg = "MessageGroupId"
                raise MissingParameter(msg)

            if (
                queue.attributes.get("ContentBasedDeduplication") == "false"
                and group_id
                and not deduplication_id
            ):
                msg = (
                    "The queue should either have ContentBasedDeduplication enabled or "
                    "MessageDeduplicationId provided explicitly"
                )
                raise InvalidParameterValue(msg)

            if delay_seconds > 0:
                raise InvalidParameterValue(
                    f"Value {delay_seconds} for parameter DelaySeconds is invalid. Reason: "
                    "The request include parameter that is not valid for this queue type."
                )

        if len(message_body) > queue.maximum_message_size:  # type: ignore
            msg = f"One or more parameters are invalid. Reason: Message must be shorter than {queue.maximum_message_size} bytes."  # type: ignore
            raise InvalidParameterValue(msg)

        if group_id is None:
            if queue.fifo_queue:
                # MessageGroupId is a mandatory parameter for all
                # messages in a fifo queue
                raise MissingParameter("MessageGroupId")
        else:
            if not queue.fifo_queue:
                msg = (
                    f"Value {group_id} for parameter MessageGroupId is invalid. "
                    "Reason: The request include parameter that is not valid for this queue type."
                )
                raise InvalidParameterValue(msg)

    def send_message(
        self,
        queue_name: str,
        message_body: str,
        message_attributes: Optional[Dict[str, Any]] = None,
        delay_seconds: Optional[int] = None,
        deduplication_id: Optional[str] = None,
        group_id: Optional[str] = None,
        system_attributes: Optional[Dict[str, Any]] = None,
    ) -> Message:
        queue = self.get_queue(queue_name)

        self._validate_message(
            queue, message_body, int(delay_seconds or 0), deduplication_id, group_id
        )

        if delay_seconds is not None:
            delay_seconds = int(delay_seconds)
        else:
            delay_seconds = queue.delay_seconds  # type: ignore

        message_id = str(random.uuid4())
        message = Message(message_id, message_body, system_attributes)

        # if content based deduplication is set then set sha256 hash of the message
        # as the deduplication_id
        if queue.attributes.get("ContentBasedDeduplication") == "true":
            sha256 = hashlib.sha256()
            sha256.update(message_body.encode("utf-8"))
            message.deduplication_id = sha256.hexdigest()

        # Attributes, but not *message* attributes
        if deduplication_id is not None:
            message.deduplication_id = deduplication_id
            message.sequence_number = "".join(
                random.choice(string.digits) for _ in range(20)
            )

        if group_id is not None:
            message.group_id = group_id

        if message_attributes:
            message.message_attributes = message_attributes

        if delay_seconds > MAXIMUM_MESSAGE_DELAY:
            msg = (
                f"Value {delay_seconds} for parameter DelaySeconds is invalid. "
                "Reason: DelaySeconds must be >= 0 and <= 900."
            )
            raise InvalidParameterValue(msg)

        message.mark_sent(delay_seconds=delay_seconds)

        queue.add_message(message)

        return message

    def send_message_batch(
        self, queue_name: str, entries: Dict[str, Dict[str, Any]]
    ) -> Tuple[List[Message], List[Dict[str, Any]]]:
        queue = self.get_queue(queue_name)

        if any(
            not re.match(r"^[\w-]{1,80}$", entry["Id"]) for entry in entries.values()
        ):
            raise InvalidBatchEntryId()

        body_length = sum(len(entry["MessageBody"]) for entry in entries.values())
        if body_length > MAXIMUM_MESSAGE_LENGTH:
            raise BatchRequestTooLong(body_length)

        duplicate_id = self._get_first_duplicate_id(
            [entry["Id"] for entry in entries.values()]
        )
        if duplicate_id:
            raise BatchEntryIdsNotDistinct(duplicate_id)

        if len(entries) > 10:
            raise TooManyEntriesInBatchRequest(len(entries))

        messages = []
        failedInvalidDelay = []

        for entry in entries.values():
            # validate ALL messages before trying to send any
            self._validate_message(
                queue,
                entry["MessageBody"],
                int(entry.get("DelaySeconds") or 0),
                entry.get("MessageDeduplicationId"),
                entry.get("MessageGroupId"),
            )

        for entry in entries.values():
            try:
                # Loop through looking for messages
                message = self.send_message(
                    queue_name,
                    entry["MessageBody"],
                    message_attributes=entry["MessageAttributes"],
                    delay_seconds=entry["DelaySeconds"],
                    group_id=entry.get("MessageGroupId"),
                    deduplication_id=entry.get("MessageDeduplicationId"),
                )
                message.user_id = entry["Id"]  # type: ignore[attr-defined]
                messages.append(message)
            except InvalidParameterValue as err:
                if "DelaySeconds is invalid" in str(err):
                    failedInvalidDelay.append(entry)
                else:
                    raise err

        return messages, failedInvalidDelay

    def _get_first_duplicate_id(self, ids: List[str]) -> Optional[str]:
        unique_ids = set()
        for _id in ids:
            if _id in unique_ids:
                return _id
            unique_ids.add(_id)
        return None

    def receive_message(
        self,
        queue_name: str,
        count: int,
        wait_seconds_timeout: int,
        visibility_timeout: int,
        message_attribute_names: Optional[List[str]] = None,
    ) -> List[Message]:
        # Attempt to retrieve visible messages from a queue.

        # If a message was read by client and not deleted it is considered to be
        # "inflight" and cannot be read. We make attempts to obtain ``count``
        # messages but we may return less if messages are in-flight or there
        # are simple not enough messages in the queue.

        if message_attribute_names is None:
            message_attribute_names = []
        queue = self.get_queue(queue_name)
        result: List[Message] = []
        previous_result_count = len(result)

        polling_end = unix_time() + wait_seconds_timeout
        currently_pending_groups = deepcopy(queue.pending_message_groups)

        # queue.messages only contains visible messages
        while True:
            if result or (wait_seconds_timeout and unix_time() > polling_end):
                break

            messages_to_dlq: List[Message] = []

            for message in queue.messages:
                if not message.visible:
                    continue

                if message in queue.pending_messages:
                    # The message is pending but is visible again, so the
                    # consumer must have timed out.
                    queue.pending_messages.remove(message)
                    currently_pending_groups = deepcopy(queue.pending_message_groups)

                if message.group_id and queue.fifo_queue:
                    if message.group_id in currently_pending_groups:
                        # A previous call is still processing messages in this group, so we cannot deliver this one.
                        continue

                if (
                    queue.dead_letter_queue is not None
                    and queue.redrive_policy
                    and message.approximate_receive_count
                    >= queue.redrive_policy["maxReceiveCount"]
                ):
                    messages_to_dlq.append(message)
                    continue

                queue.pending_messages.add(message)
                message.mark_received(visibility_timeout=visibility_timeout)
                # Create deepcopy to not mutate the message state when filtering for attributes
                message_copy = deepcopy(message)
                _filter_message_attributes(message_copy, message_attribute_names)
                if not self.is_message_valid_based_on_retention_period(
                    queue_name, message
                ):
                    break
                result.append(message_copy)
                if len(result) >= count:
                    break

            for message in messages_to_dlq:
                queue._messages.remove(message)
                queue.dead_letter_queue.add_message(message)  # type: ignore

            if previous_result_count == len(result):
                if wait_seconds_timeout == 0:
                    # There is no timeout and no additional results,
                    # so break to avoid an infinite loop.
                    break

                queue.wait_for_messages(1)
                continue

            previous_result_count = len(result)

        return result

    def delete_message(self, queue_name: str, receipt_handle: str) -> None:
        queue = self.get_queue(queue_name)

        queue.delete_message(receipt_handle)

    def delete_message_batch(
        self, queue_name: str, receipts: List[Dict[str, Any]]
    ) -> Tuple[List[str], List[Dict[str, str]]]:
        success = []
        errors = []
        for receipt_and_id in receipts:
            try:
                self.delete_message(queue_name, receipt_and_id["receipt_handle"])
                success.append(receipt_and_id["msg_user_id"])
            except ReceiptHandleIsInvalid:
                errors.append(
                    {
                        "Id": receipt_and_id["msg_user_id"],
                        "SenderFault": True,
                        "Code": "ReceiptHandleIsInvalid",
                        "Message": f'The input receipt handle "{receipt_and_id["receipt_handle"]}" is not a valid receipt handle.',
                    }
                )
        return success, errors

    def change_message_visibility(
        self, queue_name: str, receipt_handle: str, visibility_timeout: int
    ) -> None:
        queue = self.get_queue(queue_name)
        for message in queue._messages:
            if message.had_receipt_handle(receipt_handle):
                visibility_timeout_msec = int(visibility_timeout) * 1000
                given_visibility_timeout = unix_time_millis() + visibility_timeout_msec
                if given_visibility_timeout - message.sent_timestamp > 43200 * 1000:  # type: ignore
                    raise InvalidParameterValue(
                        f"Value {visibility_timeout} for parameter VisibilityTimeout is invalid. Reason: Total "
                        "VisibilityTimeout for the message is beyond the limit [43200 seconds]"
                    )

                message.change_visibility(visibility_timeout)
                if message.visible and message in queue.pending_messages:
                    # If the message is visible again, remove it from pending
                    # messages.
                    queue.pending_messages.remove(message)
                return
        raise ReceiptHandleIsInvalid

    def change_message_visibility_batch(
        self, queue_name: str, entries: List[Dict[str, Any]]
    ) -> Tuple[List[str], List[Dict[str, str]]]:
        success = []
        error = []
        for entry in entries:
            try:
                visibility_timeout = int(entry["visibility_timeout"])
                assert visibility_timeout <= MAXIMUM_VISIBILITY_TIMEOUT
            except:  # noqa: E722 Do not use bare except
                error.append(
                    {
                        "Id": entry["id"],
                        "SenderFault": "true",
                        "Code": "InvalidParameterValue",
                        "Message": "Visibility timeout invalid",
                    }
                )
                continue

            try:
                self.change_message_visibility(
                    queue_name=queue_name,
                    receipt_handle=entry["receipt_handle"],
                    visibility_timeout=visibility_timeout,
                )
                success.append(entry["id"])
            except ReceiptHandleIsInvalid as e:
                error.append(
                    {
                        "Id": entry["id"],
                        "SenderFault": "true",
                        "Code": "ReceiptHandleIsInvalid",
                        "Message": e.description,
                    }
                )
        return success, error

    def purge_queue(self, queue_name: str) -> None:
        queue = self.get_queue(queue_name)
        queue._messages = []
        queue._pending_messages = set()

    def list_dead_letter_source_queues(self, queue_name: str) -> List[Queue]:
        dlq = self.get_queue(queue_name)

        queues: List[Queue] = []
        for queue in self.queues.values():
            if queue.dead_letter_queue is dlq:
                queues.append(queue)

        return queues

    def add_permission(
        self,
        region_name: str,
        queue_name: str,
        actions: List[str],
        account_ids: List[str],
        label: str,
    ) -> None:
        queue = self.get_queue(queue_name)

        if not actions:
            raise MissingParameter("Actions")

        if not account_ids:
            raise InvalidParameterValue(
                "Value [] for parameter PrincipalId is invalid. Reason: Unable to verify."
            )

        count = len(actions)
        if count > 7:
            raise OverLimit(count)

        invalid_action = next(
            (action for action in actions if action not in Queue.ALLOWED_PERMISSIONS),
            None,
        )
        if invalid_action:
            raise InvalidParameterValue(
                f"Value SQS:{invalid_action} for parameter ActionName is invalid. "
                "Reason: Only the queue owner is allowed to invoke this action."
            )

        policy = queue._policy_json
        statement = next(
            (
                statement
                for statement in policy["Statement"]
                if statement["Sid"] == label
            ),
            None,
        )
        if statement:
            raise InvalidParameterValue(
                f"Value {label} for parameter Label is invalid. Reason: Already exists."
            )

        principals = [
            f"arn:{get_partition(region_name)}:iam::{account_id}:root"
            for account_id in account_ids
        ]
        actions = [f"SQS:{action}" for action in actions]

        statement = {
            "Sid": label,
            "Effect": "Allow",
            "Principal": {"AWS": principals[0] if len(principals) == 1 else principals},
            "Action": actions[0] if len(actions) == 1 else actions,
            "Resource": queue.queue_arn,
        }

        queue._policy_json["Statement"].append(statement)

    def remove_permission(self, queue_name: str, label: str) -> None:
        queue = self.get_queue(queue_name)

        statements = queue._policy_json["Statement"]
        statements_new = [
            statement for statement in statements if statement["Sid"] != label
        ]

        if len(statements) == len(statements_new):
            raise InvalidParameterValue(
                f"Value {label} for parameter Label is invalid. "
                "Reason: can't find label on existing policy."
            )

        queue._policy_json["Statement"] = statements_new

    def tag_queue(self, queue_name: str, tags: Dict[str, str]) -> None:
        queue = self.get_queue(queue_name)

        if not len(tags):
            raise MissingParameter("Tags")

        if len(tags) > 50:
            raise InvalidParameterValue(f"Too many tags added for queue {queue_name}.")

        queue.tags.update(tags)

    def untag_queue(self, queue_name: str, tag_keys: List[str]) -> None:
        queue = self.get_queue(queue_name)

        if not len(tag_keys):
            raise RESTError(
                "InvalidParameterValue",
                "Tag keys must be between 1 and 128 characters in length.",
            )

        for key in tag_keys:
            try:
                del queue.tags[key]
            except KeyError:
                pass

    def list_queue_tags(self, queue_name: str) -> Queue:
        return self.get_queue(queue_name)

    def is_message_valid_based_on_retention_period(
        self, queue_name: str, message: Message
    ) -> bool:
        retention_period = self.get_queue_attributes(
            queue_name, ["MessageRetentionPeriod"]
        )["MessageRetentionPeriod"]
        retain_until = retention_period + message.sent_timestamp / 1000  # type: ignore
        if retain_until <= unix_time():
            return False
        return True


sqs_backends = BackendDict(SQSBackend, "sqs")
