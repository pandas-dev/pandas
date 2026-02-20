import re
from typing import Any, Optional
from urllib.parse import urlparse

from moto.core.responses import ActionResult, BaseResponse, EmptyResult

from .constants import (
    DEFAULT_RECEIVED_MESSAGES,
    MAXIMUM_VISIBILITY_TIMEOUT,
)
from .exceptions import (
    BatchEntryIdsNotDistinct,
    EmptyBatchRequest,
    InvalidAttributeName,
    MaxVisibilityTimeout,
    QueueDoesNotExist,
    SQSException,
)
from .models import SQSBackend, sqs_backends
from .utils import validate_message_attributes


class SQSResponse(BaseResponse):
    region_regex = re.compile(r"://(.+?)\.queue\.amazonaws\.com")

    def __init__(self) -> None:
        super().__init__(service_name="sqs")
        self.automated_parameter_parsing = True

    @property
    def sqs_backend(self) -> SQSBackend:
        return sqs_backends[self.current_account][self.region]

    def _determine_resource(self) -> str:
        queue_name = self._get_queue_name()
        try:
            return self.sqs_backend.get_queue(queue_name).queue_arn
        except QueueDoesNotExist:
            return "*"

    def _get_queue_name(self) -> str:
        try:
            queue_url = self._get_param("QueueUrl")
            if queue_url.startswith("http://") or queue_url.startswith("https://"):
                return queue_url.split("/")[-1]
            else:
                # The parameter could be the name itself, which AWS also accepts
                return queue_url
        except (AttributeError, TypeError):
            pass
        # Fallback to reading from the URL for botocore
        return self.path.split("/")[-1]

    def _get_validated_visibility_timeout(self, timeout: Optional[str] = None) -> int:
        """
        :raises ValueError: If specified visibility timeout exceeds MAXIMUM_VISIBILITY_TIMEOUT
        :raises TypeError: If visibility timeout was not specified
        """
        if timeout is not None:
            visibility_timeout = int(timeout)
        else:
            visibility_timeout = self._get_param("VisibilityTimeout")
        if visibility_timeout > MAXIMUM_VISIBILITY_TIMEOUT:
            raise ValueError
        return visibility_timeout

    def create_queue(self) -> ActionResult:
        request_url = urlparse(self.uri)
        queue_name = self._get_param("QueueName")
        attributes = self._get_param("Attributes", {})
        tags = self._get_param("tags", {})
        queue = self.sqs_backend.create_queue(queue_name, tags, **attributes)
        result = {"QueueUrl": queue.url(request_url)}
        return ActionResult(result)

    def get_queue_url(self) -> ActionResult:
        request_url = urlparse(self.uri)
        queue_name = self._get_param("QueueName")
        queue = self.sqs_backend.get_queue_url(queue_name)
        result = {"QueueUrl": queue.url(request_url)}
        return ActionResult(result)

    def list_queues(self) -> ActionResult:
        request_url = urlparse(self.uri)
        queue_name_prefix = self._get_param("QueueNamePrefix")
        queues = self.sqs_backend.list_queues(queue_name_prefix)
        result = {}
        if queues:
            result["QueueUrls"] = [queue.url(request_url) for queue in queues]
        return ActionResult(result)

    def change_message_visibility(self) -> ActionResult:
        queue_name = self._get_queue_name()
        receipt_handle = self._get_param("ReceiptHandle")
        try:
            visibility_timeout = self._get_validated_visibility_timeout()
        except ValueError:
            raise MaxVisibilityTimeout()
        self.sqs_backend.change_message_visibility(
            queue_name=queue_name,
            receipt_handle=receipt_handle,
            visibility_timeout=visibility_timeout,
        )
        return EmptyResult()

    def change_message_visibility_batch(self) -> ActionResult:
        queue_name = self._get_queue_name()
        entries = self._get_param("Entries", [])
        success, failed = self.sqs_backend.change_message_visibility_batch(
            queue_name, entries
        )
        result = {"Successful": [{"Id": _id} for _id in success], "Failed": failed}
        return ActionResult(result)

    def get_queue_attributes(self) -> ActionResult:
        queue_name = self._get_queue_name()
        attribute_names = self._get_param("AttributeNames", [])
        if attribute_names and "" in attribute_names:
            raise InvalidAttributeName("")
        attributes = self.sqs_backend.get_queue_attributes(queue_name, attribute_names)
        result = {
            "Attributes": {
                key: str(value)
                for key, value in attributes.items()
                if value is not None
            }
        }
        if not result["Attributes"]:
            return EmptyResult()
        return ActionResult(result)

    def set_queue_attributes(self) -> ActionResult:
        # TODO validate self.get_param('QueueUrl')
        attributes = self._get_param("Attributes", {})
        queue_name = self._get_queue_name()
        self.sqs_backend.set_queue_attributes(queue_name, attributes)
        return EmptyResult()

    def delete_queue(self) -> ActionResult:
        # TODO validate self.get_param('QueueUrl')
        queue_name = self._get_queue_name()
        self.sqs_backend.delete_queue(queue_name)
        return EmptyResult()

    def send_message(self) -> ActionResult:
        message = self._get_param("MessageBody")
        delay_seconds = self._get_param("DelaySeconds")
        message_group_id = self._get_param("MessageGroupId")
        message_dedupe_id = self._get_param("MessageDeduplicationId")
        message_attributes = self._get_param("MessageAttributes", {})
        validate_message_attributes(message_attributes)
        system_message_attributes = self._get_param("MessageSystemAttributes")
        validate_message_attributes(system_message_attributes)
        queue_name = self._get_queue_name()
        message = self.sqs_backend.send_message(
            queue_name,
            message,
            message_attributes=message_attributes,
            delay_seconds=delay_seconds,
            deduplication_id=message_dedupe_id,
            group_id=message_group_id,
            system_attributes=system_message_attributes,
        )
        resp = {
            "MD5OfMessageBody": message.body_md5,
            "MessageId": message.id,
        }
        if len(message.message_attributes) > 0:
            resp["MD5OfMessageAttributes"] = message.attribute_md5
        return ActionResult(resp)

    def send_message_batch(self) -> ActionResult:
        queue_name = self._get_queue_name()
        self.sqs_backend.get_queue(queue_name)
        entries = self._get_param("Entries", [])
        entries = {str(idx): entry for idx, entry in enumerate(entries)}
        for entry in entries.values():
            if "MessageAttributes" in entry:
                validate_message_attributes(entry["MessageAttributes"])
            else:
                entry["MessageAttributes"] = {}
            if "DelaySeconds" not in entry:
                entry["DelaySeconds"] = None

        if entries == {}:
            raise EmptyBatchRequest()

        messages, failed_invalid_delay = self.sqs_backend.send_message_batch(
            queue_name, entries
        )

        errors = []
        for entry in failed_invalid_delay:
            errors.append(
                {
                    "Id": entry["Id"],
                    "SenderFault": "true",
                    "Code": "InvalidParameterValue",
                    "Message": "Value 1800 for parameter DelaySeconds is invalid. Reason: DelaySeconds must be &gt;= 0 and &lt;= 900.",
                }
            )

        resp: dict[str, Any] = {"Successful": [], "Failed": errors}
        for msg in messages:
            msg_dict = {
                "Id": msg.user_id,  # type: ignore
                "MessageId": msg.id,
                "MD5OfMessageBody": msg.body_md5,
            }
            if len(msg.message_attributes) > 0:
                msg_dict["MD5OfMessageAttributes"] = msg.attribute_md5
            resp["Successful"].append(msg_dict)
        return ActionResult(resp)

    def delete_message(self) -> ActionResult:
        queue_name = self._get_queue_name()
        receipt_handle = self._get_param("ReceiptHandle")
        self.sqs_backend.delete_message(queue_name, receipt_handle)
        return EmptyResult()

    def delete_message_batch(self) -> ActionResult:
        queue_name = self._get_queue_name()
        receipts = self._get_param("Entries", [])
        if not receipts:
            raise EmptyBatchRequest(action="Delete")
        receipt_seen = set()
        for receipt_and_id in receipts:
            receipt = receipt_and_id["ReceiptHandle"]
            if receipt in receipt_seen:
                raise BatchEntryIdsNotDistinct(receipt_and_id["Id"])
            receipt_seen.add(receipt)
        success, failed = self.sqs_backend.delete_message_batch(queue_name, receipts)
        result = {"Successful": [{"Id": _id} for _id in success], "Failed": failed}
        return ActionResult(result)

    def purge_queue(self) -> ActionResult:
        queue_name = self._get_queue_name()
        self.sqs_backend.purge_queue(queue_name)
        return EmptyResult()

    def receive_message(self) -> ActionResult:
        queue_name = self._get_queue_name()
        attribute_names = self._get_param("AttributeNames", [])
        message_system_attribute_names = self._get_param(
            "MessageSystemAttributeNames", []
        )
        # AttributeNames has been deprecated in favor of MessageSystemAttributeNames.
        # We combine both here (with no duplicates) in order to support either parameter.
        message_system_attribute_names = list(
            set(attribute_names + message_system_attribute_names)
        )
        message_attribute_names = self._get_param("MessageAttributeNames", [])
        queue = self.sqs_backend.get_queue(queue_name)
        message_count = self._get_param(
            "MaxNumberOfMessages", DEFAULT_RECEIVED_MESSAGES
        )
        if message_count < 1 or message_count > 10:
            raise SQSException(
                "InvalidParameterValue",
                "An error occurred (InvalidParameterValue) when calling "
                f"the ReceiveMessage operation: Value {message_count} for parameter "
                "MaxNumberOfMessages is invalid. Reason: must be between "
                "1 and 10, if provided.",
            )
        try:
            wait_time = int(self._get_param("WaitTimeSeconds"))
        except TypeError:
            wait_time = int(queue.receive_message_wait_time_seconds)  # type: ignore
        if wait_time < 0 or wait_time > 20:
            raise SQSException(
                "InvalidParameterValue",
                "An error occurred (InvalidParameterValue) when calling "
                f"the ReceiveMessage operation: Value {wait_time} for parameter "
                "WaitTimeSeconds is invalid. Reason: must be &lt;= 0 and "
                "&gt;= 20 if provided.",
            )
        try:
            visibility_timeout = self._get_validated_visibility_timeout()
        except TypeError:
            visibility_timeout = queue.visibility_timeout  # type: ignore
        except ValueError:
            raise MaxVisibilityTimeout()

        messages = self.sqs_backend.receive_message(
            queue_name,
            message_count,
            wait_time,
            visibility_timeout,
            message_attribute_names,
        )
        SUPPORTED_SYSTEM_MESSAGE_ATTRIBUTE_NAMES = [
            "AWSTraceHeader",
            "ApproximateFirstReceiveTimestamp",
            "ApproximateReceiveCount",
            "MessageDeduplicationId",
            "MessageGroupId",
            "SenderId",
            "SentTimestamp",
            "SequenceNumber",
        ]
        attributes_to_include = []
        include_all = "All" in message_system_attribute_names
        for attribute_name in SUPPORTED_SYSTEM_MESSAGE_ATTRIBUTE_NAMES:
            if include_all or attribute_name in message_system_attribute_names:
                attributes_to_include.append(attribute_name)

        msgs = []
        for message in messages:
            msg: dict[str, Any] = {
                "MessageId": message.id,
                "ReceiptHandle": message.receipt_handle,
                "MD5OfBody": message.body_md5,
                "Body": message.body,
            }
            if len(message.message_attributes):
                msg["MessageAttributes"] = message.message_attributes
                msg["MD5OfMessageAttributes"] = message.attribute_md5
            attributes = {}
            if "SenderId" in attributes_to_include:
                attributes["SenderId"] = message.sender_id
            if "SentTimestamp" in attributes_to_include:
                attributes["SentTimestamp"] = str(message.sent_timestamp)
            if "ApproximateReceiveCount" in attributes_to_include:
                attributes["ApproximateReceiveCount"] = str(
                    message.approximate_receive_count
                )
            if "ApproximateFirstReceiveTimestamp" in attributes_to_include:
                attributes["ApproximateFirstReceiveTimestamp"] = str(
                    message.approximate_first_receive_timestamp
                )
            if "MessageDeduplicationId" in attributes_to_include:
                if message.deduplication_id is not None:
                    attributes["MessageDeduplicationId"] = message.deduplication_id
            if "MessageGroupId" in attributes_to_include:
                if message.group_id is not None:
                    attributes["MessageGroupId"] = message.group_id
            if "AWSTraceHeader" in attributes_to_include:
                if message.system_attributes and message.system_attributes.get(
                    "AWSTraceHeader"
                ):
                    attributes["AWSTraceHeader"] = message.system_attributes[
                        "AWSTraceHeader"
                    ].get("StringValue")
            if "SequenceNumber" in attributes_to_include:
                if message.sequence_number is not None:
                    attributes["SequenceNumber"] = message.sequence_number
            if attributes_to_include:
                msg["Attributes"] = attributes
            msgs.append(msg)

        result = {"Messages": msgs} if msgs else {}
        return ActionResult(result)

    def list_dead_letter_source_queues(self) -> ActionResult:
        request_url = urlparse(self.uri)
        queue_name = self._get_queue_name()
        queues = self.sqs_backend.list_dead_letter_source_queues(queue_name)
        result = {"queueUrls": [queue.url(request_url) for queue in queues]}
        return ActionResult(result)

    def add_permission(self) -> ActionResult:
        queue_name = self._get_queue_name()
        actions = self._get_param("Actions", [])
        account_ids = self._get_param("AWSAccountIds", [])
        label = self._get_param("Label")
        self.sqs_backend.add_permission(
            region_name=self.region,
            queue_name=queue_name,
            actions=actions,
            account_ids=account_ids,
            label=label,
        )
        return EmptyResult()

    def remove_permission(self) -> ActionResult:
        queue_name = self._get_queue_name()
        label = self._get_param("Label")
        self.sqs_backend.remove_permission(queue_name, label)
        return EmptyResult()

    def tag_queue(self) -> ActionResult:
        queue_name = self._get_queue_name()
        tags = self._get_param("Tags", [])
        self.sqs_backend.tag_queue(queue_name, tags)
        return EmptyResult()

    def untag_queue(self) -> ActionResult:
        queue_name = self._get_queue_name()
        tag_keys = self._get_param("TagKeys", [])
        self.sqs_backend.untag_queue(queue_name, tag_keys)
        return EmptyResult()

    def list_queue_tags(self) -> ActionResult:
        queue_name = self._get_queue_name()
        queue = self.sqs_backend.list_queue_tags(queue_name)
        result = {"Tags": queue.tags}
        return ActionResult(result)
