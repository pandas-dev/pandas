import string
from typing import Any

from moto.moto_api._internal import mock_random as random

from .exceptions import MessageAttributesInvalid


def generate_receipt_handle() -> str:
    # http://docs.aws.amazon.com/AWSSimpleQueueService/latest/SQSDeveloperGuide/ImportantIdentifiers.html#ImportantIdentifiers-receipt-handles
    length = 185
    return "".join(random.choice(string.ascii_lowercase) for x in range(length))


def validate_message_attributes(message_attributes: dict[str, Any]) -> None:
    for name, value in (message_attributes or {}).items():
        data_type = value["DataType"]

        if not data_type:
            raise MessageAttributesInvalid(
                f"The message attribute '{name}' must contain non-empty message attribute value."
            )

        data_type_parts = data_type.split(".")[0]
        if data_type_parts not in [
            "String",
            "Binary",
            "Number",
        ]:
            raise MessageAttributesInvalid(
                f"The message attribute '{name}' has an invalid message attribute type, the set of supported type prefixes is Binary, Number, and String."
            )

        possible_value_fields = ["StringValue", "BinaryValue"]
        for field in possible_value_fields:
            if field in value and value[field] is None:
                raise MessageAttributesInvalid(
                    f"The message attribute '{name}' must contain non-empty message attribute value for message attribute type '{data_type}'."
                )
