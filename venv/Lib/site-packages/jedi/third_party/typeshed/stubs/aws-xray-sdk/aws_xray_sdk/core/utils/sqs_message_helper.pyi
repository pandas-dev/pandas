from _typeshed import Incomplete
from collections.abc import Mapping
from typing import Final

SQS_XRAY_HEADER: Final = "AWSTraceHeader"

class SqsMessageHelper:
    @staticmethod
    def isSampled(sqs_message: Mapping[str, Incomplete]) -> bool: ...
