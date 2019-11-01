from vortexa_utils.aws.utils.dataclasses import nested_dataclass
from .mail import Mail
from .receipt import Receipt


@nested_dataclass
class SESRecord:
    receipt: Receipt
    mail: Mail


@nested_dataclass
class Record:
    """
    """
    eventSource: str   # "aws:ses",
    eventVersion: str  # "1.0",
    ses: SESRecord
