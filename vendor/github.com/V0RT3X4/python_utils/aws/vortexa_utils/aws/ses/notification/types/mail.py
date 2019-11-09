from typing import List, Dict, Any
from dataclasses import dataclass


@dataclass
class Mail:
    """Mail Object.

    Attributes
    ----------
    destination: List[str]
        A complete list of all recipient addresses (including To: and CC:)
        from the MIME headers of the incoming email.
    messageId: str
        String that contains the unique ID assigned to the email by Amazon SES.
        If the email was delivered to Amazon S3, the message ID is also the
        Amazon S3 object key that was used to write the message to your Amazon
        S3 bucket.
    source: str
        String that contains the email address (the envelope MAIL FROM address)
        that the email was sent from.
    timestamp:
        String that contains the time at which the email was received,
        in ISO8601 format.
    headers: List[List[str]]
        A list of Amazon SES headers and your custom headers.
        Each header in the list has a name field and a value field.
    commonHeaders: List[List[str]]
        A list of headers common to all emails.
        Each header in the list is composed of a name and a value.
    headersTruncated: str
        String that specifies whether the headers were truncated,
        which will happen if the headers are larger than 10 KB.
        Possible values are true and false.

    """

    destination: List[str]
    messageId: str
    source: str
    timestamp: str
    headers: List[Dict[str, str]]
    commonHeaders: Dict[str, Any]
    headersTruncated: str
