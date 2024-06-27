from typing import List, Optional

from moto.utilities.utils import load_resource_as_str


class ReservedKeywords:
    """
    DynamoDB has an extensive list of keywords. Keywords are considered when validating the expression Tree.
    Not earlier since an update expression like "SET path = VALUE 1" fails with:
        'Invalid UpdateExpression: Syntax error; token: "1", near: "VALUE 1"'
    """

    KEYWORDS: Optional[List[str]] = None

    @classmethod
    def get_reserved_keywords(cls) -> List[str]:
        if cls.KEYWORDS is None:
            cls.KEYWORDS = cls._get_reserved_keywords()
        return cls.KEYWORDS

    @classmethod
    def _get_reserved_keywords(cls) -> List[str]:
        """
        Get a list of reserved keywords of DynamoDB
        """
        reserved_keywords = load_resource_as_str(__name__, "reserved_keywords.txt")
        return reserved_keywords.split()
