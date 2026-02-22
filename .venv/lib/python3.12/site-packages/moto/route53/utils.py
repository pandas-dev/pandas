"""Pagination control model for Route53."""

from .exceptions import InvalidPaginationToken, UnsupportedCharacter

PAGINATION_MODEL = {
    "list_hosted_zones": {
        "input_token": "marker",
        "limit_key": "max_size",
        "limit_default": 100,
        "unique_attribute": "id",
    },
    "list_query_logging_configs": {
        "input_token": "next_token",
        "limit_key": "max_results",
        "limit_default": 100,
        "unique_attribute": "hosted_zone_id",
        "fail_on_invalid_token": InvalidPaginationToken,
    },
}


def validate_domain_name(domain_name: str, code: str = "InvalidInput") -> None:
    # https://docs.aws.amazon.com/Route53/latest/DeveloperGuide/DomainNameFormat.html
    # Domain names can get wacky with international alphabets/emojis
    # (https://i❤️.ws is valid!)
    # We only validate that some reserved characters (0-40, 177-377) are prefaced with a backslash
    prev = ""
    for char in domain_name:
        if (ord(char) <= 32 or (127 < ord(char) < 255)) and prev != "\\":
            raise UnsupportedCharacter(code=code, char=char)
        prev = char
