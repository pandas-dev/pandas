from typing import TYPE_CHECKING, List

if TYPE_CHECKING:
    from ..models import FakeSecret


def name_filter(secret: "FakeSecret", names: List[str]) -> bool:
    return _matcher(names, [secret.name])


def description_filter(secret: "FakeSecret", descriptions: List[str]) -> bool:
    return _matcher(descriptions, [secret.description], match_prefix=False)  # type: ignore


def tag_key(secret: "FakeSecret", tag_keys: List[str]) -> bool:
    if not secret.tags:
        return False
    return _matcher(tag_keys, [tag["Key"] for tag in secret.tags])


def tag_value(secret: "FakeSecret", tag_values: List[str]) -> bool:
    if not secret.tags:
        return False
    return _matcher(tag_values, [tag["Value"] for tag in secret.tags])


def filter_all(secret: "FakeSecret", values: List[str]) -> bool:
    attributes = [secret.name, secret.description]
    if secret.tags:
        attributes += [tag["Key"] for tag in secret.tags] + [
            tag["Value"] for tag in secret.tags
        ]

    return _matcher(values, attributes)  # type: ignore


def _matcher(
    patterns: List[str], strings: List[str], match_prefix: bool = True
) -> bool:
    for pattern in [p for p in patterns if p.startswith("!")]:
        for string in strings:
            if not _match_pattern(pattern[1:], string, match_prefix):
                return True

    for pattern in [p for p in patterns if not p.startswith("!")]:
        for string in strings:
            if _match_pattern(pattern, string, match_prefix):
                return True
    return False


def _match_pattern(pattern: str, value: str, match_prefix: bool = True) -> bool:
    if match_prefix:
        return value.startswith(pattern)
    else:
        pattern_words = pattern.split(" ")
        value_words = value.split(" ")
        for pattern_word in pattern_words:
            # all words in value must start with pattern_word
            if not any(
                value_word.startswith(pattern_word) for value_word in value_words
            ):
                return False
    return True
