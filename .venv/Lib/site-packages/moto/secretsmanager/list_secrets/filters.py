import re
from typing import TYPE_CHECKING, Iterator, List, Union

if TYPE_CHECKING:
    from ..models import FakeSecret, ReplicaSecret


def name_filter(secret: "FakeSecret", names: List[str]) -> bool:
    return _matcher(names, [secret.name])


def description_filter(secret: "FakeSecret", descriptions: List[str]) -> bool:
    if not secret.description:
        return False
    # The documentation states that this search uses `Prefix match`
    # But actual testing determines that it uses the same approach to the `all_filter`:
    # 'Breaks the filter value string into words and then searches all attributes for matches.'
    return _matcher(
        descriptions, [secret.description], match_prefix=False, case_sensitive=False
    )


def owning_service_filter(secret: "FakeSecret", owning_services: List[str]) -> bool:
    return _matcher(
        owning_services,
        [secret.owning_service] if secret.owning_service else [],
        match_prefix=False,
        case_sensitive=False,
    )


def tag_key(secret: Union["FakeSecret", "ReplicaSecret"], tag_keys: List[str]) -> bool:
    if not secret.tags:
        return False
    return _matcher(tag_keys, [tag["Key"] for tag in secret.tags])


def tag_value(
    secret: Union["FakeSecret", "ReplicaSecret"], tag_values: List[str]
) -> bool:
    if not secret.tags:
        return False
    return _matcher(tag_values, [tag["Value"] for tag in secret.tags])


def filter_all(secret: Union["FakeSecret", "ReplicaSecret"], values: List[str]) -> bool:
    attributes = [secret.name]
    if secret.description:
        attributes.append(secret.description)
    if secret.tags:
        attributes += [tag["Key"] for tag in secret.tags] + [
            tag["Value"] for tag in secret.tags
        ]

    return _matcher(values, attributes, match_prefix=False, case_sensitive=False)


def _matcher(
    patterns: List[str],
    strings: List[str],
    match_prefix: bool = True,
    case_sensitive: bool = True,
) -> bool:
    for pattern in [p for p in patterns if p.startswith("!")]:
        for string in strings:
            if not _match_pattern(
                pattern[1:], string, match_prefix, case_sensitive=case_sensitive
            ):
                return True

    for pattern in [p for p in patterns if not p.startswith("!")]:
        for string in strings:
            if _match_pattern(
                pattern, string, match_prefix, case_sensitive=case_sensitive
            ):
                return True
    return False


def _match_pattern(
    pattern: str, value: str, match_prefix: bool = True, case_sensitive: bool = True
) -> bool:
    if match_prefix:
        if not case_sensitive:
            return value.lower().startswith(pattern.lower())
        else:
            return value.startswith(pattern)
    else:
        pattern_words = split_words(pattern)
        if not pattern_words:
            return False
        value_words = split_words(value)
        if not case_sensitive:
            pattern_words = [p.lower() for p in pattern_words]
            value_words = [v.lower() for v in value_words]
        for pattern_word in pattern_words:
            # all words in value must start with pattern_word
            if not any(
                value_word.startswith(pattern_word) for value_word in value_words
            ):
                return False
    return True


def split_words(s: str) -> List[str]:
    """
    Secrets are split by special characters first (/, +, _, etc)
    Partial results are then split again by UpperCasing
    """
    special_chars = ["/", "-", "_", "+", "=", ".", "@"]

    if s in special_chars:
        # Special case: this does not return any values
        return []

    for char in special_chars:
        if char in s:
            others = special_chars.copy()
            others.remove(char)
            contains_other = any([c in s for c in others])
            if contains_other:
                # Secret contains two different characters, i.e. my/secret+value
                # Values like this will not be split
                return [s]
            else:
                return list(split_by_uppercase(s.split(char)))
    return list(split_by_uppercase(s))


def split_by_uppercase(s: Union[str, List[str]]) -> Iterator[str]:
    """
    Split a string into words. Words are recognized by upper case letters, i.e.:
    test   -> [test]
    MyTest -> [My, Test]
    """
    if isinstance(s, str):
        for x in re.split(r"([^a-z][a-z]+)", s):
            if x:
                yield x.strip()
    else:
        for word in s:
            yield from split_by_uppercase(word)
