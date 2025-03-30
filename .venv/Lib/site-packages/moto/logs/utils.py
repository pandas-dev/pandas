from typing import Type

PAGINATION_MODEL = {
    "describe_log_groups": {
        "input_token": "next_token",
        "limit_key": "limit",
        "limit_default": 50,
        "unique_attribute": "arn",
        "fail_on_invalid_token": False,
    },
    "describe_log_streams": {
        "input_token": "next_token",
        "limit_key": "limit",
        "limit_default": 50,
        "unique_attribute": "arn",
    },
}


class FilterPattern:
    def __init__(self, term: str):
        self.term = term


class QuotedTermFilterPattern(FilterPattern):
    def matches(self, message: str) -> bool:
        # We still have the quotes around the term - we should remove those in the parser
        return self.term[1:-1] in message


class SingleTermFilterPattern(FilterPattern):
    def matches(self, message: str) -> bool:
        required_words = self.term.split(" ")
        return all([word in message for word in required_words])


class UnsupportedFilterPattern(FilterPattern):
    def matches(self, message: str) -> bool:  # pylint: disable=unused-argument
        return True


class EventMessageFilter:
    def __init__(self, pattern: str):
        current_phrase = ""
        current_type: Type[FilterPattern] = None  # type: ignore
        if pattern:
            for char in pattern:
                if not current_type:
                    if char.isalpha():
                        current_type = SingleTermFilterPattern
                    elif char == '"':
                        current_type = QuotedTermFilterPattern
                    else:
                        current_type = UnsupportedFilterPattern
                current_phrase += char
        else:
            current_type = UnsupportedFilterPattern
        self.filter_type = current_type(current_phrase)

    def matches(self, message: str) -> bool:
        return self.filter_type.matches(message)  # type: ignore
