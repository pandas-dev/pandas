"""Exceptions raised by the comprehend service."""

from typing import List

from moto.core.exceptions import JsonRESTError


class ResourceNotFound(JsonRESTError):
    def __init__(self) -> None:
        super().__init__(
            "ResourceNotFoundException",
            "RESOURCE_NOT_FOUND: Could not find specified resource.",
        )


class DetectPIIValidationException(JsonRESTError):
    def __init__(self, language: str, all_languages: List[str]) -> None:
        all_languages_str = str(all_languages).replace("'", "")
        super().__init__(
            "ValidationException",
            f"Value '{language}' at 'languageCode'failed to satisfy constraint: "
            f"Member must satisfy enum value set: {all_languages_str}",
        )


class TextSizeLimitExceededException(JsonRESTError):
    def __init__(self, size: int) -> None:
        super().__init__(
            "TextSizeLimitExceededException",
            "Input text size exceeds limit. Max length of request text allowed is 100000 bytes while in "
            f"this request the text size is {size} bytes",
        )
