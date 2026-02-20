from __future__ import annotations

import base64
import binascii


class Base64EncodedString:
    """
    A type to encapsulate a Base64 encoded string.
    """

    def __init__(self, encoded_string: str):
        if not isinstance(encoded_string, str):
            raise TypeError("Base64EncodedString must be initialized with a string.")
        try:
            base64.b64decode(encoded_string, validate=True)
        except binascii.Error:
            raise ValueError("Invalid Base64 encoded string provided.")
        self._encoded_string = encoded_string

    def as_bytes(self, encoding: str = "utf-8") -> bytes:
        return self._encoded_string.encode(encoding)

    def decode(self, encoding: str = "utf-8") -> str:
        decoded_bytes = base64.b64decode(self._encoded_string)
        return decoded_bytes.decode(encoding)

    @classmethod
    def from_encoded_bytes(cls, raw_bytes: bytes) -> Base64EncodedString:
        return cls(raw_bytes.decode("ascii"))

    @classmethod
    def from_raw_string(
        cls, raw_string: str, encoding: str = "utf-8"
    ) -> Base64EncodedString:
        encoded_bytes = base64.b64encode(raw_string.encode(encoding))
        return cls.from_encoded_bytes(encoded_bytes)

    def __str__(self) -> str:
        return self._encoded_string

    def __repr__(self) -> str:
        return f"Base64EncodedString('{self._encoded_string}')"
