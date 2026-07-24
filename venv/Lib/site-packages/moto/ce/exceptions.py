"""Exceptions raised by the ce service."""

from moto.core.exceptions import JsonRESTError


class CostCategoryNotFound(JsonRESTError):
    def __init__(self, ccd_id: str):
        super().__init__(
            "ResourceNotFoundException", f"No Cost Categories found with ID {ccd_id}"
        )
