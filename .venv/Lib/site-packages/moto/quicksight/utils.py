from typing import Any, Dict, List, Type, Union

from moto.core.common_models import BaseModel

from .data_models import QuicksightGroup
from .exceptions import ParamValidationError, SchemaException, ValidationException

PAGINATION_MODEL = {
    "list_users": {
        "input_token": "next_token",
        "limit_key": "max_results",
        "limit_default": 100,  # This should be the sum of the directory limits
        "unique_attribute": "arn",
    },
    "list_user_groups": {
        "input_token": "next_token",
        "limit_key": "max_results",
        "limit_default": 100,  # This should be the sum of the directory limits
        "unique_attribute": "arn",
    },
    "list_groups": {
        "input_token": "next_token",
        "limit_key": "max_results",
        "limit_default": 100,  # This should be the sum of the directory limits
        "unique_attribute": "arn",
    },
    "list_group_memberships": {
        "input_token": "next_token",
        "limit_key": "max_results",
        "limit_default": 100,  # This should be the sum of the directory limits
        "unique_attribute": "arn",
    },
    "search_groups": {
        "input_token": "next_token",
        "limit_key": "max_results",
        "limit_default": 100,  # This should be the sum of the directory limits
        "unique_attribute": "arn",
    },
}


class QuicksightBaseSearchFilter(BaseModel):
    """Base Search Filter."""

    schema: Dict[str, Any] = {
        "$schema": "https://json-schema.org/draft/2020-12/schema",
    }

    def __init__(self, operator: str, name: str, value: str):
        self.operator = operator
        self.name = name
        self.value = value

    def match(self, input: BaseModel) -> bool:
        raise NotImplementedError


class QuicksightGroupSearchFilter(QuicksightBaseSearchFilter):
    """Group Search Filter."""

    schema = {
        "$schema": "https://json-schema.org/draft/2020-12/schema",
        "type": "array",
        "items": {
            "type": "object",
            "properties": {
                "Operator": {"type": "string", "enum": ["StringEquals", "StartsWith"]},
                "Name": {"type": "string", "enum": ["GROUP_NAME"]},
                "Value": {"type": "string"},
            },
            "required": ["Operator", "Name", "Value"],
        },
        "minItems": 1,
        "maxItems": 1,
    }

    def __init__(self, operator: str, name: str, value: str):
        super().__init__(operator, name, value)

    def match(self, input: BaseModel) -> bool:
        if isinstance(input, QuicksightGroup):
            if self.name == "GROUP_NAME" and self.operator == "StringEquals":
                return input.group_name.lower() == self.value.lower()
            if self.name == "GROUP_NAME" and self.operator == "StartsWith":
                return input.group_name.lower().startswith(self.value.lower())
        return False

    @classmethod
    def parse_filter(cls, filter: Dict[str, str]) -> QuicksightBaseSearchFilter:
        return QuicksightGroupSearchFilter(
            operator=filter.get("Operator", ""),
            name=filter.get("Name", ""),
            value=filter.get("Value", ""),
        )


class QuicksightSearchFilterList:
    """Generic QuickSight Search Filter List."""

    def __init__(self, filters: List[QuicksightBaseSearchFilter]):
        self.filters: List[QuicksightBaseSearchFilter] = filters

    def match(self, input: BaseModel) -> bool:
        return any([filter.match(input) for filter in self.filters])


class QuicksightSearchFilterFactory:
    """Creates the appropriate search filter."""

    @classmethod
    def validate_and_create_filter(
        cls, model_type: Type[BaseModel], input: Union[List[Dict[str, str]], None]
    ) -> QuicksightSearchFilterList:
        if issubclass(model_type, QuicksightGroup):
            if input is None:
                # Should never happen as there are other validations before but just in case.
                raise ParamValidationError(
                    'Missing required parameter in input: "Filters"'
                )
            from jsonschema import validate
            from jsonschema.exceptions import SchemaError, ValidationError

            try:
                validate(instance=input, schema=QuicksightGroupSearchFilter.schema)
                return QuicksightSearchFilterList(
                    filters=[
                        QuicksightGroupSearchFilter.parse_filter(filter)
                        for filter in input
                    ]
                )
            except ValidationError as err:
                raise ValidationException(err.message)

            except SchemaError as err:
                # This exception can happen only, when the schema defined in moto is wrong.
                # The schema is not presented by the user.
                raise SchemaException(f"Schema definition wrong in moto: {err.message}")

        raise NotImplementedError(
            f"Filter for '{model_type.__name__}' not implemented."
        )
