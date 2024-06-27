import base64
import copy
from decimal import Decimal
from typing import Any, Dict, List, Optional, Union

from boto3.dynamodb.types import TypeDeserializer, TypeSerializer
from botocore.utils import merge_dicts

from moto.core.common_models import BaseModel
from moto.dynamodb.exceptions import (
    EmptyKeyAttributeException,
    IncorrectDataType,
    ItemSizeTooLarge,
)

from .utilities import bytesize, find_nested_key

deserializer = TypeDeserializer()
serializer = TypeSerializer()


class DDBType:
    """
    Official documentation at https://docs.aws.amazon.com/amazondynamodb/latest/APIReference/API_AttributeValue.html
    """

    BINARY_SET = "BS"
    NUMBER_SET = "NS"
    STRING_SET = "SS"
    STRING = "S"
    NUMBER = "N"
    MAP = "M"
    LIST = "L"
    BOOLEAN = "BOOL"
    BINARY = "B"
    NULL = "NULL"


class DDBTypeConversion:
    _human_type_mapping = {
        val: key.replace("_", " ")
        for key, val in DDBType.__dict__.items()
        if key.upper() == key
    }

    @classmethod
    def get_human_type(cls, abbreviated_type: str) -> str:
        """
        Args:
            abbreviated_type(str): An attribute of DDBType

        Returns:
            str: The human-readable form of the DDBType.
        """
        return cls._human_type_mapping.get(abbreviated_type, abbreviated_type)


class DynamoType(object):
    """
    http://docs.aws.amazon.com/amazondynamodb/latest/developerguide/DataModel.html#DataModelDataTypes
    """

    def __init__(self, type_as_dict: Union["DynamoType", Dict[str, Any]]):
        if type(type_as_dict) == DynamoType:
            self.type: str = type_as_dict.type
            self.value: Any = type_as_dict.value
        else:
            self.type = list(type_as_dict)[0]  # type: ignore[arg-type]
            self.value = list(type_as_dict.values())[0]  # type: ignore[union-attr]
        if self.is_list():
            self.value = [DynamoType(val) for val in self.value]
        elif self.is_map():
            self.value = dict((k, DynamoType(v)) for k, v in self.value.items())

    def __hash__(self) -> int:
        return hash((self.type, self.value))

    def __eq__(self, other: "DynamoType") -> bool:  # type: ignore[override]
        return self.type == other.type and self.value == other.value

    def __ne__(self, other: "DynamoType") -> bool:  # type: ignore[override]
        return self.type != other.type or self.value != other.value

    def __lt__(self, other: "DynamoType") -> bool:
        return self.cast_value < other.cast_value

    def __le__(self, other: "DynamoType") -> bool:
        return self.cast_value <= other.cast_value

    def __gt__(self, other: "DynamoType") -> bool:
        return self.cast_value > other.cast_value

    def __ge__(self, other: "DynamoType") -> bool:
        return self.cast_value >= other.cast_value

    def __repr__(self) -> str:
        return f"DynamoType: {self.to_json()}"

    def __add__(self, other: "DynamoType") -> "DynamoType":
        if self.type != other.type:
            raise TypeError("Different types of operandi is not allowed.")
        if self.is_number():
            self_value: Union[Decimal, int] = (
                Decimal(self.value) if "." in self.value else int(self.value)
            )
            other_value: Union[Decimal, int] = (
                Decimal(other.value) if "." in other.value else int(other.value)
            )
            total = self_value + other_value
            return DynamoType({DDBType.NUMBER: f"{total}"})
        else:
            raise IncorrectDataType()

    def __sub__(self, other: "DynamoType") -> "DynamoType":
        if self.type != other.type:
            raise TypeError("Different types of operandi is not allowed.")
        if self.type == DDBType.NUMBER:
            self_value = float(self.value) if "." in self.value else int(self.value)
            other_value = float(other.value) if "." in other.value else int(other.value)
            return DynamoType({DDBType.NUMBER: f"{self_value - other_value}"})
        else:
            raise TypeError("Sum only supported for Numbers.")

    def __getitem__(self, item: "DynamoType") -> "DynamoType":
        if isinstance(item, str):
            # If our DynamoType is a map it should be subscriptable with a key
            if self.type == DDBType.MAP:
                return self.value[item]
        elif isinstance(item, int):
            # If our DynamoType is a list is should be subscriptable with an index
            if self.type == DDBType.LIST:
                return self.value[item]
        raise TypeError(
            f"This DynamoType {self.type} is not subscriptable by a {type(item)}"
        )

    def __setitem__(self, key: Any, value: Any) -> None:
        if isinstance(key, int):
            if self.is_list():
                if key >= len(self.value):
                    # DynamoDB doesn't care you are out of box just add it to the end.
                    self.value.append(value)
                else:
                    self.value[key] = value
        elif isinstance(key, str):
            if self.is_map():
                self.value[key] = value
        else:
            raise NotImplementedError(f"No set_item for {type(key)}")

    @property
    def cast_value(self) -> Any:  # type: ignore[misc]
        if self.is_number():
            try:
                return int(self.value)
            except ValueError:
                return float(self.value)
        elif self.is_set():
            sub_type = self.type[0]
            return set([DynamoType({sub_type: v}).cast_value for v in self.value])
        elif self.is_list():
            return [DynamoType(v).cast_value for v in self.value]
        elif self.is_map():
            return dict([(k, DynamoType(v).cast_value) for k, v in self.value.items()])
        else:
            return self.value

    def child_attr(self, key: Union[int, str, None]) -> Optional["DynamoType"]:
        """
        Get Map or List children by key. str for Map, int for List.

        Returns DynamoType or None.
        """
        if isinstance(key, str) and self.is_map():
            if key in self.value:
                return DynamoType(self.value[key])

        if isinstance(key, int) and self.is_list():
            idx = key
            if 0 <= idx < len(self.value):
                return DynamoType(self.value[idx])

        return None

    def size(self) -> int:
        if self.is_number():
            value_size = len(str(self.value))
        elif self.is_set():
            sub_type = self.type[0]
            value_size = sum([DynamoType({sub_type: v}).size() for v in self.value])
        elif self.is_list():
            value_size = sum([v.size() for v in self.value])
        elif self.is_map():
            value_size = sum(
                [bytesize(k) + DynamoType(v).size() for k, v in self.value.items()]
            )
        elif isinstance(self.value, bool):
            value_size = 1
        else:
            value_size = bytesize(self.value)
        return value_size

    def to_json(self) -> Dict[str, Any]:
        # Returns a regular JSON object where the value can still be/contain a DynamoType
        if self.is_binary() and isinstance(self.value, bytes):
            # Binary data cannot be represented in JSON
            # AWS returns a base64-encoded value - the SDK's then convert that back
            return {self.type: base64.b64encode(self.value).decode("utf-8")}
        return {self.type: self.value}

    def to_regular_json(self) -> Dict[str, Any]:
        # Returns a regular JSON object in full
        value = copy.deepcopy(self.value)
        if isinstance(value, dict):
            for key, nested_value in value.items():
                value[key] = (
                    nested_value.to_regular_json()
                    if isinstance(nested_value, DynamoType)
                    else nested_value
                )
        if isinstance(value, list):
            value = [
                val.to_regular_json() if isinstance(val, DynamoType) else val
                for val in value
            ]
        if self.is_binary():
            value = base64.b64decode(value)
        return {self.type: value}

    def compare(self, range_comparison: str, range_objs: List[Any]) -> bool:
        """
        Compares this type against comparison filters
        """
        from moto.dynamodb.comparisons import get_comparison_func

        range_values = [obj.cast_value for obj in range_objs]
        comparison_func = get_comparison_func(range_comparison)
        return comparison_func(self.cast_value, *range_values)

    def is_number(self) -> bool:
        return self.type == DDBType.NUMBER

    def is_set(self) -> bool:
        return self.type in (DDBType.STRING_SET, DDBType.NUMBER_SET, DDBType.BINARY_SET)

    def is_list(self) -> bool:
        return self.type == DDBType.LIST

    def is_map(self) -> bool:
        return self.type == DDBType.MAP

    def is_binary(self) -> bool:
        return self.type == DDBType.BINARY

    def same_type(self, other: "DynamoType") -> bool:
        return self.type == other.type

    def pop(self, key: str, *args: Any, **kwargs: Any) -> None:
        if self.is_map() or self.is_list():
            self.value.pop(key, *args, **kwargs)
        else:
            raise TypeError(f"pop not supported for DynamoType {self.type}")


# https://github.com/getmoto/moto/issues/1874
# Ensure that the total size of an item does not exceed 400kb
class LimitedSizeDict(Dict[str, Any]):
    def __init__(self, *args: Any, **kwargs: Any):
        self.update(*args, **kwargs)

    def __setitem__(self, key: str, value: Any) -> None:
        current_item_size = sum(
            [
                item.size() if type(item) == DynamoType else bytesize(str(item))
                for item in (list(self.keys()) + list(self.values()))
            ]
        )
        new_item_size = bytesize(key) + (
            value.size() if type(value) == DynamoType else bytesize(str(value))
        )
        # Official limit is set to 400000 (400KB)
        # Manual testing confirms that the actual limit is between 409 and 410KB
        # We'll set the limit to something in between to be safe
        if (current_item_size + new_item_size) > 405000:
            raise ItemSizeTooLarge
        super().__setitem__(key, value)


class Item(BaseModel):
    def __init__(
        self,
        hash_key: DynamoType,
        range_key: Optional[DynamoType],
        attrs: Dict[str, Any],
    ):
        self.hash_key = hash_key
        self.range_key = range_key

        self.attrs = LimitedSizeDict()
        for key, value in attrs.items():
            self.attrs[key] = DynamoType(value)

    def __eq__(self, other: "Item") -> bool:  # type: ignore[override]
        return all(
            [
                self.hash_key == other.hash_key,
                self.range_key == other.range_key,  # type: ignore[operator]
                self.attrs == other.attrs,
            ]
        )

    def __repr__(self) -> str:
        return f"Item: {self.to_json()}"

    def size(self) -> int:
        return sum(bytesize(key) + value.size() for key, value in self.attrs.items())

    def to_json(self) -> Dict[str, Any]:
        attributes: Dict[str, Any] = {}
        for attribute_key, attribute in self.attrs.items():
            if isinstance(attribute.value, dict):
                attr_dict_value = {
                    key: value.to_regular_json()
                    for key, value in attribute.value.items()
                }
                attributes[attribute_key] = {attribute.type: attr_dict_value}
            elif isinstance(attribute.value, list):
                attr_list_value = [
                    value.to_regular_json() if isinstance(value, DynamoType) else value
                    for value in attribute.value
                ]
                attributes[attribute_key] = {attribute.type: attr_list_value}
            else:
                attributes[attribute_key] = {attribute.type: attribute.value}

        return {"Attributes": attributes}

    def to_regular_json(self) -> Dict[str, Any]:
        attributes = {}
        for key, attribute in self.attrs.items():
            attributes[key] = deserializer.deserialize(attribute.to_regular_json())
        return attributes

    def describe_attrs(
        self, attributes: Optional[Dict[str, Any]] = None
    ) -> Dict[str, Dict[str, Any]]:
        if attributes:
            included = {}
            for key, value in self.attrs.items():
                if key in attributes:
                    included[key] = value
        else:
            included = self.attrs
        return {"Item": included}

    def validate_no_empty_key_values(
        self, attribute_updates: Dict[str, Any], key_attributes: List[str]
    ) -> None:
        for attribute_name, update_action in attribute_updates.items():
            action = update_action.get("Action") or "PUT"  # PUT is default
            if action == "DELETE":
                continue
            new_value = next(iter(update_action["Value"].values()))
            if action == "PUT" and new_value == "" and attribute_name in key_attributes:
                raise EmptyKeyAttributeException

    def update_with_attribute_updates(self, attribute_updates: Dict[str, Any]) -> None:
        for attribute_name, update_action in attribute_updates.items():
            # Use default Action value, if no explicit Action is passed.
            # Default value is 'Put', according to
            # Boto3 DynamoDB.Client.update_item documentation.
            action = update_action.get("Action", "PUT")
            if action == "DELETE" and "Value" not in update_action:
                if attribute_name in self.attrs:
                    del self.attrs[attribute_name]
                continue
            new_value = list(update_action["Value"].values())[0]
            if action == "PUT":
                # TODO deal with other types
                if set(update_action["Value"].keys()) == set(["SS"]):
                    self.attrs[attribute_name] = DynamoType({"SS": new_value})
                elif set(update_action["Value"].keys()) == set(["NS"]):
                    self.attrs[attribute_name] = DynamoType({"NS": new_value})
                elif isinstance(new_value, list):
                    self.attrs[attribute_name] = DynamoType({"L": new_value})
                elif isinstance(new_value, dict):
                    self.attrs[attribute_name] = DynamoType({"M": new_value})
                elif set(update_action["Value"].keys()) == set(["N"]):
                    self.attrs[attribute_name] = DynamoType({"N": new_value})
                elif set(update_action["Value"].keys()) == set(["NULL"]):
                    if attribute_name in self.attrs:
                        del self.attrs[attribute_name]
                else:
                    self.attrs[attribute_name] = DynamoType({"S": new_value})
            elif action == "ADD":
                if set(update_action["Value"].keys()) == set(["N"]):
                    existing = self.attrs.get(attribute_name, DynamoType({"N": "0"}))
                    self.attrs[attribute_name] = DynamoType(
                        {"N": str(Decimal(existing.value) + Decimal(new_value))}
                    )
                elif set(update_action["Value"].keys()) == set(["SS"]):
                    existing = self.attrs.get(attribute_name, DynamoType({"SS": {}}))
                    new_set = set(existing.value).union(set(new_value))
                    self.attrs[attribute_name] = DynamoType({"SS": list(new_set)})
                elif set(update_action["Value"].keys()) == set(["NS"]):
                    existing = self.attrs.get(attribute_name, DynamoType({"NS": {}}))
                    new_set = set(existing.value).union(set(new_value))
                    self.attrs[attribute_name] = DynamoType({"NS": list(new_set)})
                elif set(update_action["Value"].keys()) == {"L"}:
                    existing = self.attrs.get(attribute_name, DynamoType({"L": []}))
                    new_list = existing.value + new_value
                    self.attrs[attribute_name] = DynamoType({"L": new_list})
                else:
                    # TODO: implement other data types
                    raise NotImplementedError(
                        "ADD not supported for %s"
                        % ", ".join(update_action["Value"].keys())
                    )
            elif action == "DELETE":
                if set(update_action["Value"].keys()) == set(["SS"]):
                    existing = self.attrs.get(attribute_name, DynamoType({"SS": {}}))
                    new_set = set(existing.value).difference(set(new_value))
                    self.attrs[attribute_name] = DynamoType({"SS": list(new_set)})
                elif set(update_action["Value"].keys()) == set(["NS"]):
                    existing = self.attrs.get(attribute_name, DynamoType({"NS": {}}))
                    new_set = set(existing.value).difference(set(new_value))
                    self.attrs[attribute_name] = DynamoType({"NS": list(new_set)})
                else:
                    raise NotImplementedError(
                        "ADD not supported for %s"
                        % ", ".join(update_action["Value"].keys())
                    )
            else:
                raise NotImplementedError(
                    f"{action} action not support for update_with_attribute_updates"
                )

    def project(self, projection_expressions: List[List[str]]) -> "Item":
        # Returns a new Item with only the dictionary-keys that match the provided projection_expression
        # Will return an empty Item if the expression does not match anything
        result: Dict[str, Any] = dict()
        for expr in projection_expressions:
            x = find_nested_key(expr, self.to_regular_json())
            merge_dicts(result, x)

        return Item(
            hash_key=self.hash_key,
            range_key=self.range_key,
            # 'result' is a normal Python dictionary ({'key': 'value'}
            # We need to convert that into DynamoDB dictionary ({'M': {'key': {'S': 'value'}}})
            attrs=serializer.serialize(result)["M"],
        )
