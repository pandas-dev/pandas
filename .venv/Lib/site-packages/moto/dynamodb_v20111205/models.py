import json
from collections import OrderedDict, defaultdict
from typing import Any, Dict, Iterable, List, Optional, Tuple, Union

from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.common_models import BaseModel
from moto.core.utils import unix_time, utcnow
from moto.utilities.utils import PARTITION_NAMES

from .comparisons import get_comparison_func


class DynamoJsonEncoder(json.JSONEncoder):
    def default(self, o: Any) -> Optional[str]:  # type: ignore[return]
        if hasattr(o, "to_json"):
            return o.to_json()


def dynamo_json_dump(dynamo_object: Any) -> str:
    return json.dumps(dynamo_object, cls=DynamoJsonEncoder)


class DynamoType(object):
    """
    http://docs.aws.amazon.com/amazondynamodb/latest/developerguide/DataModel.html#DataModelDataTypes
    """

    def __init__(self, type_as_dict: Dict[str, Any]):
        self.type = list(type_as_dict.keys())[0]
        self.value = list(type_as_dict.values())[0]

    def __hash__(self) -> int:
        return hash((self.type, self.value))

    def __eq__(self, other: Any) -> bool:
        return self.type == other.type and self.value == other.value

    def __repr__(self) -> str:
        return f"DynamoType: {self.to_json()}"

    def add(self, dyn_type: "DynamoType") -> None:
        if self.type == "SS":
            self.value.append(dyn_type.value)
        if self.type == "N":
            self.value = str(int(self.value) + int(dyn_type.value))

    def to_json(self) -> Dict[str, Any]:
        return {self.type: self.value}

    def compare(self, range_comparison: str, range_objs: List["DynamoType"]) -> Any:
        """
        Compares this type against comparison filters
        """
        range_values = [obj.value for obj in range_objs]
        comparison_func = get_comparison_func(range_comparison)
        return comparison_func(self.value, *range_values)


class Item(BaseModel):
    def __init__(
        self,
        hash_key: DynamoType,
        hash_key_type: str,
        range_key: Optional[DynamoType],
        range_key_type: Optional[str],
        attrs: Dict[str, Any],
    ):
        self.hash_key = hash_key
        self.hash_key_type = hash_key_type
        self.range_key = range_key
        self.range_key_type = range_key_type

        self.attrs = {}
        for key, value in attrs.items():
            self.attrs[key] = DynamoType(value)

    def __repr__(self) -> str:
        return f"Item: {self.to_json()}"

    def to_json(self) -> Dict[str, Any]:
        attributes = {}
        for attribute_key, attribute in self.attrs.items():
            attributes[attribute_key] = attribute.value

        return {"Attributes": attributes}

    def describe_attrs(self, attributes: List[str]) -> Dict[str, Any]:
        if attributes:
            included = {}
            for key, value in self.attrs.items():
                if key in attributes:
                    included[key] = value
        else:
            included = self.attrs
        return {"Item": included}


class Table(BaseModel):
    def __init__(
        self,
        account_id: str,
        name: str,
        hash_key_attr: str,
        hash_key_type: str,
        range_key_attr: Optional[str] = None,
        range_key_type: Optional[str] = None,
        read_capacity: Optional[str] = None,
        write_capacity: Optional[str] = None,
    ):
        self.account_id = account_id
        self.name = name
        self.hash_key_attr = hash_key_attr
        self.hash_key_type = hash_key_type
        self.range_key_attr = range_key_attr
        self.range_key_type = range_key_type
        self.read_capacity = read_capacity
        self.write_capacity = write_capacity
        self.created_at = utcnow()
        self.items: Dict[DynamoType, Union[Item, Dict[DynamoType, Item]]] = defaultdict(
            dict
        )

    @property
    def has_range_key(self) -> bool:
        return self.range_key_attr is not None

    @property
    def describe(self) -> Dict[str, Any]:  # type: ignore[misc]
        results: Dict[str, Any] = {
            "Table": {
                "CreationDateTime": unix_time(self.created_at),
                "KeySchema": {
                    "HashKeyElement": {
                        "AttributeName": self.hash_key_attr,
                        "AttributeType": self.hash_key_type,
                    }
                },
                "ProvisionedThroughput": {
                    "ReadCapacityUnits": self.read_capacity,
                    "WriteCapacityUnits": self.write_capacity,
                },
                "TableName": self.name,
                "TableStatus": "ACTIVE",
                "ItemCount": len(self),
                "TableSizeBytes": 0,
            }
        }
        if self.has_range_key:
            results["Table"]["KeySchema"]["RangeKeyElement"] = {
                "AttributeName": self.range_key_attr,
                "AttributeType": self.range_key_type,
            }
        return results

    def __len__(self) -> int:
        return sum(
            [(len(value) if self.has_range_key else 1) for value in self.items.values()]  # type: ignore
        )

    def __nonzero__(self) -> bool:
        return True

    def __bool__(self) -> bool:
        return self.__nonzero__()

    def put_item(self, item_attrs: Dict[str, Any]) -> Item:
        hash_value = DynamoType(item_attrs.get(self.hash_key_attr))  # type: ignore[arg-type]
        if self.has_range_key:
            range_value: Optional[DynamoType] = DynamoType(
                item_attrs.get(self.range_key_attr)  # type: ignore[arg-type]
            )
        else:
            range_value = None

        item = Item(
            hash_value, self.hash_key_type, range_value, self.range_key_type, item_attrs
        )

        if range_value:
            self.items[hash_value][range_value] = item  # type: ignore[index]
        else:
            self.items[hash_value] = item
        return item

    def get_item(
        self, hash_key: DynamoType, range_key: Optional[DynamoType]
    ) -> Optional[Item]:
        if self.has_range_key and not range_key:
            raise ValueError(
                "Table has a range key, but no range key was passed into get_item"
            )
        try:
            if range_key:
                return self.items[hash_key][range_key]  # type: ignore
            else:
                return self.items[hash_key]  # type: ignore
        except KeyError:
            return None

    def query(
        self, hash_key: DynamoType, range_comparison: str, range_objs: Any
    ) -> Tuple[Iterable[Item], bool]:
        results = []
        last_page = True  # Once pagination is implemented, change this

        if self.range_key_attr:
            possible_results = self.items[hash_key].values()  # type: ignore[union-attr]
        else:
            possible_results = list(self.all_items())

        if range_comparison:
            for result in possible_results:
                if result.range_key.compare(range_comparison, range_objs):  # type: ignore[union-attr]
                    results.append(result)
        else:
            # If we're not filtering on range key, return all values
            results = possible_results  # type: ignore[assignment]
        return results, last_page

    def all_items(self) -> Iterable[Item]:
        for hash_set in self.items.values():
            if self.range_key_attr:
                for item in hash_set.values():  # type: ignore
                    yield item
            else:
                yield hash_set  # type: ignore[misc]

    def scan(self, filters: Dict[str, Any]) -> Tuple[List[Item], int, bool]:
        results = []
        scanned_count = 0
        last_page = True  # Once pagination is implemented, change this

        for result in self.all_items():
            scanned_count += 1
            passes_all_conditions = True
            for (
                attribute_name,
                (comparison_operator, comparison_objs),
            ) in filters.items():
                attribute = result.attrs.get(attribute_name)

                if attribute:
                    # Attribute found
                    if not attribute.compare(comparison_operator, comparison_objs):
                        passes_all_conditions = False
                        break
                elif comparison_operator == "NULL":
                    # Comparison is NULL and we don't have the attribute
                    continue
                else:
                    # No attribute found and comparison is no NULL. This item
                    # fails
                    passes_all_conditions = False
                    break

            if passes_all_conditions:
                results.append(result)

        return results, scanned_count, last_page

    def delete_item(
        self, hash_key: DynamoType, range_key: Optional[DynamoType]
    ) -> Optional[Item]:
        try:
            if range_key:
                return self.items[hash_key].pop(range_key)  # type: ignore
            else:
                return self.items.pop(hash_key)  # type: ignore
        except KeyError:
            return None

    def update_item(
        self,
        hash_key: DynamoType,
        range_key: Optional[DynamoType],
        attr_updates: Dict[str, Any],
    ) -> Optional[Item]:
        item = self.get_item(hash_key, range_key)
        if not item:
            return None

        for attr, update in attr_updates.items():
            if update["Action"] == "PUT":
                item.attrs[attr] = DynamoType(update["Value"])
            if update["Action"] == "DELETE":
                item.attrs.pop(attr)
            if update["Action"] == "ADD":
                item.attrs[attr].add(DynamoType(update["Value"]))
        return item


class DynamoDBBackend(BaseBackend):
    def __init__(self, region_name: str, account_id: str):
        super().__init__(region_name, account_id)
        self.tables: Dict[str, Table] = OrderedDict()

    def create_table(self, name: str, **params: Any) -> Table:
        table = Table(self.account_id, name, **params)
        self.tables[name] = table
        return table

    def delete_table(self, name: str) -> Optional[Table]:
        return self.tables.pop(name, None)

    def update_table_throughput(
        self, name: str, new_read_units: str, new_write_units: str
    ) -> Table:
        table = self.tables[name]
        table.read_capacity = new_read_units
        table.write_capacity = new_write_units
        return table

    def put_item(self, table_name: str, item_attrs: Dict[str, Any]) -> Optional[Item]:
        table = self.tables.get(table_name)
        if not table:
            return None

        return table.put_item(item_attrs)

    def get_item(
        self,
        table_name: str,
        hash_key_dict: Dict[str, Any],
        range_key_dict: Optional[Dict[str, Any]],
    ) -> Optional[Item]:
        table = self.tables.get(table_name)
        if not table:
            return None

        hash_key = DynamoType(hash_key_dict)
        range_key = DynamoType(range_key_dict) if range_key_dict else None

        return table.get_item(hash_key, range_key)

    def query(
        self,
        table_name: str,
        hash_key_dict: Dict[str, Any],
        range_comparison: str,
        range_value_dicts: List[Dict[str, Any]],
    ) -> Tuple[Optional[Iterable[Item]], Optional[bool]]:
        table = self.tables.get(table_name)
        if not table:
            return None, None

        hash_key = DynamoType(hash_key_dict)
        range_values = [DynamoType(range_value) for range_value in range_value_dicts]

        return table.query(hash_key, range_comparison, range_values)

    def scan(
        self, table_name: str, filters: Dict[str, Any]
    ) -> Tuple[Optional[List[Item]], Optional[int], Optional[bool]]:
        table = self.tables.get(table_name)
        if not table:
            return None, None, None

        scan_filters = {}
        for key, (comparison_operator, comparison_values) in filters.items():
            dynamo_types = [DynamoType(value) for value in comparison_values]
            scan_filters[key] = (comparison_operator, dynamo_types)

        return table.scan(scan_filters)

    def delete_item(
        self,
        table_name: str,
        hash_key_dict: Dict[str, Any],
        range_key_dict: Optional[Dict[str, Any]],
    ) -> Optional[Item]:
        table = self.tables.get(table_name)
        if not table:
            return None

        hash_key = DynamoType(hash_key_dict)
        range_key = DynamoType(range_key_dict) if range_key_dict else None

        return table.delete_item(hash_key, range_key)

    def update_item(
        self,
        table_name: str,
        hash_key_dict: Dict[str, Any],
        range_key_dict: Optional[Dict[str, Any]],
        attr_updates: Dict[str, Any],
    ) -> Optional[Item]:
        table = self.tables.get(table_name)
        if not table:
            return None

        hash_key = DynamoType(hash_key_dict)
        range_key = DynamoType(range_key_dict) if range_key_dict else None

        return table.update_item(hash_key, range_key, attr_updates)


dynamodb_backends = BackendDict(
    DynamoDBBackend,
    "dynamodb_v20111205",
    use_boto3_regions=False,
    additional_regions=PARTITION_NAMES,
)
