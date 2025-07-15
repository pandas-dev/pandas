import json
from typing import Any, Dict, Optional, Union

from moto.core.common_types import TYPE_RESPONSE
from moto.core.responses import BaseResponse
from moto.core.utils import camelcase_to_underscores

from .models import DynamoDBBackend, dynamo_json_dump, dynamodb_backends


class DynamoHandler(BaseResponse):
    def __init__(self) -> None:
        super().__init__(service_name="dynamodb")

    def get_endpoint_name(self, headers: Dict[str, str]) -> Optional[str]:  # type: ignore[return]
        """Parses request headers and extracts part od the X-Amz-Target
        that corresponds to a method of DynamoHandler

        ie: X-Amz-Target: DynamoDB_20111205.ListTables -> ListTables
        """
        # Headers are case-insensitive. Probably a better way to do this.
        match = headers.get("x-amz-target") or headers.get("X-Amz-Target")
        if match:
            return match.split(".")[1]

    def error(self, type_: str, status: int = 400) -> TYPE_RESPONSE:
        return status, self.response_headers, dynamo_json_dump({"__type": type_})

    def call_action(self) -> TYPE_RESPONSE:
        self.body = json.loads(self.body or "{}")
        endpoint = self.get_endpoint_name(self.headers)
        if endpoint:
            endpoint = camelcase_to_underscores(endpoint)
            response = getattr(self, endpoint)()
            if isinstance(response, str):
                return 200, self.response_headers, response

            else:
                status_code, new_headers, response_content = response
                self.response_headers.update(new_headers)
                return status_code, self.response_headers, response_content
        else:
            return 404, self.response_headers, ""

    @property
    def backend(self) -> DynamoDBBackend:
        return dynamodb_backends[self.current_account][self.partition]

    def list_tables(self) -> str:
        body = self.body
        limit = body.get("Limit")
        if body.get("ExclusiveStartTableName"):
            last = body.get("ExclusiveStartTableName")
            start = list(self.backend.tables.keys()).index(last) + 1
        else:
            start = 0
        all_tables = list(self.backend.tables.keys())
        if limit:
            tables = all_tables[start : start + limit]
        else:
            tables = all_tables[start:]
        response: Dict[str, Any] = {"TableNames": tables}
        if limit and len(all_tables) > start + limit:
            response["LastEvaluatedTableName"] = tables[-1]
        return dynamo_json_dump(response)

    def create_table(self) -> str:
        body = self.body
        name = body["TableName"]

        key_schema = body["KeySchema"]
        hash_key = key_schema["HashKeyElement"]
        hash_key_attr = hash_key["AttributeName"]
        hash_key_type = hash_key["AttributeType"]

        range_key = key_schema.get("RangeKeyElement", {})
        range_key_attr = range_key.get("AttributeName")
        range_key_type = range_key.get("AttributeType")

        throughput = body["ProvisionedThroughput"]
        read_units = throughput["ReadCapacityUnits"]
        write_units = throughput["WriteCapacityUnits"]

        table = self.backend.create_table(
            name,
            hash_key_attr=hash_key_attr,
            hash_key_type=hash_key_type,
            range_key_attr=range_key_attr,
            range_key_type=range_key_type,
            read_capacity=int(read_units),
            write_capacity=int(write_units),
        )
        return dynamo_json_dump(table.describe)

    def delete_table(self) -> Union[str, TYPE_RESPONSE]:
        name = self.body["TableName"]
        table = self.backend.delete_table(name)
        if table:
            return dynamo_json_dump(table.describe)
        else:
            er = "com.amazonaws.dynamodb.v20111205#ResourceNotFoundException"
            return self.error(er)

    def update_table(self) -> str:
        name = self.body["TableName"]
        throughput = self.body["ProvisionedThroughput"]
        new_read_units = throughput["ReadCapacityUnits"]
        new_write_units = throughput["WriteCapacityUnits"]
        table = self.backend.update_table_throughput(
            name, new_read_units, new_write_units
        )
        return dynamo_json_dump(table.describe)

    def describe_table(self) -> Union[TYPE_RESPONSE, str]:
        name = self.body["TableName"]
        try:
            table = self.backend.tables[name]
        except KeyError:
            er = "com.amazonaws.dynamodb.v20111205#ResourceNotFoundException"
            return self.error(er)
        return dynamo_json_dump(table.describe)

    def put_item(self) -> Union[TYPE_RESPONSE, str]:
        name = self.body["TableName"]
        item = self.body["Item"]
        result = self.backend.put_item(name, item)
        if result:
            item_dict = result.to_json()
            item_dict["ConsumedCapacityUnits"] = 1
            return dynamo_json_dump(item_dict)
        else:
            er = "com.amazonaws.dynamodb.v20111205#ResourceNotFoundException"
            return self.error(er)

    def batch_write_item(self) -> str:
        table_batches = self.body["RequestItems"]

        for table_name, table_requests in table_batches.items():
            for table_request in table_requests:
                request_type = list(table_request)[0]
                request = list(table_request.values())[0]

                if request_type == "PutRequest":
                    item = request["Item"]
                    self.backend.put_item(table_name, item)
                elif request_type == "DeleteRequest":
                    key = request["Key"]
                    hash_key = key["HashKeyElement"]
                    range_key = key.get("RangeKeyElement")
                    self.backend.delete_item(table_name, hash_key, range_key)

        response = {
            "Responses": {
                "Thread": {"ConsumedCapacityUnits": 1.0},
                "Reply": {"ConsumedCapacityUnits": 1.0},
            },
            "UnprocessedItems": {},
        }

        return dynamo_json_dump(response)

    def get_item(self) -> Union[TYPE_RESPONSE, str]:
        name = self.body["TableName"]
        key = self.body["Key"]
        hash_key = key["HashKeyElement"]
        range_key = key.get("RangeKeyElement")
        attrs_to_get = self.body.get("AttributesToGet")
        try:
            item = self.backend.get_item(name, hash_key, range_key)
        except ValueError:
            er = "com.amazon.coral.validate#ValidationException"
            return self.error(er, status=400)
        if item:
            item_dict = item.describe_attrs(attrs_to_get)
            item_dict["ConsumedCapacityUnits"] = 0.5
            return dynamo_json_dump(item_dict)
        else:
            # Item not found
            er = "com.amazonaws.dynamodb.v20111205#ResourceNotFoundException"
            return self.error(er, status=404)

    def batch_get_item(self) -> str:
        table_batches = self.body["RequestItems"]

        results: Dict[str, Any] = {"Responses": {"UnprocessedKeys": {}}}

        for table_name, table_request in table_batches.items():
            items = []
            keys = table_request["Keys"]
            attributes_to_get = table_request.get("AttributesToGet")
            for key in keys:
                hash_key = key["HashKeyElement"]
                range_key = key.get("RangeKeyElement")
                item = self.backend.get_item(table_name, hash_key, range_key)
                if item:
                    item_describe = item.describe_attrs(attributes_to_get)
                    items.append(item_describe)
            results["Responses"][table_name] = {
                "Items": items,
                "ConsumedCapacityUnits": 1,
            }
        return dynamo_json_dump(results)

    def query(self) -> Union[TYPE_RESPONSE, str]:
        name = self.body["TableName"]
        hash_key = self.body["HashKeyValue"]
        range_condition = self.body.get("RangeKeyCondition")
        if range_condition:
            range_comparison = range_condition["ComparisonOperator"]
            range_values = range_condition["AttributeValueList"]
        else:
            range_comparison = None
            range_values = []

        items, _ = self.backend.query(name, hash_key, range_comparison, range_values)

        if items is None:
            er = "com.amazonaws.dynamodb.v20111205#ResourceNotFoundException"
            return self.error(er)

        result = {
            "Count": len(items),  # type: ignore[arg-type]
            "Items": [item.attrs for item in items],
            "ConsumedCapacityUnits": 1,
        }

        # Implement this when we do pagination
        # if not last_page:
        #     result["LastEvaluatedKey"] = {
        #         "HashKeyElement": items[-1].hash_key,
        #         "RangeKeyElement": items[-1].range_key,
        #     }
        return dynamo_json_dump(result)

    def scan(self) -> Union[TYPE_RESPONSE, str]:
        name = self.body["TableName"]

        filters = {}
        scan_filters = self.body.get("ScanFilter", {})
        for attribute_name, scan_filter in scan_filters.items():
            # Keys are attribute names. Values are tuples of (comparison,
            # comparison_value)
            comparison_operator = scan_filter["ComparisonOperator"]
            comparison_values = scan_filter.get("AttributeValueList", [])
            filters[attribute_name] = (comparison_operator, comparison_values)

        items, scanned_count, _ = self.backend.scan(name, filters)

        if items is None:
            er = "com.amazonaws.dynamodb.v20111205#ResourceNotFoundException"
            return self.error(er)

        result = {
            "Count": len(items),
            "Items": [item.attrs for item in items if item],
            "ConsumedCapacityUnits": 1,
            "ScannedCount": scanned_count,
        }

        # Implement this when we do pagination
        # if not last_page:
        #     result["LastEvaluatedKey"] = {
        #         "HashKeyElement": items[-1].hash_key,
        #         "RangeKeyElement": items[-1].range_key,
        #     }
        return dynamo_json_dump(result)

    def delete_item(self) -> Union[TYPE_RESPONSE, str]:
        name = self.body["TableName"]
        key = self.body["Key"]
        hash_key = key["HashKeyElement"]
        range_key = key.get("RangeKeyElement")
        return_values = self.body.get("ReturnValues", "")
        item = self.backend.delete_item(name, hash_key, range_key)
        if item:
            if return_values == "ALL_OLD":
                item_dict = item.to_json()
            else:
                item_dict = {"Attributes": []}
            item_dict["ConsumedCapacityUnits"] = 0.5
            return dynamo_json_dump(item_dict)
        else:
            er = "com.amazonaws.dynamodb.v20111205#ResourceNotFoundException"
            return self.error(er)

    def update_item(self) -> Union[TYPE_RESPONSE, str]:
        name = self.body["TableName"]
        key = self.body["Key"]
        hash_key = key["HashKeyElement"]
        range_key = key.get("RangeKeyElement")
        updates = self.body["AttributeUpdates"]

        item = self.backend.update_item(name, hash_key, range_key, updates)

        if item:
            item_dict = item.to_json()
            item_dict["ConsumedCapacityUnits"] = 0.5

            return dynamo_json_dump(item_dict)
        else:
            er = "com.amazonaws.dynamodb.v20111205#ResourceNotFoundException"
            return self.error(er)
