import copy
import itertools
import json
from functools import wraps
from typing import Any, Callable, Dict, List, Optional, Union

from moto.core.common_types import TYPE_RESPONSE
from moto.core.responses import BaseResponse
from moto.dynamodb.models import DynamoDBBackend, Table, dynamodb_backends
from moto.dynamodb.models.utilities import dynamo_json_dump
from moto.dynamodb.parsing.key_condition_expression import parse_expression
from moto.dynamodb.parsing.reserved_keywords import ReservedKeywords
from moto.utilities.aws_headers import amz_crc32

from .exceptions import (
    KeyIsEmptyStringException,
    MockValidationException,
    ResourceNotFoundException,
    UnknownKeyType,
)

TRANSACTION_MAX_ITEMS = 25


def include_consumed_capacity(
    val: float = 1.0,
) -> Callable[
    [Callable[["DynamoHandler"], str]],
    Callable[["DynamoHandler"], Union[str, TYPE_RESPONSE]],
]:
    def _inner(
        f: Callable[["DynamoHandler"], str],
    ) -> Callable[["DynamoHandler"], Union[str, TYPE_RESPONSE]]:
        @wraps(f)
        def _wrapper(
            *args: "DynamoHandler", **kwargs: None
        ) -> Union[str, TYPE_RESPONSE]:
            (handler,) = args
            expected_capacity = handler.body.get("ReturnConsumedCapacity", "NONE")
            if expected_capacity not in ["NONE", "TOTAL", "INDEXES"]:
                type_ = "ValidationException"
                headers = handler.response_headers.copy()
                headers["status"] = "400"
                message = f"1 validation error detected: Value '{expected_capacity}' at 'returnConsumedCapacity' failed to satisfy constraint: Member must satisfy enum value set: [INDEXES, TOTAL, NONE]"
                return (
                    400,
                    headers,
                    dynamo_json_dump({"__type": type_, "message": message}),
                )
            table_name = handler.body.get("TableName", "")
            index_name = handler.body.get("IndexName", None)

            response = f(*args, **kwargs)

            if isinstance(response, str):
                body = json.loads(response)

                if expected_capacity == "TOTAL":
                    body["ConsumedCapacity"] = {
                        "TableName": table_name,
                        "CapacityUnits": val,
                    }
                elif expected_capacity == "INDEXES":
                    body["ConsumedCapacity"] = {
                        "TableName": table_name,
                        "CapacityUnits": val,
                        "Table": {"CapacityUnits": val},
                    }
                    if index_name:
                        body["ConsumedCapacity"]["LocalSecondaryIndexes"] = {
                            index_name: {"CapacityUnits": val}
                        }

                return dynamo_json_dump(body)

            return response

        return _wrapper

    return _inner


def validate_put_has_empty_keys(
    field_updates: Dict[str, Any], table: Table, custom_error_msg: Optional[str] = None
) -> None:
    """
    Error if any keys have an empty value. Checks Global index attributes as well
    """
    if table:
        key_names = table.attribute_keys
        gsi_key_names = list(
            itertools.chain(*[gsi.schema_key_attrs for gsi in table.global_indexes])
        )

        # string/binary fields with empty string as value
        empty_str_fields = [
            key
            for (key, val) in field_updates.items()
            if next(iter(val.keys())) in ["S", "B"] and next(iter(val.values())) == ""
        ]

        # First validate that all of the GSI-keys are set
        empty_gsi_key = next(
            (kn for kn in gsi_key_names if kn in empty_str_fields), None
        )
        if empty_gsi_key:
            gsi_name = table.global_indexes[0].name
            raise MockValidationException(
                f"One or more parameter values are not valid. A value specified for a secondary index key is not supported. The AttributeValue for a key attribute cannot contain an empty string value. IndexName: {gsi_name}, IndexKey: {empty_gsi_key}"
            )

        # Then validate that all of the regular keys are set
        empty_key = next(
            (keyname for keyname in key_names if keyname in empty_str_fields), None
        )
        if empty_key:
            msg = (
                custom_error_msg
                or "One or more parameter values were invalid: An AttributeValue may not contain an empty string. Key: {}"
            )
            raise MockValidationException(msg.format(empty_key))


def put_has_empty_attrs(field_updates: Dict[str, Any], table: Table) -> bool:
    # Example invalid attribute: [{'M': {'SS': {'NS': []}}}]
    def _validate_attr(attr: Dict[str, Any]) -> bool:
        if "NS" in attr and attr["NS"] == []:
            return True
        else:
            return any(
                [_validate_attr(val) for val in attr.values() if isinstance(val, dict)]
            )

    if table:
        key_names = table.attribute_keys
        attrs_to_check = [
            val for attr, val in field_updates.items() if attr not in key_names
        ]
        return any([_validate_attr(attr) for attr in attrs_to_check])
    return False


def validate_put_has_gsi_keys_set_to_none(item: Dict[str, Any], table: Table) -> None:
    for gsi in table.global_indexes:
        for attr in gsi.schema:
            attr_name = attr["AttributeName"]
            if attr_name in item and item[attr_name] == {"NULL": True}:
                raise MockValidationException(
                    f"One or more parameter values were invalid: Type mismatch for Index Key {attr_name} Expected: S Actual: NULL IndexName: {gsi.name}"
                )


def check_projection_expression(expression: str) -> None:
    if expression.upper() in ReservedKeywords.get_reserved_keywords():
        raise MockValidationException(
            f"ProjectionExpression: Attribute name is a reserved keyword; reserved keyword: {expression}"
        )
    if expression[0].isnumeric():
        raise MockValidationException(
            "ProjectionExpression: Attribute name starts with a number"
        )
    if " " in expression:
        raise MockValidationException(
            "ProjectionExpression: Attribute name contains white space"
        )


class DynamoHandler(BaseResponse):
    def __init__(self) -> None:
        super().__init__(service_name="dynamodb")

    def get_endpoint_name(self, headers: Any) -> Optional[str]:
        """Parses request headers and extracts part od the X-Amz-Target
        that corresponds to a method of DynamoHandler

        ie: X-Amz-Target: DynamoDB_20111205.ListTables -> ListTables
        """
        # Headers are case-insensitive. Probably a better way to do this.
        match = headers.get("x-amz-target") or headers.get("X-Amz-Target")
        if match:
            return match.split(".")[1]
        return None

    @property
    def dynamodb_backend(self) -> DynamoDBBackend:
        """
        :return: DynamoDB Backend
        :rtype: moto.dynamodb.models.DynamoDBBackend
        """
        return dynamodb_backends[self.current_account][self.region]

    @amz_crc32
    def call_action(self) -> TYPE_RESPONSE:
        self.body = json.loads(self.body or "{}")
        return super().call_action()

    def list_tables(self) -> str:
        body = self.body
        limit = body.get("Limit", 100)
        exclusive_start_table_name = body.get("ExclusiveStartTableName")
        tables, last_eval = self.dynamodb_backend.list_tables(
            limit, exclusive_start_table_name
        )

        response: Dict[str, Any] = {"TableNames": tables}
        if last_eval:
            response["LastEvaluatedTableName"] = last_eval

        return dynamo_json_dump(response)

    def create_table(self) -> str:
        body = self.body
        table_name = body["TableName"]
        # check billing mode and get the throughput
        if "BillingMode" in body.keys() and body["BillingMode"] == "PAY_PER_REQUEST":
            billing_mode = "PAY_PER_REQUEST"
        else:
            billing_mode = "PROVISIONED"  # Default
        throughput = body.get("ProvisionedThroughput")
        sse_spec = body.get("SSESpecification")
        key_schema = body["KeySchema"]
        attr = body["AttributeDefinitions"]
        global_indexes = body.get("GlobalSecondaryIndexes")
        local_secondary_indexes = body.get("LocalSecondaryIndexes")
        streams = body.get("StreamSpecification")
        tags = body.get("Tags", [])
        deletion_protection_enabled = body.get("DeletionProtectionEnabled", False)

        self._validate_table_creation(
            billing_mode=billing_mode,
            throughput=throughput,
            key_schema=key_schema,
            global_indexes=global_indexes,
            local_secondary_indexes=local_secondary_indexes,
            attr=attr,
        )

        table = self.dynamodb_backend.create_table(
            table_name,
            schema=key_schema,
            throughput=throughput,
            attr=attr,
            global_indexes=global_indexes or [],
            indexes=local_secondary_indexes or [],
            streams=streams,
            billing_mode=billing_mode,
            sse_specification=sse_spec,
            tags=tags,
            deletion_protection_enabled=deletion_protection_enabled,
        )
        return dynamo_json_dump(table.describe())

    def _validate_table_creation(
        self,
        billing_mode: str,
        throughput: Optional[Dict[str, Any]],
        key_schema: List[Dict[str, str]],
        global_indexes: Optional[List[Dict[str, Any]]],
        local_secondary_indexes: Optional[List[Dict[str, Any]]],
        attr: List[Dict[str, str]],
    ) -> None:
        # Validate Throughput
        if billing_mode == "PAY_PER_REQUEST" and throughput:
            raise MockValidationException(
                "ProvisionedThroughput cannot be specified when BillingMode is PAY_PER_REQUEST"
            )
        if billing_mode == "PROVISIONED" and throughput is None:
            raise MockValidationException(
                "One or more parameter values were invalid: ReadCapacityUnits and WriteCapacityUnits must both be specified when BillingMode is PROVISIONED"
            )

        # Validate KeySchema
        for idx, _key in enumerate(key_schema, start=1):
            key_type = _key["KeyType"]
            if key_type not in ["HASH", "RANGE"]:
                raise UnknownKeyType(
                    key_type=key_type, position=f"keySchema.{idx}.member.keyType"
                )
        if len(key_schema) > 2:
            key_elements = [
                f"KeySchemaElement(attributeName={key.get('AttributeName')}, keyType={key['KeyType']})"
                for key in key_schema
            ]
            provided_keys = ", ".join(key_elements)
            err = f"1 validation error detected: Value '[{provided_keys}]' at 'keySchema' failed to satisfy constraint: Member must have length less than or equal to 2"
            raise MockValidationException(err)

        # Validate Global Indexes
        if global_indexes == []:
            raise MockValidationException(
                "One or more parameter values were invalid: List of GlobalSecondaryIndexes is empty"
            )
        for idx, g_idx in enumerate(global_indexes or [], start=1):
            for idx2, _key in enumerate(g_idx["KeySchema"], start=1):
                key_type = _key["KeyType"]
                if key_type not in ["HASH", "RANGE"]:
                    position = f"globalSecondaryIndexes.{idx}.member.keySchema.{idx2}.member.keyType"
                    raise UnknownKeyType(key_type=key_type, position=position)

        # Validate Local Indexes
        if local_secondary_indexes == []:
            raise MockValidationException(
                "One or more parameter values were invalid: List of LocalSecondaryIndexes is empty"
            )
        for idx, g_idx in enumerate(local_secondary_indexes or [], start=1):
            for idx2, _key in enumerate(g_idx["KeySchema"], start=1):
                key_type = _key["KeyType"]
                if key_type not in ["HASH", "RANGE"]:
                    position = f"localSecondaryIndexes.{idx}.member.keySchema.{idx2}.member.keyType"
                    raise UnknownKeyType(key_type=key_type, position=position)

        # Validate Attributes
        expected_attrs = []
        expected_attrs.extend([key["AttributeName"] for key in key_schema])
        local_key_schemas = itertools.chain(
            *list(idx["KeySchema"] for idx in (local_secondary_indexes or []))
        )
        expected_attrs.extend(schema["AttributeName"] for schema in local_key_schemas)

        global_key_schemas = itertools.chain(
            *list(idx["KeySchema"] for idx in (global_indexes or []))
        )
        expected_attrs.extend(schema["AttributeName"] for schema in global_key_schemas)
        expected_attrs = list(set(expected_attrs))
        expected_attrs.sort()
        actual_attrs = [item["AttributeName"] for item in attr]
        actual_attrs.sort()
        has_index = global_indexes is not None or local_secondary_indexes is not None
        if actual_attrs != expected_attrs:
            self._throw_attr_error(actual_attrs, expected_attrs, has_index)

    def _throw_attr_error(
        self, actual_attrs: List[str], expected_attrs: List[str], indexes: bool
    ) -> None:
        def dump_list(list_: List[str]) -> str:
            return str(list_).replace("'", "")

        err_head = "One or more parameter values were invalid: "
        if len(actual_attrs) > len(expected_attrs):
            if indexes:
                raise MockValidationException(
                    err_head
                    + "Some AttributeDefinitions are not used. AttributeDefinitions: "
                    + dump_list(actual_attrs)
                    + ", keys used: "
                    + dump_list(expected_attrs)
                )
            else:
                raise MockValidationException(
                    err_head
                    + "Number of attributes in KeySchema does not exactly match number of attributes defined in AttributeDefinitions"
                )
        elif len(actual_attrs) < len(expected_attrs):
            if indexes:
                raise MockValidationException(
                    err_head
                    + "Some index key attributes are not defined in AttributeDefinitions. Keys: "
                    + dump_list(list(set(expected_attrs) - set(actual_attrs)))
                    + ", AttributeDefinitions: "
                    + dump_list(actual_attrs)
                )
            else:
                raise MockValidationException(
                    "Invalid KeySchema: Some index key attribute have no definition"
                )
        else:
            if indexes:
                raise MockValidationException(
                    err_head
                    + "Some index key attributes are not defined in AttributeDefinitions. Keys: "
                    + dump_list(list(set(expected_attrs) - set(actual_attrs)))
                    + ", AttributeDefinitions: "
                    + dump_list(actual_attrs)
                )
            else:
                raise MockValidationException(
                    err_head
                    + "Some index key attributes are not defined in AttributeDefinitions. Keys: "
                    + dump_list(expected_attrs)
                    + ", AttributeDefinitions: "
                    + dump_list(actual_attrs)
                )

    def _get_filter_expression(self) -> Optional[str]:
        filter_expression = self.body.get("FilterExpression")
        if filter_expression == "":
            raise MockValidationException(
                "Invalid FilterExpression: The expression can not be empty;"
            )
        return filter_expression

    def _get_projection_expression(self) -> Optional[str]:
        expression = self.body.get("ProjectionExpression")
        if expression == "":
            raise MockValidationException(
                "Invalid ProjectionExpression: The expression can not be empty;"
            )
        return expression

    def delete_table(self) -> str:
        name = self.body["TableName"]
        table = self.dynamodb_backend.delete_table(name)
        return dynamo_json_dump(table.describe())

    def describe_endpoints(self) -> str:
        response = {"Endpoints": self.dynamodb_backend.describe_endpoints()}
        return dynamo_json_dump(response)

    def tag_resource(self) -> str:
        table_arn = self.body["ResourceArn"]
        tags = self.body["Tags"]
        self.dynamodb_backend.tag_resource(table_arn, tags)
        return ""

    def untag_resource(self) -> str:
        table_arn = self.body["ResourceArn"]
        tags = self.body["TagKeys"]
        self.dynamodb_backend.untag_resource(table_arn, tags)
        return ""

    def list_tags_of_resource(self) -> str:
        table_arn = self.body["ResourceArn"]
        all_tags = self.dynamodb_backend.list_tags_of_resource(table_arn)
        all_tag_keys = [tag["Key"] for tag in all_tags]
        marker = self.body.get("NextToken")
        if marker:
            start = all_tag_keys.index(marker) + 1
        else:
            start = 0
        max_items = 10  # there is no default, but using 10 to make testing easier
        tags_resp = all_tags[start : start + max_items]
        next_marker = None
        if len(all_tags) > start + max_items:
            next_marker = tags_resp[-1]["Key"]
        if next_marker:
            return json.dumps({"Tags": tags_resp, "NextToken": next_marker})
        return json.dumps({"Tags": tags_resp})

    def update_table(self) -> str:
        name = self.body["TableName"]
        attr_definitions = self.body.get("AttributeDefinitions", None)
        global_index = self.body.get("GlobalSecondaryIndexUpdates", None)
        throughput = self.body.get("ProvisionedThroughput", None)
        billing_mode = self.body.get("BillingMode", None)
        stream_spec = self.body.get("StreamSpecification", None)
        deletion_protection_enabled = self.body.get("DeletionProtectionEnabled")
        table = self.dynamodb_backend.update_table(
            name=name,
            attr_definitions=attr_definitions,
            global_index=global_index,
            throughput=throughput,
            billing_mode=billing_mode,
            stream_spec=stream_spec,
            deletion_protection_enabled=deletion_protection_enabled,
        )
        return dynamo_json_dump(table.describe())

    def describe_table(self) -> str:
        name = self.body["TableName"]
        table = self.dynamodb_backend.describe_table(name)
        return dynamo_json_dump(table)

    @include_consumed_capacity()
    def put_item(self) -> str:
        name = self.body["TableName"]
        item = self.body["Item"]
        return_values = self.body.get("ReturnValues", "NONE")
        return_values_on_condition_check_failure = self.body.get(
            "ReturnValuesOnConditionCheckFailure"
        )

        if return_values not in ("ALL_OLD", "NONE"):
            raise MockValidationException("Return values set to invalid value")

        table = self.dynamodb_backend.get_table(name)
        validate_put_has_empty_keys(item, table)
        if put_has_empty_attrs(item, table):
            raise MockValidationException(
                "One or more parameter values were invalid: An number set  may not be empty"
            )
        validate_put_has_gsi_keys_set_to_none(item, table)

        overwrite = "Expected" not in self.body
        if not overwrite:
            expected = self.body["Expected"]
        else:
            expected = None

        if return_values == "ALL_OLD":
            existing_item = self.dynamodb_backend.get_item(name, item)
            if existing_item:
                existing_attributes = existing_item.to_json()["Attributes"]
            else:
                existing_attributes = {}

        # Attempt to parse simple ConditionExpressions into an Expected
        # expression
        condition_expression = self.body.get("ConditionExpression")
        expression_attribute_names = self.body.get("ExpressionAttributeNames", {})
        expression_attribute_values = self._get_expr_attr_values()

        if condition_expression:
            overwrite = False

        result = self.dynamodb_backend.put_item(
            name,
            item,
            expected,
            condition_expression,
            expression_attribute_names,
            expression_attribute_values,
            overwrite,
            return_values_on_condition_check_failure,
        )

        item_dict = result.to_json()
        if return_values == "ALL_OLD":
            item_dict["Attributes"] = existing_attributes
        else:
            item_dict.pop("Attributes", None)
        return dynamo_json_dump(item_dict)

    def batch_write_item(self) -> str:
        table_batches = self.body["RequestItems"]
        put_requests = []
        delete_requests = []
        for table_name, table_requests in table_batches.items():
            table = self.dynamodb_backend.get_table(table_name)
            for table_request in table_requests:
                request_type = list(table_request.keys())[0]
                request = list(table_request.values())[0]
                if request_type == "PutRequest":
                    item = request["Item"]
                    validate_put_has_empty_keys(
                        item,
                        table,
                        custom_error_msg="One or more parameter values are not valid. The AttributeValue for a key attribute cannot contain an empty string value. Key: {}",
                    )
                    put_requests.append((table_name, item))
                elif request_type == "DeleteRequest":
                    keys = request["Key"]
                    delete_requests.append((table_name, keys))
        if self._contains_duplicates(
            [json.dumps(k[1]) for k in delete_requests]
        ) or self._contains_duplicates([json.dumps(k[1]) for k in put_requests]):
            raise MockValidationException(
                "Provided list of item keys contains duplicates"
            )
        for table_name, item in put_requests:
            self.dynamodb_backend.put_item(table_name, item)
        for table_name, keys in delete_requests:
            self.dynamodb_backend.delete_item(table_name, keys)

        response = {
            "ConsumedCapacity": [
                {
                    "TableName": table_name,
                    "CapacityUnits": 1.0,
                    "Table": {"CapacityUnits": 1.0},
                }
                for table_name, table_requests in table_batches.items()
            ],
            "ItemCollectionMetrics": {},
            "UnprocessedItems": {},
        }

        return dynamo_json_dump(response)

    @include_consumed_capacity(0.5)
    def get_item(self) -> str:
        name = self.body["TableName"]
        self.dynamodb_backend.get_table(name)
        key = self.body["Key"]
        empty_keys = [k for k, v in key.items() if not next(iter(v.values()))]
        if empty_keys:
            raise KeyIsEmptyStringException(empty_keys[0])

        projection_expression = self._get_projection_expression()
        attributes_to_get = self.body.get("AttributesToGet")
        if projection_expression and attributes_to_get:
            raise MockValidationException(
                "Can not use both expression and non-expression parameters in the same request: Non-expression parameters: {AttributesToGet} Expression parameters: {ProjectionExpression}"
            )

        expression_attribute_names = self.body.get("ExpressionAttributeNames")
        if expression_attribute_names == {}:
            if projection_expression is None:
                raise MockValidationException(
                    "ExpressionAttributeNames can only be specified when using expressions"
                )
            else:
                raise MockValidationException(
                    "ExpressionAttributeNames must not be empty"
                )

        expression_attribute_names = expression_attribute_names or {}
        projection_expressions = self._adjust_projection_expression(
            projection_expression, expression_attribute_names
        )

        item = self.dynamodb_backend.get_item(name, key, projection_expressions)
        if item:
            item_dict = item.describe_attrs(attributes=None)
            return dynamo_json_dump(item_dict)
        else:
            # Item not found
            return dynamo_json_dump({})

    def batch_get_item(self) -> str:
        table_batches = self.body["RequestItems"]

        results: Dict[str, Any] = {
            "ConsumedCapacity": [],
            "Responses": {},
            "UnprocessedKeys": {},
        }

        # Validation: Can only request up to 100 items at the same time
        # Scenario 1: We're requesting more than a 100 keys from a single table
        for table_name, table_request in table_batches.items():
            if len(table_request["Keys"]) > 100:
                raise MockValidationException(
                    "1 validation error detected: Value at 'requestItems."
                    + table_name
                    + ".member.keys' failed to satisfy constraint: Member must have length less than or equal to 100"
                )
        # Scenario 2: We're requesting more than a 100 keys across all tables
        nr_of_keys_across_all_tables = sum(
            [len(req["Keys"]) for _, req in table_batches.items()]
        )
        if nr_of_keys_across_all_tables > 100:
            raise MockValidationException(
                "Too many items requested for the BatchGetItem call"
            )

        result_size: int = 0
        for table_name, table_request in table_batches.items():
            keys = table_request["Keys"]
            if self._contains_duplicates(keys):
                raise MockValidationException(
                    "Provided list of item keys contains duplicates"
                )
            attributes_to_get = table_request.get("AttributesToGet")
            projection_expression = table_request.get("ProjectionExpression")
            expression_attribute_names = table_request.get(
                "ExpressionAttributeNames", {}
            )

            projection_expressions = self._adjust_projection_expression(
                projection_expression, expression_attribute_names
            )

            results["Responses"][table_name] = []
            for key in keys:
                item = self.dynamodb_backend.get_item(
                    table_name, key, projection_expressions
                )
                if item:
                    # A single operation can retrieve up to 16 MB of data [and] returns a partial result if the response size limit is exceeded
                    if result_size + item.size() > (16 * 1024 * 1024):
                        # Result is already getting too big - next results should be part of UnprocessedKeys
                        if table_name not in results["UnprocessedKeys"]:
                            results["UnprocessedKeys"][table_name] = {"Keys": []}
                        results["UnprocessedKeys"][table_name]["Keys"].append(key)
                    else:
                        item_describe = item.describe_attrs(attributes_to_get)
                        results["Responses"][table_name].append(item_describe["Item"])
                        result_size += item.size()

            results["ConsumedCapacity"].append(
                {"CapacityUnits": len(keys), "TableName": table_name}
            )
        return dynamo_json_dump(results)

    def _contains_duplicates(self, keys: List[str]) -> bool:
        unique_keys = []
        for k in keys:
            if k in unique_keys:
                return True
            else:
                unique_keys.append(k)
        return False

    @include_consumed_capacity()
    def query(self) -> str:
        name = self.body["TableName"]
        key_condition_expression = self.body.get("KeyConditionExpression")
        projection_expression = self._get_projection_expression()
        expression_attribute_names = self.body.get("ExpressionAttributeNames", {})
        filter_expression = self._get_filter_expression()
        expression_attribute_values = self._get_expr_attr_values()

        projection_expressions = self._adjust_projection_expression(
            projection_expression, expression_attribute_names
        )

        filter_kwargs = {}

        if key_condition_expression:
            index_name = self.body.get("IndexName")
            schema = self.dynamodb_backend.get_schema(
                table_name=name, index_name=index_name
            )
            hash_key, range_comparison, range_values = parse_expression(
                key_condition_expression=key_condition_expression,
                expression_attribute_names=expression_attribute_names,
                expression_attribute_values=expression_attribute_values,
                schema=schema,
            )
        else:
            # 'KeyConditions': {u'forum_name': {u'ComparisonOperator': u'EQ', u'AttributeValueList': [{u'S': u'the-key'}]}}
            key_conditions = self.body.get("KeyConditions")
            query_filters = self.body.get("QueryFilter")

            if not (key_conditions or query_filters):
                raise MockValidationException(
                    "Either KeyConditions or QueryFilter should be present"
                )

            if key_conditions:
                (
                    hash_key_name,
                    range_key_name,
                ) = self.dynamodb_backend.get_table_keys_name(
                    name, key_conditions.keys()
                )
                for key, value in key_conditions.items():
                    if key not in (hash_key_name, range_key_name):
                        filter_kwargs[key] = value
                if hash_key_name is None:
                    raise ResourceNotFoundException
                hash_key = key_conditions[hash_key_name]["AttributeValueList"][0]
                if len(key_conditions) == 1:
                    range_comparison = None
                    range_values = []
                else:
                    if range_key_name is None and not filter_kwargs:
                        raise MockValidationException("Validation Exception")
                    else:
                        range_condition = key_conditions.get(range_key_name)
                        if range_condition:
                            range_comparison = range_condition["ComparisonOperator"]
                            range_values = range_condition["AttributeValueList"]
                        else:
                            range_comparison = None
                            range_values = []
            if query_filters:
                filter_kwargs.update(query_filters)
        index_name = self.body.get("IndexName")
        exclusive_start_key = self.body.get("ExclusiveStartKey")
        limit = self.body.get("Limit")
        scan_index_forward = self.body.get("ScanIndexForward", True)
        consistent_read = self.body.get("ConsistentRead", False)

        items, scanned_count, last_evaluated_key = self.dynamodb_backend.query(
            name,
            hash_key,
            range_comparison,
            range_values,
            limit,
            exclusive_start_key,
            scan_index_forward,
            projection_expressions,
            index_name=index_name,
            consistent_read=consistent_read,
            expr_names=expression_attribute_names,
            expr_values=expression_attribute_values,
            filter_expression=filter_expression,
            **filter_kwargs,
        )

        result: Dict[str, Any] = {
            "Count": len(items),
            "ScannedCount": scanned_count,
        }

        if self.body.get("Select", "").upper() != "COUNT":
            result["Items"] = [item.attrs for item in items]

        if last_evaluated_key is not None:
            result["LastEvaluatedKey"] = last_evaluated_key

        return dynamo_json_dump(result)

    def _adjust_projection_expression(
        self, projection_expression: Optional[str], expr_attr_names: Dict[str, str]
    ) -> List[List[str]]:
        """
        lvl1.lvl2.attr1,lvl1.attr2 --> [["lvl1", "lvl2", "attr1"], ["lvl1", "attr2]]
        """

        def _adjust(expression: str) -> str:
            return (expr_attr_names or {}).get(expression, expression)

        if projection_expression:
            expressions = [x.strip() for x in projection_expression.split(",")]
            for expression in expressions:
                check_projection_expression(expression)
            return [
                [_adjust(expr) for expr in nested_expr.split(".")]
                for nested_expr in expressions
            ]

        return []

    @include_consumed_capacity()
    def scan(self) -> str:
        name = self.body["TableName"]

        filters = {}
        scan_filters = self.body.get("ScanFilter", {})
        for attribute_name, scan_filter in scan_filters.items():
            # Keys are attribute names. Values are tuples of (comparison,
            # comparison_value)
            comparison_operator = scan_filter["ComparisonOperator"]
            comparison_values = scan_filter.get("AttributeValueList", [])
            filters[attribute_name] = (comparison_operator, comparison_values)

        filter_expression = self._get_filter_expression()
        expression_attribute_values = self._get_expr_attr_values()
        expression_attribute_names = self.body.get("ExpressionAttributeNames", {})
        projection_expression = self._get_projection_expression()
        exclusive_start_key = self.body.get("ExclusiveStartKey")
        limit = self.body.get("Limit")
        index_name = self.body.get("IndexName")
        consistent_read = self.body.get("ConsistentRead", False)

        projection_expressions = self._adjust_projection_expression(
            projection_expression, expression_attribute_names
        )

        try:
            items, scanned_count, last_evaluated_key = self.dynamodb_backend.scan(
                name,
                filters,
                limit,
                exclusive_start_key,
                filter_expression,
                expression_attribute_names,
                expression_attribute_values,
                index_name,
                consistent_read,
                projection_expressions,
            )
        except ValueError as err:
            raise MockValidationException(f"Bad Filter Expression: {err}")

        result = {
            "Count": len(items),
            "Items": [item.attrs for item in items],
            "ScannedCount": scanned_count,
        }
        if last_evaluated_key is not None:
            result["LastEvaluatedKey"] = last_evaluated_key
        return dynamo_json_dump(result)

    def delete_item(self) -> str:
        name = self.body["TableName"]
        key = self.body["Key"]
        return_values = self.body.get("ReturnValues", "NONE")
        if return_values not in ("ALL_OLD", "NONE"):
            raise MockValidationException("Return values set to invalid value")

        self.dynamodb_backend.get_table(name)

        # Attempt to parse simple ConditionExpressions into an Expected
        # expression
        condition_expression = self.body.get("ConditionExpression")
        expression_attribute_names = self.body.get("ExpressionAttributeNames", {})
        expression_attribute_values = self._get_expr_attr_values()

        item = self.dynamodb_backend.delete_item(
            name,
            key,
            expression_attribute_names,
            expression_attribute_values,
            condition_expression,
        )

        if item and return_values == "ALL_OLD":
            item_dict = item.to_json()
        else:
            item_dict = {"Attributes": {}}
        item_dict["ConsumedCapacityUnits"] = 0.5
        return dynamo_json_dump(item_dict)

    def update_item(self) -> str:
        name = self.body["TableName"]
        key = self.body["Key"]
        return_values = self.body.get("ReturnValues", "NONE")
        update_expression = self.body.get("UpdateExpression", "").strip()
        attribute_updates = self.body.get("AttributeUpdates")
        if update_expression and attribute_updates:
            raise MockValidationException(
                "Can not use both expression and non-expression parameters in the same request: Non-expression parameters: {AttributeUpdates} Expression parameters: {UpdateExpression}"
            )
        return_values_on_condition_check_failure = self.body.get(
            "ReturnValuesOnConditionCheckFailure"
        )
        # We need to copy the item in order to avoid it being modified by the update_item operation
        existing_item = copy.deepcopy(self.dynamodb_backend.get_item(name, key))
        if existing_item:
            existing_attributes = existing_item.to_json()["Attributes"]
        else:
            existing_attributes = {}

        if return_values not in (
            "NONE",
            "ALL_OLD",
            "ALL_NEW",
            "UPDATED_OLD",
            "UPDATED_NEW",
        ):
            raise MockValidationException("Return values set to invalid value")

        if "Expected" in self.body:
            expected = self.body["Expected"]
        else:
            expected = None

        # Attempt to parse simple ConditionExpressions into an Expected
        # expression
        condition_expression = self.body.get("ConditionExpression")
        expression_attribute_names = self.body.get("ExpressionAttributeNames", {})
        expression_attribute_values = self._get_expr_attr_values()

        item = self.dynamodb_backend.update_item(
            name,
            key,
            update_expression=update_expression,
            attribute_updates=attribute_updates,
            expression_attribute_names=expression_attribute_names,
            expression_attribute_values=expression_attribute_values,
            expected=expected,
            condition_expression=condition_expression,
            return_values_on_condition_check_failure=return_values_on_condition_check_failure,
        )

        item_dict = item.to_json()
        item_dict["ConsumedCapacity"] = {"TableName": name, "CapacityUnits": 0.5}
        unchanged_attributes = {
            k
            for k in existing_attributes.keys()
            if existing_attributes[k] == item_dict["Attributes"].get(k)
        }
        changed_attributes = (
            set(existing_attributes.keys())
            .union(item_dict["Attributes"].keys())
            .difference(unchanged_attributes)
        )

        if return_values == "NONE":
            item_dict["Attributes"] = {}
        elif return_values == "ALL_OLD":
            item_dict["Attributes"] = existing_attributes
        elif return_values == "UPDATED_OLD":
            item_dict["Attributes"] = {
                k: v for k, v in existing_attributes.items() if k in changed_attributes
            }
        elif return_values == "UPDATED_NEW":
            item_dict["Attributes"] = self._build_updated_new_attributes(
                existing_attributes, item_dict["Attributes"]
            )
        return dynamo_json_dump(item_dict)

    def _get_expr_attr_values(self) -> Dict[str, Dict[str, str]]:
        values = self.body.get("ExpressionAttributeValues", {})
        for key in values.keys():
            if not key.startswith(":"):
                raise MockValidationException(
                    f'ExpressionAttributeValues contains invalid key: Syntax error; key: "{key}"'
                )
        return values

    def _build_updated_new_attributes(self, original: Any, changed: Any) -> Any:
        if type(changed) != type(original):
            return changed
        else:
            if isinstance(changed, dict):
                return {
                    key: self._build_updated_new_attributes(
                        original.get(key, None), changed[key]
                    )
                    for key in changed.keys()
                    if key not in original or changed[key] != original[key]
                }
            elif type(changed) in (set, list):
                if len(changed) != len(original):
                    return changed
                else:
                    any_element_has_changed = any(
                        changed[index] != original[index]
                        for index in range(len(changed))
                    )
                    return changed if any_element_has_changed else original
            else:
                return changed

    def describe_limits(self) -> str:
        return json.dumps(
            {
                "AccountMaxReadCapacityUnits": 20000,
                "TableMaxWriteCapacityUnits": 10000,
                "AccountMaxWriteCapacityUnits": 20000,
                "TableMaxReadCapacityUnits": 10000,
            }
        )

    def update_time_to_live(self) -> str:
        name = self.body["TableName"]
        ttl_spec = self.body["TimeToLiveSpecification"]

        self.dynamodb_backend.update_time_to_live(name, ttl_spec)

        return json.dumps({"TimeToLiveSpecification": ttl_spec})

    def describe_time_to_live(self) -> str:
        name = self.body["TableName"]

        ttl_spec = self.dynamodb_backend.describe_time_to_live(name)

        return json.dumps({"TimeToLiveDescription": ttl_spec})

    def transact_get_items(self) -> str:
        transact_items = self.body["TransactItems"]
        responses: List[Dict[str, Any]] = list()

        if len(transact_items) > TRANSACTION_MAX_ITEMS:
            msg = "1 validation error detected: Value '["
            err_list = list()
            request_id = 268435456
            for _ in transact_items:
                request_id += 1
                hex_request_id = format(request_id, "x")
                err_list.append(
                    "com.amazonaws.dynamodb.v20120810.TransactGetItem@%s"
                    % hex_request_id
                )
            msg += ", ".join(err_list)
            msg += (
                "'] at 'transactItems' failed to satisfy constraint: "
                "Member must have length less than or equal to %s"
                % TRANSACTION_MAX_ITEMS
            )

            raise MockValidationException(msg)

        ret_consumed_capacity = self.body.get("ReturnConsumedCapacity", "NONE")
        consumed_capacity: Dict[str, Any] = dict()

        for transact_item in transact_items:
            table_name = transact_item["Get"]["TableName"]
            key = transact_item["Get"]["Key"]
            item = self.dynamodb_backend.get_item(table_name, key)

            if not item:
                responses.append({})
                continue

            item_describe = item.describe_attrs(attributes=None)
            responses.append(item_describe)

            table_capacity = consumed_capacity.get(table_name, {})
            table_capacity["TableName"] = table_name
            capacity_units = table_capacity.get("CapacityUnits", 0) + 2.0
            table_capacity["CapacityUnits"] = capacity_units
            read_capacity_units = table_capacity.get("ReadCapacityUnits", 0) + 2.0
            table_capacity["ReadCapacityUnits"] = read_capacity_units
            consumed_capacity[table_name] = table_capacity

            if ret_consumed_capacity == "INDEXES":
                table_capacity["Table"] = {
                    "CapacityUnits": capacity_units,
                    "ReadCapacityUnits": read_capacity_units,
                }

        result = dict()
        result.update({"Responses": responses})
        if ret_consumed_capacity != "NONE":
            result.update({"ConsumedCapacity": [v for v in consumed_capacity.values()]})

        return dynamo_json_dump(result)

    def transact_write_items(self) -> str:
        transact_items = self.body["TransactItems"]
        # Validate first - we should error before we start the transaction
        for item in transact_items:
            if "Put" in item:
                item_attrs = item["Put"]["Item"]
                table = self.dynamodb_backend.get_table(item["Put"]["TableName"])
                validate_put_has_empty_keys(item_attrs, table)
            if "Update" in item:
                item_attrs = item["Update"]["Key"]
                table = self.dynamodb_backend.get_table(item["Update"]["TableName"])
                validate_put_has_empty_keys(
                    item_attrs,
                    table,
                    custom_error_msg="One or more parameter values are not valid. The AttributeValue for a key attribute cannot contain an empty string value. Key: {}",
                )
        self.dynamodb_backend.transact_write_items(transact_items)
        response: Dict[str, Any] = {"ConsumedCapacity": [], "ItemCollectionMetrics": {}}
        return dynamo_json_dump(response)

    def describe_continuous_backups(self) -> str:
        name = self.body["TableName"]

        response = self.dynamodb_backend.describe_continuous_backups(name)

        return json.dumps({"ContinuousBackupsDescription": response})

    def update_continuous_backups(self) -> str:
        name = self.body["TableName"]
        point_in_time_spec = self.body["PointInTimeRecoverySpecification"]

        response = self.dynamodb_backend.update_continuous_backups(
            name, point_in_time_spec
        )

        return json.dumps({"ContinuousBackupsDescription": response})

    def list_backups(self) -> str:
        body = self.body
        table_name = body.get("TableName")
        backups = self.dynamodb_backend.list_backups(table_name)
        response = {"BackupSummaries": [backup.summary for backup in backups]}
        return dynamo_json_dump(response)

    def create_backup(self) -> str:
        body = self.body
        table_name = body.get("TableName")
        backup_name = body.get("BackupName")
        backup = self.dynamodb_backend.create_backup(table_name, backup_name)
        response = {"BackupDetails": backup.details}
        return dynamo_json_dump(response)

    def delete_backup(self) -> str:
        body = self.body
        backup_arn = body.get("BackupArn")
        backup = self.dynamodb_backend.delete_backup(backup_arn)
        response = {"BackupDescription": backup.description}
        return dynamo_json_dump(response)

    def describe_backup(self) -> str:
        body = self.body
        backup_arn = body.get("BackupArn")
        backup = self.dynamodb_backend.describe_backup(backup_arn)
        response = {"BackupDescription": backup.description}
        return dynamo_json_dump(response)

    def restore_table_from_backup(self) -> str:
        body = self.body
        target_table_name = body.get("TargetTableName")
        backup_arn = body.get("BackupArn")
        restored_table = self.dynamodb_backend.restore_table_from_backup(
            target_table_name, backup_arn
        )
        return dynamo_json_dump(restored_table.describe())

    def restore_table_to_point_in_time(self) -> str:
        body = self.body
        target_table_name = body.get("TargetTableName")
        source_table_name = body.get("SourceTableName")
        restored_table = self.dynamodb_backend.restore_table_to_point_in_time(
            target_table_name, source_table_name
        )
        return dynamo_json_dump(restored_table.describe())

    def execute_statement(self) -> str:
        stmt = self.body.get("Statement", "")
        parameters = self.body.get("Parameters", [])
        items = self.dynamodb_backend.execute_statement(
            statement=stmt, parameters=parameters
        )
        return dynamo_json_dump({"Items": items})

    def execute_transaction(self) -> str:
        stmts = self.body.get("TransactStatements", [])
        items = self.dynamodb_backend.execute_transaction(stmts)
        return dynamo_json_dump({"Responses": items})

    def batch_execute_statement(self) -> str:
        stmts = self.body.get("Statements", [])
        items = self.dynamodb_backend.batch_execute_statement(stmts)
        return dynamo_json_dump({"Responses": items})

    def import_table(self) -> str:
        params = self.body
        s3_source = params.get("S3BucketSource")
        input_format = params.get("InputFormat") or "DYNAMODB_JSON"
        compression_type = params.get("InputCompressionType") or "NONE"
        table_parameters = params.get("TableCreationParameters")
        table_name = table_parameters["TableName"]
        table_attrs = table_parameters["AttributeDefinitions"]
        table_schema = table_parameters["KeySchema"]
        table_billing = table_parameters.get("BillingMode")
        table_throughput = table_parameters.get("ProvisionedThroughput")
        global_indexes = table_parameters.get("GlobalSecondaryIndexes")

        self._validate_table_creation(
            billing_mode=table_billing,
            throughput=table_throughput,
            key_schema=table_schema,
            global_indexes=global_indexes,
            local_secondary_indexes=None,
            attr=table_attrs,
        )

        import_table = self.dynamodb_backend.import_table(
            s3_source=s3_source,
            input_format=input_format,
            compression_type=compression_type,
            table_name=table_name,
            billing_mode=table_billing or "PROVISIONED",
            throughput=table_throughput,
            key_schema=table_schema,
            global_indexes=global_indexes,
            attrs=table_attrs,
        )
        return json.dumps({"ImportTableDescription": import_table.response()})

    def describe_import(self) -> str:
        import_arn = self.body["ImportArn"]
        import_table = self.dynamodb_backend.describe_import(import_arn)
        return json.dumps({"ImportTableDescription": import_table.response()})
