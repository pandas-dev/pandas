import json
from typing import Any, Dict, List, Optional

from moto.core.exceptions import JsonRESTError
from moto.dynamodb.limits import HASH_KEY_MAX_LENGTH, RANGE_KEY_MAX_LENGTH

ERROR_TYPE_PREFIX = "com.amazonaws.dynamodb.v20120810#"


class DynamodbException(JsonRESTError):
    pass


class MockValidationException(DynamodbException):
    error_type = ERROR_TYPE_PREFIX + "ValidationException"

    def __init__(self, message: str):
        super().__init__(MockValidationException.error_type, message=message)
        self.exception_msg = message


class KeyIsEmptyStringException(MockValidationException):
    def __init__(self, empty_key: str):
        super().__init__(
            message=f"One or more parameter values are not valid. The AttributeValue for a key attribute cannot contain an empty string value. Key: {empty_key}"
        )


class InvalidIndexNameError(MockValidationException):
    pass


class InvalidUpdateExpressionInvalidDocumentPath(MockValidationException):
    invalid_update_expression_msg = (
        "The document path provided in the update expression is invalid for update"
    )

    def __init__(self) -> None:
        super().__init__(self.invalid_update_expression_msg)


class InvalidUpdateExpression(MockValidationException):
    invalid_update_expr_msg = "Invalid UpdateExpression: {update_expression_error}"

    def __init__(self, update_expression_error: str):
        self.update_expression_error = update_expression_error
        super().__init__(
            self.invalid_update_expr_msg.format(
                update_expression_error=update_expression_error
            )
        )


class InvalidConditionExpression(MockValidationException):
    invalid_condition_expr_msg = (
        "Invalid ConditionExpression: {condition_expression_error}"
    )

    def __init__(self, condition_expression_error: str):
        self.condition_expression_error = condition_expression_error
        super().__init__(
            self.invalid_condition_expr_msg.format(
                condition_expression_error=condition_expression_error
            )
        )


class ConditionAttributeIsReservedKeyword(InvalidConditionExpression):
    attribute_is_keyword_msg = (
        "Attribute name is a reserved keyword; reserved keyword: {keyword}"
    )

    def __init__(self, keyword: str):
        self.keyword = keyword
        super().__init__(self.attribute_is_keyword_msg.format(keyword=keyword))


class AttributeDoesNotExist(MockValidationException):
    attr_does_not_exist_msg = (
        "The provided expression refers to an attribute that does not exist in the item"
    )

    def __init__(self) -> None:
        super().__init__(self.attr_does_not_exist_msg)


class ProvidedKeyDoesNotExist(MockValidationException):
    provided_key_does_not_exist_msg = (
        "The provided key element does not match the schema"
    )

    def __init__(self) -> None:
        super().__init__(self.provided_key_does_not_exist_msg)


class ExpressionAttributeNameNotDefined(InvalidUpdateExpression):
    name_not_defined_msg = "An expression attribute name used in the document path is not defined; attribute name: {n}"

    def __init__(self, attribute_name: str):
        self.not_defined_attribute_name = attribute_name
        super().__init__(self.name_not_defined_msg.format(n=attribute_name))


class AttributeIsReservedKeyword(InvalidUpdateExpression):
    attribute_is_keyword_msg = (
        "Attribute name is a reserved keyword; reserved keyword: {keyword}"
    )

    def __init__(self, keyword: str):
        self.keyword = keyword
        super().__init__(self.attribute_is_keyword_msg.format(keyword=keyword))


class ExpressionAttributeValueNotDefined(InvalidUpdateExpression):
    attr_value_not_defined_msg = "An expression attribute value used in expression is not defined; attribute value: {attribute_value}"

    def __init__(self, attribute_value: str):
        self.attribute_value = attribute_value
        super().__init__(
            self.attr_value_not_defined_msg.format(attribute_value=attribute_value)
        )


class UpdateExprSyntaxError(InvalidUpdateExpression):
    update_expr_syntax_error_msg = "Syntax error; {error_detail}"

    def __init__(self, error_detail: str):
        self.error_detail = error_detail
        super().__init__(
            self.update_expr_syntax_error_msg.format(error_detail=error_detail)
        )


class InvalidTokenException(UpdateExprSyntaxError):
    token_detail_msg = 'token: "{token}", near: "{near}"'

    def __init__(self, token: str, near: str):
        self.token = token
        self.near = near
        super().__init__(self.token_detail_msg.format(token=token, near=near))


class InvalidExpressionAttributeNameKey(MockValidationException):
    invalid_expr_attr_name_msg = (
        'ExpressionAttributeNames contains invalid key: Syntax error; key: "{key}"'
    )

    def __init__(self, key: str):
        self.key = key
        super().__init__(self.invalid_expr_attr_name_msg.format(key=key))


class ItemSizeTooLarge(MockValidationException):
    item_size_too_large_msg = "Item size has exceeded the maximum allowed size"

    def __init__(self) -> None:
        super().__init__(self.item_size_too_large_msg)


class ItemSizeToUpdateTooLarge(MockValidationException):
    item_size_to_update_too_large_msg = (
        "Item size to update has exceeded the maximum allowed size"
    )

    def __init__(self) -> None:
        super().__init__(self.item_size_to_update_too_large_msg)


class HashKeyTooLong(MockValidationException):
    # deliberately no space between of and {lim}
    key_too_large_msg = f"One or more parameter values were invalid: Size of hashkey has exceeded the maximum size limit of{HASH_KEY_MAX_LENGTH} bytes"

    def __init__(self) -> None:
        super().__init__(self.key_too_large_msg)


class RangeKeyTooLong(MockValidationException):
    key_too_large_msg = f"One or more parameter values were invalid: Aggregated size of all range keys has exceeded the size limit of {RANGE_KEY_MAX_LENGTH} bytes"

    def __init__(self) -> None:
        super().__init__(self.key_too_large_msg)


class IncorrectOperandType(InvalidUpdateExpression):
    inv_operand_msg = "Incorrect operand type for operator or function; operator or function: {f}, operand type: {t}"

    def __init__(self, operator_or_function: str, operand_type: str):
        self.operator_or_function = operator_or_function
        self.operand_type = operand_type
        super().__init__(
            self.inv_operand_msg.format(f=operator_or_function, t=operand_type)
        )


class IncorrectDataType(MockValidationException):
    inc_data_type_msg = "An operand in the update expression has an incorrect data type"

    def __init__(self) -> None:
        super().__init__(self.inc_data_type_msg)


class ConditionalCheckFailed(DynamodbException):
    error_type = ERROR_TYPE_PREFIX + "ConditionalCheckFailedException"

    def __init__(
        self, msg: Optional[str] = None, item: Optional[Dict[str, Any]] = None
    ):
        _msg = msg or "The conditional request failed"
        super().__init__(ConditionalCheckFailed.error_type, _msg)
        if item:
            self.description = json.dumps(
                {
                    "__type": ConditionalCheckFailed.error_type,
                    # Note the uppercase Message
                    # This ensures the message is only part of the 'error': {'message': .., 'code': ..}
                    "Message": _msg,
                    "Item": item,
                }
            )
        else:
            self.description = json.dumps(
                {
                    "__type": ConditionalCheckFailed.error_type,
                    # Note that lowercase 'message'
                    # This ensures that 'message' is a top-level field in the response
                    # (in addition to being part of the 'error': {'message': .., 'code': ..}
                    "message": _msg,
                }
            )


class TransactionCanceledException(DynamodbException):
    cancel_reason_msg = "Transaction cancelled, please refer cancellation reasons for specific reasons [{}]"
    error_type = "com.amazonaws.dynamodb.v20120810#TransactionCanceledException"

    def __init__(self, errors: List[Any]):
        msg = self.cancel_reason_msg.format(
            ", ".join([str(code) for code, _, _ in errors])
        )
        super().__init__(
            error_type=TransactionCanceledException.error_type, message=msg
        )
        reasons = []
        for code, message, item in errors:
            r = {"Code": code, "Message": message} if code else {"Code": "None"}
            if item:
                r["Item"] = item
            reasons.append(r)

        self.description = json.dumps(
            {
                "__type": TransactionCanceledException.error_type,
                "CancellationReasons": reasons,
                "Message": msg,
            }
        )


class MultipleTransactionsException(MockValidationException):
    msg = "Transaction request cannot include multiple operations on one item"

    def __init__(self) -> None:
        super().__init__(self.msg)


class TooManyTransactionsException(MockValidationException):
    msg = (
        "1 validation error detected at 'transactItems' failed to satisfy constraint: "
        "Member must have length less than or equal to 100."
    )

    def __init__(self) -> None:
        super().__init__(self.msg)


class EmptyKeyAttributeException(MockValidationException):
    empty_str_msg = "One or more parameter values were invalid: An AttributeValue may not contain an empty string"
    # AWS has a different message for empty index keys
    empty_index_msg = "One or more parameter values are not valid. The update expression attempted to update a secondary index key to a value that is not supported. The AttributeValue for a key attribute cannot contain an empty string value."

    def __init__(self, key_in_index: bool = False):
        super().__init__(self.empty_index_msg if key_in_index else self.empty_str_msg)


class UpdateHashRangeKeyException(MockValidationException):
    msg = "One or more parameter values were invalid: Cannot update attribute {}. This attribute is part of the key"

    def __init__(self, key_name: str):
        super().__init__(self.msg.format(key_name))


class InvalidAttributeTypeError(MockValidationException):
    msg = "One or more parameter values were invalid: Type mismatch for key {} expected: {} actual: {}"

    def __init__(
        self, name: Optional[str], expected_type: Optional[str], actual_type: str
    ):
        super().__init__(self.msg.format(name, expected_type, actual_type))


class DuplicateUpdateExpression(InvalidUpdateExpression):
    def __init__(self, names: List[str]):
        super().__init__(
            f"Two document paths overlap with each other; must remove or rewrite one of these paths; path one: [{names[0]}], path two: [{names[1]}]"
        )


class TooManyAddClauses(InvalidUpdateExpression):
    msg = 'The "ADD" section can only be used once in an update expression;'

    def __init__(self) -> None:
        super().__init__(self.msg)


class ResourceNotFoundException(JsonRESTError):
    def __init__(self, msg: Optional[str] = None, table_name: Optional[str] = None):
        err = ERROR_TYPE_PREFIX + "ResourceNotFoundException"
        default_msg = "Requested resource not found"
        if table_name is not None:
            default_msg += f": Table: {table_name} not found"
        super().__init__(err, msg or default_msg)


class TableNotFoundException(JsonRESTError):
    def __init__(self, name: str):
        err = ERROR_TYPE_PREFIX + "TableNotFoundException"
        super().__init__(err, f"Table not found: {name}")


class SourceTableNotFoundException(JsonRESTError):
    def __init__(self, source_table_name: str):
        er = ERROR_TYPE_PREFIX + "SourceTableNotFoundException"
        super().__init__(er, f"Source table not found: {source_table_name}")


class BackupNotFoundException(JsonRESTError):
    def __init__(self, backup_arn: str):
        er = ERROR_TYPE_PREFIX + "BackupNotFoundException"
        super().__init__(er, f"Backup not found: {backup_arn}")


class TableAlreadyExistsException(JsonRESTError):
    def __init__(self, target_table_name: str):
        er = ERROR_TYPE_PREFIX + "TableAlreadyExistsException"
        super().__init__(er, f"Table already exists: {target_table_name}")


class ResourceInUseException(JsonRESTError):
    def __init__(self, msg: Optional[str] = None) -> None:
        er = ERROR_TYPE_PREFIX + "ResourceInUseException"
        super().__init__(er, msg or "Resource in use")


class StreamAlreadyEnabledException(JsonRESTError):
    def __init__(self) -> None:
        er = ERROR_TYPE_PREFIX + "ResourceInUseException"
        super().__init__(er, "Cannot enable stream")


class InvalidConversion(JsonRESTError):
    def __init__(self) -> None:
        er = "SerializationException"
        super().__init__(er, "NUMBER_VALUE cannot be converted to String")


class TransactWriteSingleOpException(MockValidationException):
    there_can_be_only_one = (
        "TransactItems can only contain one of Check, Put, Update or Delete"
    )

    def __init__(self) -> None:
        super().__init__(self.there_can_be_only_one)


class SerializationException(DynamodbException):
    def __init__(self, msg: str):
        super().__init__(error_type="SerializationException", message=msg)


class UnknownKeyType(MockValidationException):
    def __init__(self, key_type: str, position: str):
        msg = f"1 validation error detected: Value '{key_type}' at '{position}' failed to satisfy constraint: Member must satisfy enum value set: [HASH, RANGE]"
        super().__init__(msg)


class DeletionProtectedException(MockValidationException):
    def __init__(self, table_name: str):
        msg = f"1 validation error detected: Table '{table_name}' can't be deleted while DeletionProtectionEnabled is set to True"
        super().__init__(msg)
