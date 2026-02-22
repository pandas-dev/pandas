from typing import TYPE_CHECKING, Any, Optional

from moto.dynamodb.exceptions import MockValidationException, ResourceNotFoundException

if TYPE_CHECKING:
    from py_partiql_parser import QueryMetadata


def query(
    statement: str, source_data: dict[str, list[Any]], parameters: list[dict[str, Any]]
) -> tuple[
    list[dict[str, Any]],
    dict[str, list[tuple[Optional[dict[str, Any]], Optional[dict[str, Any]]]]],
]:
    from py_partiql_parser import DynamoDBStatementParser
    from py_partiql_parser.exceptions import DocumentNotFoundException

    try:
        return DynamoDBStatementParser(source_data).parse(statement, parameters)
    except DocumentNotFoundException as dnfe:
        if "." in dnfe.name:
            table_name = dnfe.name.split(".")[0]
            if table_name in source_data:
                raise MockValidationException(
                    message="The table does not have the specified index"
                )
        raise ResourceNotFoundException()


def get_query_metadata(statement: str) -> "QueryMetadata":
    from py_partiql_parser import DynamoDBStatementParser

    return DynamoDBStatementParser.get_query_metadata(query=statement)
