from typing import TYPE_CHECKING, Any, Dict, List, Optional, Tuple

if TYPE_CHECKING:
    from py_partiql_parser import QueryMetadata


def query(
    statement: str, source_data: Dict[str, str], parameters: List[Dict[str, Any]]
) -> Tuple[
    List[Dict[str, Any]],
    Dict[str, List[Tuple[Optional[Dict[str, Any]], Optional[Dict[str, Any]]]]],
]:
    from py_partiql_parser import DynamoDBStatementParser

    return DynamoDBStatementParser(source_data).parse(statement, parameters)


def get_query_metadata(statement: str) -> "QueryMetadata":
    from py_partiql_parser import DynamoDBStatementParser

    return DynamoDBStatementParser.get_query_metadata(query=statement)
