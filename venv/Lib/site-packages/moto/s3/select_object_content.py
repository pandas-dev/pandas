from typing import Any


def parse_query(text_input: str, query: str) -> tuple[list[dict[str, Any]], int]:
    from py_partiql_parser import S3SelectParser

    parser = S3SelectParser(source_data=text_input)
    result = parser.parse(query)
    return result, parser.bytes_scanned
