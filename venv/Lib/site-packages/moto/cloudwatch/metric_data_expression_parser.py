from datetime import datetime
from typing import Any, SupportsFloat


def parse_expression(
    expression: str, results: list[dict[str, Any]]
) -> tuple[list[SupportsFloat], list[datetime]]:
    values: list[SupportsFloat] = []
    timestamps: list[datetime] = []
    for result in results:
        if result.get("id") == expression:
            values.extend(result["values"])
            timestamps.extend(result["timestamps"])
    return values, timestamps
