from datetime import datetime
from typing import Any, Dict, List, SupportsFloat, Tuple


def parse_expression(
    expression: str, results: List[Dict[str, Any]]
) -> Tuple[List[SupportsFloat], List[datetime]]:
    values: List[SupportsFloat] = []
    timestamps: List[datetime] = []
    for result in results:
        if result.get("id") == expression:
            values.extend(result["values"])
            timestamps.extend(result["timestamps"])
    return values, timestamps
