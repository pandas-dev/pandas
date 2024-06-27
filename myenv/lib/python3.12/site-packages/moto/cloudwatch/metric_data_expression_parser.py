from typing import Any, Dict, List, SupportsFloat, Tuple


def parse_expression(
    expression: str, results: List[Dict[str, Any]]
) -> Tuple[List[SupportsFloat], List[str]]:
    values: List[SupportsFloat] = []
    timestamps: List[str] = []
    for result in results:
        if result.get("id") == expression:
            values.extend(result["vals"])
            timestamps.extend(result["timestamps"])
    return values, timestamps
