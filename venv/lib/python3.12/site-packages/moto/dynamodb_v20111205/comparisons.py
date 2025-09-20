from typing import Any, Callable

# TODO add tests for all of these
COMPARISON_FUNCS = {
    "EQ": lambda item_value, test_value: item_value == test_value,
    "NE": lambda item_value, test_value: item_value != test_value,
    "LE": lambda item_value, test_value: item_value <= test_value,
    "LT": lambda item_value, test_value: item_value < test_value,
    "GE": lambda item_value, test_value: item_value >= test_value,
    "GT": lambda item_value, test_value: item_value > test_value,
    "NULL": lambda item_value: item_value is None,
    "NOT_NULL": lambda item_value: item_value is not None,
    "CONTAINS": lambda item_value, test_value: test_value in item_value,
    "NOT_CONTAINS": lambda item_value, test_value: test_value not in item_value,
    "BEGINS_WITH": lambda item_value, test_value: item_value.startswith(test_value),
    "IN": lambda item_value, *test_values: item_value in test_values,
    "BETWEEN": lambda item_value, lower_test_value, upper_test_value: lower_test_value
    <= item_value
    <= upper_test_value,
}


def get_comparison_func(range_comparison: str) -> Callable[..., Any]:
    return COMPARISON_FUNCS.get(range_comparison)  # type: ignore[return-value]
