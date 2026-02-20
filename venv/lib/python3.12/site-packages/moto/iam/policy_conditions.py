from collections.abc import Iterator
from typing import Union


def string_equals_operation(
    expected_value: Union[list[str], str], actual_value: str
) -> bool:
    if isinstance(expected_value, str):
        return actual_value == expected_value

    if isinstance(expected_value, list):
        return actual_value in expected_value

    return False


CONDITION_OPERATIONS = {"StringEquals": string_equals_operation}


class TrustCondition:
    """Single condition object"""

    def __init__(
        self,
        condition: str,
        expected_value: list[Union[list[str], str]],
        expected_value_source: list[str],
    ) -> None:
        self._condition = condition

        self._expected_value_source = expected_value_source
        self._expected_value = expected_value

    @property
    def context_keys(self) -> list[str]:
        return self._expected_value_source

    def verify_condition(self, actual_values: list[str]) -> bool:
        verify_action = CONDITION_OPERATIONS.get(self._condition)

        # If the evaluation operation is not supported, ignore condition evaluation
        if not verify_action:
            return True

        for index, actual_value in enumerate(actual_values):
            if not verify_action(self._expected_value[index], actual_value):
                return False

        return True


ConditionData = dict[str, dict[str, Union[str, list[str]]]]


class TrustRelationShipConditions:
    """Trust relationship conditions segment model"""

    def __init__(self, conditions: ConditionData) -> None:
        self._conditions: list[TrustCondition] = []
        for conditions_operation, condition_values in conditions.items():
            expected_value_source: list[str] = list(condition_values.keys())
            expected_value_data: list[Union[str, list[str]]] = list(
                condition_values.values()
            )

            self._conditions.append(
                TrustCondition(
                    conditions_operation, expected_value_data, expected_value_source
                )
            )

    def __iter__(self) -> Iterator[TrustCondition]:
        return self._conditions.__iter__()
