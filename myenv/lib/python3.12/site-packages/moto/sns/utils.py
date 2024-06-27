import json
import re
from typing import Any, Callable, Dict, Iterable, List, Optional, Tuple, Union

from moto.moto_api._internal import mock_random
from moto.utilities.utils import get_partition

E164_REGEX = re.compile(r"^\+?[1-9]\d{1,14}$")


def make_arn_for_topic(account_id: str, name: str, region_name: str) -> str:
    return f"arn:{get_partition(region_name)}:sns:{region_name}:{account_id}:{name}"


def make_arn_for_subscription(topic_arn: str) -> str:
    subscription_id = mock_random.uuid4()
    return f"{topic_arn}:{subscription_id}"


def is_e164(number: str) -> bool:
    return E164_REGEX.match(number) is not None


class FilterPolicyMatcher:
    class CheckException(Exception):
        pass

    def __init__(self, filter_policy: Dict[str, Any], filter_policy_scope: str):
        self.filter_policy = filter_policy
        self.filter_policy_scope = (
            filter_policy_scope
            if filter_policy_scope is not None
            else "MessageAttributes"
        )

        if self.filter_policy_scope not in ("MessageAttributes", "MessageBody"):
            raise FilterPolicyMatcher.CheckException(
                f"Unsupported filter_policy_scope: {filter_policy_scope}"
            )

    def matches(
        self, message_attributes: Optional[Dict[str, Any]], message: str
    ) -> bool:
        if not self.filter_policy:
            return True

        if self.filter_policy_scope == "MessageAttributes":
            if message_attributes is None:
                message_attributes = {}

            return FilterPolicyMatcher._attributes_based_match(
                message_attributes, source=self.filter_policy
            )
        else:
            try:
                message_dict = json.loads(message)
            except ValueError:
                return False
            return self._body_based_match(message_dict)

    @staticmethod
    def _attributes_based_match(  # type: ignore[misc]
        message_attributes: Dict[str, Any], source: Dict[str, Any]
    ) -> bool:
        return all(
            FilterPolicyMatcher._field_match(field, rules, message_attributes)
            for field, rules in source.items()
        )

    def _body_based_match(self, message_dict: Dict[str, Any]) -> bool:
        try:
            checks = self._compute_body_checks(self.filter_policy, message_dict)
        except FilterPolicyMatcher.CheckException:
            return False

        return self._perform_body_checks(checks)

    def _perform_body_checks(self, check: Any) -> bool:
        # If the checks are a list, only a single elem has to pass
        # otherwise all the entries have to pass

        if isinstance(check, tuple):
            if len(check) == 2:
                # (any|all, checks)
                aggregate_func, checks = check
                return aggregate_func(
                    self._perform_body_checks(single_check) for single_check in checks
                )
            elif len(check) == 3:
                field, rules, dict_body = check
                return FilterPolicyMatcher._field_match(field, rules, dict_body, False)

        raise FilterPolicyMatcher.CheckException(f"Check is not a tuple: {str(check)}")

    def _compute_body_checks(
        self,
        filter_policy: Dict[str, Union[Dict[str, Any], List[Any]]],
        message_body: Union[Dict[str, Any], List[Any]],
    ) -> Tuple[Callable[[Iterable[Any]], bool], Any]:
        """
        Generate (possibly nested) list of checks to be performed based on the filter policy
        Returned list is of format (any|all, checks), where first elem defines what aggregation should be used in checking
        and the second argument is a list containing sublists of the same format or concrete checks (field, rule, body): Tuple[str, List[Any], Dict[str, Any]]

        All the checks returned by this function will only require one-level-deep entry into dict in _field_match function
        This is done this way to simplify the actual check logic and keep it as close as possible between MessageAttributes and MessageBody

        Given message_body:
        {"Records": [
            {
                "eventName": "ObjectCreated:Put",
            },
            {
                "eventName": "ObjectCreated:Delete",
            },
        ]}

        and filter policy:
        {"Records": {
            "eventName": [{"prefix": "ObjectCreated:"}],
        }}

        the following check list would be computed:
        (<built-in function all>, (
            (<built-in function all>, (
                (<built-in function any>, (
                    ('eventName', [{'prefix': 'ObjectCreated:'}], {'eventName': 'ObjectCreated:Put'}),
                    ('eventName', [{'prefix': 'ObjectCreated:'}], {'eventName': 'ObjectCreated:Delete'}))
                ),
            )
        ),))
        """
        rules = []
        for filter_key, filter_value in filter_policy.items():
            if isinstance(filter_value, dict):
                if isinstance(message_body, dict):
                    message_value = message_body.get(filter_key)
                    if message_value is not None:
                        rules.append(
                            self._compute_body_checks(filter_value, message_value)
                        )
                    else:
                        raise FilterPolicyMatcher.CheckException
                elif isinstance(message_body, list):
                    subchecks = []
                    for entry in message_body:
                        subchecks.append(
                            self._compute_body_checks(filter_policy, entry)
                        )
                    rules.append((any, tuple(subchecks)))
                else:
                    raise FilterPolicyMatcher.CheckException

            elif isinstance(filter_value, list):
                # These are the real rules, same as in MessageAttributes case

                concrete_checks = []
                if isinstance(message_body, dict):
                    if message_body is not None:
                        concrete_checks.append((filter_key, filter_value, message_body))
                    else:
                        raise FilterPolicyMatcher.CheckException
                elif isinstance(message_body, list):
                    # Apply policy to each element of the list, pass if at list one element matches
                    for list_elem in message_body:
                        concrete_checks.append((filter_key, filter_value, list_elem))
                else:
                    raise FilterPolicyMatcher.CheckException
                rules.append((any, tuple(concrete_checks)))
            else:
                raise FilterPolicyMatcher.CheckException

        return (all, tuple(rules))

    @staticmethod
    def _field_match(  # type: ignore # decorated function contains type Any
        field: str,
        rules: List[Any],
        dict_body: Dict[str, Any],
        attributes_based_check: bool = True,
    ) -> bool:
        # dict_body = MessageAttributes if attributes_based_check is True
        # otherwise it's the cut-out part of the MessageBody (so only single-level nesting must be supported)

        # Iterate over every rule from the list of rules
        # At least one rule has to match the field for the function to return a match

        def _str_exact_match(value: str, rule: Union[str, List[str]]) -> bool:
            if value == rule:
                return True

            if isinstance(value, list):
                if rule in value:
                    return True

            try:
                json_data = json.loads(value)
                if rule in json_data:
                    return True
            except (ValueError, TypeError):
                pass

            return False

        def _number_match(values: List[float], rule: float) -> bool:
            for value in values:
                # Even the official documentation states a 5 digits of accuracy after the decimal point for numerics, in reality it is 6
                # https://docs.aws.amazon.com/sns/latest/dg/sns-subscription-filter-policies.html#subscription-filter-policy-constraints
                if int(value * 1000000) == int(rule * 1000000):
                    return True

            return False

        def _exists_match(
            should_exist: bool, field: str, dict_body: Dict[str, Any]
        ) -> bool:
            if should_exist and field in dict_body:
                return True
            elif not should_exist and field not in dict_body:
                return True

            return False

        def _prefix_match(prefix: str, value: str) -> bool:
            return value.startswith(prefix)

        def _suffix_match(prefix: str, value: str) -> bool:
            return value.endswith(prefix)

        def _anything_but_match(
            filter_value: Union[Dict[str, Any], List[str], str],
            actual_values: List[str],
        ) -> bool:
            if isinstance(filter_value, dict):
                # We can combine anything-but with the prefix-filter
                anything_but_key = list(filter_value.keys())[0]
                anything_but_val = list(filter_value.values())[0]
                if anything_but_key != "prefix":
                    return False
                if all([not v.startswith(anything_but_val) for v in actual_values]):
                    return True
            else:
                undesired_values = (
                    [filter_value] if isinstance(filter_value, str) else filter_value
                )
                if all([v not in undesired_values for v in actual_values]):
                    return True

            return False

        def _numeric_match(
            numeric_ranges: Iterable[Tuple[str, float]], numeric_value: float
        ) -> bool:
            # numeric_ranges' format:
            # [(< x), (=, y), (>=, z)]
            msg_value = numeric_value
            matches = []
            for operator, test_value in numeric_ranges:
                if operator == ">":
                    matches.append((msg_value > test_value))
                if operator == ">=":
                    matches.append((msg_value >= test_value))
                if operator == "=":
                    matches.append((msg_value == test_value))
                if operator == "<":
                    matches.append((msg_value < test_value))
                if operator == "<=":
                    matches.append((msg_value <= test_value))
            return all(matches)

        for rule in rules:
            #  TODO: boolean value matching is not supported, SNS behavior unknown
            if isinstance(rule, str):
                if attributes_based_check:
                    if field not in dict_body:
                        return False
                    if _str_exact_match(dict_body[field]["Value"], rule):
                        return True
                else:
                    if field not in dict_body:
                        return False
                    if _str_exact_match(dict_body[field], rule):
                        return True

            if isinstance(rule, (int, float)):
                if attributes_based_check:
                    if field not in dict_body:
                        return False
                    if dict_body[field]["Type"] == "Number":
                        attribute_values = [dict_body[field]["Value"]]
                    elif dict_body[field]["Type"] == "String.Array":
                        try:
                            attribute_values = json.loads(dict_body[field]["Value"])
                            if not isinstance(attribute_values, list):
                                attribute_values = [attribute_values]
                        except (ValueError, TypeError):
                            return False
                    else:
                        return False

                    values = [float(value) for value in attribute_values]
                    if _number_match(values, rule):
                        return True
                else:
                    if field not in dict_body:
                        return False

                    if isinstance(dict_body[field], (int, float)):
                        values = [dict_body[field]]
                    elif isinstance(dict_body[field], list):
                        values = [float(value) for value in dict_body[field]]
                    else:
                        return False

                    if _number_match(values, rule):
                        return True

            if isinstance(rule, dict):
                keyword = list(rule.keys())[0]
                value = list(rule.values())[0]
                if keyword == "exists":
                    if attributes_based_check:
                        if _exists_match(value, field, dict_body):
                            return True
                    else:
                        if _exists_match(value, field, dict_body):
                            return True

                elif keyword == "equals-ignore-case" and isinstance(value, str):
                    if attributes_based_check:
                        if field not in dict_body:
                            return False
                        if _str_exact_match(dict_body[field]["Value"].lower(), value):
                            return True
                    else:
                        if field not in dict_body:
                            return False
                        if _str_exact_match(dict_body[field].lower(), value):
                            return True

                elif keyword == "prefix" and isinstance(value, str):
                    if attributes_based_check:
                        if field in dict_body:
                            attr = dict_body[field]
                            if attr["Type"] == "String":
                                if _prefix_match(value, attr["Value"]):
                                    return True
                    else:
                        if field in dict_body:
                            if _prefix_match(value, dict_body[field]):
                                return True
                elif keyword == "suffix" and isinstance(value, str):
                    if attributes_based_check:
                        if field in dict_body:
                            attr = dict_body[field]
                            if attr["Type"] == "String":
                                if _suffix_match(value, attr["Value"]):
                                    return True
                    else:
                        if field in dict_body:
                            if _suffix_match(value, dict_body[field]):
                                return True

                elif keyword == "anything-but":
                    if attributes_based_check:
                        if field not in dict_body:
                            return False
                        attr = dict_body[field]
                        if isinstance(value, dict):
                            # We can combine anything-but with the prefix-filter
                            if attr["Type"] == "String":
                                actual_values = [attr["Value"]]
                            else:
                                actual_values = [v for v in attr["Value"]]
                        else:
                            if attr["Type"] == "Number":
                                actual_values = [float(attr["Value"])]
                            elif attr["Type"] == "String":
                                actual_values = [attr["Value"]]
                            else:
                                actual_values = [v for v in attr["Value"]]

                        if _anything_but_match(value, actual_values):
                            return True
                    else:
                        if field not in dict_body:
                            return False
                        attr = dict_body[field]
                        if isinstance(value, dict):
                            if isinstance(attr, str):
                                actual_values = [attr]
                            elif isinstance(attr, list):
                                actual_values = attr
                            else:
                                return False
                        else:
                            if isinstance(attr, (int, float, str)):
                                actual_values = [attr]
                            elif isinstance(attr, list):
                                actual_values = attr
                            else:
                                return False

                        if _anything_but_match(value, actual_values):
                            return True

                elif keyword == "numeric" and isinstance(value, list):
                    if attributes_based_check:
                        if dict_body.get(field, {}).get("Type", "") == "Number":
                            checks = value
                            numeric_ranges = zip(checks[0::2], checks[1::2])
                            if _numeric_match(
                                numeric_ranges, float(dict_body[field]["Value"])
                            ):
                                return True
                    else:
                        if field not in dict_body:
                            return False

                        checks = value
                        numeric_ranges = zip(checks[0::2], checks[1::2])
                        if _numeric_match(numeric_ranges, dict_body[field]):
                            return True

            if field == "$or" and isinstance(rules, list):
                return any(
                    [
                        FilterPolicyMatcher._attributes_based_match(dict_body, rule)
                        for rule in rules
                    ]
                )

        return False
