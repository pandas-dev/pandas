import json
from typing import Any, Dict, List

from moto.ecr.exceptions import InvalidParameterException

REQUIRED_RULE_PROPERTIES = {"rulePriority", "selection", "action"}
VALID_RULE_PROPERTIES = {"description", *REQUIRED_RULE_PROPERTIES}

REQUIRED_ACTION_PROPERTIES = {"type"}
VALID_ACTION_PROPERTIES = REQUIRED_ACTION_PROPERTIES

VALID_ACTION_TYPE_VALUES = {"expire"}

REQUIRED_SELECTION_PROPERTIES = {"tagStatus", "countType", "countNumber"}
VALID_SELECTION_PROPERTIES = {
    "tagPrefixList",
    "countUnit",
    *REQUIRED_SELECTION_PROPERTIES,
}

VALID_SELECTION_TAG_STATUS_VALUES = {"tagged", "untagged", "any"}
VALID_SELECTION_COUNT_TYPE_VALUES = {"imageCountMoreThan", "sinceImagePushed"}
VALID_SELECTION_COUNT_UNIT_VALUES = {"days"}


class EcrLifecyclePolicyValidator:
    INVALID_PARAMETER_ERROR_MESSAGE = (
        "Invalid parameter at 'LifecyclePolicyText' failed to satisfy constraint: "
        "'Lifecycle policy validation failure: "
    )

    def __init__(self, policy_text: str):
        self._policy_text = policy_text
        self._policy_json: Dict[str, Any] = {}
        self._rules: List[Any] = []

    def validate(self) -> None:
        try:
            self._parse_policy()
        except Exception:
            raise InvalidParameterException(
                "".join(
                    [
                        self.INVALID_PARAMETER_ERROR_MESSAGE,
                        "Could not map policyString into LifecyclePolicy.'",
                    ]
                )
            )

        try:
            self._extract_rules()
        except Exception:
            raise InvalidParameterException(
                "".join(
                    [
                        self.INVALID_PARAMETER_ERROR_MESSAGE,
                        'object has missing required properties (["rules"])\'',
                    ]
                )
            )

        self._validate_rule_type()
        self._validate_rule_top_properties()

    def _parse_policy(self) -> None:
        self._policy_json = json.loads(self._policy_text)
        assert isinstance(self._policy_json, dict)

    def _extract_rules(self) -> None:
        assert "rules" in self._policy_json
        assert isinstance(self._policy_json["rules"], list)

        self._rules = self._policy_json["rules"]

    def _validate_rule_type(self) -> None:
        for rule in self._rules:
            if not isinstance(rule, dict):
                raise InvalidParameterException(
                    "".join(
                        [
                            self.INVALID_PARAMETER_ERROR_MESSAGE,
                            f'instance type ({type(rule)}) does not match any allowed primitive type (allowed: ["object"])\'',
                        ]
                    )
                )

    def _validate_rule_top_properties(self) -> None:
        for rule in self._rules:
            rule_properties = set(rule.keys())
            missing_properties = REQUIRED_RULE_PROPERTIES - rule_properties
            if missing_properties:
                raise InvalidParameterException(
                    "".join(
                        [
                            self.INVALID_PARAMETER_ERROR_MESSAGE,
                            f"object has missing required properties ({json.dumps(sorted(missing_properties))})'",
                        ]
                    )
                )

            for rule_property in rule_properties:
                if rule_property not in VALID_RULE_PROPERTIES:
                    raise InvalidParameterException(
                        "".join(
                            [
                                self.INVALID_PARAMETER_ERROR_MESSAGE,
                                f'object instance has properties which are not allowed by the schema: (["{rule_property}"])\'',
                            ]
                        )
                    )

            self._validate_action(rule["action"])
            self._validate_selection(rule["selection"])

    def _validate_action(self, action: Any) -> None:
        given_properties = set(action.keys())
        missing_properties = REQUIRED_ACTION_PROPERTIES - given_properties

        if missing_properties:
            raise InvalidParameterException(
                "".join(
                    [
                        self.INVALID_PARAMETER_ERROR_MESSAGE,
                        f"object has missing required properties ({json.dumps(sorted(missing_properties))})'",
                    ]
                )
            )

        for given_property in given_properties:
            if given_property not in VALID_ACTION_PROPERTIES:
                raise InvalidParameterException(
                    "".join(
                        [
                            self.INVALID_PARAMETER_ERROR_MESSAGE,
                            "object instance has properties "
                            f'which are not allowed by the schema: (["{given_property}"])\'',
                        ]
                    )
                )

            self._validate_action_type(action["type"])

    def _validate_action_type(self, action_type: str) -> None:
        if action_type not in VALID_ACTION_TYPE_VALUES:
            raise InvalidParameterException(
                "".join(
                    [
                        self.INVALID_PARAMETER_ERROR_MESSAGE,
                        f"instance value ({action_type}) not found in enum "
                        f":(possible values: {json.dumps(sorted(VALID_ACTION_TYPE_VALUES))})'",
                    ]
                )
            )

    def _validate_selection(self, selection: Any) -> None:
        given_properties = set(selection.keys())
        missing_properties = REQUIRED_SELECTION_PROPERTIES - given_properties

        if missing_properties:
            raise InvalidParameterException(
                "".join(
                    [
                        self.INVALID_PARAMETER_ERROR_MESSAGE,
                        f"object has missing required properties ({json.dumps(sorted(missing_properties))})'",
                    ]
                )
            )

        for given_property in given_properties:
            if given_property not in VALID_SELECTION_PROPERTIES:
                raise InvalidParameterException(
                    "".join(
                        [
                            self.INVALID_PARAMETER_ERROR_MESSAGE,
                            "object instance has properties "
                            f'which are not allowed by the schema: (["{given_property}"])\'',
                        ]
                    )
                )

            self._validate_selection_tag_status(selection["tagStatus"])
            self._validate_selection_count_type(selection["countType"])
            self._validate_selection_count_unit(selection.get("countUnit"))
            self._validate_selection_count_number(selection["countNumber"])

    def _validate_selection_tag_status(self, tag_status: Any) -> None:
        if tag_status not in VALID_SELECTION_TAG_STATUS_VALUES:
            raise InvalidParameterException(
                "".join(
                    [
                        self.INVALID_PARAMETER_ERROR_MESSAGE,
                        f"instance value ({tag_status}) not found in enum "
                        f":(possible values: {json.dumps(sorted(VALID_SELECTION_TAG_STATUS_VALUES))})'",
                    ]
                )
            )

    def _validate_selection_count_type(self, count_type: Any) -> None:
        if count_type not in VALID_SELECTION_COUNT_TYPE_VALUES:
            raise InvalidParameterException(
                "".join(
                    [
                        self.INVALID_PARAMETER_ERROR_MESSAGE,
                        "instance failed to match exactly one schema (matched 0 out of 2)",
                    ]
                )
            )

    def _validate_selection_count_unit(self, count_unit: Any) -> None:
        if not count_unit:
            return None

        if count_unit not in VALID_SELECTION_COUNT_UNIT_VALUES:
            raise InvalidParameterException(
                "".join(
                    [
                        self.INVALID_PARAMETER_ERROR_MESSAGE,
                        f"instance value ({count_unit}) not found in enum "
                        f":(possible values: {json.dumps(sorted(VALID_SELECTION_COUNT_UNIT_VALUES))})'",
                    ]
                )
            )

    def _validate_selection_count_number(self, count_number: int) -> None:
        if count_number < 1:
            raise InvalidParameterException(
                "".join(
                    [
                        self.INVALID_PARAMETER_ERROR_MESSAGE,
                        "numeric instance is lower than the required minimum "
                        f"(minimum: 1, found: {count_number})",
                    ]
                )
            )
