import json
from typing import TYPE_CHECKING, Any, Dict, List, TypedDict

if TYPE_CHECKING:
    from typing_extensions import Any, Dict, Required, Union


# NOTE: Typing is based on the following document https://docs.aws.amazon.com/eventbridge/latest/userguide/eb-event-patterns.html
EventMessageType = TypedDict(
    "EventMessageType",
    {
        "version": str,
        "id": str,
        "detail-type": "Required[Union[str, List[str]]]",
        "source": "Required[Union[str, List[str]]]",
        "account": str,
        # support float type for internal use of moto.
        "time": "Union[str, float]",
        "replay-name": str,
        "region": str,
        "resources": List[str],
        "detail": "Required[Dict[str, Any]]",
        "ingestion-time": float,
    },
    total=False,
)

PAGINATION_MODEL = {
    "list_rules": {
        "input_token": "next_token",
        "limit_key": "limit",
        "limit_default": 50,
        "unique_attribute": "arn",
        "fail_on_invalid_token": False,
    },
    "list_rule_names_by_target": {
        "input_token": "next_token",
        "limit_key": "limit",
        "limit_default": 50,
        "unique_attribute": "arn",
        "fail_on_invalid_token": False,
    },
    "list_targets_by_rule": {
        "input_token": "next_token",
        "limit_key": "limit",
        "limit_default": 50,
        "unique_attribute": "Arn",
        "fail_on_invalid_token": False,
    },
}

_BASE_EVENT_MESSAGE: EventMessageType = {
    "version": "0",
    "id": "17793124-05d4-b198-2fde-7ededc63b103",
    "detail-type": "",
    "source": "",
    "account": "",
    "time": "",
    "region": "",
    "resources": [],
    "detail": {},
}


class EventTemplateParser:
    DEFAULT_EVENT_INPUT_TEMPLATE = '{"id": "<id>", "time": <time>, "version": "0", "detail-type": "<detail-type>", "source": "<source>", "region": "<region>", "resources": <resources>, "detail": <detail>}'
    DEFAULT_EVENT_INPUT_PATHS_MAP = {
        "account": "$.account",
        "detail": "$.detail",
        "detail-type": "$.detail-type",
        "source": "$.source",
        "id": "$.id",
        "region": "$.region",
        "resources": "$.resources",
        "time": "$.time",
    }

    @staticmethod
    def _stringify(result: Any) -> str:  # type: ignore[misc]
        if isinstance(result, dict):
            result = json.dumps(result)
        elif isinstance(result, list):
            result = json.dumps([EventTemplateParser._stringify(x) for x in result])
        elif isinstance(result, (int, float)):
            result = str(result)
        return result

    @staticmethod
    def parse(  # type: ignore[misc]
        input_template: str, input_paths_map: Dict[str, Any], event: EventMessageType
    ) -> Dict[str, Any]:
        from jsonpath_ng.ext import parse

        template_to_use = (
            input_template or EventTemplateParser.DEFAULT_EVENT_INPUT_TEMPLATE
        )
        input_paths_map = (
            input_paths_map or EventTemplateParser.DEFAULT_EVENT_INPUT_PATHS_MAP
        )
        for input_path in input_paths_map:
            input_expr = parse(input_paths_map[input_path])
            matches = input_expr.find(event)
            result = (
                EventTemplateParser._stringify(matches[0].value) if matches else None
            )
            if result:
                template_to_use = template_to_use.replace(f"<{input_path}>", result)

        default_inputs_map = {
            "aws.events.event.json": event,
            "aws.events.event.ingestion-time": event["ingestion-time"],
        }
        for input_path in default_inputs_map:
            result = EventTemplateParser._stringify(default_inputs_map[input_path])
            template_to_use = template_to_use.replace(f"<{input_path}>", result)

        return json.loads(template_to_use)
