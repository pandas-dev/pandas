from typing import TYPE_CHECKING, List, TypedDict

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
