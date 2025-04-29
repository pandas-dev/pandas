from typing import Dict, List

PAGINATION_MODEL = {
    "list_executions": {
        "input_token": "next_token",
        "limit_key": "max_results",
        "limit_default": 100,
        "unique_attribute": ["start_date", "execution_arn"],
    },
    "list_state_machines": {
        "input_token": "next_token",
        "limit_key": "max_results",
        "limit_default": 100,
        "unique_attribute": ["creation_date", "arn"],
    },
}


def cfn_to_api_tags(cfn_tags_entry: List[Dict[str, str]]) -> List[Dict[str, str]]:
    return [{k.lower(): v for k, v in d.items()} for d in cfn_tags_entry]


def api_to_cfn_tags(api_tags: List[Dict[str, str]]) -> List[Dict[str, str]]:
    return [{k.capitalize(): v for k, v in d.items()} for d in api_tags]
