from typing import Any

PAGINATION_MODEL = {
    "list_job_executions_for_job": {
        "input_token": "token",
        "limit_key": "max_results",
        "limit_default": 100,
        "unique_attribute": ["job_id", "thing_arn"],
    },
    "list_job_executions_for_thing": {
        "input_token": "next_token",
        "limit_key": "max_results",
        "limit_default": 100,
        "unique_attribute": "job_id",
    },
    "list_job_templates": {
        "input_token": "next_token",
        "limit_key": "max_results",
        "limit_default": 100,
        "unique_attribute": "jobTemplateId",
    },
    "list_jobs": {
        "input_token": "next_token",
        "limit_key": "max_results",
        "limit_default": 100,
        "unique_attribute": "jobId",
    },
    "list_things": {
        "input_token": "token",
        "limit_key": "max_results",
        "limit_default": 100,
        "unique_attribute": "thingArn",
    },
    "list_billing_groups": {
        "input_token": "token",
        "limit_key": "max_results",
        "limit_default": 100,
        "unique_attribute": "groupName",
    },
    "list_things_in_billing_group": {
        "input_token": "token",
        "limit_key": "max_results",
        "limit_default": 100,
        "unique_attribute": "thing_name",
    },
}


def decapitalize_str(obj: str) -> str:
    return obj[0].lower() + obj[1:]


def decapitalize_dict(obj: Any) -> Any:
    if isinstance(obj, dict):
        return {
            decapitalize_str(key): decapitalize_dict(value)
            for key, value in obj.items()
        }
    else:
        return obj
