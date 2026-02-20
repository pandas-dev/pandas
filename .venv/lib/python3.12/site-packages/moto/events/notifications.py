import json
from typing import Any

from moto.utilities.utils import get_partition

from .utils import EventMessageType

_EVENT_S3_OBJECT_CREATED: EventMessageType = {
    "version": "0",
    "id": "17793124-05d4-b198-2fde-7ededc63b103",
    "detail-type": "Object Created",
    "source": "aws.s3",
    "account": "123456789012",
    "time": "2021-11-12T00:00:00Z",
    "region": "us-west-2",
    "resources": [],
    "detail": {},
}


def send_notification(
    source: str, event_name: str, region: str, resources: Any, detail: Any
) -> None:
    try:
        _send_safe_notification(source, event_name, region, resources, detail)
    except:  # noqa
        # If anything goes wrong, we should never fail
        pass


def _send_safe_notification(
    source: str, event_name: str, region: str, resources: Any, detail: Any
) -> None:
    from .models import events_backends

    event = None
    if source == "aws.s3" and event_name == "CreateBucket":
        event = _EVENT_S3_OBJECT_CREATED.copy()
        event["region"] = region
        event["resources"] = resources
        event["detail"] = detail

    if event is None:
        return

    partition = get_partition(region)

    for account_id, account in events_backends.items():
        for backend in account.values():
            applicable_targets = []
            for event_bus in backend.event_buses.values():
                for rule in event_bus.rules.values():
                    if rule.state != "ENABLED":
                        continue
                    pattern = rule.event_pattern.get_pattern()
                    if source in pattern.get("source", []):
                        if event_name in pattern.get("detail", {}).get("eventName", []):
                            applicable_targets.extend(rule.targets)

            for target in applicable_targets:
                target_arn = target.get("Arn", "")
                if target_arn.startswith(f"arn:{partition}:lambda"):
                    _invoke_lambda(account_id, fn_arn=target_arn, event=event)  # type: ignore[arg-type]
                if target_arn.startswith(f"arn:{partition}:states"):
                    _start_sfn_execution(
                        account_id, region=region, sfn_arn=target_arn, event=event
                    )


def _invoke_lambda(account_id: str, fn_arn: str, event: Any) -> None:
    from moto.awslambda.utils import get_backend

    lambda_region = fn_arn.split(":")[3]

    body = json.dumps(event)
    get_backend(account_id, lambda_region).invoke(
        function_name=fn_arn,
        qualifier=None,
        body=body,
        headers={},
        response_headers={},
    )


def _start_sfn_execution(
    account_id: str, region: str, sfn_arn: str, event: Any
) -> None:
    from moto.stepfunctions.models import StepFunctionBackend
    from moto.stepfunctions.responses import get_backend

    backend: StepFunctionBackend = get_backend(account_id, region)

    backend.start_execution(
        state_machine_arn=sfn_arn,
        name="",
        execution_input=json.dumps(event),
    )
