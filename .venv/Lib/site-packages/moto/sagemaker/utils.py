import json
import typing
from datetime import datetime
from typing import Any, DefaultDict, Dict, List, Optional

from dateutil.tz import tzutc

from moto.s3.models import s3_backends
from moto.utilities.utils import get_partition

from .exceptions import ValidationError

if typing.TYPE_CHECKING:
    from .models import FakeModelCard, FakePipeline, FakePipelineExecution


def get_pipeline_from_name(
    pipelines: Dict[str, "FakePipeline"], pipeline_name: str
) -> "FakePipeline":
    try:
        return pipelines[pipeline_name]
    except KeyError:
        raise ValidationError(
            message=f"Could not find pipeline with PipelineName {pipeline_name}."
        )


def get_pipeline_name_from_execution_arn(pipeline_execution_arn: str) -> str:
    return pipeline_execution_arn.split("/")[1].split(":")[-1]


def get_pipeline_execution_from_arn(
    pipelines: Dict[str, "FakePipeline"], pipeline_execution_arn: str
) -> "FakePipelineExecution":
    try:
        pipeline_name = get_pipeline_name_from_execution_arn(pipeline_execution_arn)
        pipeline = get_pipeline_from_name(pipelines, pipeline_name)
        return pipeline.pipeline_executions[pipeline_execution_arn]
    except KeyError:
        raise ValidationError(
            message=f"Could not find pipeline execution with PipelineExecutionArn {pipeline_execution_arn}."
        )


def load_pipeline_definition_from_s3(
    pipeline_definition_s3_location: Dict[str, Any], account_id: str, partition: str
) -> Dict[str, Any]:
    s3_backend = s3_backends[account_id][partition]
    result = s3_backend.get_object(
        bucket_name=pipeline_definition_s3_location["Bucket"],
        key_name=pipeline_definition_s3_location["ObjectKey"],
    )
    return json.loads(result.value)  # type: ignore[union-attr]


def arn_formatter(_type: str, _id: str, account_id: str, region_name: str) -> str:
    return f"arn:{get_partition(region_name)}:sagemaker:{region_name}:{account_id}:{_type}/{_id}"


def validate_model_approval_status(model_approval_status: typing.Optional[str]) -> None:
    if model_approval_status is not None and model_approval_status not in [
        "Approved",
        "Rejected",
        "PendingManualApproval",
    ]:
        raise ValidationError(
            f"Value '{model_approval_status}' at 'modelApprovalStatus' failed to satisfy constraint: "
            "Member must satisfy enum value set: [PendingManualApproval, Approved, Rejected]"
        )


def filter_model_cards(
    model_cards: DefaultDict[str, List["FakeModelCard"]],
    creation_time_after: Optional[datetime],
    creation_time_before: Optional[datetime],
    name_contains: Optional[str],
    model_card_status: Optional[str],
    sort_by: Optional[str],
    sort_order: Optional[str],
) -> List["FakeModelCard"]:
    reverse = sort_order == "Descending"

    if name_contains:
        lowercase_name = name_contains.lower()
        filtered_cards = {
            k: v for k, v in model_cards.items() if lowercase_name in k.lower()
        }
    else:
        filtered_cards = {k: v for k, v in model_cards.items()}

    result: List[FakeModelCard] = []
    for _, versions in filtered_cards.items():
        filtered_versions = versions

        if creation_time_after:
            if isinstance(creation_time_after, int):
                creation_time_after = datetime.fromtimestamp(
                    creation_time_after, tz=tzutc()
                )
            filtered_versions = [
                v
                for v in filtered_versions
                if v.last_modified_time > str(creation_time_after)
            ]

        if creation_time_before:
            if isinstance(creation_time_before, int):
                creation_time_before = datetime.fromtimestamp(
                    creation_time_before, tz=tzutc()
                )
            filtered_versions = [
                v
                for v in filtered_versions
                if v.last_modified_time < str(creation_time_before)
            ]

        if model_card_status:
            filtered_versions = [
                v for v in filtered_versions if v.model_card_status == model_card_status
            ]

        if filtered_versions:
            latest_version = max(filtered_versions, key=lambda v: v.last_modified_time)
            result.append(latest_version)

    if not result:
        return []

    def sort_key(x: "FakeModelCard") -> str:
        if sort_by == "Name":
            return x.model_card_name
        return x.creation_time

    return sorted(result, key=sort_key, reverse=reverse)
