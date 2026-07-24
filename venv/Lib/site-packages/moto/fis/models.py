"""FISBackend class with methods for supported APIs."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Optional

from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.common_models import BaseModel
from moto.core.utils import unix_time
from moto.moto_api._internal import mock_random
from moto.utilities.tagging_service import TaggingService
from moto.utilities.utils import get_partition

from .exceptions import ResourceNotFoundException


@dataclass
class ExperimentTemplate(BaseModel):
    account_id: str
    region_name: str
    id: str
    client_token: str
    description: str
    stop_conditions: list[dict[str, Any]]
    targets: dict[str, Any]
    actions: dict[str, Any]
    role_arn: str
    log_configuration: Optional[dict[str, Any]]
    experiment_options: Optional[dict[str, Any]]
    experiment_report_configuration: Optional[dict[str, Any]]
    creation_time: float
    last_update_time: float

    @property
    def arn(self) -> str:
        return f"arn:{get_partition(self.region_name)}:fis:{self.region_name}:{self.account_id}:experiment-template/{self.id}"

    def to_dict(self, tags: Optional[dict[str, str]] = None) -> dict[str, Any]:
        dct: dict[str, Any] = {
            "id": self.id,
            "arn": self.arn,
            "description": self.description,
            "targets": self.targets or {},
            "actions": self.actions or {},
            "stopConditions": self.stop_conditions or [],
            "creationTime": self.creation_time,
            "lastUpdateTime": self.last_update_time,
            "roleArn": self.role_arn,
            "tags": tags or {},
            "logConfiguration": self.log_configuration,
            "experimentOptions": self.experiment_options,
            "targetAccountConfigurationsCount": 0,
            "experimentReportConfiguration": self.experiment_report_configuration,
        }
        return {k: v for k, v in dct.items() if v is not None}


class FISBackend(BaseBackend):
    def __init__(self, region_name: str, account_id: str):
        super().__init__(region_name, account_id)
        self.experiment_templates: dict[str, ExperimentTemplate] = {}
        self._client_token_to_template_id: dict[str, str] = {}
        self.tagger = TaggingService()

    def create_experiment_template(
        self,
        client_token: Optional[str],
        description: str,
        stop_conditions: list[dict[str, Any]],
        targets: Optional[dict[str, Any]],
        actions: dict[str, Any],
        role_arn: str,
        tags: Optional[dict[str, str]] = None,
        log_configuration: Optional[dict[str, Any]] = None,
        experiment_options: Optional[dict[str, Any]] = None,
        experiment_report_configuration: Optional[dict[str, Any]] = None,
    ) -> dict[str, Any]:
        # clientToken is required per the AWS API; auto-populate to match boto3's behaviour.
        token = client_token or str(mock_random.uuid4())
        existing_id = self._client_token_to_template_id.get(token)
        if existing_id and existing_id in self.experiment_templates:
            tmpl = self.experiment_templates[existing_id]
            return tmpl.to_dict(tags=self.tagger.get_tag_dict_for_resource(tmpl.arn))

        template_id = mock_random.uuid4().hex
        now = unix_time()
        template = ExperimentTemplate(
            account_id=self.account_id,
            region_name=self.region_name,
            id=template_id,
            client_token=token,
            description=description,
            stop_conditions=stop_conditions,
            targets=targets or {},
            actions=actions,
            role_arn=role_arn,
            log_configuration=log_configuration,
            experiment_options=experiment_options,
            experiment_report_configuration=experiment_report_configuration,
            creation_time=now,
            last_update_time=now,
        )
        self.experiment_templates[template_id] = template
        self._client_token_to_template_id[token] = template_id
        if tags:
            self.tagger.tag_resource(
                template.arn, TaggingService.convert_dict_to_tags_input(tags)
            )
        return template.to_dict(tags=tags or {})

    def delete_experiment_template(self, id: str) -> dict[str, Any]:
        if id not in self.experiment_templates:
            raise ResourceNotFoundException(f"Experiment template {id} does not exist")
        template = self.experiment_templates.pop(id)
        tags = self.tagger.get_tag_dict_for_resource(template.arn)
        self._client_token_to_template_id.pop(template.client_token, None)
        return template.to_dict(tags=tags)

    def tag_resource(self, resource_arn: str, tags: dict[str, str]) -> None:
        self.tagger.tag_resource(
            resource_arn, TaggingService.convert_dict_to_tags_input(tags or {})
        )
        for tmpl in self.experiment_templates.values():
            if tmpl.arn == resource_arn:
                tmpl.last_update_time = unix_time()
                break

    def untag_resource(self, resource_arn: str, tag_keys: list[str]) -> None:
        self.tagger.untag_resource_using_names(resource_arn, tag_keys)
        for tmpl in self.experiment_templates.values():
            if tmpl.arn == resource_arn:
                tmpl.last_update_time = unix_time()
                break

    def list_tags_for_resource(self, resource_arn: str) -> dict[str, str]:
        return self.tagger.get_tag_dict_for_resource(resource_arn)


fis_backends = BackendDict(
    FISBackend,
    "fis",
    use_boto3_regions=False,
    additional_regions=[
        "us-east-1",
        "us-east-2",
        "us-west-1",
        "us-west-2",
        "af-south-1",
        "ap-east-1",
        "ap-south-1",
        "ap-northeast-1",
        "ap-northeast-2",
        "ap-northeast-3",
        "ap-southeast-1",
        "ap-southeast-2",
        "ca-central-1",
        "eu-central-1",
        "eu-west-1",
        "eu-west-2",
        "eu-west-3",
        "eu-north-1",
        "eu-south-1",
        "me-south-1",
        "sa-east-1",
    ],
)
