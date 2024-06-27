import json
from typing import Any, Dict, List, Tuple

from moto.codepipeline.exceptions import (
    InvalidStructureException,
    InvalidTagsException,
    PipelineNotFoundException,
    ResourceNotFoundException,
    TooManyTagsException,
)
from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.common_models import BaseModel
from moto.core.utils import iso_8601_datetime_with_milliseconds, utcnow
from moto.iam.exceptions import IAMNotFoundException
from moto.iam.models import IAMBackend, iam_backends
from moto.utilities.utils import get_partition


class CodePipeline(BaseModel):
    def __init__(self, account_id: str, region: str, pipeline: Dict[str, Any]):
        # the version number for a new pipeline is always 1
        pipeline["version"] = 1

        self.pipeline = self.add_default_values(pipeline)
        self.tags: Dict[str, str] = {}

        self._arn = f"arn:{get_partition(region)}:codepipeline:{region}:{account_id}:{pipeline['name']}"
        self._created = utcnow()
        self._updated = utcnow()

    @property
    def metadata(self) -> Dict[str, str]:
        return {
            "pipelineArn": self._arn,
            "created": iso_8601_datetime_with_milliseconds(self._created),
            "updated": iso_8601_datetime_with_milliseconds(self._updated),
        }

    def add_default_values(self, pipeline: Dict[str, Any]) -> Dict[str, Any]:
        for stage in pipeline["stages"]:
            for action in stage["actions"]:
                if "runOrder" not in action:
                    action["runOrder"] = 1
                if "configuration" not in action:
                    action["configuration"] = {}
                if "outputArtifacts" not in action:
                    action["outputArtifacts"] = []
                if "inputArtifacts" not in action:
                    action["inputArtifacts"] = []

        return pipeline

    def validate_tags(self, tags: List[Dict[str, str]]) -> None:
        for tag in tags:
            if tag["key"].startswith("aws:"):
                raise InvalidTagsException(
                    "Not allowed to modify system tags. "
                    "System tags start with 'aws:'. "
                    "msg=[Caller is an end user and not allowed to mutate system tags]"
                )

        if (len(self.tags) + len(tags)) > 50:
            raise TooManyTagsException(self._arn)


class CodePipelineBackend(BaseBackend):
    def __init__(self, region_name: str, account_id: str):
        super().__init__(region_name, account_id)
        self.pipelines: Dict[str, CodePipeline] = {}

    @staticmethod
    def default_vpc_endpoint_service(
        service_region: str, zones: List[str]
    ) -> List[Dict[str, str]]:
        """Default VPC endpoint service."""
        return BaseBackend.default_vpc_endpoint_service_factory(
            service_region, zones, "codepipeline", policy_supported=False
        )

    @property
    def iam_backend(self) -> IAMBackend:
        return iam_backends[self.account_id][self.partition]

    def create_pipeline(
        self, pipeline: Dict[str, Any], tags: List[Dict[str, str]]
    ) -> Tuple[Dict[str, Any], List[Dict[str, str]]]:
        name = pipeline["name"]
        if name in self.pipelines:
            raise InvalidStructureException(
                f"A pipeline with the name '{name}' already exists in account '{self.account_id}'"
            )

        try:
            role = self.iam_backend.get_role_by_arn(pipeline["roleArn"])
            trust_policy_statements = json.loads(role.assume_role_policy_document)[
                "Statement"
            ]
            trusted_service_principals = [
                i["Principal"]["Service"] for i in trust_policy_statements
            ]
            if "codepipeline.amazonaws.com" not in trusted_service_principals:
                raise IAMNotFoundException("")
        except IAMNotFoundException:
            raise InvalidStructureException(
                f"CodePipeline is not authorized to perform AssumeRole on role {pipeline['roleArn']}"
            )

        if len(pipeline["stages"]) < 2:
            raise InvalidStructureException(
                "Pipeline has only 1 stage(s). There should be a minimum of 2 stages in a pipeline"
            )

        self.pipelines[pipeline["name"]] = CodePipeline(
            self.account_id, self.region_name, pipeline
        )

        if tags is not None:
            self.pipelines[pipeline["name"]].validate_tags(tags)

            new_tags = {tag["key"]: tag["value"] for tag in tags}
            self.pipelines[pipeline["name"]].tags.update(new_tags)
        else:
            tags = []

        return pipeline, sorted(tags, key=lambda i: i["key"])

    def get_pipeline(self, name: str) -> Tuple[Dict[str, Any], Dict[str, str]]:
        codepipeline = self.pipelines.get(name)

        if not codepipeline:
            raise PipelineNotFoundException(
                f"Account '{self.account_id}' does not have a pipeline with name '{name}'"
            )

        return codepipeline.pipeline, codepipeline.metadata

    def update_pipeline(self, pipeline: Dict[str, Any]) -> Dict[str, Any]:
        codepipeline = self.pipelines.get(pipeline["name"])

        if not codepipeline:
            raise ResourceNotFoundException(
                f"The account with id '{self.account_id}' does not include a pipeline with the name '{pipeline['name']}'"
            )

        # version number is auto incremented
        pipeline["version"] = codepipeline.pipeline["version"] + 1
        codepipeline._updated = utcnow()
        codepipeline.pipeline = codepipeline.add_default_values(pipeline)

        return codepipeline.pipeline

    def list_pipelines(self) -> List[Dict[str, str]]:
        pipelines = []

        for name, codepipeline in self.pipelines.items():
            pipelines.append(
                {
                    "name": name,
                    "version": codepipeline.pipeline["version"],
                    "created": codepipeline.metadata["created"],
                    "updated": codepipeline.metadata["updated"],
                }
            )

        return sorted(pipelines, key=lambda i: i["name"])

    def delete_pipeline(self, name: str) -> None:
        self.pipelines.pop(name, None)

    def list_tags_for_resource(self, arn: str) -> List[Dict[str, str]]:
        name = arn.split(":")[-1]
        pipeline = self.pipelines.get(name)

        if not pipeline:
            raise ResourceNotFoundException(
                f"The account with id '{self.account_id}' does not include a pipeline with the name '{name}'"
            )

        tags = [{"key": key, "value": value} for key, value in pipeline.tags.items()]

        return sorted(tags, key=lambda i: i["key"])

    def tag_resource(self, arn: str, tags: List[Dict[str, str]]) -> None:
        name = arn.split(":")[-1]
        pipeline = self.pipelines.get(name)

        if not pipeline:
            raise ResourceNotFoundException(
                f"The account with id '{self.account_id}' does not include a pipeline with the name '{name}'"
            )

        pipeline.validate_tags(tags)

        for tag in tags:
            pipeline.tags.update({tag["key"]: tag["value"]})

    def untag_resource(self, arn: str, tag_keys: List[str]) -> None:
        name = arn.split(":")[-1]
        pipeline = self.pipelines.get(name)

        if not pipeline:
            raise ResourceNotFoundException(
                f"The account with id '{self.account_id}' does not include a pipeline with the name '{name}'"
            )

        for key in tag_keys:
            pipeline.tags.pop(key, None)


codepipeline_backends = BackendDict(CodePipelineBackend, "codepipeline")
