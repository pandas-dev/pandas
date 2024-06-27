import datetime
from collections import OrderedDict
from typing import Any, Dict, Iterable, List

from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.common_models import BaseModel, CloudFormationModel
from moto.core.utils import utcnow

from .utils import get_random_pipeline_id, remove_capitalization_of_dict_keys


class PipelineObject(BaseModel):
    def __init__(self, object_id: str, name: str, fields: Any):
        self.object_id = object_id
        self.name = name
        self.fields = fields

    def to_json(self) -> Dict[str, Any]:
        return {"fields": self.fields, "id": self.object_id, "name": self.name}


class Pipeline(CloudFormationModel):
    def __init__(self, name: str, unique_id: str, **kwargs: Any):
        self.name = name
        self.unique_id = unique_id
        self.description = kwargs.get("description", "")
        self.pipeline_id = get_random_pipeline_id()
        self.creation_time = utcnow()
        self.objects: List[Any] = []
        self.status = "PENDING"
        self.tags = kwargs.get("tags", [])

    @property
    def physical_resource_id(self) -> str:
        return self.pipeline_id

    def to_meta_json(self) -> Dict[str, str]:
        return {"id": self.pipeline_id, "name": self.name}

    def to_json(self) -> Dict[str, Any]:
        return {
            "description": self.description,
            "fields": [
                {"key": "@pipelineState", "stringValue": self.status},
                {"key": "description", "stringValue": self.description},
                {"key": "name", "stringValue": self.name},
                {
                    "key": "@creationTime",
                    "stringValue": datetime.datetime.strftime(
                        self.creation_time, "%Y-%m-%dT%H-%M-%S"
                    ),
                },
                {"key": "@id", "stringValue": self.pipeline_id},
                {"key": "@sphere", "stringValue": "PIPELINE"},
                {"key": "@version", "stringValue": "1"},
                {"key": "@userId", "stringValue": "924374875933"},
                {"key": "@accountId", "stringValue": "924374875933"},
                {"key": "uniqueId", "stringValue": self.unique_id},
            ],
            "name": self.name,
            "pipelineId": self.pipeline_id,
            "tags": self.tags,
        }

    def set_pipeline_objects(self, pipeline_objects: Any) -> None:
        self.objects = [
            PipelineObject(
                pipeline_object["id"],
                pipeline_object["name"],
                pipeline_object["fields"],
            )
            for pipeline_object in remove_capitalization_of_dict_keys(pipeline_objects)
        ]

    def activate(self) -> None:
        self.status = "SCHEDULED"

    @staticmethod
    def cloudformation_name_type() -> str:
        return "Name"

    @staticmethod
    def cloudformation_type() -> str:
        # https://docs.aws.amazon.com/AWSCloudFormation/latest/UserGuide/aws-resource-datapipeline-pipeline.html
        return "AWS::DataPipeline::Pipeline"

    @classmethod
    def create_from_cloudformation_json(  # type: ignore[misc]
        cls,
        resource_name: str,
        cloudformation_json: Dict[str, Any],
        account_id: str,
        region_name: str,
        **kwargs: Any,
    ) -> "Pipeline":
        datapipeline_backend = datapipeline_backends[account_id][region_name]
        properties = cloudformation_json["Properties"]

        cloudformation_unique_id = "cf-" + resource_name
        pipeline = datapipeline_backend.create_pipeline(
            resource_name, cloudformation_unique_id
        )
        datapipeline_backend.put_pipeline_definition(
            pipeline.pipeline_id, properties["PipelineObjects"]
        )

        if properties["Activate"]:
            pipeline.activate()
        return pipeline


class DataPipelineBackend(BaseBackend):
    def __init__(self, region_name: str, account_id: str):
        super().__init__(region_name, account_id)
        self.pipelines: Dict[str, Pipeline] = OrderedDict()

    def create_pipeline(self, name: str, unique_id: str, **kwargs: Any) -> Pipeline:
        pipeline = Pipeline(name, unique_id, **kwargs)
        self.pipelines[pipeline.pipeline_id] = pipeline
        return pipeline

    def list_pipelines(self) -> Iterable[Pipeline]:
        return self.pipelines.values()

    def describe_pipelines(self, pipeline_ids: List[str]) -> List[Pipeline]:
        pipelines = [
            pipeline
            for pipeline in self.pipelines.values()
            if pipeline.pipeline_id in pipeline_ids
        ]
        return pipelines

    def get_pipeline(self, pipeline_id: str) -> Pipeline:
        return self.pipelines[pipeline_id]

    def delete_pipeline(self, pipeline_id: str) -> None:
        self.pipelines.pop(pipeline_id, None)

    def put_pipeline_definition(self, pipeline_id: str, pipeline_objects: Any) -> None:
        pipeline = self.get_pipeline(pipeline_id)
        pipeline.set_pipeline_objects(pipeline_objects)

    def get_pipeline_definition(self, pipeline_id: str) -> Any:
        pipeline = self.get_pipeline(pipeline_id)
        return pipeline.objects

    def describe_objects(self, object_ids: List[str], pipeline_id: str) -> List[Any]:
        pipeline = self.get_pipeline(pipeline_id)
        pipeline_objects = [
            pipeline_object
            for pipeline_object in pipeline.objects
            if pipeline_object.object_id in object_ids
        ]
        return pipeline_objects

    def activate_pipeline(self, pipeline_id: str) -> None:
        pipeline = self.get_pipeline(pipeline_id)
        pipeline.activate()


datapipeline_backends = BackendDict(DataPipelineBackend, "datapipeline")
