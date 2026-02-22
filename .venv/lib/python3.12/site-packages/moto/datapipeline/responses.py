import json

from moto.core.responses import BaseResponse

from .models import DataPipelineBackend, datapipeline_backends


class DataPipelineResponse(BaseResponse):
    def __init__(self) -> None:
        super().__init__(service_name="datapipeline")

    @property
    def datapipeline_backend(self) -> DataPipelineBackend:
        return datapipeline_backends[self.current_account][self.region]

    def create_pipeline(self) -> str:
        name = self._get_param("name")
        unique_id = self._get_param("uniqueId")
        description = self._get_param("description", "")
        tags = self._get_param("tags", [])
        pipeline = self.datapipeline_backend.create_pipeline(
            name, unique_id, description=description, tags=tags
        )
        return json.dumps({"pipelineId": pipeline.pipeline_id})

    def list_pipelines(self) -> str:
        pipelines = list(self.datapipeline_backend.list_pipelines())
        pipeline_ids = [pipeline.pipeline_id for pipeline in pipelines]
        max_pipelines = 50
        marker = self._get_param("marker")
        if marker:
            start = pipeline_ids.index(marker) + 1
        else:
            start = 0
        pipelines_resp = pipelines[start : start + max_pipelines]
        has_more_results = False
        marker = None
        if start + max_pipelines < len(pipeline_ids) - 1:
            has_more_results = True
            marker = pipelines_resp[-1].pipeline_id
        return json.dumps(
            {
                "hasMoreResults": has_more_results,
                "marker": marker,
                "pipelineIdList": [
                    pipeline.to_meta_json() for pipeline in pipelines_resp
                ],
            }
        )

    def describe_pipelines(self) -> str:
        pipeline_ids = self._get_param("pipelineIds")
        pipelines = self.datapipeline_backend.describe_pipelines(pipeline_ids)

        return json.dumps(
            {"pipelineDescriptionList": [pipeline.to_json() for pipeline in pipelines]}
        )

    def delete_pipeline(self) -> str:
        pipeline_id = self._get_param("pipelineId")
        self.datapipeline_backend.delete_pipeline(pipeline_id)
        return json.dumps({})

    def put_pipeline_definition(self) -> str:
        pipeline_id = self._get_param("pipelineId")
        pipeline_objects = self._get_param("pipelineObjects")

        self.datapipeline_backend.put_pipeline_definition(pipeline_id, pipeline_objects)
        return json.dumps({"errored": False})

    def get_pipeline_definition(self) -> str:
        pipeline_id = self._get_param("pipelineId")
        pipeline_definition = self.datapipeline_backend.get_pipeline_definition(
            pipeline_id
        )
        return json.dumps(
            {
                "pipelineObjects": [
                    pipeline_object.to_json() for pipeline_object in pipeline_definition
                ]
            }
        )

    def describe_objects(self) -> str:
        pipeline_id = self._get_param("pipelineId")
        object_ids = self._get_param("objectIds")

        pipeline_objects = self.datapipeline_backend.describe_objects(
            object_ids, pipeline_id
        )
        return json.dumps(
            {
                "hasMoreResults": False,
                "marker": None,
                "pipelineObjects": [
                    pipeline_object.to_json() for pipeline_object in pipeline_objects
                ],
            }
        )

    def activate_pipeline(self) -> str:
        pipeline_id = self._get_param("pipelineId")
        self.datapipeline_backend.activate_pipeline(pipeline_id)
        return json.dumps({})
