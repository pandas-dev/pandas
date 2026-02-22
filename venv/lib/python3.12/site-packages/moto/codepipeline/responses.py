import json

from moto.core.responses import BaseResponse

from .models import CodePipelineBackend, codepipeline_backends


class CodePipelineResponse(BaseResponse):
    def __init__(self) -> None:
        super().__init__(service_name="codepipeline")

    @property
    def codepipeline_backend(self) -> CodePipelineBackend:
        return codepipeline_backends[self.current_account][self.region]

    def create_pipeline(self) -> str:
        pipeline, tags = self.codepipeline_backend.create_pipeline(
            self._get_param("pipeline"), self._get_param("tags")
        )

        return json.dumps({"pipeline": pipeline, "tags": tags})

    def get_pipeline(self) -> str:
        pipeline, metadata = self.codepipeline_backend.get_pipeline(
            self._get_param("name")
        )

        return json.dumps({"pipeline": pipeline, "metadata": metadata})

    def update_pipeline(self) -> str:
        pipeline = self.codepipeline_backend.update_pipeline(
            self._get_param("pipeline")
        )

        return json.dumps({"pipeline": pipeline})

    def list_pipelines(self) -> str:
        pipelines = self.codepipeline_backend.list_pipelines()

        return json.dumps({"pipelines": pipelines})

    def delete_pipeline(self) -> str:
        self.codepipeline_backend.delete_pipeline(self._get_param("name"))

        return ""

    def list_tags_for_resource(self) -> str:
        tags = self.codepipeline_backend.list_tags_for_resource(
            self._get_param("resourceArn")
        )

        return json.dumps({"tags": tags})

    def tag_resource(self) -> str:
        self.codepipeline_backend.tag_resource(
            self._get_param("resourceArn"), self._get_param("tags")
        )

        return ""

    def untag_resource(self) -> str:
        self.codepipeline_backend.untag_resource(
            self._get_param("resourceArn"), self._get_param("tagKeys")
        )

        return ""
