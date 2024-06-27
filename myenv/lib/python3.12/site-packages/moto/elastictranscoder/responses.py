import json
import re
from typing import Any, Optional

from moto.core.common_types import TYPE_RESPONSE
from moto.core.responses import BaseResponse

from .models import ElasticTranscoderBackend, elastictranscoder_backends


class ElasticTranscoderResponse(BaseResponse):
    def __init__(self) -> None:
        super().__init__(service_name="elastictranscoder")

    @property
    def elastictranscoder_backend(self) -> ElasticTranscoderBackend:
        return elastictranscoder_backends[self.current_account][self.region]

    def pipelines(self, request: Any, full_url: str, headers: Any) -> TYPE_RESPONSE:
        self.setup_class(request, full_url, headers)
        if request.method == "POST":
            return self.create_pipeline()
        elif request.method == "GET":
            return self.list_pipelines()
        else:
            return self.update_pipeline()

    def individual_pipeline(
        self, request: Any, full_url: str, headers: Any
    ) -> TYPE_RESPONSE:
        self.setup_class(request, full_url, headers)
        if request.method == "GET":
            return self.read_pipeline()
        elif request.method == "DELETE":
            return self.delete_pipeline()
        else:
            return self.update_pipeline()

    def create_pipeline(self) -> TYPE_RESPONSE:
        name = self._get_param("Name")
        input_bucket = self._get_param("InputBucket")
        output_bucket = self._get_param("OutputBucket")
        role = self._get_param("Role")
        content_config = self._get_param("ContentConfig")
        thumbnail_config = self._get_param("ThumbnailConfig")
        if not role:
            return self.err("Role cannot be blank")
        if not re.match("^arn:aws:iam::[0-9]+:role/.+", role):
            return self.err(f"Role ARN is invalid: {role}")
        if not output_bucket and not content_config:
            return self.err(
                "[OutputBucket and ContentConfig:Bucket are not allowed to both be null.]"
            )
        if output_bucket and content_config:
            return self.err("[OutputBucket and ContentConfig are mutually exclusive.]")
        if content_config and not thumbnail_config:
            return self.err(
                "[ThumbnailConfig:Bucket is not allowed to be null if ContentConfig is specified.]"
            )
        pipeline, warnings = self.elastictranscoder_backend.create_pipeline(
            name=name,
            input_bucket=input_bucket,
            output_bucket=output_bucket,
            role=role,
            content_config=content_config,
            thumbnail_config=thumbnail_config,
        )
        return (
            201,
            {"status": 201},
            json.dumps({"Pipeline": pipeline.to_dict(), "Warnings": warnings}),
        )

    def list_pipelines(self) -> TYPE_RESPONSE:
        return (
            200,
            {},
            json.dumps({"Pipelines": self.elastictranscoder_backend.list_pipelines()}),
        )

    def validate_pipeline_id(self, pipeline_id: str) -> Optional[TYPE_RESPONSE]:
        r = "^\\d{13}-\\w{6}$"
        if not re.match(r, pipeline_id):
            msg = f"1 validation error detected: Value '{pipeline_id}' at 'id' failed to satisfy constraint: Member must satisfy regular expression pattern: {r}"
            return self.err(msg)
        try:
            self.elastictranscoder_backend.read_pipeline(pipeline_id)
            return None
        except KeyError:
            msg = f"The specified pipeline was not found: account={self.current_account}, pipelineId={pipeline_id}."
            return (
                404,
                {"status": 404, "x-amzn-ErrorType": "ResourceNotFoundException"},
                json.dumps({"message": msg}),
            )

    def read_pipeline(self) -> TYPE_RESPONSE:
        _id = self.path.rsplit("/", 1)[-1]
        err = self.validate_pipeline_id(_id)
        if err:
            return err
        pipeline = self.elastictranscoder_backend.read_pipeline(_id)
        return 200, {}, json.dumps({"Pipeline": pipeline.to_dict()})

    def update_pipeline(self) -> TYPE_RESPONSE:
        _id = self.path.rsplit("/", 1)[-1]
        name = self._get_param("Name")
        input_bucket = self._get_param("InputBucket")
        role = self._get_param("Role")
        err = self.validate_pipeline_id(_id)
        if err:
            return err
        pipeline, warnings = self.elastictranscoder_backend.update_pipeline(
            pipeline_id=_id, name=name, input_bucket=input_bucket, role=role
        )

        return (
            200,
            {},
            json.dumps({"Pipeline": pipeline.to_dict(), "Warnings": warnings}),
        )

    def delete_pipeline(self) -> TYPE_RESPONSE:
        _id = self.path.rsplit("/", 1)[-1]
        self.elastictranscoder_backend.delete_pipeline(pipeline_id=_id)
        return 200, {}, "{}"

    def err(self, msg: str) -> TYPE_RESPONSE:
        return (
            400,
            {"status": 400, "x-amzn-ErrorType": "ValidationException"},
            json.dumps({"message": msg}),
        )
