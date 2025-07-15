import string
from typing import Any, Dict, List, Tuple

from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.common_models import BaseModel
from moto.moto_api._internal import mock_random as random
from moto.utilities.utils import get_partition


class Pipeline(BaseModel):
    def __init__(
        self,
        account_id: str,
        region: str,
        name: str,
        input_bucket: str,
        output_bucket: str,
        role: str,
        content_config: Dict[str, Any],
        thumbnail_config: Dict[str, Any],
    ):
        a = "".join(random.choice(string.digits) for _ in range(13))
        b = "".join(random.choice(string.ascii_lowercase) for _ in range(6))
        self.id = f"{a}-{b}"
        self.name = name
        self.arn = f"arn:{get_partition(region)}:elastictranscoder:{region}:{account_id}:pipeline/{self.id}"
        self.status = "Active"
        self.input_bucket = input_bucket
        self.output_bucket = output_bucket or content_config["Bucket"]
        self.role = role
        self.content_config = content_config or {"Bucket": self.output_bucket}
        if "Permissions" not in self.content_config:
            self.content_config["Permissions"] = []
        self.thumbnail_config = thumbnail_config or {"Bucket": self.output_bucket}
        if "Permissions" not in self.thumbnail_config:
            self.thumbnail_config["Permissions"] = []

    def update(self, name: str, input_bucket: str, role: str) -> None:
        if name:
            self.name = name
        if input_bucket:
            self.input_bucket = input_bucket
        if role:
            self.role = role

    def to_dict(self) -> Dict[str, Any]:
        return {
            "Id": self.id,
            "Name": self.name,
            "Arn": self.arn,
            "Status": self.status,
            "InputBucket": self.input_bucket,
            "OutputBucket": self.output_bucket,
            "Role": self.role,
            "Notifications": {
                "Progressing": "",
                "Completed": "",
                "Warning": "",
                "Error": "",
            },
            "ContentConfig": self.content_config,
            "ThumbnailConfig": self.thumbnail_config,
        }


class ElasticTranscoderBackend(BaseBackend):
    def __init__(self, region_name: str, account_id: str):
        super().__init__(region_name, account_id)
        self.pipelines: Dict[str, Pipeline] = {}

    def create_pipeline(
        self,
        name: str,
        input_bucket: str,
        output_bucket: str,
        role: str,
        content_config: Dict[str, Any],
        thumbnail_config: Dict[str, Any],
    ) -> Tuple[Pipeline, List[str]]:
        """
        The following parameters are not yet implemented:
        AWSKMSKeyArn, Notifications
        """
        pipeline = Pipeline(
            self.account_id,
            self.region_name,
            name,
            input_bucket,
            output_bucket,
            role,
            content_config,
            thumbnail_config,
        )
        self.pipelines[pipeline.id] = pipeline
        warnings: List[str] = []
        return pipeline, warnings

    def list_pipelines(self) -> List[Dict[str, Any]]:
        return [p.to_dict() for _, p in self.pipelines.items()]

    def read_pipeline(self, pipeline_id: str) -> Pipeline:
        return self.pipelines[pipeline_id]

    def update_pipeline(
        self, pipeline_id: str, name: str, input_bucket: str, role: str
    ) -> Tuple[Pipeline, List[str]]:
        """
        The following parameters are not yet implemented:
        AWSKMSKeyArn, Notifications, ContentConfig, ThumbnailConfig
        """
        pipeline = self.read_pipeline(pipeline_id)
        pipeline.update(name, input_bucket, role)
        warnings: List[str] = []
        return pipeline, warnings

    def delete_pipeline(self, pipeline_id: str) -> None:
        self.pipelines.pop(pipeline_id)


elastictranscoder_backends = BackendDict(ElasticTranscoderBackend, "elastictranscoder")
