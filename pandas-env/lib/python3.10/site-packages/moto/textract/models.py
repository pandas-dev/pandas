from collections import defaultdict
from typing import Any, Dict, List

from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.common_models import BaseModel
from moto.moto_api._internal import mock_random

from .exceptions import InvalidJobIdException, InvalidParameterException


class TextractJobStatus:
    in_progress = "IN_PROGRESS"
    succeeded = "SUCCEEDED"
    failed = "FAILED"
    partial_success = "PARTIAL_SUCCESS"


class TextractJob(BaseModel):
    def __init__(self, job: Dict[str, Any]):
        self.job = job

    def to_dict(self) -> Dict[str, Any]:
        return self.job


class TextractBackend(BaseBackend):
    """Implementation of Textract APIs."""

    JOB_STATUS = TextractJobStatus.succeeded
    PAGES = {"Pages": mock_random.randint(5, 500)}
    BLOCKS: List[Dict[str, Any]] = []

    def __init__(self, region_name: str, account_id: str):
        super().__init__(region_name, account_id)
        self.async_text_detection_jobs: Dict[str, TextractJob] = defaultdict()

    def get_document_text_detection(self, job_id: str) -> TextractJob:
        """
        Pagination has not yet been implemented
        """
        job = self.async_text_detection_jobs.get(job_id)
        if not job:
            raise InvalidJobIdException()
        return job

    def detect_document_text(self) -> Dict[str, Any]:
        return {
            "Blocks": TextractBackend.BLOCKS,
            "DetectDocumentTextModelVersion": "1.0",
            "DocumentMetadata": {"Pages": TextractBackend.PAGES},
        }

    def start_document_text_detection(self, document_location: str) -> str:
        """
        The following parameters have not yet been implemented: ClientRequestToken, JobTag, NotificationChannel, OutputConfig, KmsKeyID
        """
        if not document_location:
            raise InvalidParameterException()
        job_id = str(mock_random.uuid4())
        self.async_text_detection_jobs[job_id] = TextractJob(
            {
                "Blocks": TextractBackend.BLOCKS,
                "DetectDocumentTextModelVersion": "1.0",
                "DocumentMetadata": {"Pages": TextractBackend.PAGES},
                "JobStatus": TextractBackend.JOB_STATUS,
            }
        )
        return job_id


textract_backends = BackendDict(TextractBackend, "textract")
