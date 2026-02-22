import json
import time
from collections import defaultdict
from typing import Any, Optional

from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.common_models import BaseModel
from moto.moto_api._internal import mock_random
from moto.sns import sns_backends
from moto.sns.exceptions import TopicNotFound

from .exceptions import InvalidJobIdException, InvalidParameterException


class TextractJobStatus:
    in_progress = "IN_PROGRESS"
    succeeded = "SUCCEEDED"
    failed = "FAILED"
    partial_success = "PARTIAL_SUCCESS"


class TextractJob(BaseModel):
    def __init__(
        self, job: dict[str, Any], notification_channel: Optional[dict[str, str]] = None
    ):
        self.job = job
        self.notification_channel = notification_channel
        self.job_id = str(mock_random.uuid4())

    def to_dict(self) -> dict[str, Any]:
        return self.job

    def send_completion_notification(
        self, account_id: str, region_name: str, document_location: dict[str, Any]
    ) -> None:
        if not self.notification_channel:
            return

        topic_arn = self.notification_channel.get("SNSTopicArn")
        if not topic_arn:
            return

        # Convert document_location from {'S3Object': {'Bucket': '...', 'Name': '...'}} format
        # to {'S3Bucket': '...', 'S3ObjectName': '...'} format as per AWS docs
        s3_object = document_location.get("S3Object", {})
        doc_location = {
            "S3Bucket": s3_object.get("Bucket", ""),
            "S3ObjectName": s3_object.get("Name", ""),
        }

        notification = {
            "JobId": self.job_id,
            "Status": self.job["JobStatus"],
            "API": "StartDocumentTextDetection",
            "JobTag": "",  # Not implemented yet
            "Timestamp": int(time.time() * 1000),  # Convert to milliseconds
            "DocumentLocation": doc_location,
        }

        sns_backend = sns_backends[account_id][region_name]
        try:
            sns_backend.publish(
                message=json.dumps(notification),  # SNS requires message to be a string
                arn=topic_arn,
                subject="Amazon Textract Job Completion",
            )
        except TopicNotFound:
            pass


class TextractBackend(BaseBackend):
    """Implementation of Textract APIs."""

    JOB_STATUS = TextractJobStatus.succeeded
    PAGES = {"Pages": mock_random.randint(5, 500)}
    BLOCKS: list[dict[str, Any]] = []

    def __init__(self, region_name: str, account_id: str):
        super().__init__(region_name, account_id)
        self.async_text_detection_jobs: dict[str, TextractJob] = defaultdict()
        self.async_document_analysis_jobs: dict[str, TextractJob] = defaultdict()

    def get_document_text_detection(self, job_id: str) -> TextractJob:
        """
        Pagination has not yet been implemented
        """
        job = self.async_text_detection_jobs.get(job_id)
        if not job:
            raise InvalidJobIdException()
        return job

    def detect_document_text(self) -> dict[str, Any]:
        return {
            "Blocks": TextractBackend.BLOCKS,
            "DetectDocumentTextModelVersion": "1.0",
            "DocumentMetadata": {"Pages": TextractBackend.PAGES},
        }

    def start_document_text_detection(
        self,
        document_location: dict[str, Any],
        notification_channel: Optional[dict[str, str]] = None,
    ) -> str:
        """
        The following parameters have not yet been implemented: ClientRequestToken, JobTag, OutputConfig, KmsKeyID
        """
        if not document_location:
            raise InvalidParameterException()

        job = TextractJob(
            {
                "Blocks": TextractBackend.BLOCKS,
                "DetectDocumentTextModelVersion": "1.0",
                "DocumentMetadata": {"Pages": TextractBackend.PAGES},
                "JobStatus": TextractBackend.JOB_STATUS,
            },
            notification_channel=notification_channel,
        )

        self.async_text_detection_jobs[job.job_id] = job

        # Send completion notification since we're mocking an immediate completion
        job.send_completion_notification(
            self.account_id, self.region_name, document_location
        )

        return job.job_id

    def start_document_analysis(
        self, document_location: dict[str, Any], feature_types: list[str]
    ) -> str:
        """
        The following parameters have not yet been implemented: ClientRequestToken, JobTag, NotificationChannel, OutputConfig, KmsKeyID
        """
        if not document_location or not feature_types:
            raise InvalidParameterException()
        job_id = str(mock_random.uuid4())
        self.async_document_analysis_jobs[job_id] = TextractJob(
            {
                "Blocks": TextractBackend.BLOCKS,
                "DetectDocumentTextModelVersion": "1.0",
                "DocumentMetadata": {"Pages": TextractBackend.PAGES},
                "JobStatus": TextractBackend.JOB_STATUS,
            }
        )
        return job_id

    def get_document_analysis(
        self, job_id: str, max_results: Optional[int], next_token: Optional[str] = None
    ) -> TextractJob:
        job = self.async_document_analysis_jobs.get(job_id)
        if not job:
            raise InvalidJobIdException()
        return job


textract_backends = BackendDict(TextractBackend, "textract")
