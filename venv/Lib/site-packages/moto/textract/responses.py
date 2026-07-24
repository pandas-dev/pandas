"""Handles incoming textract requests, invokes methods, returns responses."""

import json

from moto.core.responses import BaseResponse

from .models import TextractBackend, textract_backends


class TextractResponse(BaseResponse):
    """Handler for Textract requests and responses."""

    def __init__(self) -> None:
        super().__init__(service_name="textract")

    @property
    def textract_backend(self) -> TextractBackend:
        """Return backend instance specific for this region."""
        return textract_backends[self.current_account][self.region]

    def get_document_text_detection(self) -> str:
        params = json.loads(self.body)
        job_id = params.get("JobId")
        job = self.textract_backend.get_document_text_detection(job_id=job_id).to_dict()
        return json.dumps(job)

    def detect_document_text(self) -> str:
        result = self.textract_backend.detect_document_text()
        return json.dumps(result)

    def start_document_text_detection(self) -> str:
        params = json.loads(self.body)
        document_location = params.get("DocumentLocation")
        notification_channel = params.get("NotificationChannel")
        job_id = self.textract_backend.start_document_text_detection(
            document_location=document_location,
            notification_channel=notification_channel,
        )
        return json.dumps({"JobId": job_id})

    def start_document_analysis(self) -> str:
        params = json.loads(self.body)
        document_location = params.get("DocumentLocation")
        feature_types = params.get("FeatureTypes")
        job_id = self.textract_backend.start_document_analysis(
            document_location=document_location, feature_types=feature_types
        )
        return json.dumps({"JobId": job_id})

    def get_document_analysis(self) -> str:
        params = json.loads(self.body)
        job_id = params.get("JobId")
        max_results = params.get("MaxResults")
        next_token = params.get("NextToken")
        job = self.textract_backend.get_document_analysis(
            job_id=job_id, max_results=max_results, next_token=next_token
        ).to_dict()
        return json.dumps(job)
