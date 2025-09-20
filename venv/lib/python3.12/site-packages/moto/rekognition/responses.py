import json

from moto.core.common_types import TYPE_RESPONSE
from moto.core.responses import BaseResponse

from .models import RekognitionBackend, rekognition_backends


class RekognitionResponse(BaseResponse):
    """Handler for Rekognition requests and responses."""

    def __init__(self) -> None:
        super().__init__(service_name="rekognition")

    @property
    def rekognition_backend(self) -> RekognitionBackend:
        return rekognition_backends[self.current_account][self.region]

    def get_face_search(self) -> str:
        (
            job_status,
            status_message,
            video_metadata,
            persons,
            next_token,
            text_model_version,
        ) = self.rekognition_backend.get_face_search()

        return json.dumps(
            dict(
                JobStatus=job_status,
                StatusMessage=status_message,
                VideoMetadata=video_metadata,
                Persons=persons,
                NextToken=next_token,
                TextModelVersion=text_model_version,
            )
        )

    def get_text_detection(self) -> str:
        (
            job_status,
            status_message,
            video_metadata,
            text_detections,
            next_token,
            text_model_version,
        ) = self.rekognition_backend.get_text_detection()

        return json.dumps(
            dict(
                JobStatus=job_status,
                StatusMessage=status_message,
                VideoMetadata=video_metadata,
                TextDetections=text_detections,
                NextToken=next_token,
                TextModelVersion=text_model_version,
            )
        )

    def compare_faces(self) -> str:
        (
            face_matches,
            source_image_orientation_correction,
            target_image_orientation_correction,
            unmatched_faces,
            source_image_face,
        ) = self.rekognition_backend.compare_faces()

        return json.dumps(
            dict(
                FaceMatches=face_matches,
                SourceImageOrientationCorrection=source_image_orientation_correction,
                TargetImageOrientationCorrection=target_image_orientation_correction,
                UnmatchedFaces=unmatched_faces,
                SourceImageFace=source_image_face,
            )
        )

    def detect_labels(self) -> str:
        (
            labels,
            image_properties,
            label_model_version,
        ) = self.rekognition_backend.detect_labels()
        return json.dumps(
            dict(
                Labels=labels,
                ImageProperties=image_properties,
                LabelModelVersion=label_model_version,
            )
        )

    def detect_text(self) -> str:
        (
            text_detections,
            text_model_version,
        ) = self.rekognition_backend.detect_text()
        return json.dumps(
            dict(
                TextDetections=text_detections,
                TextModelVersion=text_model_version,
            )
        )

    def detect_custom_labels(self) -> str:
        (custom_labels,) = self.rekognition_backend.detect_custom_labels()
        return json.dumps(
            dict(
                CustomLabels=custom_labels,
            )
        )

    def start_face_search(self) -> TYPE_RESPONSE:
        headers = {"Content-Type": "application/x-amz-json-1.1"}
        job_id = self.rekognition_backend.start_face_search()
        response = ('{"JobId":"' + job_id + '"}').encode()

        return 200, headers, response

    def start_text_detection(self) -> TYPE_RESPONSE:
        headers = {"Content-Type": "application/x-amz-json-1.1"}
        job_id = self.rekognition_backend.start_text_detection()
        response = ('{"JobId":"' + job_id + '"}').encode()

        return 200, headers, response
