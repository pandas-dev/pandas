"""ComprehendBackend class with methods for supported APIs."""

from typing import Any, Dict, Iterable, List

from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.common_models import BaseModel
from moto.utilities.tagging_service import TaggingService
from moto.utilities.utils import get_partition

from .exceptions import (
    DetectPIIValidationException,
    ResourceNotFound,
    TextSizeLimitExceededException,
)

CANNED_DETECT_RESPONSE = [
    {
        "Score": 0.9999890923500061,
        "Type": "NAME",
        "BeginOffset": 50,
        "EndOffset": 58,
    },
    {
        "Score": 0.9999966621398926,
        "Type": "EMAIL",
        "BeginOffset": 230,
        "EndOffset": 259,
    },
    {
        "Score": 0.9999954700469971,
        "Type": "BANK_ACCOUNT_NUMBER",
        "BeginOffset": 334,
        "EndOffset": 349,
    },
]

CANNED_PHRASES_RESPONSE = [
    {
        "Score": 0.9999890923500061,
        "BeginOffset": 50,
        "EndOffset": 58,
    },
    {
        "Score": 0.9999966621398926,
        "BeginOffset": 230,
        "EndOffset": 259,
    },
    {
        "Score": 0.9999954700469971,
        "BeginOffset": 334,
        "EndOffset": 349,
    },
]

CANNED_SENTIMENT_RESPONSE = {
    "Sentiment": "NEUTRAL",
    "SentimentScore": {
        "Positive": 0.008101312443614006,
        "Negative": 0.0002824589901138097,
        "Neutral": 0.9916020035743713,
        "Mixed": 1.4156351426208857e-05,
    },
}


class EntityRecognizer(BaseModel):
    def __init__(
        self,
        region_name: str,
        account_id: str,
        language_code: str,
        input_data_config: Dict[str, Any],
        data_access_role_arn: str,
        version_name: str,
        recognizer_name: str,
        volume_kms_key_id: str,
        vpc_config: Dict[str, List[str]],
        model_kms_key_id: str,
        model_policy: str,
    ):
        self.name = recognizer_name
        self.arn = f"arn:{get_partition(region_name)}:comprehend:{region_name}:{account_id}:entity-recognizer/{recognizer_name}"
        if version_name:
            self.arn += f"/version/{version_name}"
        self.language_code = language_code
        self.input_data_config = input_data_config
        self.data_access_role_arn = data_access_role_arn
        self.version_name = version_name
        self.volume_kms_key_id = volume_kms_key_id
        self.vpc_config = vpc_config
        self.model_kms_key_id = model_kms_key_id
        self.model_policy = model_policy
        self.status = "TRAINED"

    def to_dict(self) -> Dict[str, Any]:
        return {
            "EntityRecognizerArn": self.arn,
            "LanguageCode": self.language_code,
            "Status": self.status,
            "InputDataConfig": self.input_data_config,
            "DataAccessRoleArn": self.data_access_role_arn,
            "VersionName": self.version_name,
            "VolumeKmsKeyId": self.volume_kms_key_id,
            "VpcConfig": self.vpc_config,
            "ModelKmsKeyId": self.model_kms_key_id,
            "ModelPolicy": self.model_policy,
        }


class ComprehendBackend(BaseBackend):
    """Implementation of Comprehend APIs."""

    # https://boto3.amazonaws.com/v1/documentation/api/latest/reference/services/comprehend/client/detect_key_phrases.html
    detect_key_phrases_languages = [
        "ar",
        "hi",
        "ko",
        "zh-TW",
        "ja",
        "zh",
        "de",
        "pt",
        "en",
        "it",
        "fr",
        "es",
    ]
    # https://boto3.amazonaws.com/v1/documentation/api/latest/reference/services/comprehend/client/detect_pii_entities.html
    detect_pii_entities_languages = ["en"]

    def __init__(self, region_name: str, account_id: str):
        super().__init__(region_name, account_id)
        self.recognizers: Dict[str, EntityRecognizer] = dict()
        self.tagger = TaggingService()

    def list_entity_recognizers(
        self, _filter: Dict[str, Any]
    ) -> Iterable[EntityRecognizer]:
        """
        Pagination is not yet implemented.
        The following filters are not yet implemented: Status, SubmitTimeBefore, SubmitTimeAfter
        """
        if "RecognizerName" in _filter:
            return [
                entity
                for entity in self.recognizers.values()
                if entity.name == _filter["RecognizerName"]
            ]
        return self.recognizers.values()

    def create_entity_recognizer(
        self,
        recognizer_name: str,
        version_name: str,
        data_access_role_arn: str,
        tags: List[Dict[str, str]],
        input_data_config: Dict[str, Any],
        language_code: str,
        volume_kms_key_id: str,
        vpc_config: Dict[str, List[str]],
        model_kms_key_id: str,
        model_policy: str,
    ) -> str:
        """
        The ClientRequestToken-parameter is not yet implemented
        """
        recognizer = EntityRecognizer(
            region_name=self.region_name,
            account_id=self.account_id,
            language_code=language_code,
            input_data_config=input_data_config,
            data_access_role_arn=data_access_role_arn,
            version_name=version_name,
            recognizer_name=recognizer_name,
            volume_kms_key_id=volume_kms_key_id,
            vpc_config=vpc_config,
            model_kms_key_id=model_kms_key_id,
            model_policy=model_policy,
        )
        self.recognizers[recognizer.arn] = recognizer
        self.tagger.tag_resource(recognizer.arn, tags)
        return recognizer.arn

    def describe_entity_recognizer(
        self, entity_recognizer_arn: str
    ) -> EntityRecognizer:
        if entity_recognizer_arn not in self.recognizers:
            raise ResourceNotFound
        return self.recognizers[entity_recognizer_arn]

    def stop_training_entity_recognizer(self, entity_recognizer_arn: str) -> None:
        recognizer = self.describe_entity_recognizer(entity_recognizer_arn)
        if recognizer.status == "TRAINING":
            recognizer.status = "STOP_REQUESTED"

    def list_tags_for_resource(self, resource_arn: str) -> List[Dict[str, str]]:
        return self.tagger.list_tags_for_resource(resource_arn)["Tags"]

    def delete_entity_recognizer(self, entity_recognizer_arn: str) -> None:
        self.recognizers.pop(entity_recognizer_arn, None)

    def tag_resource(self, resource_arn: str, tags: List[Dict[str, str]]) -> None:
        self.tagger.tag_resource(resource_arn, tags)

    def untag_resource(self, resource_arn: str, tag_keys: List[str]) -> None:
        self.tagger.untag_resource_using_names(resource_arn, tag_keys)

    def detect_pii_entities(self, text: str, language: str) -> List[Dict[str, Any]]:
        if language not in self.detect_pii_entities_languages:
            raise DetectPIIValidationException(
                language, self.detect_pii_entities_languages
            )
        text_size = len(text)
        if text_size > 100000:
            raise TextSizeLimitExceededException(text_size)
        return CANNED_DETECT_RESPONSE

    def detect_key_phrases(self, text: str, language: str) -> List[Dict[str, Any]]:
        if language not in self.detect_key_phrases_languages:
            raise DetectPIIValidationException(
                language, self.detect_key_phrases_languages
            )
        text_size = len(text)
        if text_size > 100000:
            raise TextSizeLimitExceededException(text_size)
        return CANNED_PHRASES_RESPONSE

    def detect_sentiment(self, text: str, language: str) -> Dict[str, Any]:
        if language not in self.detect_key_phrases_languages:
            raise DetectPIIValidationException(
                language, self.detect_key_phrases_languages
            )
        text_size = len(text)
        if text_size > 5000:
            raise TextSizeLimitExceededException(text_size)
        return CANNED_SENTIMENT_RESPONSE


comprehend_backends = BackendDict(ComprehendBackend, "comprehend")
