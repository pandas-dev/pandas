import re
from datetime import datetime, timedelta
from typing import Any, Dict, List, Optional

from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.common_models import BaseModel
from moto.moto_api._internal import mock_random
from moto.moto_api._internal.managed_state_model import ManagedState

from .exceptions import BadRequestException, ConflictException


class BaseObject(BaseModel):
    def camelCase(self, key: str) -> str:
        words = []
        for word in key.split("_"):
            words.append(word.title())
        return "".join(words)

    def gen_response_object(self) -> Dict[str, Any]:
        response_object: Dict[str, Any] = dict()
        for key, value in self.__dict__.items():
            if "_" in key:
                response_object[self.camelCase(key)] = value
            else:
                response_object[key[0].upper() + key[1:]] = value
        return response_object

    @property
    def response_object(self) -> Dict[str, Any]:  # type: ignore[misc]
        return self.gen_response_object()


class FakeTranscriptionJob(BaseObject, ManagedState):
    def __init__(
        self,
        account_id: str,
        region_name: str,
        transcription_job_name: str,
        language_code: Optional[str],
        media_sample_rate_hertz: Optional[int],
        media_format: Optional[str],
        media: Dict[str, str],
        output_bucket_name: Optional[str],
        output_key: Optional[str],
        output_encryption_kms_key_id: Optional[str],
        settings: Optional[Dict[str, Any]],
        model_settings: Optional[Dict[str, Optional[str]]],
        job_execution_settings: Optional[Dict[str, Any]],
        content_redaction: Optional[Dict[str, Any]],
        identify_language: Optional[bool],
        identify_multiple_languages: Optional[bool],
        language_options: Optional[List[str]],
        subtitles: Optional[Dict[str, Any]],
    ):
        ManagedState.__init__(
            self,
            "transcribe::transcriptionjob",
            transitions=[
                (None, "QUEUED"),
                ("QUEUED", "IN_PROGRESS"),
                ("IN_PROGRESS", "COMPLETED"),
            ],
        )
        self._account_id = account_id
        self._region_name = region_name
        self.transcription_job_name = transcription_job_name
        self.language_code = language_code
        self.language_codes: Optional[List[Dict[str, Any]]] = None
        self.media_sample_rate_hertz = media_sample_rate_hertz
        self.media_format = media_format
        self.media = media
        self.transcript: Optional[Dict[str, str]] = None
        self.start_time: Optional[str] = None
        self.completion_time: Optional[str] = None
        self.creation_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        self.failure_reason = None
        self.settings = settings or {
            "ChannelIdentification": False,
            "ShowAlternatives": False,
            "ShowSpeakerLabels": False,
        }
        self.model_settings = model_settings or {"LanguageModelName": None}
        self.job_execution_settings = job_execution_settings or {
            "AllowDeferredExecution": False,
            "DataAccessRoleArn": None,
        }
        self.content_redaction = content_redaction or {
            "RedactionType": None,
            "RedactionOutput": None,
        }
        self.identify_language = identify_language
        self.identify_multiple_languages = identify_multiple_languages
        self.language_options = language_options
        self.identified_language_score: Optional[float] = None
        self._output_bucket_name = output_bucket_name
        self.output_key = output_key
        self._output_encryption_kms_key_id = output_encryption_kms_key_id
        self.output_location_type = (
            "CUSTOMER_BUCKET" if self._output_bucket_name else "SERVICE_BUCKET"
        )
        self.subtitles = subtitles or {"Formats": [], "OutputStartIndex": 0}

    def response_object(self, response_type: str) -> Dict[str, Any]:  # type: ignore
        response_field_dict = {
            "CREATE": [
                "TranscriptionJobName",
                "TranscriptionJobStatus",
                "LanguageCode",
                "LanguageCodes",
                "MediaFormat",
                "Media",
                "Settings",
                "StartTime",
                "CreationTime",
                "IdentifyLanguage",
                "IdentifyMultipleLanguages",
                "LanguageOptions",
                "JobExecutionSettings",
                "Subtitles",
            ],
            "GET": [
                "TranscriptionJobName",
                "TranscriptionJobStatus",
                "LanguageCode",
                "LanguageCodes",
                "MediaSampleRateHertz",
                "MediaFormat",
                "Media",
                "Settings",
                "Transcript",
                "StartTime",
                "CreationTime",
                "CompletionTime",
                "IdentifyLanguage",
                "IdentifyMultipleLanguages",
                "LanguageOptions",
                "IdentifiedLanguageScore",
                "Subtitles",
            ],
            "LIST": [
                "TranscriptionJobName",
                "CreationTime",
                "StartTime",
                "CompletionTime",
                "LanguageCode",
                "LanguageCodes",
                "TranscriptionJobStatus",
                "FailureReason",
                "IdentifyLanguage",
                "IdentifyMultipleLanguages",
                "IdentifiedLanguageScore",
                "OutputLocationType",
            ],
        }
        response_fields = response_field_dict[response_type]
        response_object = self.gen_response_object()
        response_object["TranscriptionJobStatus"] = self.status
        if response_type != "LIST":
            return {
                "TranscriptionJob": {
                    k: v
                    for k, v in response_object.items()
                    if k in response_fields and v is not None and v != [None]
                }
            }
        else:
            return {
                k: v
                for k, v in response_object.items()
                if k in response_fields and v is not None and v != [None]
            }

    def advance(self) -> None:
        old_status = self.status
        super().advance()
        new_status = self.status

        if old_status == new_status:
            return

        if new_status == "IN_PROGRESS":
            self.start_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            if not self.media_sample_rate_hertz:
                self.media_sample_rate_hertz = 44100
            if not self.media_format:
                file_ext = self.media["MediaFileUri"].split(".")[-1].lower()
                self.media_format = (
                    file_ext if file_ext in ["mp3", "mp4", "wav", "flac"] else "mp3"
                )
            if self.identify_language:
                self.identified_language_score = 0.999645948
                # Simply identify first language passed in language_options
                # If none is set, default to "en-US"
                if self.language_options is not None and len(self.language_options) > 0:
                    self.language_code = self.language_options[0]
                else:
                    self.language_code = "en-US"
            if self.identify_multiple_languages:
                self.identified_language_score = 0.999645948
                # Identify first two languages passed in language_options
                # If none is set, default to "en-US"
                self.language_codes: List[Dict[str, Any]] = []  # type: ignore[no-redef]
                if self.language_options is None or len(self.language_options) == 0:
                    self.language_codes.append(
                        {"LanguageCode": "en-US", "DurationInSeconds": 123.0}
                    )
                else:
                    self.language_codes.append(
                        {
                            "LanguageCode": self.language_options[0],
                            "DurationInSeconds": 123.0,
                        }
                    )
                    if len(self.language_options) > 1:
                        self.language_codes.append(
                            {
                                "LanguageCode": self.language_options[1],
                                "DurationInSeconds": 321.0,
                            }
                        )
        elif new_status == "COMPLETED":
            self.completion_time = (datetime.now() + timedelta(seconds=10)).strftime(
                "%Y-%m-%d %H:%M:%S"
            )
            if self._output_bucket_name:
                remove_json_extension = re.compile("\\.json$")
                transcript_file_prefix = (
                    f"https://s3.{self._region_name}.amazonaws.com/"
                    f"{self._output_bucket_name}/"
                    f"{remove_json_extension.sub('', self.output_key or self.transcription_job_name)}"
                )
                self.output_location_type = "CUSTOMER_BUCKET"
            else:
                transcript_file_prefix = (
                    f"https://s3.{self._region_name}.amazonaws.com/"
                    f"aws-transcribe-{self._region_name}-prod/"
                    f"{self._account_id}/"
                    f"{self.transcription_job_name}/"
                    f"{mock_random.uuid4()}/"
                    "asrOutput"
                )
                self.output_location_type = "SERVICE_BUCKET"
            self.transcript = {"TranscriptFileUri": f"{transcript_file_prefix}.json"}
            self.subtitles["SubtitleFileUris"] = [
                f"{transcript_file_prefix}.{format}"
                for format in self.subtitles["Formats"]
            ]


class FakeVocabulary(BaseObject, ManagedState):
    def __init__(
        self,
        account_id: str,
        region_name: str,
        vocabulary_name: str,
        language_code: str,
        phrases: Optional[List[str]],
        vocabulary_file_uri: Optional[str],
    ):
        # Configured ManagedState
        super().__init__(
            "transcribe::vocabulary",
            transitions=[(None, "PENDING"), ("PENDING", "READY")],
        )
        # Configure internal properties
        self._region_name = region_name
        self.vocabulary_name = vocabulary_name
        self.language_code = language_code
        self.phrases = phrases
        self.vocabulary_file_uri = vocabulary_file_uri
        self.last_modified_time: Optional[str] = None
        self.failure_reason = None
        self.download_uri = f"https://s3.{region_name}.amazonaws.com/aws-transcribe-dictionary-model-{region_name}-prod/{account_id}/{vocabulary_name}/{mock_random.uuid4()}/input.txt"

    def response_object(self, response_type: str) -> Dict[str, Any]:  # type: ignore
        response_field_dict = {
            "CREATE": [
                "VocabularyName",
                "LanguageCode",
                "VocabularyState",
                "LastModifiedTime",
                "FailureReason",
            ],
            "GET": [
                "VocabularyName",
                "LanguageCode",
                "VocabularyState",
                "LastModifiedTime",
                "FailureReason",
                "DownloadUri",
            ],
            "LIST": [
                "VocabularyName",
                "LanguageCode",
                "LastModifiedTime",
                "VocabularyState",
            ],
        }
        response_fields = response_field_dict[response_type]
        response_object = self.gen_response_object()
        response_object["VocabularyState"] = self.status
        return {
            k: v
            for k, v in response_object.items()
            if k in response_fields and v is not None and v != [None]
        }

    def advance(self) -> None:
        old_status = self.status
        super().advance()
        new_status = self.status

        if old_status != new_status:
            self.last_modified_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")


class FakeMedicalTranscriptionJob(BaseObject, ManagedState):
    def __init__(
        self,
        region_name: str,
        medical_transcription_job_name: str,
        language_code: str,
        media_sample_rate_hertz: Optional[int],
        media_format: Optional[str],
        media: Dict[str, str],
        output_bucket_name: str,
        output_encryption_kms_key_id: Optional[str],
        settings: Optional[Dict[str, Any]],
        specialty: str,
        job_type: str,
    ):
        ManagedState.__init__(
            self,
            "transcribe::medicaltranscriptionjob",
            transitions=[
                (None, "QUEUED"),
                ("QUEUED", "IN_PROGRESS"),
                ("IN_PROGRESS", "COMPLETED"),
            ],
        )
        self._region_name = region_name
        self.medical_transcription_job_name = medical_transcription_job_name
        self.language_code = language_code
        self.media_sample_rate_hertz = media_sample_rate_hertz
        self.media_format = media_format
        self.media = media
        self.transcript: Optional[Dict[str, str]] = None
        self.start_time: Optional[str] = None
        self.completion_time: Optional[str] = None
        self.creation_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        self.failure_reason = None
        self.settings = settings or {
            "ChannelIdentification": False,
            "ShowAlternatives": False,
        }
        self.specialty = specialty
        self.type = job_type
        self._output_bucket_name = output_bucket_name
        self._output_encryption_kms_key_id = output_encryption_kms_key_id
        self.output_location_type = "CUSTOMER_BUCKET"

    def response_object(self, response_type: str) -> Dict[str, Any]:  # type: ignore
        response_field_dict = {
            "CREATE": [
                "MedicalTranscriptionJobName",
                "TranscriptionJobStatus",
                "LanguageCode",
                "MediaFormat",
                "Media",
                "StartTime",
                "CreationTime",
                "Specialty",
                "Type",
            ],
            "GET": [
                "MedicalTranscriptionJobName",
                "TranscriptionJobStatus",
                "LanguageCode",
                "MediaSampleRateHertz",
                "MediaFormat",
                "Media",
                "Transcript",
                "StartTime",
                "CreationTime",
                "CompletionTime",
                "Settings",
                "Specialty",
                "Type",
            ],
            "LIST": [
                "MedicalTranscriptionJobName",
                "CreationTime",
                "StartTime",
                "CompletionTime",
                "LanguageCode",
                "TranscriptionJobStatus",
                "FailureReason",
                "OutputLocationType",
                "Specialty",
                "Type",
            ],
        }
        response_fields = response_field_dict[response_type]
        response_object = self.gen_response_object()
        response_object["TranscriptionJobStatus"] = self.status
        if response_type != "LIST":
            return {
                "MedicalTranscriptionJob": {
                    k: v
                    for k, v in response_object.items()
                    if k in response_fields and v is not None and v != [None]
                }
            }
        else:
            return {
                k: v
                for k, v in response_object.items()
                if k in response_fields and v is not None and v != [None]
            }

    def advance(self) -> None:
        old_status = self.status
        super().advance()
        new_status = self.status

        if old_status == new_status:
            return

        if new_status == "IN_PROGRESS":
            self.start_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            if not self.media_sample_rate_hertz:
                self.media_sample_rate_hertz = 44100
            if not self.media_format:
                file_ext = self.media["MediaFileUri"].split(".")[-1].lower()
                self.media_format = (
                    file_ext if file_ext in ["mp3", "mp4", "wav", "flac"] else "mp3"
                )
        elif new_status == "COMPLETED":
            self.completion_time = (datetime.now() + timedelta(seconds=10)).strftime(
                "%Y-%m-%d %H:%M:%S"
            )
            self.transcript = {
                "TranscriptFileUri": f"https://s3.{self._region_name}.amazonaws.com/{self._output_bucket_name}/medical/{self.medical_transcription_job_name}.json"
            }


class FakeMedicalVocabulary(FakeVocabulary):
    def __init__(
        self,
        account_id: str,
        region_name: str,
        vocabulary_name: str,
        language_code: str,
        vocabulary_file_uri: Optional[str],
    ):
        super().__init__(
            account_id,
            region_name,
            vocabulary_name,
            language_code=language_code,
            phrases=None,
            vocabulary_file_uri=vocabulary_file_uri,
        )
        self.model_name = "transcribe::medicalvocabulary"
        self._region_name = region_name
        self.vocabulary_name = vocabulary_name
        self.language_code = language_code
        self.vocabulary_file_uri = vocabulary_file_uri
        self.last_modified_time = None
        self.failure_reason = None
        self.download_uri = f"https://s3.us-east-1.amazonaws.com/aws-transcribe-dictionary-model-{region_name}-prod/{account_id}/medical/{self.vocabulary_name}/{mock_random.uuid4()}/input.txt"


class TranscribeBackend(BaseBackend):
    def __init__(self, region_name: str, account_id: str):
        super().__init__(region_name, account_id)
        self.medical_transcriptions: Dict[str, FakeMedicalTranscriptionJob] = {}
        self.transcriptions: Dict[str, FakeTranscriptionJob] = {}
        self.medical_vocabularies: Dict[str, FakeMedicalVocabulary] = {}
        self.vocabularies: Dict[str, FakeVocabulary] = {}

    def start_transcription_job(
        self,
        transcription_job_name: str,
        language_code: Optional[str],
        media_sample_rate_hertz: Optional[int],
        media_format: Optional[str],
        media: Dict[str, str],
        output_bucket_name: Optional[str],
        output_key: Optional[str],
        output_encryption_kms_key_id: Optional[str],
        settings: Optional[Dict[str, Any]],
        model_settings: Optional[Dict[str, Optional[str]]],
        job_execution_settings: Optional[Dict[str, Any]],
        content_redaction: Optional[Dict[str, Any]],
        identify_language: Optional[bool],
        identify_multiple_languages: Optional[bool],
        language_options: Optional[List[str]],
        subtitles: Optional[Dict[str, Any]],
    ) -> Dict[str, Any]:
        if transcription_job_name in self.transcriptions:
            raise ConflictException(
                message="The requested job name already exists. Use a different job name."
            )

        vocabulary_name = settings.get("VocabularyName") if settings else None
        if vocabulary_name and vocabulary_name not in self.vocabularies:
            raise BadRequestException(
                message="The requested vocabulary couldn't be found. "
                "Check the vocabulary name and try your request again."
            )

        transcription_job_object = FakeTranscriptionJob(
            account_id=self.account_id,
            region_name=self.region_name,
            transcription_job_name=transcription_job_name,
            language_code=language_code,
            media_sample_rate_hertz=media_sample_rate_hertz,
            media_format=media_format,
            media=media,
            output_bucket_name=output_bucket_name,
            output_key=output_key,
            output_encryption_kms_key_id=output_encryption_kms_key_id,
            settings=settings,
            model_settings=model_settings,
            job_execution_settings=job_execution_settings,
            content_redaction=content_redaction,
            identify_language=identify_language,
            identify_multiple_languages=identify_multiple_languages,
            language_options=language_options,
            subtitles=subtitles,
        )
        self.transcriptions[transcription_job_name] = transcription_job_object

        return transcription_job_object.response_object("CREATE")

    def start_medical_transcription_job(
        self,
        medical_transcription_job_name: str,
        language_code: str,
        media_sample_rate_hertz: Optional[int],
        media_format: Optional[str],
        media: Dict[str, str],
        output_bucket_name: str,
        output_encryption_kms_key_id: Optional[str],
        settings: Optional[Dict[str, Any]],
        specialty: str,
        type_: str,
    ) -> Dict[str, Any]:
        if medical_transcription_job_name in self.medical_transcriptions:
            raise ConflictException(
                message="The requested job name already exists. Use a different job name."
            )

        vocabulary_name = settings.get("VocabularyName") if settings else None
        if vocabulary_name and vocabulary_name not in self.medical_vocabularies:
            raise BadRequestException(
                message="The requested vocabulary couldn't be found. "
                "Check the vocabulary name and try your request again."
            )

        transcription_job_object = FakeMedicalTranscriptionJob(
            region_name=self.region_name,
            medical_transcription_job_name=medical_transcription_job_name,
            language_code=language_code,
            media_sample_rate_hertz=media_sample_rate_hertz,
            media_format=media_format,
            media=media,
            output_bucket_name=output_bucket_name,
            output_encryption_kms_key_id=output_encryption_kms_key_id,
            settings=settings,
            specialty=specialty,
            job_type=type_,
        )

        self.medical_transcriptions[medical_transcription_job_name] = (
            transcription_job_object
        )

        return transcription_job_object.response_object("CREATE")

    def get_transcription_job(self, transcription_job_name: str) -> Dict[str, Any]:
        try:
            job = self.transcriptions[transcription_job_name]
            job.advance()  # Fakes advancement through statuses.
            return job.response_object("GET")
        except KeyError:
            raise BadRequestException(
                message="The requested job couldn't be found. "
                "Check the job name and try your request again."
            )

    def get_medical_transcription_job(
        self, medical_transcription_job_name: str
    ) -> Dict[str, Any]:
        try:
            job = self.medical_transcriptions[medical_transcription_job_name]
            job.advance()  # Fakes advancement through statuses.
            return job.response_object("GET")
        except KeyError:
            raise BadRequestException(
                message="The requested job couldn't be found. "
                "Check the job name and try your request again."
            )

    def delete_transcription_job(self, transcription_job_name: str) -> None:
        try:
            del self.transcriptions[transcription_job_name]
        except KeyError:
            raise BadRequestException(
                message="The requested job couldn't be found. "
                "Check the job name and try your request again."
            )

    def delete_medical_transcription_job(
        self, medical_transcription_job_name: str
    ) -> None:
        try:
            del self.medical_transcriptions[medical_transcription_job_name]
        except KeyError:
            raise BadRequestException(
                message="The requested job couldn't be found. "
                "Check the job name and try your request again."
            )

    def list_transcription_jobs(
        self,
        state_equals: str,
        job_name_contains: str,
        next_token: str,
        max_results: int,
    ) -> Dict[str, Any]:
        jobs = list(self.transcriptions.values())

        if state_equals:
            jobs = [job for job in jobs if job.status == state_equals]

        if job_name_contains:
            jobs = [
                job for job in jobs if job_name_contains in job.transcription_job_name
            ]

        start_offset = int(next_token) if next_token else 0
        end_offset = start_offset + (
            max_results if max_results else 100
        )  # Arbitrarily selected...
        jobs_paginated = jobs[start_offset:end_offset]

        response: Dict[str, Any] = {
            "TranscriptionJobSummaries": [
                job.response_object("LIST") for job in jobs_paginated
            ]
        }
        if end_offset < len(jobs):
            response["NextToken"] = str(end_offset)
        if state_equals:
            response["Status"] = state_equals
        return response

    def list_medical_transcription_jobs(
        self, status: str, job_name_contains: str, next_token: str, max_results: int
    ) -> Dict[str, Any]:
        jobs = list(self.medical_transcriptions.values())

        if status:
            jobs = [job for job in jobs if job.status == status]

        if job_name_contains:
            jobs = [
                job
                for job in jobs
                if job_name_contains in job.medical_transcription_job_name
            ]

        start_offset = int(next_token) if next_token else 0
        end_offset = start_offset + (
            max_results if max_results else 100
        )  # Arbitrarily selected...
        jobs_paginated = jobs[start_offset:end_offset]

        response: Dict[str, Any] = {
            "MedicalTranscriptionJobSummaries": [
                job.response_object("LIST") for job in jobs_paginated
            ]
        }
        if end_offset < len(jobs):
            response["NextToken"] = str(end_offset)
        if status:
            response["Status"] = status
        return response

    def create_vocabulary(
        self,
        vocabulary_name: str,
        language_code: str,
        phrases: Optional[List[str]],
        vocabulary_file_uri: Optional[str],
    ) -> Dict[str, Any]:
        if (
            phrases is not None
            and vocabulary_file_uri is not None
            or phrases is None
            and vocabulary_file_uri is None
        ):
            raise BadRequestException(
                message="Either Phrases or VocabularyFileUri field should be provided."
            )
        if phrases is not None and len(phrases) < 1:
            raise BadRequestException(
                message="1 validation error detected: Value '[]' at 'phrases' failed to "
                "satisfy constraint: Member must have length greater than or "
                "equal to 1"
            )
        if vocabulary_name in self.vocabularies:
            raise ConflictException(
                message="The requested vocabulary name already exists. "
                "Use a different vocabulary name."
            )

        vocabulary_object = FakeVocabulary(
            account_id=self.account_id,
            region_name=self.region_name,
            vocabulary_name=vocabulary_name,
            language_code=language_code,
            phrases=phrases,
            vocabulary_file_uri=vocabulary_file_uri,
        )

        self.vocabularies[vocabulary_name] = vocabulary_object

        return vocabulary_object.response_object("CREATE")

    def create_medical_vocabulary(
        self,
        vocabulary_name: str,
        language_code: str,
        vocabulary_file_uri: Optional[str],
    ) -> Dict[str, Any]:
        if vocabulary_name in self.medical_vocabularies:
            raise ConflictException(
                message="The requested vocabulary name already exists. "
                "Use a different vocabulary name."
            )

        medical_vocabulary_object = FakeMedicalVocabulary(
            account_id=self.account_id,
            region_name=self.region_name,
            vocabulary_name=vocabulary_name,
            language_code=language_code,
            vocabulary_file_uri=vocabulary_file_uri,
        )

        self.medical_vocabularies[vocabulary_name] = medical_vocabulary_object

        return medical_vocabulary_object.response_object("CREATE")

    def get_vocabulary(self, vocabulary_name: str) -> Dict[str, Any]:
        try:
            job = self.vocabularies[vocabulary_name]
            job.advance()  # Fakes advancement through statuses.
            return job.response_object("GET")
        except KeyError:
            raise BadRequestException(
                message="The requested vocabulary couldn't be found. "
                "Check the vocabulary name and try your request again."
            )

    def get_medical_vocabulary(self, vocabulary_name: str) -> Dict[str, Any]:
        try:
            job = self.medical_vocabularies[vocabulary_name]
            job.advance()  # Fakes advancement through statuses.
            return job.response_object("GET")
        except KeyError:
            raise BadRequestException(
                message="The requested vocabulary couldn't be found. "
                "Check the vocabulary name and try your request again."
            )

    def delete_vocabulary(self, vocabulary_name: str) -> None:
        try:
            del self.vocabularies[vocabulary_name]
        except KeyError:
            raise BadRequestException(
                message="The requested vocabulary couldn't be found. Check the vocabulary name and try your request again."
            )

    def delete_medical_vocabulary(self, vocabulary_name: str) -> None:
        try:
            del self.medical_vocabularies[vocabulary_name]
        except KeyError:
            raise BadRequestException(
                message="The requested vocabulary couldn't be found. Check the vocabulary name and try your request again."
            )

    def list_vocabularies(
        self, state_equals: str, name_contains: str, next_token: str, max_results: int
    ) -> Dict[str, Any]:
        vocabularies = list(self.vocabularies.values())

        if state_equals:
            vocabularies = [
                vocabulary
                for vocabulary in vocabularies
                if vocabulary.status == state_equals
            ]

        if name_contains:
            vocabularies = [
                vocabulary
                for vocabulary in vocabularies
                if name_contains in vocabulary.vocabulary_name
            ]

        start_offset = int(next_token) if next_token else 0
        end_offset = start_offset + (
            max_results if max_results else 100
        )  # Arbitrarily selected...
        vocabularies_paginated = vocabularies[start_offset:end_offset]

        response: Dict[str, Any] = {
            "Vocabularies": [
                vocabulary.response_object("LIST")
                for vocabulary in vocabularies_paginated
            ]
        }
        if end_offset < len(vocabularies):
            response["NextToken"] = str(end_offset)
        if state_equals:
            response["Status"] = state_equals
        return response

    def list_medical_vocabularies(
        self, state_equals: str, name_contains: str, next_token: str, max_results: int
    ) -> Dict[str, Any]:
        vocabularies = list(self.medical_vocabularies.values())

        if state_equals:
            vocabularies = [
                vocabulary
                for vocabulary in vocabularies
                if vocabulary.status == state_equals
            ]

        if name_contains:
            vocabularies = [
                vocabulary
                for vocabulary in vocabularies
                if name_contains in vocabulary.vocabulary_name
            ]

        start_offset = int(next_token) if next_token else 0
        end_offset = start_offset + (
            max_results if max_results else 100
        )  # Arbitrarily selected...
        vocabularies_paginated = vocabularies[start_offset:end_offset]

        response: Dict[str, Any] = {
            "Vocabularies": [
                vocabulary.response_object("LIST")
                for vocabulary in vocabularies_paginated
            ]
        }
        if end_offset < len(vocabularies):
            response["NextToken"] = str(end_offset)
        if state_equals:
            response["Status"] = state_equals
        return response


transcribe_backends = BackendDict(TranscribeBackend, "transcribe")
