import json

from moto.core.responses import BaseResponse

from .models import TranscribeBackend, transcribe_backends


class TranscribeResponse(BaseResponse):
    def __init__(self) -> None:
        super().__init__(service_name="transcribe")

    @property
    def transcribe_backend(self) -> TranscribeBackend:
        return transcribe_backends[self.current_account][self.region]

    def start_transcription_job(self) -> str:
        name = self._get_param("TranscriptionJobName")
        response = self.transcribe_backend.start_transcription_job(
            transcription_job_name=name,
            language_code=self._get_param("LanguageCode"),
            media_sample_rate_hertz=self._get_param("MediaSampleRateHertz"),
            media_format=self._get_param("MediaFormat"),
            media=self._get_param("Media"),
            output_bucket_name=self._get_param("OutputBucketName"),
            output_key=self._get_param("OutputKey"),
            output_encryption_kms_key_id=self._get_param("OutputEncryptionKMSKeyId"),
            settings=self._get_param("Settings"),
            model_settings=self._get_param("ModelSettings"),
            job_execution_settings=self._get_param("JobExecutionSettings"),
            content_redaction=self._get_param("ContentRedaction"),
            identify_language=self._get_param("IdentifyLanguage"),
            identify_multiple_languages=self._get_param("IdentifyMultipleLanguages"),
            language_options=self._get_param("LanguageOptions"),
            subtitles=self._get_param("Subtitles"),
        )
        return json.dumps(response)

    def start_medical_transcription_job(self) -> str:
        name = self._get_param("MedicalTranscriptionJobName")
        response = self.transcribe_backend.start_medical_transcription_job(
            medical_transcription_job_name=name,
            language_code=self._get_param("LanguageCode"),
            media_sample_rate_hertz=self._get_param("MediaSampleRateHertz"),
            media_format=self._get_param("MediaFormat"),
            media=self._get_param("Media"),
            output_bucket_name=self._get_param("OutputBucketName"),
            output_encryption_kms_key_id=self._get_param("OutputEncryptionKMSKeyId"),
            settings=self._get_param("Settings"),
            specialty=self._get_param("Specialty"),
            type_=self._get_param("Type"),
        )
        return json.dumps(response)

    def list_transcription_jobs(self) -> str:
        state_equals = self._get_param("Status")
        job_name_contains = self._get_param("JobNameContains")
        next_token = self._get_param("NextToken")
        max_results = self._get_param("MaxResults")

        response = self.transcribe_backend.list_transcription_jobs(
            state_equals=state_equals,
            job_name_contains=job_name_contains,
            next_token=next_token,
            max_results=max_results,
        )
        return json.dumps(response)

    def list_medical_transcription_jobs(self) -> str:
        status = self._get_param("Status")
        job_name_contains = self._get_param("JobNameContains")
        next_token = self._get_param("NextToken")
        max_results = self._get_param("MaxResults")

        response = self.transcribe_backend.list_medical_transcription_jobs(
            status=status,
            job_name_contains=job_name_contains,
            next_token=next_token,
            max_results=max_results,
        )
        return json.dumps(response)

    def get_transcription_job(self) -> str:
        transcription_job_name = self._get_param("TranscriptionJobName")
        response = self.transcribe_backend.get_transcription_job(
            transcription_job_name=transcription_job_name
        )
        return json.dumps(response)

    def get_medical_transcription_job(self) -> str:
        medical_transcription_job_name = self._get_param("MedicalTranscriptionJobName")
        response = self.transcribe_backend.get_medical_transcription_job(
            medical_transcription_job_name=medical_transcription_job_name
        )
        return json.dumps(response)

    def delete_transcription_job(self) -> str:
        transcription_job_name = self._get_param("TranscriptionJobName")
        self.transcribe_backend.delete_transcription_job(
            transcription_job_name=transcription_job_name
        )
        return "{}"

    def delete_medical_transcription_job(self) -> str:
        medical_transcription_job_name = self._get_param("MedicalTranscriptionJobName")
        self.transcribe_backend.delete_medical_transcription_job(
            medical_transcription_job_name=medical_transcription_job_name
        )
        return "{}"

    def create_vocabulary(self) -> str:
        vocabulary_name = self._get_param("VocabularyName")
        language_code = self._get_param("LanguageCode")
        phrases = self._get_param("Phrases")
        vocabulary_file_uri = self._get_param("VocabularyFileUri")
        response = self.transcribe_backend.create_vocabulary(
            vocabulary_name=vocabulary_name,
            language_code=language_code,
            phrases=phrases,
            vocabulary_file_uri=vocabulary_file_uri,
        )
        return json.dumps(response)

    def create_medical_vocabulary(self) -> str:
        vocabulary_name = self._get_param("VocabularyName")
        language_code = self._get_param("LanguageCode")
        vocabulary_file_uri = self._get_param("VocabularyFileUri")
        response = self.transcribe_backend.create_medical_vocabulary(
            vocabulary_name=vocabulary_name,
            language_code=language_code,
            vocabulary_file_uri=vocabulary_file_uri,
        )
        return json.dumps(response)

    def get_vocabulary(self) -> str:
        vocabulary_name = self._get_param("VocabularyName")
        response = self.transcribe_backend.get_vocabulary(
            vocabulary_name=vocabulary_name
        )
        return json.dumps(response)

    def get_medical_vocabulary(self) -> str:
        vocabulary_name = self._get_param("VocabularyName")
        response = self.transcribe_backend.get_medical_vocabulary(
            vocabulary_name=vocabulary_name
        )
        return json.dumps(response)

    def list_vocabularies(self) -> str:
        state_equals = self._get_param("StateEquals")
        name_contains = self._get_param("NameContains")
        next_token = self._get_param("NextToken")
        max_results = self._get_param("MaxResults")

        response = self.transcribe_backend.list_vocabularies(
            state_equals=state_equals,
            name_contains=name_contains,
            next_token=next_token,
            max_results=max_results,
        )
        return json.dumps(response)

    def list_medical_vocabularies(self) -> str:
        state_equals = self._get_param("StateEquals")
        name_contains = self._get_param("NameContains")
        next_token = self._get_param("NextToken")
        max_results = self._get_param("MaxResults")

        response = self.transcribe_backend.list_medical_vocabularies(
            state_equals=state_equals,
            name_contains=name_contains,
            next_token=next_token,
            max_results=max_results,
        )
        return json.dumps(response)

    def delete_vocabulary(self) -> str:
        vocabulary_name = self._get_param("VocabularyName")
        self.transcribe_backend.delete_vocabulary(vocabulary_name=vocabulary_name)
        return "{}"

    def delete_medical_vocabulary(self) -> str:
        vocabulary_name = self._get_param("VocabularyName")
        self.transcribe_backend.delete_medical_vocabulary(
            vocabulary_name=vocabulary_name
        )
        return "{}"
