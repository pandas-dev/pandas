"""Handles incoming comprehend requests, invokes methods, returns responses."""

import json
from typing import Any, Dict, List

from moto.core.responses import BaseResponse

from .models import ComprehendBackend, comprehend_backends


class ComprehendResponse(BaseResponse):
    """Handler for Comprehend requests and responses."""

    def __init__(self) -> None:
        super().__init__(service_name="comprehend")

    @property
    def comprehend_backend(self) -> ComprehendBackend:
        """Return backend instance specific for this region."""
        return comprehend_backends[self.current_account][self.region]

    def list_entity_recognizers(self) -> str:
        params = json.loads(self.body)
        _filter = params.get("Filter", {})
        recognizers = self.comprehend_backend.list_entity_recognizers(_filter=_filter)
        return json.dumps(
            dict(EntityRecognizerPropertiesList=[r.to_dict() for r in recognizers])
        )

    def create_entity_recognizer(self) -> str:
        params = json.loads(self.body)
        recognizer_name = params.get("RecognizerName")
        version_name = params.get("VersionName")
        data_access_role_arn = params.get("DataAccessRoleArn")
        tags = params.get("Tags")
        input_data_config = params.get("InputDataConfig")
        language_code = params.get("LanguageCode")
        volume_kms_key_id = params.get("VolumeKmsKeyId")
        vpc_config = params.get("VpcConfig")
        model_kms_key_id = params.get("ModelKmsKeyId")
        model_policy = params.get("ModelPolicy")
        entity_recognizer_arn = self.comprehend_backend.create_entity_recognizer(
            recognizer_name=recognizer_name,
            version_name=version_name,
            data_access_role_arn=data_access_role_arn,
            tags=tags,
            input_data_config=input_data_config,
            language_code=language_code,
            volume_kms_key_id=volume_kms_key_id,
            vpc_config=vpc_config,
            model_kms_key_id=model_kms_key_id,
            model_policy=model_policy,
        )
        return json.dumps(dict(EntityRecognizerArn=entity_recognizer_arn))

    def describe_entity_recognizer(self) -> str:
        params = json.loads(self.body)
        entity_recognizer_arn = params.get("EntityRecognizerArn")
        recognizer = self.comprehend_backend.describe_entity_recognizer(
            entity_recognizer_arn=entity_recognizer_arn,
        )
        return json.dumps(dict(EntityRecognizerProperties=recognizer.to_dict()))

    def stop_training_entity_recognizer(self) -> str:
        params = json.loads(self.body)
        entity_recognizer_arn = params.get("EntityRecognizerArn")
        self.comprehend_backend.stop_training_entity_recognizer(
            entity_recognizer_arn=entity_recognizer_arn,
        )
        return json.dumps(dict())

    def list_tags_for_resource(self) -> str:
        params = json.loads(self.body)
        resource_arn = params.get("ResourceArn")
        tags = self.comprehend_backend.list_tags_for_resource(
            resource_arn=resource_arn,
        )
        return json.dumps(dict(ResourceArn=resource_arn, Tags=tags))

    def delete_entity_recognizer(self) -> str:
        params = json.loads(self.body)
        entity_recognizer_arn = params.get("EntityRecognizerArn")
        self.comprehend_backend.delete_entity_recognizer(
            entity_recognizer_arn=entity_recognizer_arn,
        )
        return "{}"

    def tag_resource(self) -> str:
        params = json.loads(self.body)
        resource_arn = params.get("ResourceArn")
        tags = params.get("Tags")
        self.comprehend_backend.tag_resource(resource_arn, tags)
        return "{}"

    def untag_resource(self) -> str:
        params = json.loads(self.body)
        resource_arn = params.get("ResourceArn")
        tag_keys = params.get("TagKeys")
        self.comprehend_backend.untag_resource(resource_arn, tag_keys)
        return "{}"

    def detect_pii_entities(self) -> str:
        params = json.loads(self.body)
        text = params.get("Text")
        language = params.get("LanguageCode")
        resp = self.comprehend_backend.detect_pii_entities(text, language)
        return json.dumps(dict(Entities=resp))

    def detect_key_phrases(self) -> str:
        params = json.loads(self.body)
        text = params.get("Text")
        language = params.get("LanguageCode")
        resp = self.comprehend_backend.detect_key_phrases(text, language)
        return json.dumps(dict(KeyPhrases=resp))

    def detect_sentiment(self) -> str:
        params = json.loads(self.body)
        text = params.get("Text")
        language = params.get("LanguageCode")
        resp = self.comprehend_backend.detect_sentiment(text, language)
        return json.dumps(resp)

    def create_document_classifier(self) -> str:
        params = json.loads(self.body)
        document_classifier_name = params.get("DocumentClassifierName")
        version_name = params.get("VersionName")
        data_access_role_arn = params.get("DataAccessRoleArn")
        tags = params.get("Tags")
        input_data_config = params.get("InputDataConfig")
        output_data_config = params.get("OutputDataConfig")
        client_request_token = params.get("ClientRequestToken")
        language_code = params.get("LanguageCode")
        volume_kms_key_id = params.get("VolumeKmsKeyId")
        vpc_config = params.get("VpcConfig")
        mode = params.get("Mode")
        model_kms_key_id = params.get("ModelKmsKeyId")
        model_policy = params.get("ModelPolicy")

        document_classifier_arn = self.comprehend_backend.create_document_classifier(
            document_classifier_name=document_classifier_name,
            version_name=version_name,
            data_access_role_arn=data_access_role_arn,
            tags=tags,
            input_data_config=input_data_config,
            output_data_config=output_data_config,
            client_request_token=client_request_token,
            language_code=language_code,
            volume_kms_key_id=volume_kms_key_id,
            vpc_config=vpc_config,
            mode=mode,
            model_kms_key_id=model_kms_key_id,
            model_policy=model_policy,
        )

        return json.dumps(dict(DocumentClassifierArn=document_classifier_arn))

    def create_endpoint(self) -> str:
        params = json.loads(self.body)
        endpoint_name = params.get("EndpointName")
        model_arn = params.get("ModelArn")
        desired_inference_units = params.get("DesiredInferenceUnits")
        client_request_token = params.get("ClientRequestToken")
        tags = params.get("Tags")
        data_access_role_arn = params.get("DataAccessRoleArn")
        flywheel_arn = params.get("FlywheelArn")
        endpoint_arn, model_arn = self.comprehend_backend.create_endpoint(
            endpoint_name=endpoint_name,
            model_arn=model_arn,
            desired_inference_units=desired_inference_units,
            client_request_token=client_request_token,
            tags=tags,
            data_access_role_arn=data_access_role_arn,
            flywheel_arn=flywheel_arn,
        )

        return json.dumps(dict(EndpointArn=endpoint_arn, ModelArn=model_arn))

    def create_flywheel(self) -> str:
        params = json.loads(self.body)
        flywheel_name = params.get("FlywheelName")
        active_model_arn = params.get("ActiveModelArn")
        data_access_role_arn = params.get("DataAccessRoleArn")
        task_config = params.get("TaskConfig")
        model_type = params.get("ModelType")
        data_lake_s3_uri = params.get("DataLakeS3Uri")
        data_security_config = params.get("DataSecurityConfig")
        client_request_token = params.get("ClientRequestToken")
        tags = params.get("Tags")
        flywheel_arn, active_model_arn = self.comprehend_backend.create_flywheel(
            flywheel_name=flywheel_name,
            active_model_arn=active_model_arn,
            data_access_role_arn=data_access_role_arn,
            task_config=task_config,
            model_type=model_type,
            data_lake_s3_uri=data_lake_s3_uri,
            data_security_config=data_security_config,
            client_request_token=client_request_token,
            tags=tags,
        )

        return json.dumps(
            dict(FlywheelArn=flywheel_arn, activeModelArn=active_model_arn)
        )

    def describe_document_classifier(self) -> str:
        params = json.loads(self.body)
        document_classifier_arn = params.get("DocumentClassifierArn")
        document_classifier = self.comprehend_backend.describe_document_classifier(
            document_classifier_arn=document_classifier_arn,
        )

        return json.dumps(
            dict(DocumentClassifierProperties=document_classifier.to_dict())
        )

    def describe_endpoint(self) -> str:
        params = json.loads(self.body)
        endpoint_arn = params.get("EndpointArn")
        endpoint_properties = self.comprehend_backend.describe_endpoint(
            endpoint_arn=endpoint_arn,
        )

        return json.dumps(dict(EndpointProperties=endpoint_properties.to_dict()))

    def describe_flywheel(self) -> str:
        params = json.loads(self.body)
        flywheel_arn = params.get("FlywheelArn")
        flywheel_properties = self.comprehend_backend.describe_flywheel(
            flywheel_arn=flywheel_arn,
        )

        return json.dumps(dict(FlywheelProperties=flywheel_properties.to_dict()))

    def delete_document_classifier(self) -> str:
        params = json.loads(self.body)
        document_classifier_arn = params.get("DocumentClassifierArn")
        self.comprehend_backend.delete_document_classifier(
            document_classifier_arn=document_classifier_arn,
        )

        return json.dumps(dict())

    def delete_endpoint(self) -> str:
        params = json.loads(self.body)
        endpoint_arn = params.get("EndpointArn")
        self.comprehend_backend.delete_endpoint(
            endpoint_arn=endpoint_arn,
        )

        return json.dumps(dict())

    def delete_flywheel(self) -> str:
        params = json.loads(self.body)
        flywheel_arn = params.get("FlywheelArn")
        self.comprehend_backend.delete_flywheel(
            flywheel_arn=flywheel_arn,
        )

        return json.dumps(dict())

    def list_document_classifiers(self) -> str:
        params = json.loads(self.body)
        filter = params.get("Filter")
        next_token = params.get("NextToken")
        max_results = params.get("MaxResults")
        document_classifier_properties_list, next_token = (
            self.comprehend_backend.list_document_classifiers(
                filter=filter,
                next_token=next_token,
                max_results=max_results,
            )
        )

        return json.dumps(
            dict(
                DocumentClassifierPropertiesList=document_classifier_properties_list,
                NextToken=next_token,
            )
        )

    def list_endpoints(self) -> str:
        params = json.loads(self.body)
        filter = params.get("Filter")
        next_token = params.get("NextToken")
        max_results = params.get("MaxResults")
        endpoint_properties_list, next_token = self.comprehend_backend.list_endpoints(
            filter=filter,
            next_token=next_token,
            max_results=max_results,
        )

        return json.dumps(
            dict(EndpointPropertiesList=endpoint_properties_list, NextToken=next_token)
        )

    def list_flywheels(self) -> str:
        params = json.loads(self.body)
        filter = params.get("Filter")
        next_token = params.get("NextToken")
        max_results = params.get("MaxResults")
        flywheel_summary_list, next_token = self.comprehend_backend.list_flywheels(
            filter=filter,
            next_token=next_token,
            max_results=max_results,
        )

        return json.dumps(
            dict(FlywheelSummaryList=flywheel_summary_list, nextToken=next_token)
        )

    def stop_training_document_classifier(self) -> str:
        params = json.loads(self.body)
        document_classifier_arn = params.get("DocumentClassifierArn")
        self.comprehend_backend.stop_training_document_classifier(
            document_classifier_arn=document_classifier_arn,
        )

        return json.dumps(dict())

    def start_flywheel_iteration(self) -> str:
        params = json.loads(self.body)
        flywheel_arn = params.get("FlywheelArn")
        client_request_token = params.get("ClientRequestToken")
        flywheel_arn, flywheel_iteration_id = (
            self.comprehend_backend.start_flywheel_iteration(
                flywheel_arn=flywheel_arn,
                client_request_token=client_request_token,
            )
        )

        return json.dumps(
            dict(FlywheelArn=flywheel_arn, FlywheelIterationId=flywheel_iteration_id)
        )

    def update_endpoint(self) -> str:
        params = json.loads(self.body)
        endpoint_arn = params.get("EndpointArn")
        desired_model_arn = params.get("DesiredModelArn")
        desired_inference_units = params.get("DesiredInferenceUnits")
        desired_data_access_role_arn = params.get("DesiredDataAccessRoleArn")
        flywheel_arn = params.get("FlywheelArn")
        desired_model_arn = self.comprehend_backend.update_endpoint(
            endpoint_arn=endpoint_arn,
            desired_model_arn=desired_model_arn,
            desired_inference_units=desired_inference_units,
            desired_data_access_role_arn=desired_data_access_role_arn,
            flywheel_arn=flywheel_arn,
        )

        return json.dumps(dict(DesiredModelArn=desired_model_arn))

    def _job_to_dict_resp(self, job_properties: Dict[str, Any]) -> str:
        job_type_key = job_properties.pop("job_type")

        if job_properties.get("SubmitTime"):
            job_properties["SubmitTime"] = job_properties["SubmitTime"].isoformat()
        if job_properties.get("EndTime"):
            job_properties["EndTime"] = job_properties["EndTime"].isoformat()

        key_name = f"{job_type_key}JobProperties"
        return json.dumps({key_name: job_properties})

    def _list_jobs_to_dict_resp(
        self, job_list: List[Dict[str, Any]], job_type: str
    ) -> str:
        for job_properties in job_list:
            job_properties.pop("job_type")
            if job_properties.get("SubmitTime"):
                job_properties["SubmitTime"] = job_properties["SubmitTime"].isoformat()
            if job_properties.get("EndTime"):
                job_properties["EndTime"] = job_properties["EndTime"].isoformat()

        key_name = f"{job_type}JobPropertiesList"
        return json.dumps({key_name: job_list, "NextToken": None})

    def start_pii_entities_detection_job(self) -> str:
        params = json.loads(self.body)
        job = self.comprehend_backend.start_pii_entities_detection_job(**params)
        return json.dumps(
            {"JobId": job.job_id, "JobArn": job.job_arn, "JobStatus": job.job_status}
        )

    def describe_pii_entities_detection_job(self) -> str:
        params = json.loads(self.body)
        job = self.comprehend_backend.describe_pii_entities_detection_job(
            job_id=params["JobId"]
        )
        return self._job_to_dict_resp(job.to_dict())

    def stop_pii_entities_detection_job(self) -> str:
        params = json.loads(self.body)
        self.comprehend_backend.stop_pii_entities_detection_job(job_id=params["JobId"])
        return json.dumps({"JobId": params["JobId"], "JobStatus": "STOP_REQUESTED"})

    def list_pii_entities_detection_jobs(self) -> str:
        params = json.loads(self.body)
        job_filter = params.get("Filter")
        jobs = self.comprehend_backend.list_pii_entities_detection_jobs(
            filter=job_filter
        )
        job_list = [job.to_dict() for job in jobs]
        return self._list_jobs_to_dict_resp(job_list, "PiiEntitiesDetection")

    def start_key_phrases_detection_job(self) -> str:
        params = json.loads(self.body)
        job = self.comprehend_backend.start_key_phrases_detection_job(**params)
        return json.dumps(
            {"JobId": job.job_id, "JobArn": job.job_arn, "JobStatus": job.job_status}
        )

    def describe_key_phrases_detection_job(self) -> str:
        params = json.loads(self.body)
        job = self.comprehend_backend.describe_key_phrases_detection_job(
            job_id=params["JobId"]
        )
        return self._job_to_dict_resp(job.to_dict())

    def stop_key_phrases_detection_job(self) -> str:
        params = json.loads(self.body)
        self.comprehend_backend.stop_key_phrases_detection_job(job_id=params["JobId"])
        return json.dumps({"JobId": params["JobId"], "JobStatus": "STOP_REQUESTED"})

    def list_key_phrases_detection_jobs(self) -> str:
        params = json.loads(self.body)
        job_filter = params.get("Filter")
        jobs = self.comprehend_backend.list_key_phrases_detection_jobs(
            filter=job_filter
        )
        job_list = [job.to_dict() for job in jobs]
        return self._list_jobs_to_dict_resp(job_list, "KeyPhrasesDetection")

    def start_sentiment_detection_job(self) -> str:
        params = json.loads(self.body)
        job = self.comprehend_backend.start_sentiment_detection_job(**params)
        return json.dumps(
            {"JobId": job.job_id, "JobArn": job.job_arn, "JobStatus": job.job_status}
        )

    def describe_sentiment_detection_job(self) -> str:
        params = json.loads(self.body)
        job = self.comprehend_backend.describe_sentiment_detection_job(
            job_id=params["JobId"]
        )
        return self._job_to_dict_resp(job.to_dict())

    def stop_sentiment_detection_job(self) -> str:
        params = json.loads(self.body)
        self.comprehend_backend.stop_sentiment_detection_job(job_id=params["JobId"])
        return json.dumps({"JobId": params["JobId"], "JobStatus": "STOP_REQUESTED"})

    def list_sentiment_detection_jobs(self) -> str:
        params = json.loads(self.body)
        job_filter = params.get("Filter")
        jobs = self.comprehend_backend.list_sentiment_detection_jobs(filter=job_filter)
        job_list = [job.to_dict() for job in jobs]
        return self._list_jobs_to_dict_resp(job_list, "SentimentDetection")

    def put_resource_policy(self) -> str:
        params = json.loads(self.body)
        resource_arn = params.get("ResourceArn")
        resource_policy = params.get("ResourcePolicy")
        policy_revision_id = params.get("PolicyRevisionId")

        revision_id = self.comprehend_backend.put_resource_policy(
            resource_arn=resource_arn,
            resource_policy=resource_policy,
            policy_revision_id=policy_revision_id,
        )

        return json.dumps({"PolicyRevisionId": revision_id})

    def describe_resource_policy(self) -> str:
        params = json.loads(self.body)
        resource_arn = params.get("ResourceArn")

        policy_details = self.comprehend_backend.describe_resource_policy(
            resource_arn=resource_arn
        )

        response_payload = {
            "ResourcePolicy": policy_details["ResourcePolicy"],
            "CreationTime": policy_details["CreationTime"].isoformat(),
            "LastModifiedTime": policy_details["LastModifiedTime"].isoformat(),
            "PolicyRevisionId": policy_details["PolicyRevisionId"],
        }

        return json.dumps(response_payload)

    def delete_resource_policy(self) -> str:
        params = json.loads(self.body)
        resource_arn = params.get("ResourceArn")
        policy_revision_id = params.get("PolicyRevisionId")

        self.comprehend_backend.delete_resource_policy(
            resource_arn=resource_arn,
            policy_revision_id=policy_revision_id,
        )

        return "{}"

    def start_targeted_sentiment_detection_job(self) -> str:
        params = json.loads(self.body)
        job = self.comprehend_backend.start_targeted_sentiment_detection_job(**params)
        return json.dumps(
            {"JobId": job.job_id, "JobArn": job.job_arn, "JobStatus": job.job_status}
        )

    def describe_targeted_sentiment_detection_job(self) -> str:
        params = json.loads(self.body)
        job = self.comprehend_backend.describe_targeted_sentiment_detection_job(
            job_id=params["JobId"]
        )
        return self._job_to_dict_resp(job.to_dict())

    def stop_targeted_sentiment_detection_job(self) -> str:
        params = json.loads(self.body)
        self.comprehend_backend.stop_targeted_sentiment_detection_job(
            job_id=params["JobId"]
        )
        return json.dumps({"JobId": params["JobId"], "JobStatus": "STOP_REQUESTED"})

    def list_targeted_sentiment_detection_jobs(self) -> str:
        params = json.loads(self.body)
        job_filter = params.get("Filter")
        jobs = self.comprehend_backend.list_targeted_sentiment_detection_jobs(
            filter=job_filter
        )
        job_list = [job.to_dict() for job in jobs]
        return self._list_jobs_to_dict_resp(job_list, "TargetedSentimentDetection")

    def start_dominant_language_detection_job(self) -> str:
        params = json.loads(self.body)
        job = self.comprehend_backend.start_dominant_language_detection_job(**params)
        return json.dumps(
            {"JobId": job.job_id, "JobArn": job.job_arn, "JobStatus": job.job_status}
        )

    def describe_dominant_language_detection_job(self) -> str:
        params = json.loads(self.body)
        job = self.comprehend_backend.describe_dominant_language_detection_job(
            job_id=params["JobId"]
        )
        return self._job_to_dict_resp(job.to_dict())

    def stop_dominant_language_detection_job(self) -> str:
        params = json.loads(self.body)
        self.comprehend_backend.stop_dominant_language_detection_job(
            job_id=params["JobId"]
        )
        return json.dumps({"JobId": params["JobId"], "JobStatus": "STOP_REQUESTED"})

    def list_dominant_language_detection_jobs(self) -> str:
        params = json.loads(self.body)
        job_filter = params.get("Filter")
        jobs = self.comprehend_backend.list_dominant_language_detection_jobs(
            filter=job_filter
        )
        job_list = [job.to_dict() for job in jobs]
        return self._list_jobs_to_dict_resp(job_list, "DominantLanguageDetection")

    def start_entities_detection_job(self) -> str:
        params = json.loads(self.body)
        job = self.comprehend_backend.start_entities_detection_job(**params)
        return json.dumps(
            {"JobId": job.job_id, "JobArn": job.job_arn, "JobStatus": job.job_status}
        )

    def describe_entities_detection_job(self) -> str:
        params = json.loads(self.body)
        job = self.comprehend_backend.describe_entities_detection_job(
            job_id=params["JobId"]
        )
        return self._job_to_dict_resp(job.to_dict())

    def stop_entities_detection_job(self) -> str:
        params = json.loads(self.body)
        self.comprehend_backend.stop_entities_detection_job(job_id=params["JobId"])
        return json.dumps({"JobId": params["JobId"], "JobStatus": "STOP_REQUESTED"})

    def list_entities_detection_jobs(self) -> str:
        params = json.loads(self.body)
        job_filter = params.get("Filter")
        jobs = self.comprehend_backend.list_entities_detection_jobs(filter=job_filter)
        job_list = [job.to_dict() for job in jobs]
        return self._list_jobs_to_dict_resp(job_list, "EntitiesDetection")

    def start_topics_detection_job(self) -> str:
        params = json.loads(self.body)
        job = self.comprehend_backend.start_topics_detection_job(**params)
        return json.dumps(
            {"JobId": job.job_id, "JobArn": job.job_arn, "JobStatus": job.job_status}
        )

    def describe_topics_detection_job(self) -> str:
        params = json.loads(self.body)
        job = self.comprehend_backend.describe_topics_detection_job(
            job_id=params["JobId"]
        )
        return self._job_to_dict_resp(job.to_dict())

    def list_topics_detection_jobs(self) -> str:
        params = json.loads(self.body)
        job_filter = params.get("Filter")
        jobs = self.comprehend_backend.list_topics_detection_jobs(filter=job_filter)
        job_list = [job.to_dict() for job in jobs]
        return self._list_jobs_to_dict_resp(job_list, "TopicsDetection")

    def start_document_classification_job(self) -> str:
        params = json.loads(self.body)
        job = self.comprehend_backend.start_document_classification_job(**params)
        return json.dumps(
            {"JobId": job.job_id, "JobArn": job.job_arn, "JobStatus": job.job_status}
        )

    def describe_document_classification_job(self) -> str:
        params = json.loads(self.body)
        job = self.comprehend_backend.describe_document_classification_job(
            job_id=params["JobId"]
        )
        return self._job_to_dict_resp(job.to_dict())

    def list_document_classification_jobs(self) -> str:
        params = json.loads(self.body)
        job_filter = params.get("Filter")
        jobs = self.comprehend_backend.list_document_classification_jobs(
            filter=job_filter
        )
        job_list = [job.to_dict() for job in jobs]
        return self._list_jobs_to_dict_resp(job_list, "DocumentClassification")

    def start_events_detection_job(self) -> str:
        params = json.loads(self.body)
        job = self.comprehend_backend.start_events_detection_job(**params)
        return json.dumps(
            {"JobId": job.job_id, "JobArn": job.job_arn, "JobStatus": job.job_status}
        )

    def describe_events_detection_job(self) -> str:
        params = json.loads(self.body)
        job = self.comprehend_backend.describe_events_detection_job(
            job_id=params["JobId"]
        )
        return self._job_to_dict_resp(job.to_dict())

    def stop_events_detection_job(self) -> str:
        params = json.loads(self.body)
        self.comprehend_backend.stop_events_detection_job(job_id=params["JobId"])
        return json.dumps({"JobId": params["JobId"], "JobStatus": "STOP_REQUESTED"})

    def list_events_detection_jobs(self) -> str:
        params = json.loads(self.body)
        job_filter = params.get("Filter")
        jobs = self.comprehend_backend.list_events_detection_jobs(filter=job_filter)
        job_list = [job.to_dict() for job in jobs]
        return self._list_jobs_to_dict_resp(job_list, "EventsDetection")
