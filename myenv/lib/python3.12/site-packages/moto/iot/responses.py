import json
from typing import Any
from urllib.parse import unquote

from moto.core.common_types import TYPE_RESPONSE
from moto.core.responses import BaseResponse

from .models import IoTBackend, iot_backends


class IoTResponse(BaseResponse):
    def __init__(self) -> None:
        super().__init__(service_name="iot")

    @property
    def iot_backend(self) -> IoTBackend:
        return iot_backends[self.current_account][self.region]

    def create_certificate_from_csr(self) -> str:
        certificate_signing_request = self._get_param("certificateSigningRequest")
        set_as_active = self._get_param("setAsActive")
        cert = self.iot_backend.create_certificate_from_csr(
            certificate_signing_request, set_as_active=set_as_active
        )
        return json.dumps(
            {
                "certificateId": cert.certificate_id,
                "certificateArn": cert.arn,
                "certificatePem": cert.certificate_pem,
            }
        )

    def create_thing(self) -> str:
        thing_name = self._get_param("thingName")
        thing_type_name = self._get_param("thingTypeName")
        attribute_payload = self._get_param("attributePayload")
        thing = self.iot_backend.create_thing(
            thing_name=thing_name,
            thing_type_name=thing_type_name,
            attribute_payload=attribute_payload,
        )
        return json.dumps(
            dict(thingName=thing.thing_name, thingArn=thing.arn, thingId=thing.thing_id)
        )

    def create_thing_type(self) -> str:
        thing_type_name = self._get_param("thingTypeName")
        thing_type_properties = self._get_param("thingTypeProperties")
        thing_type_name, thing_type_arn = self.iot_backend.create_thing_type(
            thing_type_name=thing_type_name, thing_type_properties=thing_type_properties
        )
        return json.dumps(
            dict(thingTypeName=thing_type_name, thingTypeArn=thing_type_arn)
        )

    def list_thing_types(self) -> str:
        previous_next_token = self._get_param("nextToken")
        max_results = self._get_int_param(
            "maxResults", 50
        )  # not the default, but makes testing easier
        thing_type_name = self._get_param("thingTypeName")
        thing_types = self.iot_backend.list_thing_types(thing_type_name=thing_type_name)

        thing_types = [_.to_dict() for _ in thing_types]
        if previous_next_token is None:
            result = thing_types[0:max_results]
            next_token = str(max_results) if len(thing_types) > max_results else None
        else:
            token = int(previous_next_token)
            result = thing_types[token : token + max_results]
            next_token = (
                str(token + max_results)
                if len(thing_types) > token + max_results
                else None
            )

        return json.dumps(dict(thingTypes=result, nextToken=next_token))

    def list_things(self) -> str:
        previous_next_token = self._get_param("nextToken")
        max_results = self._get_int_param(
            "maxResults", 50
        )  # not the default, but makes testing easier
        attribute_name = self._get_param("attributeName")
        attribute_value = self._get_param("attributeValue")
        thing_type_name = self._get_param("thingTypeName")
        things, next_token = self.iot_backend.list_things(
            attribute_name=attribute_name,
            attribute_value=attribute_value,
            thing_type_name=thing_type_name,
            max_results=max_results,
            token=previous_next_token,
        )

        return json.dumps(dict(things=things, nextToken=next_token))

    def describe_thing(self) -> str:
        thing_name = self._get_param("thingName")
        thing = self.iot_backend.describe_thing(thing_name=thing_name)
        return json.dumps(
            thing.to_dict(include_default_client_id=True, include_thing_id=True)
        )

    def describe_thing_type(self) -> str:
        thing_type_name = self._get_param("thingTypeName")
        thing_type = self.iot_backend.describe_thing_type(
            thing_type_name=thing_type_name
        )
        return json.dumps(thing_type.to_dict())

    def describe_endpoint(self) -> str:
        endpoint_type = self._get_param("endpointType", "iot:Data-ATS")
        endpoint = self.iot_backend.describe_endpoint(endpoint_type=endpoint_type)
        return json.dumps(endpoint.to_dict())

    def delete_thing(self) -> str:
        thing_name = self._get_param("thingName")
        self.iot_backend.delete_thing(thing_name=thing_name)
        return json.dumps(dict())

    def delete_thing_type(self) -> str:
        thing_type_name = self._get_param("thingTypeName")
        self.iot_backend.delete_thing_type(thing_type_name=thing_type_name)
        return json.dumps(dict())

    def deprecate_thing_type(self) -> str:
        thing_type_name = self._get_param("thingTypeName")
        undo_deprecate = self._get_param("undoDeprecate")
        thing_type = self.iot_backend.deprecate_thing_type(
            thing_type_name=thing_type_name, undo_deprecate=undo_deprecate
        )
        return json.dumps(thing_type.to_dict())

    def update_thing(self) -> str:
        thing_name = self._get_param("thingName")
        thing_type_name = self._get_param("thingTypeName")
        attribute_payload = self._get_param("attributePayload")
        remove_thing_type = self._get_param("removeThingType")
        self.iot_backend.update_thing(
            thing_name=thing_name,
            thing_type_name=thing_type_name,
            attribute_payload=attribute_payload,
            remove_thing_type=remove_thing_type,
        )
        return json.dumps(dict())

    def create_job(self) -> str:
        job_arn, job_id, description = self.iot_backend.create_job(
            job_id=self._get_param("jobId"),
            targets=self._get_param("targets"),
            description=self._get_param("description"),
            document_source=self._get_param("documentSource"),
            document=self._get_param("document"),
            presigned_url_config=self._get_param("presignedUrlConfig"),
            target_selection=self._get_param("targetSelection"),
            job_executions_rollout_config=self._get_param("jobExecutionsRolloutConfig"),
            document_parameters=self._get_param("documentParameters"),
        )

        return json.dumps(dict(jobArn=job_arn, jobId=job_id, description=description))

    def describe_job(self) -> str:
        job = self.iot_backend.describe_job(job_id=self._get_param("jobId"))
        return json.dumps(
            dict(
                documentSource=job.document_source,
                job=dict(
                    comment=job.comment,
                    completedAt=job.completed_at,
                    createdAt=job.created_at,
                    description=job.description,
                    documentParameters=job.document_parameters,
                    forceCanceled=job.force,
                    reasonCode=job.reason_code,
                    jobArn=job.job_arn,
                    jobExecutionsRolloutConfig=job.job_executions_rollout_config,
                    jobId=job.job_id,
                    jobProcessDetails=job.job_process_details,
                    lastUpdatedAt=job.last_updated_at,
                    presignedUrlConfig=job.presigned_url_config,
                    status=job.status,
                    targets=job.targets,
                    targetSelection=job.target_selection,
                ),
            )
        )

    def delete_job(self) -> str:
        job_id = self._get_param("jobId")
        force = self._get_bool_param("force")

        self.iot_backend.delete_job(job_id=job_id, force=force)

        return json.dumps(dict())

    def cancel_job(self) -> str:
        job_id = self._get_param("jobId")
        reason_code = self._get_param("reasonCode")
        comment = self._get_param("comment")
        force = self._get_bool_param("force")

        job = self.iot_backend.cancel_job(
            job_id=job_id, reason_code=reason_code, comment=comment, force=force
        )

        return json.dumps(job.to_dict())

    def get_job_document(self) -> str:
        job = self.iot_backend.get_job_document(job_id=self._get_param("jobId"))

        if job.document is not None:
            return json.dumps({"document": job.document})
        else:
            # job.document_source is not None:
            # TODO: needs to be implemented to get document_source's content from S3
            return json.dumps({"document": ""})

    def list_jobs(self) -> str:
        # not the default, but makes testing easier
        max_results = self._get_int_param("maxResults", 50)
        previous_next_token = self._get_param("nextToken")
        jobs, next_token = self.iot_backend.list_jobs(
            max_results=max_results, token=previous_next_token
        )

        return json.dumps(dict(jobs=jobs, nextToken=next_token))

    def describe_job_execution(self) -> str:
        job_id = self._get_param("jobId")
        thing_name = self._get_param("thingName")
        execution_number = self._get_int_param("executionNumber")
        job_execution = self.iot_backend.describe_job_execution(
            job_id=job_id, thing_name=thing_name, execution_number=execution_number
        )

        return json.dumps(dict(execution=job_execution.to_get_dict()))

    def cancel_job_execution(self) -> str:
        job_id = self._get_param("jobId")
        thing_name = self._get_param("thingName")
        force = self._get_bool_param("force")

        self.iot_backend.cancel_job_execution(
            job_id=job_id, thing_name=thing_name, force=force
        )

        return json.dumps(dict())

    def delete_job_execution(self) -> str:
        job_id = self._get_param("jobId")
        thing_name = self._get_param("thingName")
        execution_number = self._get_int_param("executionNumber")
        force = self._get_bool_param("force")

        self.iot_backend.delete_job_execution(
            job_id=job_id,
            thing_name=thing_name,
            execution_number=execution_number,
            force=force,
        )

        return json.dumps(dict())

    def list_job_executions_for_job(self) -> str:
        job_id = self._get_param("jobId")
        status = self._get_param("status")
        max_results = self._get_int_param(
            "maxResults", 50
        )  # not the default, but makes testing easier
        next_token = self._get_param("nextToken")
        job_executions, next_token = self.iot_backend.list_job_executions_for_job(
            job_id=job_id, status=status, max_results=max_results, token=next_token
        )

        return json.dumps(dict(executionSummaries=job_executions, nextToken=next_token))

    def list_job_executions_for_thing(self) -> str:
        thing_name = self._get_param("thingName")
        status = self._get_param("status")
        max_results = self._get_int_param(
            "maxResults", 50
        )  # not the default, but makes testing easier
        next_token = self._get_param("nextToken")
        job_executions, next_token = self.iot_backend.list_job_executions_for_thing(
            thing_name=thing_name,
            status=status,
            max_results=max_results,
            next_token=next_token,
        )

        return json.dumps(dict(executionSummaries=job_executions, nextToken=next_token))

    def create_keys_and_certificate(self) -> str:
        set_as_active = self._get_bool_param("setAsActive")
        cert, key_pair = self.iot_backend.create_keys_and_certificate(
            set_as_active=set_as_active
        )
        return json.dumps(
            dict(
                certificateArn=cert.arn,
                certificateId=cert.certificate_id,
                certificatePem=cert.certificate_pem,
                keyPair=key_pair,
            )
        )

    def delete_ca_certificate(self) -> str:
        certificate_id = self.path.split("/")[-1]
        self.iot_backend.delete_ca_certificate(certificate_id=certificate_id)
        return json.dumps(dict())

    def delete_certificate(self) -> str:
        certificate_id = self._get_param("certificateId")
        force_delete = self._get_bool_param("forceDelete", False)
        self.iot_backend.delete_certificate(certificate_id, force_delete)
        return json.dumps(dict())

    def describe_ca_certificate(self) -> str:
        certificate_id = self.path.split("/")[-1]
        certificate = self.iot_backend.describe_ca_certificate(
            certificate_id=certificate_id
        )
        return json.dumps(
            {
                "certificateDescription": certificate.to_description_dict(),
                "registrationConfig": certificate.registration_config,
            }
        )

    def describe_certificate(self) -> str:
        certificate_id = self._get_param("certificateId")
        certificate = self.iot_backend.describe_certificate(
            certificate_id=certificate_id
        )
        return json.dumps(
            dict(certificateDescription=certificate.to_description_dict())
        )

    def get_registration_code(self) -> str:
        code = self.iot_backend.get_registration_code()
        return json.dumps(dict(registrationCode=code))

    def list_certificates(self) -> str:
        # page_size = self._get_int_param("pageSize")
        # marker = self._get_param("marker")
        # ascending_order = self._get_param("ascendingOrder")
        certificates = self.iot_backend.list_certificates()
        return json.dumps(dict(certificates=[_.to_dict() for _ in certificates]))

    def list_certificates_by_ca(self) -> str:
        ca_certificate_id = self._get_param("caCertificateId")
        certificates = self.iot_backend.list_certificates_by_ca(ca_certificate_id)
        return json.dumps(dict(certificates=[_.to_dict() for _ in certificates]))

    def register_ca_certificate(self) -> str:
        ca_certificate = self._get_param("caCertificate")
        set_as_active = self._get_bool_param("setAsActive")
        registration_config = self._get_param("registrationConfig")

        cert = self.iot_backend.register_ca_certificate(
            ca_certificate=ca_certificate,
            set_as_active=set_as_active,
            registration_config=registration_config,
        )
        return json.dumps(
            dict(certificateId=cert.certificate_id, certificateArn=cert.arn)
        )

    def register_certificate(self) -> str:
        certificate_pem = self._get_param("certificatePem")
        ca_certificate_pem = self._get_param("caCertificatePem")
        set_as_active = self._get_bool_param("setAsActive")
        status = self._get_param("status")

        cert = self.iot_backend.register_certificate(
            certificate_pem=certificate_pem,
            ca_certificate_pem=ca_certificate_pem,
            set_as_active=set_as_active,
            status=status,
        )
        return json.dumps(
            dict(certificateId=cert.certificate_id, certificateArn=cert.arn)
        )

    def register_certificate_without_ca(self) -> str:
        certificate_pem = self._get_param("certificatePem")
        status = self._get_param("status")

        cert = self.iot_backend.register_certificate_without_ca(
            certificate_pem=certificate_pem, status=status
        )
        return json.dumps(
            dict(certificateId=cert.certificate_id, certificateArn=cert.arn)
        )

    def update_ca_certificate(self) -> str:
        certificate_id = self.path.split("/")[-1]
        new_status = self._get_param("newStatus")
        config = self._get_param("registrationConfig")
        self.iot_backend.update_ca_certificate(
            certificate_id=certificate_id, new_status=new_status, config=config
        )
        return json.dumps(dict())

    def update_certificate(self) -> str:
        certificate_id = self._get_param("certificateId")
        new_status = self._get_param("newStatus")
        self.iot_backend.update_certificate(
            certificate_id=certificate_id, new_status=new_status
        )
        return json.dumps(dict())

    def create_policy(self) -> str:
        policy_name = self._get_param("policyName")
        policy_document = self._get_param("policyDocument")
        policy = self.iot_backend.create_policy(
            policy_name=policy_name, policy_document=policy_document
        )
        return json.dumps(policy.to_dict_at_creation())

    def list_policies(self) -> str:
        policies = self.iot_backend.list_policies()

        return json.dumps(dict(policies=[_.to_dict() for _ in policies]))

    def get_policy(self) -> str:
        policy_name = self._get_param("policyName")
        policy = self.iot_backend.get_policy(policy_name=policy_name)
        return json.dumps(policy.to_get_dict())

    def delete_policy(self) -> str:
        policy_name = self._get_param("policyName")
        self.iot_backend.delete_policy(policy_name=policy_name)
        return json.dumps(dict())

    def create_policy_version(self) -> str:
        policy_name = self._get_param("policyName")
        policy_document = self._get_param("policyDocument")
        set_as_default = self._get_bool_param("setAsDefault")
        policy_version = self.iot_backend.create_policy_version(
            policy_name, policy_document, set_as_default
        )

        return json.dumps(dict(policy_version.to_dict_at_creation()))

    def set_default_policy_version(self) -> str:
        policy_name = self._get_param("policyName")
        version_id = self._get_param("policyVersionId")
        self.iot_backend.set_default_policy_version(policy_name, version_id)

        return json.dumps(dict())

    def get_policy_version(self) -> str:
        policy_name = self._get_param("policyName")
        version_id = self._get_param("policyVersionId")
        policy_version = self.iot_backend.get_policy_version(policy_name, version_id)
        return json.dumps(dict(policy_version.to_get_dict()))

    def list_policy_versions(self) -> str:
        policy_name = self._get_param("policyName")
        policiy_versions = self.iot_backend.list_policy_versions(
            policy_name=policy_name
        )

        return json.dumps(dict(policyVersions=[_.to_dict() for _ in policiy_versions]))

    def delete_policy_version(self) -> str:
        policy_name = self._get_param("policyName")
        version_id = self._get_param("policyVersionId")
        self.iot_backend.delete_policy_version(policy_name, version_id)

        return json.dumps(dict())

    def attach_policy(self) -> str:
        policy_name = self._get_param("policyName")
        target = self._get_param("target")
        self.iot_backend.attach_policy(policy_name=policy_name, target=target)
        return json.dumps(dict())

    def dispatch_attached_policies(
        self, request: Any, full_url: str, headers: Any
    ) -> TYPE_RESPONSE:
        # This endpoint requires specialized handling because it has
        # a uri parameter containing forward slashes that is not
        # correctly url encoded when we're running in server mode.
        # https://github.com/pallets/flask/issues/900
        self.setup_class(request, full_url, headers)
        self.querystring["Action"] = ["ListAttachedPolicies"]
        target = self.path.partition("/attached-policies/")[-1]
        self.querystring["target"] = [unquote(target)] if "%" in target else [target]
        return self.call_action()

    def list_attached_policies(self) -> str:
        principal = self._get_param("target")
        policies = self.iot_backend.list_attached_policies(target=principal)
        return json.dumps(dict(policies=[_.to_dict() for _ in policies]))

    def attach_principal_policy(self) -> str:
        policy_name = self._get_param("policyName")
        principal = self.headers.get("x-amzn-iot-principal")
        self.iot_backend.attach_principal_policy(
            policy_name=policy_name, principal_arn=principal
        )
        return json.dumps(dict())

    def detach_policy(self) -> str:
        policy_name = self._get_param("policyName")
        target = self._get_param("target")
        self.iot_backend.detach_policy(policy_name=policy_name, target=target)
        return json.dumps(dict())

    def detach_principal_policy(self) -> str:
        policy_name = self._get_param("policyName")
        principal = self.headers.get("x-amzn-iot-principal")
        self.iot_backend.detach_principal_policy(
            policy_name=policy_name, principal_arn=principal
        )
        return json.dumps(dict())

    def list_principal_policies(self) -> str:
        principal = self.headers.get("x-amzn-iot-principal")
        policies = self.iot_backend.list_principal_policies(principal_arn=principal)
        return json.dumps(dict(policies=[_.to_dict() for _ in policies]))

    def list_policy_principals(self) -> str:
        policy_name = self.headers.get("x-amzn-iot-policy")
        principals = self.iot_backend.list_policy_principals(policy_name=policy_name)
        return json.dumps(dict(principals=principals))

    def list_targets_for_policy(self) -> str:
        """https://docs.aws.amazon.com/iot/latest/apireference/API_ListTargetsForPolicy.html"""
        policy_name = self._get_param("policyName")
        principals = self.iot_backend.list_targets_for_policy(policy_name=policy_name)
        return json.dumps(dict(targets=principals))

    def attach_thing_principal(self) -> str:
        thing_name = self._get_param("thingName")
        principal = self.headers.get("x-amzn-principal")
        self.iot_backend.attach_thing_principal(
            thing_name=thing_name, principal_arn=principal
        )
        return json.dumps(dict())

    def detach_thing_principal(self) -> str:
        thing_name = self._get_param("thingName")
        principal = self.headers.get("x-amzn-principal")
        self.iot_backend.detach_thing_principal(
            thing_name=thing_name, principal_arn=principal
        )
        return json.dumps(dict())

    def list_principal_things(self) -> str:
        next_token = self._get_param("nextToken")
        # max_results = self._get_int_param("maxResults")
        principal = self.headers.get("x-amzn-principal")
        things = self.iot_backend.list_principal_things(principal_arn=principal)
        # TODO: implement pagination in the future
        next_token = None
        return json.dumps(dict(things=things, nextToken=next_token))

    def list_thing_principals(self) -> str:
        thing_name = self._get_param("thingName")
        principals = self.iot_backend.list_thing_principals(thing_name=thing_name)
        return json.dumps(dict(principals=principals))

    def describe_thing_group(self) -> str:
        thing_group_name = unquote(self.path.split("/thing-groups/")[-1])
        thing_group = self.iot_backend.describe_thing_group(
            thing_group_name=thing_group_name
        )
        return json.dumps(thing_group.to_dict())

    def create_thing_group(self) -> str:
        thing_group_name = unquote(self.path.split("/thing-groups/")[-1])
        parent_group_name = self._get_param("parentGroupName")
        thing_group_properties = self._get_param("thingGroupProperties")
        (
            thing_group_name,
            thing_group_arn,
            thing_group_id,
        ) = self.iot_backend.create_thing_group(
            thing_group_name=thing_group_name,
            parent_group_name=parent_group_name,
            thing_group_properties=thing_group_properties,
        )
        return json.dumps(
            dict(
                thingGroupName=thing_group_name,
                thingGroupArn=thing_group_arn,
                thingGroupId=thing_group_id,
            )
        )

    def delete_thing_group(self) -> str:
        thing_group_name = unquote(self.path.split("/thing-groups/")[-1])
        self.iot_backend.delete_thing_group(thing_group_name=thing_group_name)
        return json.dumps(dict())

    def list_thing_groups(self) -> str:
        # next_token = self._get_param("nextToken")
        # max_results = self._get_int_param("maxResults")
        parent_group = self._get_param("parentGroup")
        name_prefix_filter = self._get_param("namePrefixFilter")
        recursive = self._get_bool_param("recursive")
        thing_groups = self.iot_backend.list_thing_groups(
            parent_group=parent_group,
            name_prefix_filter=name_prefix_filter,
            recursive=recursive,
        )
        next_token = None
        rets = [
            {"groupName": _.thing_group_name, "groupArn": _.arn} for _ in thing_groups
        ]
        # TODO: implement pagination in the future
        return json.dumps(dict(thingGroups=rets, nextToken=next_token))

    def update_thing_group(self) -> str:
        thing_group_name = unquote(self.path.split("/thing-groups/")[-1])
        thing_group_properties = self._get_param("thingGroupProperties")
        expected_version = self._get_param("expectedVersion")
        version = self.iot_backend.update_thing_group(
            thing_group_name=thing_group_name,
            thing_group_properties=thing_group_properties,
            expected_version=expected_version,
        )
        return json.dumps(dict(version=version))

    def add_thing_to_thing_group(self) -> str:
        thing_group_name = self._get_param("thingGroupName")
        thing_group_arn = self._get_param("thingGroupArn")
        thing_name = self._get_param("thingName")
        thing_arn = self._get_param("thingArn")
        self.iot_backend.add_thing_to_thing_group(
            thing_group_name=thing_group_name,
            thing_group_arn=thing_group_arn,
            thing_name=thing_name,
            thing_arn=thing_arn,
        )
        return json.dumps(dict())

    def remove_thing_from_thing_group(self) -> str:
        thing_group_name = self._get_param("thingGroupName")
        thing_group_arn = self._get_param("thingGroupArn")
        thing_name = self._get_param("thingName")
        thing_arn = self._get_param("thingArn")
        self.iot_backend.remove_thing_from_thing_group(
            thing_group_name=thing_group_name,
            thing_group_arn=thing_group_arn,
            thing_name=thing_name,
            thing_arn=thing_arn,
        )
        return json.dumps(dict())

    def list_things_in_thing_group(self) -> str:
        thing_group_name = self._get_param("thingGroupName")
        things = self.iot_backend.list_things_in_thing_group(
            thing_group_name=thing_group_name
        )
        next_token = None
        thing_names = [_.thing_name for _ in things]
        return json.dumps(dict(things=thing_names, nextToken=next_token))

    def list_thing_groups_for_thing(self) -> str:
        thing_name = self._get_param("thingName")
        # next_token = self._get_param("nextToken")
        # max_results = self._get_int_param("maxResults")
        thing_groups = self.iot_backend.list_thing_groups_for_thing(
            thing_name=thing_name
        )
        next_token = None
        return json.dumps(dict(thingGroups=thing_groups, nextToken=next_token))

    def update_thing_groups_for_thing(self) -> str:
        thing_name = self._get_param("thingName")
        thing_groups_to_add = self._get_param("thingGroupsToAdd") or []
        thing_groups_to_remove = self._get_param("thingGroupsToRemove") or []
        self.iot_backend.update_thing_groups_for_thing(
            thing_name=thing_name,
            thing_groups_to_add=thing_groups_to_add,
            thing_groups_to_remove=thing_groups_to_remove,
        )
        return json.dumps(dict())

    def list_topic_rules(self) -> str:
        return json.dumps(dict(rules=self.iot_backend.list_topic_rules()))

    def get_topic_rule(self) -> str:
        return json.dumps(
            self.iot_backend.get_topic_rule(rule_name=self._get_param("ruleName"))
        )

    def create_topic_rule(self) -> str:
        self.iot_backend.create_topic_rule(
            rule_name=self._get_param("ruleName"),
            description=self._get_param("description"),
            rule_disabled=self._get_param("ruleDisabled"),
            actions=self._get_param("actions"),
            error_action=self._get_param("errorAction"),
            sql=self._get_param("sql"),
            aws_iot_sql_version=self._get_param("awsIotSqlVersion"),
        )
        return json.dumps(dict())

    def replace_topic_rule(self) -> str:
        self.iot_backend.replace_topic_rule(
            rule_name=self._get_param("ruleName"),
            description=self._get_param("description"),
            rule_disabled=self._get_param("ruleDisabled"),
            actions=self._get_param("actions"),
            error_action=self._get_param("errorAction"),
            sql=self._get_param("sql"),
            aws_iot_sql_version=self._get_param("awsIotSqlVersion"),
        )
        return json.dumps(dict())

    def delete_topic_rule(self) -> str:
        self.iot_backend.delete_topic_rule(rule_name=self._get_param("ruleName"))
        return json.dumps(dict())

    def enable_topic_rule(self) -> str:
        self.iot_backend.enable_topic_rule(rule_name=self._get_param("ruleName"))
        return json.dumps(dict())

    def disable_topic_rule(self) -> str:
        self.iot_backend.disable_topic_rule(rule_name=self._get_param("ruleName"))
        return json.dumps(dict())

    def create_domain_configuration(self) -> str:
        domain_configuration = self.iot_backend.create_domain_configuration(
            domain_configuration_name=self._get_param("domainConfigurationName"),
            domain_name=self._get_param("domainName"),
            server_certificate_arns=self._get_param("serverCertificateArns"),
            authorizer_config=self._get_param("authorizerConfig"),
            service_type=self._get_param("serviceType"),
        )
        return json.dumps(domain_configuration.to_dict())

    def delete_domain_configuration(self) -> str:
        self.iot_backend.delete_domain_configuration(
            domain_configuration_name=self._get_param("domainConfigurationName")
        )
        return json.dumps(dict())

    def describe_domain_configuration(self) -> str:
        domain_configuration = self.iot_backend.describe_domain_configuration(
            domain_configuration_name=self._get_param("domainConfigurationName")
        )
        return json.dumps(domain_configuration.to_description_dict())

    def list_domain_configurations(self) -> str:
        return json.dumps(
            dict(domainConfigurations=self.iot_backend.list_domain_configurations())
        )

    def update_domain_configuration(self) -> str:
        domain_configuration = self.iot_backend.update_domain_configuration(
            domain_configuration_name=self._get_param("domainConfigurationName"),
            authorizer_config=self._get_param("authorizerConfig"),
            domain_configuration_status=self._get_param("domainConfigurationStatus"),
            remove_authorizer_config=self._get_bool_param("removeAuthorizerConfig"),
        )
        return json.dumps(domain_configuration.to_dict())

    def search_index(self) -> str:
        query = self._get_param("queryString")
        things = self.iot_backend.search_index(query)
        return json.dumps({"things": things, "thingGroups": []})
