import json
from uuid import uuid4

from moto.core.responses import BaseResponse
from moto.core.utils import unix_time

from .models import SSOAdminBackend, ssoadmin_backends


class SSOAdminResponse(BaseResponse):
    """Handler for SSOAdmin requests and responses."""

    def __init__(self) -> None:
        super().__init__(service_name="sso-admin")

    @property
    def ssoadmin_backend(self) -> SSOAdminBackend:
        """Return backend instance specific for this region."""
        return ssoadmin_backends[self.current_account][self.region]

    def create_account_assignment(self) -> str:
        params = json.loads(self.body)
        instance_arn = params.get("InstanceArn")
        target_id = params.get("TargetId")
        target_type = params.get("TargetType")
        permission_set_arn = params.get("PermissionSetArn")
        principal_type = params.get("PrincipalType")
        principal_id = params.get("PrincipalId")
        summary = self.ssoadmin_backend.create_account_assignment(
            instance_arn=instance_arn,
            target_id=target_id,
            target_type=target_type,
            permission_set_arn=permission_set_arn,
            principal_type=principal_type,
            principal_id=principal_id,
        )
        summary["Status"] = "SUCCEEDED"
        return json.dumps({"AccountAssignmentCreationStatus": summary})

    def delete_account_assignment(self) -> str:
        params = json.loads(self.body)
        instance_arn = params.get("InstanceArn")
        target_id = params.get("TargetId")
        target_type = params.get("TargetType")
        permission_set_arn = params.get("PermissionSetArn")
        principal_type = params.get("PrincipalType")
        principal_id = params.get("PrincipalId")
        summary = self.ssoadmin_backend.delete_account_assignment(
            instance_arn=instance_arn,
            target_id=target_id,
            target_type=target_type,
            permission_set_arn=permission_set_arn,
            principal_type=principal_type,
            principal_id=principal_id,
        )
        summary["Status"] = "SUCCEEDED"
        return json.dumps({"AccountAssignmentDeletionStatus": summary})

    def list_account_assignments(self) -> str:
        params = json.loads(self.body)
        instance_arn = params.get("InstanceArn")
        account_id = params.get("AccountId")
        permission_set_arn = params.get("PermissionSetArn")
        max_results = self._get_param("MaxResults")
        next_token = self._get_param("NextToken")

        assignments, next_token = self.ssoadmin_backend.list_account_assignments(
            instance_arn=instance_arn,
            account_id=account_id,
            permission_set_arn=permission_set_arn,
            next_token=next_token,
            max_results=max_results,
        )

        return json.dumps(dict(AccountAssignments=assignments, NextToken=next_token))

    def list_account_assignments_for_principal(self) -> str:
        filter_ = self._get_param("Filter", {})
        instance_arn = self._get_param("InstanceArn")
        max_results = self._get_param("MaxResults")
        next_token = self._get_param("NextToken")
        principal_id = self._get_param("PrincipalId")
        principal_type = self._get_param("PrincipalType")

        (
            assignments,
            next_token,
        ) = self.ssoadmin_backend.list_account_assignments_for_principal(
            filter_=filter_,
            instance_arn=instance_arn,
            max_results=max_results,
            next_token=next_token,
            principal_id=principal_id,
            principal_type=principal_type,
        )

        return json.dumps(dict(AccountAssignments=assignments, NextToken=next_token))

    def create_permission_set(self) -> str:
        name = self._get_param("Name")
        description = self._get_param("Description")
        instance_arn = self._get_param("InstanceArn")
        session_duration = self._get_param("SessionDuration", 3600)
        relay_state = self._get_param("RelayState")
        tags = self._get_param("Tags")

        permission_set = self.ssoadmin_backend.create_permission_set(
            name=name,
            description=description,
            instance_arn=instance_arn,
            session_duration=session_duration,
            relay_state=relay_state,
            tags=tags,
        )

        return json.dumps({"PermissionSet": permission_set})

    def delete_permission_set(self) -> str:
        params = json.loads(self.body)
        instance_arn = params.get("InstanceArn")
        permission_set_arn = params.get("PermissionSetArn")
        self.ssoadmin_backend.delete_permission_set(
            instance_arn=instance_arn,
            permission_set_arn=permission_set_arn,
        )
        return "{}"

    def update_permission_set(self) -> str:
        instance_arn = self._get_param("InstanceArn")
        permission_set_arn = self._get_param("PermissionSetArn")
        description = self._get_param("Description")
        session_duration = self._get_param("SessionDuration", 3600)
        relay_state = self._get_param("RelayState")

        self.ssoadmin_backend.update_permission_set(
            instance_arn=instance_arn,
            permission_set_arn=permission_set_arn,
            description=description,
            session_duration=session_duration,
            relay_state=relay_state,
        )
        return "{}"

    def describe_permission_set(self) -> str:
        instance_arn = self._get_param("InstanceArn")
        permission_set_arn = self._get_param("PermissionSetArn")

        permission_set = self.ssoadmin_backend.describe_permission_set(
            instance_arn=instance_arn,
            permission_set_arn=permission_set_arn,
        )
        return json.dumps({"PermissionSet": permission_set})

    def list_permission_sets(self) -> str:
        instance_arn = self._get_param("InstanceArn")
        max_results = self._get_int_param("MaxResults")
        next_token = self._get_param("NextToken")
        permission_sets, next_token = self.ssoadmin_backend.list_permission_sets(
            instance_arn=instance_arn, max_results=max_results, next_token=next_token
        )
        permission_set_ids = []
        for permission_set in permission_sets:
            permission_set_ids.append(permission_set.permission_set_arn)
        response = {"PermissionSets": permission_set_ids}
        if next_token:
            response["NextToken"] = next_token
        return json.dumps(response)

    def put_inline_policy_to_permission_set(self) -> str:
        instance_arn = self._get_param("InstanceArn")
        permission_set_arn = self._get_param("PermissionSetArn")
        inline_policy = self._get_param("InlinePolicy")
        self.ssoadmin_backend.put_inline_policy_to_permission_set(
            instance_arn=instance_arn,
            permission_set_arn=permission_set_arn,
            inline_policy=inline_policy,
        )
        return json.dumps({})

    def get_inline_policy_for_permission_set(self) -> str:
        instance_arn = self._get_param("InstanceArn")
        permission_set_arn = self._get_param("PermissionSetArn")
        inline_policy = self.ssoadmin_backend.get_inline_policy_for_permission_set(
            instance_arn=instance_arn,
            permission_set_arn=permission_set_arn,
        )
        return json.dumps({"InlinePolicy": inline_policy})

    def delete_inline_policy_from_permission_set(self) -> str:
        instance_arn = self._get_param("InstanceArn")
        permission_set_arn = self._get_param("PermissionSetArn")
        self.ssoadmin_backend.delete_inline_policy_from_permission_set(
            instance_arn=instance_arn,
            permission_set_arn=permission_set_arn,
        )
        return json.dumps({})

    def attach_managed_policy_to_permission_set(self) -> str:
        instance_arn = self._get_param("InstanceArn")
        permission_set_arn = self._get_param("PermissionSetArn")
        managed_policy_arn = self._get_param("ManagedPolicyArn")
        self.ssoadmin_backend.attach_managed_policy_to_permission_set(
            instance_arn=instance_arn,
            permission_set_arn=permission_set_arn,
            managed_policy_arn=managed_policy_arn,
        )
        return json.dumps({})

    def list_managed_policies_in_permission_set(self) -> str:
        instance_arn = self._get_param("InstanceArn")
        permission_set_arn = self._get_param("PermissionSetArn")
        max_results = self._get_int_param("MaxResults")
        next_token = self._get_param("NextToken")

        (
            managed_policies,
            next_token,
        ) = self.ssoadmin_backend.list_managed_policies_in_permission_set(
            instance_arn=instance_arn,
            permission_set_arn=permission_set_arn,
            max_results=max_results,
            next_token=next_token,
        )

        managed_policies_response = [
            {"Arn": managed_policy.arn, "Name": managed_policy.name}
            for managed_policy in managed_policies
        ]
        return json.dumps(
            {
                "AttachedManagedPolicies": managed_policies_response,
                "NextToken": next_token,
            }
        )

    def detach_managed_policy_from_permission_set(self) -> str:
        instance_arn = self._get_param("InstanceArn")
        permission_set_arn = self._get_param("PermissionSetArn")
        managed_policy_arn = self._get_param("ManagedPolicyArn")
        self.ssoadmin_backend.detach_managed_policy_from_permission_set(
            instance_arn=instance_arn,
            permission_set_arn=permission_set_arn,
            managed_policy_arn=managed_policy_arn,
        )
        return json.dumps({})

    def attach_customer_managed_policy_reference_to_permission_set(self) -> str:
        instance_arn = self._get_param("InstanceArn")
        permission_set_arn = self._get_param("PermissionSetArn")
        customer_managed_policy_reference = self._get_param(
            "CustomerManagedPolicyReference"
        )
        self.ssoadmin_backend.attach_customer_managed_policy_reference_to_permission_set(
            instance_arn=instance_arn,
            permission_set_arn=permission_set_arn,
            customer_managed_policy_reference=customer_managed_policy_reference,
        )
        return json.dumps({})

    def list_customer_managed_policy_references_in_permission_set(self) -> str:
        instance_arn = self._get_param("InstanceArn")
        permission_set_arn = self._get_param("PermissionSetArn")
        max_results = self._get_int_param("MaxResults")
        next_token = self._get_param("NextToken")

        (
            customer_managed_policy_references,
            next_token,
        ) = self.ssoadmin_backend.list_customer_managed_policy_references_in_permission_set(
            instance_arn=instance_arn,
            permission_set_arn=permission_set_arn,
            max_results=max_results,
            next_token=next_token,
        )

        customer_managed_policy_references_response = [
            {
                "Name": customer_managed_policy_reference.name,
                "Path": customer_managed_policy_reference.path,
            }
            for customer_managed_policy_reference in customer_managed_policy_references
        ]
        return json.dumps(
            {
                "CustomerManagedPolicyReferences": customer_managed_policy_references_response,
                "NextToken": next_token,
            }
        )

    def detach_customer_managed_policy_reference_from_permission_set(self) -> str:
        instance_arn = self._get_param("InstanceArn")
        permission_set_arn = self._get_param("PermissionSetArn")
        customer_managed_policy_reference = self._get_param(
            "CustomerManagedPolicyReference"
        )
        self.ssoadmin_backend.detach_customer_managed_policy_reference_from_permission_set(
            instance_arn=instance_arn,
            permission_set_arn=permission_set_arn,
            customer_managed_policy_reference=customer_managed_policy_reference,
        )
        return json.dumps({})

    def describe_account_assignment_creation_status(self) -> str:
        account_assignment_creation_request_id = self._get_param(
            "AccountAssignmentCreationRequestId"
        )
        instance_arn = self._get_param("InstanceArn")
        account_assignment_creation_status = self.ssoadmin_backend.describe_account_assignment_creation_status(
            account_assignment_creation_request_id=account_assignment_creation_request_id,
            instance_arn=instance_arn,
        )
        account_assignment_creation_status["Status"] = "SUCCEEDED"
        return json.dumps(
            dict(AccountAssignmentCreationStatus=account_assignment_creation_status)
        )

    def describe_account_assignment_deletion_status(self) -> str:
        account_assignment_deletion_request_id = self._get_param(
            "AccountAssignmentDeletionRequestId"
        )
        instance_arn = self._get_param("InstanceArn")
        account_assignment_deletion_status = self.ssoadmin_backend.describe_account_assignment_deletion_status(
            account_assignment_deletion_request_id=account_assignment_deletion_request_id,
            instance_arn=instance_arn,
        )
        account_assignment_deletion_status["Status"] = "SUCCEEDED"
        return json.dumps(
            dict(AccountAssignmentDeletionStatus=account_assignment_deletion_status)
        )

    def list_instances(self) -> str:
        instances = self.ssoadmin_backend.list_instances()

        return json.dumps({"Instances": [i.to_json() for i in instances]})

    def update_instance(self) -> str:
        instance_arn = self._get_param("InstanceArn")
        name = self._get_param("Name")

        self.ssoadmin_backend.update_instance(instance_arn=instance_arn, name=name)

        return "{}"

    def provision_permission_set(self) -> str:
        instance_arn = self._get_param("InstanceArn")
        permission_set_arn = self._get_param("PermissionSetArn")

        self.ssoadmin_backend.provision_permission_set(
            instance_arn=instance_arn,
            permission_set_arn=permission_set_arn,
        )
        return json.dumps(
            {
                "PermissionSetProvisioningStatus": {
                    "AccountId": self.current_account,
                    "CreatedDate": unix_time(),
                    "PermissionSetArn": permission_set_arn,
                    "RequestId": str(uuid4()),
                    "Status": "SUCCEEDED",
                }
            }
        )

    def list_permission_sets_provisioned_to_account(self) -> str:
        instance_arn = self._get_param("InstanceArn")

        permission_sets = (
            self.ssoadmin_backend.list_permission_sets_provisioned_to_account(
                instance_arn
            )
        )
        arns = [p.permission_set_arn for p in permission_sets]
        return json.dumps({"PermissionSets": arns})

    def list_accounts_for_provisioned_permission_set(self) -> str:
        instance_arn = self._get_param("InstanceArn")
        permission_set_arn = self._get_param("PermissionSetArn")

        account_ids = (
            self.ssoadmin_backend.list_accounts_for_provisioned_permission_set(
                instance_arn=instance_arn, permission_set_arn=permission_set_arn
            )
        )
        return json.dumps({"AccountIds": account_ids})
