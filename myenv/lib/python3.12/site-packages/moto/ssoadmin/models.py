import json
from typing import Any, Dict, List, Optional

from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.common_models import BaseModel
from moto.core.utils import unix_time
from moto.iam.aws_managed_policies import aws_managed_policies_data
from moto.moto_api._internal import mock_random as random
from moto.utilities.paginator import paginate
from moto.utilities.utils import get_partition

from .exceptions import (
    ConflictException,
    ResourceNotFoundException,
    ServiceQuotaExceededException,
)
from .utils import PAGINATION_MODEL

# https://docs.aws.amazon.com/singlesignon/latest/userguide/limits.html
MAX_MANAGED_POLICIES_PER_PERMISSION_SET = 20


class AccountAssignment(BaseModel):
    def __init__(
        self,
        instance_arn: str,
        target_id: str,
        target_type: str,
        permission_set_arn: str,
        principal_type: str,
        principal_id: str,
    ):
        self.request_id = str(random.uuid4())
        self.instance_arn = instance_arn
        self.target_id = target_id
        self.target_type = target_type
        self.permission_set_arn = permission_set_arn
        self.principal_type = principal_type
        self.principal_id = principal_id
        self.created_date = unix_time()

    def to_json(
        self, include_creation_date: bool = False, include_request_id: bool = False
    ) -> Dict[str, Any]:
        summary: Dict[str, Any] = {
            "TargetId": self.target_id,
            "TargetType": self.target_type,
            "PermissionSetArn": self.permission_set_arn,
            "PrincipalType": self.principal_type,
            "PrincipalId": self.principal_id,
        }
        if include_creation_date:
            summary["CreatedDate"] = self.created_date
        if include_request_id:
            summary["RequestId"] = self.request_id
        return summary


class PermissionSet(BaseModel):
    def __init__(
        self,
        name: str,
        description: str,
        instance_arn: str,
        session_duration: str,
        relay_state: str,
        tags: List[Dict[str, str]],
    ):
        self.name = name
        self.description = description
        self.instance_arn = instance_arn
        self.permission_set_arn = PermissionSet.generate_id(instance_arn)
        self.session_duration = session_duration
        self.relay_state = relay_state
        self.tags = tags
        self.created_date = unix_time()
        self.inline_policy = ""
        self.managed_policies: List[ManagedPolicy] = list()
        self.customer_managed_policies: List[CustomerManagedPolicy] = list()
        self.total_managed_policies_attached = (
            0  # this will also include customer managed policies
        )

    def to_json(self, include_creation_date: bool = False) -> Dict[str, Any]:
        summary: Dict[str, Any] = {
            "Name": self.name,
            "Description": self.description,
            "PermissionSetArn": self.permission_set_arn,
            "SessionDuration": self.session_duration,
            "RelayState": self.relay_state,
        }
        if include_creation_date:
            summary["CreatedDate"] = self.created_date
        return summary

    @staticmethod
    def generate_id(instance_arn: str) -> str:
        return instance_arn + "/ps-" + random.get_random_string(length=16).lower()


class ManagedPolicy(BaseModel):
    def __init__(self, arn: str, name: str):
        self.arn = arn
        self.name = name

    def __eq__(self, other: Any) -> bool:
        if not isinstance(other, ManagedPolicy):
            return False
        return self.arn == other.arn


class CustomerManagedPolicy(BaseModel):
    def __init__(self, name: str, path: str = "/"):
        self.name = name
        self.path = path

    def __eq__(self, other: Any) -> bool:
        if not isinstance(other, CustomerManagedPolicy):
            return False
        return f"{self.path}{self.name}" == f"{other.path}{other.name}"


class Instance:
    def __init__(self, account_id: str, region: str):
        self.created_date = unix_time()
        self.identity_store_id = (
            f"d-{random.get_random_string(length=10, lower_case=True)}"
        )
        self.instance_arn = f"arn:{get_partition(region)}:sso:::instance/ssoins-{random.get_random_string(length=16, lower_case=True)}"
        self.account_id = account_id
        self.status = "ACTIVE"
        self.name: Optional[str] = None

        self.provisioned_permission_sets: List[PermissionSet] = []

    def to_json(self) -> Dict[str, Any]:
        return {
            "CreatedDate": self.created_date,
            "IdentityStoreId": self.identity_store_id,
            "InstanceArn": self.instance_arn,
            "Name": self.name,
            "OwnerAccountId": self.account_id,
            "Status": self.status,
        }


class SSOAdminBackend(BaseBackend):
    """Implementation of SSOAdmin APIs."""

    def __init__(self, region_name: str, account_id: str):
        super().__init__(region_name, account_id)
        self.account_assignments: List[AccountAssignment] = list()
        self.deleted_account_assignments: List[AccountAssignment] = list()
        self.permission_sets: List[PermissionSet] = list()
        self.aws_managed_policies: Optional[Dict[str, Any]] = None
        self.instances: List[Instance] = []

        self.instances.append(Instance(self.account_id, self.region_name))

    def create_account_assignment(
        self,
        instance_arn: str,
        target_id: str,
        target_type: str,
        permission_set_arn: str,
        principal_type: str,
        principal_id: str,
    ) -> Dict[str, Any]:
        assignment = AccountAssignment(
            instance_arn,
            target_id,
            target_type,
            permission_set_arn,
            principal_type,
            principal_id,
        )
        self.account_assignments.append(assignment)
        return assignment.to_json(include_creation_date=True, include_request_id=True)

    def delete_account_assignment(
        self,
        instance_arn: str,
        target_id: str,
        target_type: str,
        permission_set_arn: str,
        principal_type: str,
        principal_id: str,
    ) -> Dict[str, Any]:
        account = self._find_account(
            instance_arn,
            target_id,
            target_type,
            permission_set_arn,
            principal_type,
            principal_id,
        )
        self.deleted_account_assignments.append(account)
        self.account_assignments.remove(account)
        return account.to_json(include_creation_date=True, include_request_id=True)

    def _find_account(
        self,
        instance_arn: str,
        target_id: str,
        target_type: str,
        permission_set_arn: str,
        principal_type: str,
        principal_id: str,
    ) -> AccountAssignment:
        for account in self.account_assignments:
            instance_arn_match = account.instance_arn == instance_arn
            target_id_match = account.target_id == target_id
            target_type_match = account.target_type == target_type
            permission_set_match = account.permission_set_arn == permission_set_arn
            principal_type_match = account.principal_type == principal_type
            principal_id_match = account.principal_id == principal_id
            if (
                instance_arn_match
                and target_id_match
                and target_type_match
                and permission_set_match
                and principal_type_match
                and principal_id_match
            ):
                return account
        raise ResourceNotFoundException

    def _find_managed_policy(self, managed_policy_arn: str) -> ManagedPolicy:
        """
        Checks to make sure the managed policy exists.
        This pulls from moto/iam/aws_managed_policies.py
        """
        # Lazy loading of aws managed policies file
        if self.aws_managed_policies is None:
            self.aws_managed_policies = json.loads(aws_managed_policies_data)

        policy_name = managed_policy_arn.split("/")[-1]
        managed_policy = self.aws_managed_policies.get(policy_name, None)
        if managed_policy is not None:
            path = managed_policy.get("path", "/")
            expected_arn = f"arn:{self.partition}:iam::aws:policy{path}{policy_name}"

            if managed_policy_arn == expected_arn:
                return ManagedPolicy(managed_policy_arn, policy_name)
        raise ResourceNotFoundException(
            f"Policy does not exist with ARN: {managed_policy_arn}"
        )

    @paginate(PAGINATION_MODEL)  # type: ignore[misc]
    def list_account_assignments(
        self, instance_arn: str, account_id: str, permission_set_arn: str
    ) -> List[Dict[str, Any]]:
        account_assignments = []
        for assignment in self.account_assignments:
            if (
                assignment.instance_arn == instance_arn
                and assignment.target_id == account_id
                and assignment.permission_set_arn == permission_set_arn
            ):
                account_assignments.append(
                    {
                        "AccountId": assignment.target_id,
                        "PermissionSetArn": assignment.permission_set_arn,
                        "PrincipalType": assignment.principal_type,
                        "PrincipalId": assignment.principal_id,
                    }
                )
        return account_assignments

    @paginate(PAGINATION_MODEL)  # type: ignore[misc]
    def list_account_assignments_for_principal(
        self,
        filter_: Dict[str, Any],
        instance_arn: str,
        principal_id: str,
        principal_type: str,
    ) -> List[Dict[str, Any]]:
        return [
            {
                "AccountId": account_assignment.target_id,
                "PermissionSetArn": account_assignment.permission_set_arn,
                "PrincipalId": account_assignment.principal_id,
                "PrincipalType": account_assignment.principal_type,
            }
            for account_assignment in self.account_assignments
            if all(
                [
                    filter_.get("AccountId", account_assignment.target_id)
                    == account_assignment.target_id,
                    principal_id == account_assignment.principal_id,
                    principal_type == account_assignment.principal_type,
                    instance_arn == account_assignment.instance_arn,
                ]
            )
        ]

    def create_permission_set(
        self,
        name: str,
        description: str,
        instance_arn: str,
        session_duration: str,
        relay_state: str,
        tags: List[Dict[str, str]],
    ) -> Dict[str, Any]:
        permission_set = PermissionSet(
            name,
            description,
            instance_arn,
            session_duration,
            relay_state,
            tags,
        )
        self.permission_sets.append(permission_set)
        return permission_set.to_json(True)

    def update_permission_set(
        self,
        instance_arn: str,
        permission_set_arn: str,
        description: str,
        session_duration: str,
        relay_state: str,
    ) -> Dict[str, Any]:
        permission_set = self._find_permission_set(
            instance_arn,
            permission_set_arn,
        )
        self.permission_sets.remove(permission_set)
        permission_set.description = description
        permission_set.session_duration = session_duration
        permission_set.relay_state = relay_state
        self.permission_sets.append(permission_set)
        return permission_set.to_json(True)

    def describe_permission_set(
        self, instance_arn: str, permission_set_arn: str
    ) -> Dict[str, Any]:
        permission_set = self._find_permission_set(
            instance_arn,
            permission_set_arn,
        )
        return permission_set.to_json(True)

    def delete_permission_set(
        self, instance_arn: str, permission_set_arn: str
    ) -> Dict[str, Any]:
        permission_set = self._find_permission_set(
            instance_arn,
            permission_set_arn,
        )
        self.permission_sets.remove(permission_set)

        for instance in self.instances:
            try:
                instance.provisioned_permission_sets.remove(permission_set)
            except ValueError:
                pass

        return permission_set.to_json(include_creation_date=True)

    def _find_permission_set(
        self, instance_arn: str, permission_set_arn: str
    ) -> PermissionSet:
        for permission_set in self.permission_sets:
            instance_arn_match = permission_set.instance_arn == instance_arn
            permission_set_match = (
                permission_set.permission_set_arn == permission_set_arn
            )
            if instance_arn_match and permission_set_match:
                return permission_set
        ps_id = permission_set_arn.split("/")[-1]
        raise ResourceNotFoundException(
            message=f"Could not find PermissionSet with id {ps_id}"
        )

    @paginate(pagination_model=PAGINATION_MODEL)
    def list_permission_sets(self, instance_arn: str) -> List[PermissionSet]:
        permission_sets = []
        for permission_set in self.permission_sets:
            if permission_set.instance_arn == instance_arn:
                permission_sets.append(permission_set)
        return permission_sets

    def put_inline_policy_to_permission_set(
        self, instance_arn: str, permission_set_arn: str, inline_policy: str
    ) -> None:
        permission_set = self._find_permission_set(
            instance_arn,
            permission_set_arn,
        )
        permission_set.inline_policy = inline_policy

    def get_inline_policy_for_permission_set(
        self, instance_arn: str, permission_set_arn: str
    ) -> str:
        permission_set = self._find_permission_set(
            instance_arn,
            permission_set_arn,
        )
        return permission_set.inline_policy

    def delete_inline_policy_from_permission_set(
        self, instance_arn: str, permission_set_arn: str
    ) -> None:
        permission_set = self._find_permission_set(
            instance_arn,
            permission_set_arn,
        )
        permission_set.inline_policy = ""

    def attach_managed_policy_to_permission_set(
        self, instance_arn: str, permission_set_arn: str, managed_policy_arn: str
    ) -> None:
        permissionset = self._find_permission_set(
            instance_arn,
            permission_set_arn,
        )
        managed_policy = self._find_managed_policy(managed_policy_arn)

        permissionset_id = permission_set_arn.split("/")[-1]
        if managed_policy in permissionset.managed_policies:
            raise ConflictException(
                f"Permission set with id {permissionset_id} already has a typed link attachment to a manged policy with {managed_policy_arn}"
            )

        if (
            permissionset.total_managed_policies_attached
            >= MAX_MANAGED_POLICIES_PER_PERMISSION_SET
        ):
            permissionset_id = permission_set_arn.split("/")[-1]
            raise ServiceQuotaExceededException(
                f"You have exceeded AWS SSO limits. Cannot create ManagedPolicy more than {MAX_MANAGED_POLICIES_PER_PERMISSION_SET} for id {permissionset_id}. Please refer to https://docs.aws.amazon.com/singlesignon/latest/userguide/limits.html"
            )

        permissionset.managed_policies.append(managed_policy)
        permissionset.total_managed_policies_attached += 1

    @paginate(pagination_model=PAGINATION_MODEL)
    def list_managed_policies_in_permission_set(
        self,
        instance_arn: str,
        permission_set_arn: str,
    ) -> List[ManagedPolicy]:
        permissionset = self._find_permission_set(
            instance_arn,
            permission_set_arn,
        )
        return permissionset.managed_policies

    def _detach_managed_policy(
        self, instance_arn: str, permission_set_arn: str, managed_policy_arn: str
    ) -> None:
        # ensure permission_set exists
        permissionset = self._find_permission_set(
            instance_arn,
            permission_set_arn,
        )

        for managed_policy in permissionset.managed_policies:
            if managed_policy.arn == managed_policy_arn:
                permissionset.managed_policies.remove(managed_policy)
                permissionset.total_managed_policies_attached -= 1
                return

        raise ResourceNotFoundException(
            f"Could not find ManagedPolicy with arn {managed_policy_arn}"
        )

    def detach_managed_policy_from_permission_set(
        self, instance_arn: str, permission_set_arn: str, managed_policy_arn: str
    ) -> None:
        self._detach_managed_policy(
            instance_arn, permission_set_arn, managed_policy_arn
        )

    def attach_customer_managed_policy_reference_to_permission_set(
        self,
        instance_arn: str,
        permission_set_arn: str,
        customer_managed_policy_reference: Dict[str, str],
    ) -> None:
        permissionset = self._find_permission_set(
            permission_set_arn=permission_set_arn, instance_arn=instance_arn
        )

        name = customer_managed_policy_reference["Name"]
        path = customer_managed_policy_reference.get("Path", "/")  # default path is "/"
        customer_managed_policy = CustomerManagedPolicy(name=name, path=path)

        if customer_managed_policy in permissionset.customer_managed_policies:
            raise ConflictException(
                f"Given customer managed policy with name: {name}  and path {path} already attached"
            )

        if (
            permissionset.total_managed_policies_attached
            >= MAX_MANAGED_POLICIES_PER_PERMISSION_SET
        ):
            raise ServiceQuotaExceededException(
                f"Cannot attach managed policy: number of attached managed policies is already at maximum {MAX_MANAGED_POLICIES_PER_PERMISSION_SET}"
            )

        permissionset.customer_managed_policies.append(customer_managed_policy)
        permissionset.total_managed_policies_attached += 1

    @paginate(pagination_model=PAGINATION_MODEL)
    def list_customer_managed_policy_references_in_permission_set(
        self, instance_arn: str, permission_set_arn: str
    ) -> List[CustomerManagedPolicy]:
        permissionset = self._find_permission_set(
            permission_set_arn=permission_set_arn, instance_arn=instance_arn
        )
        return permissionset.customer_managed_policies

    def _detach_customer_managed_policy_from_permissionset(
        self,
        instance_arn: str,
        permission_set_arn: str,
        customer_managed_policy_reference: Dict[str, str],
    ) -> None:
        permissionset = self._find_permission_set(
            permission_set_arn=permission_set_arn, instance_arn=instance_arn
        )
        path: str = customer_managed_policy_reference.get("Path", "/")
        name: str = customer_managed_policy_reference["Name"]

        for customer_managed_policy in permissionset.customer_managed_policies:
            if (
                customer_managed_policy.name == name
                and customer_managed_policy.path == path
            ):
                permissionset.customer_managed_policies.remove(customer_managed_policy)
                permissionset.total_managed_policies_attached -= 1
                return

        raise ResourceNotFoundException(
            f"Given managed policy with name: {name}  and path {path} does not exist on PermissionSet"
        )

    def detach_customer_managed_policy_reference_from_permission_set(
        self,
        instance_arn: str,
        permission_set_arn: str,
        customer_managed_policy_reference: Dict[str, str],
    ) -> None:
        self._detach_customer_managed_policy_from_permissionset(
            instance_arn=instance_arn,
            permission_set_arn=permission_set_arn,
            customer_managed_policy_reference=customer_managed_policy_reference,
        )

    def describe_account_assignment_creation_status(
        self, account_assignment_creation_request_id: str, instance_arn: str
    ) -> Dict[str, Any]:
        for account in self.account_assignments:
            if account.request_id == account_assignment_creation_request_id:
                return account.to_json(
                    include_creation_date=True, include_request_id=True
                )

        raise ResourceNotFoundException

    def describe_account_assignment_deletion_status(
        self, account_assignment_deletion_request_id: str, instance_arn: str
    ) -> Dict[str, Any]:
        for account in self.deleted_account_assignments:
            if account.request_id == account_assignment_deletion_request_id:
                return account.to_json(
                    include_creation_date=True, include_request_id=True
                )

        raise ResourceNotFoundException

    def list_instances(self) -> List[Instance]:
        return self.instances

    def update_instance(self, instance_arn: str, name: str) -> None:
        for instance in self.instances:
            if instance.instance_arn == instance_arn:
                instance.name = name

    def provision_permission_set(
        self, instance_arn: str, permission_set_arn: str
    ) -> None:
        """
        The TargetType/TargetId parameters are currently ignored - PermissionSets are simply provisioned to the caller's account
        """
        permission_set = self._find_permission_set(instance_arn, permission_set_arn)
        instance = [i for i in self.instances if i.instance_arn == instance_arn][0]
        instance.provisioned_permission_sets.append(permission_set)

    def list_permission_sets_provisioned_to_account(
        self, instance_arn: str
    ) -> List[PermissionSet]:
        """
        The following parameters are not yet implemented: AccountId, ProvisioningStatus, MaxResults, NextToken
        """
        for instance in self.instances:
            if instance.instance_arn == instance_arn:
                return instance.provisioned_permission_sets
        return []

    def list_accounts_for_provisioned_permission_set(
        self, instance_arn: str, permission_set_arn: str
    ) -> List[str]:
        """
        The following parameters are not yet implemented: MaxResults, NextToken, ProvisioningStatus
        """
        for instance in self.instances:
            if instance.instance_arn == instance_arn:
                for ps in instance.provisioned_permission_sets:
                    if ps.permission_set_arn == permission_set_arn:
                        return [self.account_id]
        return []


ssoadmin_backends = BackendDict(SSOAdminBackend, "sso-admin")
