import json
import re
from typing import Any, Dict, List, Optional

from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.common_models import BaseModel
from moto.core.exceptions import RESTError
from moto.core.utils import unix_time, utcnow
from moto.organizations import utils
from moto.organizations.exceptions import (
    AccountAlreadyRegisteredException,
    AccountNotFoundException,
    AccountNotRegisteredException,
    AWSOrganizationsNotInUseException,
    ConstraintViolationException,
    DuplicateOrganizationalUnitException,
    DuplicatePolicyException,
    InvalidInputException,
    PolicyNotFoundException,
    PolicyTypeAlreadyEnabledException,
    PolicyTypeNotEnabledException,
    RootNotFoundException,
    TargetNotFoundException,
)
from moto.utilities.paginator import paginate
from moto.utilities.utils import PARTITION_NAMES, get_partition

from .utils import PAGINATION_MODEL


class FakeOrganization(BaseModel):
    def __init__(self, account_id: str, region_name: str, feature_set: str):
        self.id = utils.make_random_org_id()
        self.root_id = utils.make_random_root_id()
        self.feature_set = feature_set
        self.master_account_id = account_id
        self.master_account_email = utils.MASTER_ACCOUNT_EMAIL
        self.available_policy_types = [
            # This policy is available, but not applied
            # User should use enable_policy_type/disable_policy_type to do anything else
            # This field is deprecated in AWS, but we'll return it for old time's sake
            {"Type": "SERVICE_CONTROL_POLICY", "Status": "ENABLED"}
        ]
        self.region = region_name

    @property
    def arn(self) -> str:
        partition = get_partition(self.region)
        return utils.ORGANIZATION_ARN_FORMAT.format(
            partition, self.master_account_id, self.id
        )

    @property
    def master_account_arn(self) -> str:
        partition = get_partition(self.region)
        return utils.MASTER_ACCOUNT_ARN_FORMAT.format(
            partition, self.master_account_id, self.id
        )

    def describe(self) -> Dict[str, Any]:
        return {
            "Organization": {
                "Id": self.id,
                "Arn": self.arn,
                "FeatureSet": self.feature_set,
                "MasterAccountArn": self.master_account_arn,
                "MasterAccountId": self.master_account_id,
                "MasterAccountEmail": self.master_account_email,
                "AvailablePolicyTypes": self.available_policy_types,
            }
        }


class FakeAccount(BaseModel):
    def __init__(self, organization: FakeOrganization, **kwargs: Any):
        self.type = "ACCOUNT"
        self.region = organization.region
        self.organization_id = organization.id
        self.master_account_id = organization.master_account_id
        self.create_account_status_id = utils.make_random_create_account_status_id()
        self.id = utils.make_random_account_id()
        self.name = kwargs["AccountName"]
        self.email = kwargs["Email"]
        self.create_time = utcnow()
        self.status = "ACTIVE"
        self.joined_method = "CREATED"
        self.parent_id = organization.root_id
        self.attached_policies: List[FakePolicy] = []
        self.tags = {tag["Key"]: tag["Value"] for tag in kwargs.get("Tags", [])}

    @property
    def arn(self) -> str:
        partition = get_partition(self.region)
        return utils.ACCOUNT_ARN_FORMAT.format(
            partition, self.master_account_id, self.organization_id, self.id
        )

    @property
    def create_account_status(self) -> Dict[str, Any]:  # type: ignore[misc]
        return {
            "CreateAccountStatus": {
                "Id": self.create_account_status_id,
                "AccountName": self.name,
                "State": "SUCCEEDED",
                "RequestedTimestamp": unix_time(self.create_time),
                "CompletedTimestamp": unix_time(self.create_time),
                "AccountId": self.id,
            }
        }

    def describe(self) -> Dict[str, Any]:
        return {
            "Id": self.id,
            "Arn": self.arn,
            "Email": self.email,
            "Name": self.name,
            "Status": self.status,
            "JoinedMethod": self.joined_method,
            "JoinedTimestamp": unix_time(self.create_time),
        }

    def close(self) -> None:
        # TODO: The CloseAccount spec allows the account to pass through a
        # "PENDING_CLOSURE" state before reaching the SUSPENDED state.
        self.status = "SUSPENDED"


class FakeOrganizationalUnit(BaseModel):
    def __init__(self, organization: FakeOrganization, **kwargs: Any):
        self.type = "ORGANIZATIONAL_UNIT"
        self.region = organization.region
        self.organization_id = organization.id
        self.master_account_id = organization.master_account_id
        self.id = utils.make_random_ou_id(organization.root_id)
        self.name = kwargs.get("Name")
        self.parent_id = kwargs.get("ParentId")
        self._arn_format = utils.OU_ARN_FORMAT
        self.attached_policies: List[FakePolicy] = []
        self.tags = {tag["Key"]: tag["Value"] for tag in kwargs.get("Tags", [])}

    @property
    def arn(self) -> str:
        partition = get_partition(self.region)
        return self._arn_format.format(
            partition, self.master_account_id, self.organization_id, self.id
        )

    def describe(self) -> Dict[str, Dict[str, Any]]:
        return {
            "OrganizationalUnit": {"Id": self.id, "Arn": self.arn, "Name": self.name}
        }


class FakeRoot(FakeOrganizationalUnit):
    SUPPORTED_POLICY_TYPES = [
        "AISERVICES_OPT_OUT_POLICY",
        "BACKUP_POLICY",
        "SERVICE_CONTROL_POLICY",
        "TAG_POLICY",
    ]

    def __init__(self, organization: FakeOrganization, **kwargs: Any):
        super().__init__(organization, **kwargs)
        self.type = "ROOT"
        self.id = organization.root_id
        self.name = "Root"
        self.policy_types: List[Dict[str, str]] = []
        self._arn_format = utils.ROOT_ARN_FORMAT
        self.attached_policies = []
        self.tags = {tag["Key"]: tag["Value"] for tag in kwargs.get("Tags", [])}

    def describe(self) -> Dict[str, Any]:
        return {
            "Id": self.id,
            "Arn": self.arn,
            "Name": self.name,
            "PolicyTypes": self.policy_types,
        }

    def add_policy_type(self, policy_type: str) -> None:
        if policy_type not in self.SUPPORTED_POLICY_TYPES:
            raise InvalidInputException("You specified an invalid value.")

        if any(type["Type"] == policy_type for type in self.policy_types):
            raise PolicyTypeAlreadyEnabledException

        self.policy_types.append({"Type": policy_type, "Status": "ENABLED"})

    def remove_policy_type(self, policy_type: str) -> None:
        if not FakePolicy.supported_policy_type(policy_type):
            raise InvalidInputException("You specified an invalid value.")

        if all(type["Type"] != policy_type for type in self.policy_types):
            raise PolicyTypeNotEnabledException

        self.policy_types.remove({"Type": policy_type, "Status": "ENABLED"})


class FakePolicy(BaseModel):
    SUPPORTED_POLICY_TYPES = [
        "AISERVICES_OPT_OUT_POLICY",
        "BACKUP_POLICY",
        "SERVICE_CONTROL_POLICY",
        "TAG_POLICY",
    ]

    def __init__(self, organization: FakeOrganization, **kwargs: Any):
        self.content = kwargs.get("Content")
        self.description = kwargs.get("Description")
        self.name = kwargs.get("Name")
        self.type = kwargs.get("Type", "")
        self.id = utils.make_random_policy_id()
        self.aws_managed = False
        self.region = organization.region
        self.organization_id = organization.id
        self.master_account_id = organization.master_account_id
        self.attachments: List[Any] = []
        self.tags = {tag["Key"]: tag["Value"] for tag in kwargs.get("Tags", [])}

        if not FakePolicy.supported_policy_type(self.type):
            raise InvalidInputException("You specified an invalid value.")
        elif self.type == "AISERVICES_OPT_OUT_POLICY":
            self._arn_format = utils.AI_POLICY_ARN_FORMAT
        elif self.type == "SERVICE_CONTROL_POLICY":
            self._arn_format = utils.SCP_ARN_FORMAT
        else:
            raise NotImplementedError(
                f"The {self.type} policy type has not been implemented"
            )

    @property
    def arn(self) -> str:
        partition = get_partition(self.region)
        return self._arn_format.format(
            partition, self.master_account_id, self.organization_id, self.id
        )

    def describe(self) -> Dict[str, Any]:
        return {
            "Policy": {
                "PolicySummary": {
                    "Id": self.id,
                    "Arn": self.arn,
                    "Name": self.name,
                    "Description": self.description,
                    "Type": self.type,
                    "AwsManaged": self.aws_managed,
                },
                "Content": self.content,
            }
        }

    @staticmethod
    def supported_policy_type(policy_type: str) -> bool:
        return policy_type in FakePolicy.SUPPORTED_POLICY_TYPES


class FakeServiceAccess(BaseModel):
    # List of trusted services, which support trusted access with Organizations
    # https://docs.aws.amazon.com/organizations/latest/userguide/orgs_integrated-services-list.html
    TRUSTED_SERVICES = [
        "access-analyzer.amazonaws.com",
        "account.amazonaws.com",  # AWS Account Management
        "auditmanager.amazonaws.com",  # AWS Audit Manager
        "aws-artifact-account-sync.amazonaws.com",
        "backup.amazonaws.com",  # AWS Backup
        "cloudtrail.amazonaws.com",  # AWS Cloudtrail
        "compute-optimizer.amazonaws.com",  # AWS Compute Optimizer
        "config.amazonaws.com",  # AWS Config
        "config-multiaccountsetup.amazonaws.com",
        "controltower.amazonaws.com",  # AWS Control Tower
        "detective.amazonaws.com",  # AWS Detective
        "devops-guru.amazonaws.com",  # Amazon DevOps Guru
        "ds.amazonaws.com",  # AWS Directory Service
        "fms.amazonaws.com",  # AWS Firewall Manager
        "guardduty.amazonaws.com",  # Amazon GuardDuty
        "health.amazonaws.com",  # Amazon Health
        "inspector2.amazonaws.com",  # Amazon Inspector
        "ipam.amazonaws.com",  # AWS VPC IP Address Manager
        "license-manager.amazonaws.com",  # AWS License Manager
        "license-manager.member-account.amazonaws.com.",  # AWS License Manager
        "license-manager-linux-subscriptions.amazonaws.com",  # AWS License Manager
        "license-management.marketplace.amazonaws.com",  # AWS Marketplace
        "macie.amazonaws.com",  # Amazon Macie
        "member.org.stacksets.cloudformation.amazonaws.com",
        "mgn.amazonaws.com",  # AWS Application Migration Service
        "ram.amazonaws.com",  # AWS Resource Access Manager
        "reporting.trustedadvisor.amazonaws.com",  # AWS Trusted Advisor
        "reachabilityanalyzer.networkinsights.amazonaws.com",  # Reachability Analyzer
        "securityhub.amazonaws.com",  # AWS Security Hub
        "storage-lens.s3.amazonaws.com",  # Amazon S3 Storage Lens
        "securitylake.amazonaws.com",  # Amazon Security Lake
        "servicecatalog.amazonaws.com",  # AWS Service Catalog
        "servicequotas.amazonaws.com",  # Service Quotas
        "stacksets.cloudformation.amazonaws.com",
        "sso.amazonaws.com",  # AWS SSO
        "ssm.amazonaws.com",  # AWS Systems Manager
        "tagpolicies.tag.amazonaws.com",  # Tag policies
        "wellarchitected.amazonaws.com",  # AWS Well Architected Tool
    ]

    def __init__(self, **kwargs: Any):
        if not self.trusted_service(kwargs["ServicePrincipal"]):
            raise InvalidInputException(
                "You specified an unrecognized service principal."
            )

        self.service_principal = kwargs["ServicePrincipal"]
        self.date_enabled = utcnow()

    def describe(self) -> Dict[str, Any]:
        return {
            "ServicePrincipal": self.service_principal,
            "DateEnabled": unix_time(self.date_enabled),
        }

    @staticmethod
    def trusted_service(service_principal: str) -> bool:
        return service_principal in FakeServiceAccess.TRUSTED_SERVICES


class FakeDelegatedAdministrator(BaseModel):
    # List of services, which support a different Account to ba a delegated administrator
    # https://docs.aws.amazon.com/organizations/latest/userguide/orgs_integrated-services-list.html
    SUPPORTED_SERVICES = [
        "config-multiaccountsetup.amazonaws.com",
        "guardduty.amazonaws.com",
        "access-analyzer.amazonaws.com",
        "macie.amazonaws.com",
        "servicecatalog.amazonaws.com",
        "ssm.amazonaws.com",
    ]

    def __init__(self, account: FakeAccount):
        self.account = account
        self.enabled_date = utcnow()
        self.services: Dict[str, Any] = {}

    def add_service_principal(self, service_principal: str) -> None:
        if service_principal in self.services:
            raise AccountAlreadyRegisteredException

        if not self.supported_service(service_principal):
            raise InvalidInputException(
                "You specified an unrecognized service principal."
            )

        self.services[service_principal] = {
            "ServicePrincipal": service_principal,
            "DelegationEnabledDate": unix_time(),
        }

    def remove_service_principal(self, service_principal: str) -> None:
        if service_principal not in self.services:
            raise InvalidInputException(
                "You specified an unrecognized service principal."
            )

        self.services.pop(service_principal)

    def describe(self) -> Dict[str, Any]:
        admin = self.account.describe()
        admin["DelegationEnabledDate"] = unix_time(self.enabled_date)

        return admin

    @staticmethod
    def supported_service(service_principal: str) -> bool:
        return service_principal in FakeDelegatedAdministrator.SUPPORTED_SERVICES


class OrganizationsBackend(BaseBackend):
    def __init__(self, region_name: str, account_id: str):
        super().__init__(region_name, account_id)
        self._reset()

    def _reset(self) -> None:
        self.org: Optional[FakeOrganization] = None
        self.accounts: List[FakeAccount] = []
        self.ou: List[FakeOrganizationalUnit] = []
        self.policies: List[FakePolicy] = []
        self.services: List[Dict[str, Any]] = []
        self.admins: List[FakeDelegatedAdministrator] = []

    def _get_root_by_id(self, root_id: str) -> FakeRoot:
        root = next((ou for ou in self.ou if ou.id == root_id), None)
        if not root:
            raise RootNotFoundException

        return root  # type: ignore[return-value]

    def create_organization(self, region: str, **kwargs: Any) -> Dict[str, Any]:
        self.org = FakeOrganization(
            self.account_id,
            region_name=region,
            feature_set=kwargs.get("FeatureSet") or "ALL",
        )
        root_ou = FakeRoot(self.org)
        self.ou.append(root_ou)
        master_account = FakeAccount(
            self.org, AccountName="master", Email=self.org.master_account_email
        )
        master_account.id = self.org.master_account_id
        self.accounts.append(master_account)
        default_policy = FakePolicy(
            self.org,
            Name="FullAWSAccess",
            Description="Allows access to every operation",
            Type="SERVICE_CONTROL_POLICY",
            Content=json.dumps(
                {
                    "Version": "2012-10-17",
                    "Statement": [{"Effect": "Allow", "Action": "*", "Resource": "*"}],
                }
            ),
        )
        default_policy.id = utils.DEFAULT_POLICY_ID
        default_policy.aws_managed = True
        self.policies.append(default_policy)
        self.attach_policy(PolicyId=default_policy.id, TargetId=root_ou.id)
        self.attach_policy(PolicyId=default_policy.id, TargetId=master_account.id)
        return self.org.describe()

    def describe_organization(self) -> Dict[str, Any]:
        if not self.org:
            raise AWSOrganizationsNotInUseException
        return self.org.describe()

    def delete_organization(self) -> None:
        if [account for account in self.accounts if account.name != "master"]:
            raise RESTError(
                "OrganizationNotEmptyException",
                "To delete an organization you must first remove all member accounts (except the master).",
            )
        self._reset()

    def list_roots(self) -> Dict[str, Any]:
        return dict(Roots=[ou.describe() for ou in self.ou if isinstance(ou, FakeRoot)])

    def create_organizational_unit(self, **kwargs: Any) -> Dict[str, Any]:
        new_ou = FakeOrganizationalUnit(self.org, **kwargs)  # type: ignore
        self.ou.append(new_ou)
        self.attach_policy(PolicyId=utils.DEFAULT_POLICY_ID, TargetId=new_ou.id)
        return new_ou.describe()

    def delete_organizational_unit(self, **kwargs: Any) -> None:
        ou_to_delete = self.get_organizational_unit_by_id(
            kwargs["OrganizationalUnitId"]
        )
        self.ou.remove(ou_to_delete)

    def update_organizational_unit(self, **kwargs: Any) -> Dict[str, Any]:
        for ou in self.ou:
            if ou.name == kwargs["Name"]:
                raise DuplicateOrganizationalUnitException
        ou = self.get_organizational_unit_by_id(kwargs["OrganizationalUnitId"])
        ou.name = kwargs["Name"]
        return ou.describe()

    def get_organizational_unit_by_id(self, ou_id: str) -> FakeOrganizationalUnit:
        ou = next((ou for ou in self.ou if ou.id == ou_id), None)
        if ou is None:
            raise RESTError(
                "OrganizationalUnitNotFoundException",
                "You specified an organizational unit that doesn't exist.",
            )
        return ou

    def validate_parent_id(self, parent_id: str) -> str:
        try:
            self.get_organizational_unit_by_id(parent_id)
        except RESTError:
            raise RESTError(
                "ParentNotFoundException", "You specified parent that doesn't exist."
            )
        return parent_id

    def describe_organizational_unit(self, **kwargs: Any) -> Dict[str, Any]:
        ou = self.get_organizational_unit_by_id(kwargs["OrganizationalUnitId"])
        return ou.describe()

    @paginate(pagination_model=PAGINATION_MODEL)  # type: ignore[misc]
    def list_organizational_units_for_parent(
        self, **kwargs: Any
    ) -> List[Dict[str, Any]]:
        parent_id = self.validate_parent_id(kwargs["parent_id"])
        return [
            {"Id": ou.id, "Arn": ou.arn, "Name": ou.name}
            for ou in self.ou
            if ou.parent_id == parent_id
        ]

    def create_account(self, **kwargs: Any) -> Dict[str, Any]:
        new_account = FakeAccount(self.org, **kwargs)  # type: ignore
        self.accounts.append(new_account)
        self.attach_policy(PolicyId=utils.DEFAULT_POLICY_ID, TargetId=new_account.id)
        return new_account.create_account_status

    def close_account(self, **kwargs: Any) -> None:
        for account in self.accounts:
            if account.id == kwargs["AccountId"]:
                account.close()
                return
        raise AccountNotFoundException

    def get_account_by_id(self, account_id: str) -> FakeAccount:
        account = next(
            (account for account in self.accounts if account.id == account_id), None
        )
        if account is None:
            raise AccountNotFoundException
        return account

    def get_account_by_attr(self, attr: str, value: Any) -> FakeAccount:
        account = next(
            (
                account
                for account in self.accounts
                if hasattr(account, attr) and getattr(account, attr) == value
            ),
            None,
        )
        if account is None:
            raise AccountNotFoundException
        return account

    def describe_account(self, **kwargs: Any) -> Dict[str, Any]:
        account = self.get_account_by_id(kwargs["AccountId"])
        return dict(Account=account.describe())

    def describe_create_account_status(self, **kwargs: Any) -> Dict[str, Any]:
        account = self.get_account_by_attr(
            "create_account_status_id", kwargs["CreateAccountRequestId"]
        )
        return account.create_account_status

    def list_create_account_status(self, **kwargs: Any) -> Dict[str, Any]:
        requested_states = kwargs.get("States")
        if not requested_states:
            requested_states = ["IN_PROGRESS", "SUCCEEDED", "FAILED"]
        accountStatuses = []
        for account in self.accounts:
            create_account_status = account.create_account_status["CreateAccountStatus"]
            if create_account_status["State"] in requested_states:
                accountStatuses.append(create_account_status)
        token = kwargs.get("NextToken")
        if token:
            start = int(token)
        else:
            start = 0
        max_results = int(kwargs.get("MaxResults", 123))
        accounts_resp = accountStatuses[start : start + max_results]
        next_token = None
        if max_results and len(accountStatuses) > (start + max_results):
            next_token = str(len(accounts_resp))
        return dict(CreateAccountStatuses=accounts_resp, NextToken=next_token)

    @paginate(pagination_model=PAGINATION_MODEL)
    def list_accounts(self) -> List[FakeAccount]:
        accounts = [account.describe() for account in self.accounts]
        return sorted(accounts, key=lambda x: x["JoinedTimestamp"])  # type: ignore

    @paginate(pagination_model=PAGINATION_MODEL)  # type: ignore[misc]
    def list_accounts_for_parent(self, **kwargs: Any) -> Any:
        parent_id = self.validate_parent_id(kwargs["parent_id"])
        accounts = [
            account.describe()
            for account in self.accounts
            if account.parent_id == parent_id
        ]
        return sorted(accounts, key=lambda x: x["JoinedTimestamp"])

    def move_account(self, **kwargs: Any) -> None:
        new_parent_id = self.validate_parent_id(kwargs["DestinationParentId"])
        self.validate_parent_id(kwargs["SourceParentId"])
        account = self.get_account_by_id(kwargs["AccountId"])
        index = self.accounts.index(account)
        self.accounts[index].parent_id = new_parent_id

    def list_parents(self, **kwargs: Any) -> Dict[str, Any]:
        if re.compile(r"[0-9]{12}").match(kwargs["ChildId"]):
            child_object: Any = self.get_account_by_id(kwargs["ChildId"])
        else:
            child_object = self.get_organizational_unit_by_id(kwargs["ChildId"])
        return dict(
            Parents=[
                {"Id": ou.id, "Type": ou.type}
                for ou in self.ou
                if ou.id == child_object.parent_id
            ]
        )

    def list_children(self, **kwargs: Any) -> Dict[str, Any]:
        parent_id = self.validate_parent_id(kwargs["ParentId"])
        if kwargs["ChildType"] == "ACCOUNT":
            obj_list: List[Any] = self.accounts
        elif kwargs["ChildType"] == "ORGANIZATIONAL_UNIT":
            obj_list = self.ou
        else:
            raise InvalidInputException("You specified an invalid value.")
        return dict(
            Children=[
                {"Id": obj.id, "Type": kwargs["ChildType"]}
                for obj in obj_list
                if obj.parent_id == parent_id
            ]
        )

    def create_policy(self, **kwargs: Any) -> Dict[str, Any]:
        new_policy = FakePolicy(self.org, **kwargs)  # type: ignore
        for policy in self.policies:
            if kwargs["Name"] == policy.name:
                raise DuplicatePolicyException
        self.policies.append(new_policy)
        return new_policy.describe()

    def describe_policy(self, **kwargs: Any) -> Dict[str, Any]:
        if re.compile(utils.POLICY_ID_REGEX).match(kwargs["PolicyId"]):
            policy = next(
                (p for p in self.policies if p.id == kwargs["PolicyId"]), None
            )
            if policy is None:
                raise PolicyNotFoundException(
                    "You specified a policy that doesn't exist.",
                )
        else:
            raise InvalidInputException("You specified an invalid value.")
        return policy.describe()

    def get_policy_by_id(self, policy_id: str) -> FakePolicy:
        policy = next(
            (policy for policy in self.policies if policy.id == policy_id), None
        )
        if policy is None:
            raise PolicyNotFoundException(
                "We can't find a policy with the PolicyId that you specified.",
            )
        return policy

    def update_policy(self, **kwargs: Any) -> Dict[str, Any]:
        policy = self.get_policy_by_id(kwargs["PolicyId"])
        policy.name = kwargs.get("Name", policy.name)
        policy.description = kwargs.get("Description", policy.description)
        policy.content = kwargs.get("Content", policy.content)
        return policy.describe()

    def attach_policy(self, **kwargs: Any) -> None:
        policy = self.get_policy_by_id(kwargs["PolicyId"])
        if re.compile(utils.ROOT_ID_REGEX).match(kwargs["TargetId"]) or re.compile(
            utils.OU_ID_REGEX
        ).match(kwargs["TargetId"]):
            ou = next((ou for ou in self.ou if ou.id == kwargs["TargetId"]), None)
            if ou is not None:
                if policy not in ou.attached_policies:
                    ou.attached_policies.append(policy)
                    policy.attachments.append(ou)
            else:
                raise RESTError(
                    "OrganizationalUnitNotFoundException",
                    "You specified an organizational unit that doesn't exist.",
                )
        elif re.compile(utils.ACCOUNT_ID_REGEX).match(kwargs["TargetId"]):
            account = next(
                (a for a in self.accounts if a.id == kwargs["TargetId"]), None
            )
            if account is not None:
                if policy not in account.attached_policies:
                    account.attached_policies.append(policy)
                    policy.attachments.append(account)
            else:
                raise AccountNotFoundException
        else:
            raise InvalidInputException("You specified an invalid value.")

    def list_policies(self) -> Dict[str, Any]:
        return dict(
            Policies=[p.describe()["Policy"]["PolicySummary"] for p in self.policies]
        )

    def delete_policy(self, **kwargs: Any) -> None:
        for idx, policy in enumerate(self.policies):
            if policy.id == kwargs["PolicyId"]:
                if self.list_targets_for_policy(PolicyId=policy.id)["Targets"]:
                    raise RESTError(
                        "PolicyInUseException",
                        "The policy is attached to one or more entities. You must detach it from all roots, OUs, and accounts before performing this operation.",
                    )
                del self.policies[idx]
                return
        raise PolicyNotFoundException(
            "We can't find a policy with the PolicyId that you specified.",
        )

    def list_policies_for_target(self, **kwargs: Any) -> Dict[str, Any]:
        _filter = kwargs["Filter"]

        if re.match(utils.ROOT_ID_REGEX, kwargs["TargetId"]):
            obj: Any = next((ou for ou in self.ou if ou.id == kwargs["TargetId"]), None)
            if obj is None:
                raise TargetNotFoundException
        elif re.compile(utils.OU_ID_REGEX).match(kwargs["TargetId"]):
            obj = next((ou for ou in self.ou if ou.id == kwargs["TargetId"]), None)
            if obj is None:
                raise RESTError(
                    "OrganizationalUnitNotFoundException",
                    "You specified an organizational unit that doesn't exist.",
                )
        elif re.compile(utils.ACCOUNT_ID_REGEX).match(kwargs["TargetId"]):
            obj = next((a for a in self.accounts if a.id == kwargs["TargetId"]), None)
            if obj is None:
                raise AccountNotFoundException
        else:
            raise InvalidInputException("You specified an invalid value.")

        if not FakePolicy.supported_policy_type(_filter):
            raise InvalidInputException("You specified an invalid value.")

        if _filter not in ["AISERVICES_OPT_OUT_POLICY", "SERVICE_CONTROL_POLICY"]:
            raise NotImplementedError(
                f"The {_filter} policy type has not been implemented"
            )

        return dict(
            Policies=[
                p.describe()["Policy"]["PolicySummary"]
                for p in obj.attached_policies
                if p.type == _filter
            ]
        )

    def _get_resource_for_tagging(self, resource_id: str) -> Any:
        if utils.fullmatch(
            re.compile(utils.OU_ID_REGEX), resource_id
        ) or utils.fullmatch(utils.ROOT_ID_REGEX, resource_id):
            resource: Any = next((a for a in self.ou if a.id == resource_id), None)
        elif utils.fullmatch(re.compile(utils.ACCOUNT_ID_REGEX), resource_id):
            resource = next((a for a in self.accounts if a.id == resource_id), None)
        elif utils.fullmatch(re.compile(utils.POLICY_ID_REGEX), resource_id):
            resource = next((a for a in self.policies if a.id == resource_id), None)
        else:
            raise InvalidInputException(
                "You provided a value that does not match the required pattern."
            )

        if resource is None:
            raise TargetNotFoundException

        return resource

    def list_targets_for_policy(self, **kwargs: Any) -> Dict[str, Any]:
        if re.compile(utils.POLICY_ID_REGEX).match(kwargs["PolicyId"]):
            policy = next(
                (p for p in self.policies if p.id == kwargs["PolicyId"]), None
            )
            if policy is None:
                raise PolicyNotFoundException(
                    "You specified a policy that doesn't exist.",
                )
        else:
            raise InvalidInputException("You specified an invalid value.")
        objects = [
            {"TargetId": obj.id, "Arn": obj.arn, "Name": obj.name, "Type": obj.type}
            for obj in policy.attachments
        ]
        return dict(Targets=objects)

    def tag_resource(self, **kwargs: Any) -> None:
        resource = self._get_resource_for_tagging(kwargs["ResourceId"])
        new_tags = {tag["Key"]: tag["Value"] for tag in kwargs["Tags"]}
        resource.tags.update(new_tags)

    def list_tags_for_resource(self, **kwargs: str) -> Dict[str, Any]:
        resource = self._get_resource_for_tagging(kwargs["ResourceId"])
        tags = [{"Key": key, "Value": value} for key, value in resource.tags.items()]
        return dict(Tags=tags)

    def untag_resource(self, **kwargs: Any) -> None:
        resource = self._get_resource_for_tagging(kwargs["ResourceId"])
        for key in kwargs["TagKeys"]:
            resource.tags.pop(key, None)

    def enable_aws_service_access(self, **kwargs: str) -> None:
        service = FakeServiceAccess(**kwargs)

        # enabling an existing service results in no changes
        if any(
            service["ServicePrincipal"] == kwargs["ServicePrincipal"]
            for service in self.services
        ):
            return

        self.services.append(service.describe())

    def list_aws_service_access_for_organization(self) -> Dict[str, Any]:
        return dict(EnabledServicePrincipals=self.services)

    def disable_aws_service_access(self, **kwargs: str) -> None:
        if not FakeServiceAccess.trusted_service(kwargs["ServicePrincipal"]):
            raise InvalidInputException(
                "You specified an unrecognized service principal."
            )

        service_principal = next(
            (
                service
                for service in self.services
                if service["ServicePrincipal"] == kwargs["ServicePrincipal"]
            ),
            None,
        )

        if service_principal:
            self.services.remove(service_principal)

    def register_delegated_administrator(self, **kwargs: str) -> None:
        account_id = kwargs["AccountId"]

        if account_id == self.account_id:
            raise ConstraintViolationException(
                "You cannot register master account/yourself as delegated administrator for your organization."
            )

        account = self.get_account_by_id(account_id)

        admin = next(
            (admin for admin in self.admins if admin.account.id == account_id), None
        )
        if admin is None:
            admin = FakeDelegatedAdministrator(account)
            self.admins.append(admin)

        admin.add_service_principal(kwargs["ServicePrincipal"])

    def list_delegated_administrators(self, **kwargs: str) -> Dict[str, Any]:
        admins = self.admins
        service = kwargs.get("ServicePrincipal")

        if service:
            if not FakeDelegatedAdministrator.supported_service(service):
                raise InvalidInputException(
                    "You specified an unrecognized service principal."
                )

            admins = [admin for admin in admins if service in admin.services]

        delegated_admins = [admin.describe() for admin in admins]

        return dict(DelegatedAdministrators=delegated_admins)

    def list_delegated_services_for_account(self, **kwargs: str) -> Dict[str, Any]:
        admin = next(
            (admin for admin in self.admins if admin.account.id == kwargs["AccountId"]),
            None,
        )
        if admin is None:
            account = next(
                (
                    account
                    for account in self.accounts
                    if account.id == kwargs["AccountId"]
                ),
                None,
            )
            if account:
                raise AccountNotRegisteredException

            raise AWSOrganizationsNotInUseException

        services = [service for service in admin.services.values()]

        return dict(DelegatedServices=services)

    def deregister_delegated_administrator(self, **kwargs: str) -> None:
        account_id = kwargs["AccountId"]
        service = kwargs["ServicePrincipal"]

        if account_id == self.account_id:
            raise ConstraintViolationException(
                "You cannot register master account/yourself as delegated administrator for your organization."
            )

        admin = next(
            (admin for admin in self.admins if admin.account.id == account_id), None
        )
        if admin is None:
            account = next(
                (
                    account
                    for account in self.accounts
                    if account.id == kwargs["AccountId"]
                ),
                None,
            )
            if account:
                raise AccountNotRegisteredException

            raise AccountNotFoundException

        admin.remove_service_principal(service)

        # remove account, when no services attached
        if not admin.services:
            self.admins.remove(admin)

    def enable_policy_type(self, **kwargs: str) -> Dict[str, Any]:
        root = self._get_root_by_id(kwargs["RootId"])

        root.add_policy_type(kwargs["PolicyType"])

        return dict(Root=root.describe())

    def disable_policy_type(self, **kwargs: str) -> Dict[str, Any]:
        root = self._get_root_by_id(kwargs["RootId"])

        root.remove_policy_type(kwargs["PolicyType"])

        return dict(Root=root.describe())

    def detach_policy(self, **kwargs: str) -> None:
        policy = self.get_policy_by_id(kwargs["PolicyId"])
        root_id_regex = utils.ROOT_ID_REGEX
        ou_id_regex = utils.OU_ID_REGEX
        account_id_regex = utils.ACCOUNT_ID_REGEX
        target_id = kwargs["TargetId"]

        if re.match(root_id_regex, target_id) or re.match(ou_id_regex, target_id):
            ou = next((ou for ou in self.ou if ou.id == target_id), None)
            if ou is not None:
                if policy in ou.attached_policies:
                    ou.attached_policies.remove(policy)
                    policy.attachments.remove(ou)
            else:
                raise RESTError(
                    "OrganizationalUnitNotFoundException",
                    "You specified an organizational unit that doesn't exist.",
                )
        elif re.match(account_id_regex, target_id):
            account = next(
                (account for account in self.accounts if account.id == target_id), None
            )
            if account is not None:
                if policy in account.attached_policies:
                    account.attached_policies.remove(policy)
                    policy.attachments.remove(account)
            else:
                raise AccountNotFoundException
        else:
            raise InvalidInputException("You specified an invalid value.")

    def remove_account_from_organization(self, **kwargs: str) -> None:
        account = self.get_account_by_id(kwargs["AccountId"])
        for policy in account.attached_policies:
            policy.attachments.remove(account)
        self.accounts.remove(account)


organizations_backends = BackendDict(
    OrganizationsBackend,
    "organizations",
    use_boto3_regions=False,
    additional_regions=PARTITION_NAMES,
)
