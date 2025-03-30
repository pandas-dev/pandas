from copy import deepcopy
from typing import Any, Dict, List, Optional, Tuple

from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.common_models import BaseModel
from moto.core.utils import unix_time
from moto.moto_api._internal import mock_random
from moto.utilities.tagging_service import TaggingService
from moto.utilities.utils import get_partition

from .exceptions import AlreadyExistsException, ResourceNotFoundException


class Plan(BaseModel):
    def __init__(
        self,
        backup_plan: Dict[str, Any],
        creator_request_id: str,
        backend: "BackupBackend",
    ):
        self.backup_plan_id = str(mock_random.uuid4())
        self.backup_plan_arn = f"arn:{get_partition(backend.region_name)}:backup:{backend.region_name}:{backend.account_id}:backup-plan:{self.backup_plan_id}"
        self.creation_date = unix_time()
        ran_str = mock_random.get_random_string(length=48)
        self.version_id = ran_str
        self.creator_request_id = creator_request_id
        self.backup_plan = backup_plan
        adv_settings = backup_plan.get("AdvancedBackupSettings")
        self.advanced_backup_settings = adv_settings or []
        self.deletion_date: Optional[float] = None
        # Deletion Date is updated when the backup_plan is deleted
        self.last_execution_date = None  # start_restore_job not yet supported
        rules = backup_plan["Rules"]
        for rule in rules:
            rule["ScheduleExpression"] = rule.get(
                "ScheduleExpression", "cron(0 5 ? * * *)"
            )  # Default CRON expression in UTC
            rule["StartWindowMinutes"] = rule.get(
                "StartWindowMinutes", 480
            )  # Default=480
            rule["CompletionWindowMinutes"] = rule.get(
                "CompletionWindowMinutes", 10080
            )  # Default=10080
            rule["ScheduleExpressionTimezone"] = rule.get(
                "ScheduleExpressionTimezone", "Etc/UTC"
            )  # set to Etc/UTc by default
            rule["RuleId"] = str(mock_random.uuid4())

    def to_dict(self) -> Dict[str, Any]:
        dct = {
            "BackupPlanId": self.backup_plan_id,
            "BackupPlanArn": self.backup_plan_arn,
            "CreationDate": self.creation_date,
            "VersionId": self.version_id,
            "AdvancedBackupSettings": self.advanced_backup_settings,
        }
        return {k: v for k, v in dct.items() if v}

    def to_get_dict(self) -> Dict[str, Any]:
        dct = self.to_dict()
        dct_options = {
            "BackupPlan": self.backup_plan,
            "CreatorRequestId": self.creator_request_id,
            "DeletionDate": self.deletion_date,
            "LastExecutionDate": self.last_execution_date,
        }
        for key, value in dct_options.items():
            if value is not None:
                dct[key] = value
        return dct

    def to_list_dict(self) -> Dict[str, Any]:
        dct = self.to_get_dict()
        dct.pop("BackupPlan")
        dct["BackupPlanName"] = self.backup_plan.get("BackupPlanName")
        return dct


class Vault(BaseModel):
    def __init__(
        self,
        backup_vault_name: str,
        encryption_key_arn: str,
        creator_request_id: str,
        backend: "BackupBackend",
    ):
        self.backup_vault_name = backup_vault_name
        self.backup_vault_arn = f"arn:{get_partition(backend.region_name)}:backup:{backend.region_name}:{backend.account_id}:backup-vault:{backup_vault_name}"
        self.creation_date = unix_time()
        self.encryption_key_arn = encryption_key_arn
        self.creator_request_id = creator_request_id
        self.num_of_recovery_points = 0  # start_backup_job not yet supported
        self.locked = False  # put_backup_vault_lock_configuration
        self.min_retention_days = 0  # put_backup_vault_lock_configuration
        self.max_retention_days = 0  # put_backup_vault_lock_configuration
        self.lock_date = None  # put_backup_vault_lock_configuration

    def to_dict(self) -> Dict[str, Any]:
        dct = {
            "BackupVaultName": self.backup_vault_name,
            "BackupVaultArn": self.backup_vault_arn,
            "CreationDate": self.creation_date,
        }
        return dct

    def to_list_dict(self) -> Dict[str, Any]:
        dct = self.to_dict()
        dct_options: Dict[str, Any] = dict()
        dct_options = {
            "EncryptionKeyArn": self.encryption_key_arn,
            "CreatorRequestId": self.creator_request_id,
            "NumberOfRecoveryPoints": self.num_of_recovery_points,
            "Locked": self.locked,
            "MinRetentionDays": self.min_retention_days,
            "MaxRetentionDays": self.max_retention_days,
            "LockDate": self.lock_date,
        }
        for key, value in dct_options.items():
            if value is not None:
                dct[key] = value
        return dct


class BackupBackend(BaseBackend):
    """Implementation of Backup APIs."""

    def __init__(self, region_name: str, account_id: str):
        super().__init__(region_name, account_id)

        self.vaults: Dict[str, Vault] = dict()
        self.plans: Dict[str, Plan] = dict()
        self.tagger = TaggingService()

    def create_backup_plan(
        self,
        backup_plan: Dict[str, Any],
        backup_plan_tags: Dict[str, str],
        creator_request_id: str,
    ) -> Plan:
        if backup_plan["BackupPlanName"] in list(
            p.backup_plan["BackupPlanName"] for p in list(self.plans.values())
        ):
            raise AlreadyExistsException(
                msg="Backup plan with the same plan document already exists"
            )
        plan = Plan(
            backup_plan=backup_plan,
            creator_request_id=creator_request_id,
            backend=self,
        )
        if backup_plan_tags:
            self.tag_resource(plan.backup_plan_arn, backup_plan_tags)
        self.plans[plan.backup_plan_id] = plan
        return plan

    def get_backup_plan(self, backup_plan_id: str, version_id: Optional[Any]) -> Plan:
        msg = "Failed reading Backup plan with provided version"
        if backup_plan_id not in self.plans:
            raise ResourceNotFoundException(msg=msg)
        plan = self.plans[backup_plan_id]
        if version_id:
            if plan.version_id == version_id:
                return plan
            else:
                raise ResourceNotFoundException(msg=msg)
        return plan

    def delete_backup_plan(self, backup_plan_id: str) -> Tuple[str, str, float, str]:
        if backup_plan_id not in self.plans:
            raise ResourceNotFoundException(
                msg="Failed reading Backup plan with provided version"
            )
        deletion_date = unix_time()
        res = self.plans[backup_plan_id]
        res.deletion_date = deletion_date
        return res.backup_plan_id, res.backup_plan_arn, deletion_date, res.version_id

    def list_backup_plans(self, include_deleted: Any) -> List[Plan]:
        """
        Pagination is not yet implemented
        """
        plans_list = deepcopy(self.plans)

        for plan in list(plans_list.values()):
            backup_plan_id = plan.backup_plan_id
            if plan.deletion_date is not None:
                plans_list.pop(backup_plan_id)
        if include_deleted:
            return list(self.plans.values())
        return list(plans_list.values())

    def create_backup_vault(
        self,
        backup_vault_name: str,
        backup_vault_tags: Dict[str, str],
        encryption_key_arn: str,
        creator_request_id: str,
    ) -> Vault:
        if backup_vault_name in self.vaults:
            raise AlreadyExistsException(
                msg="Backup vault with the same name already exists"
            )
        vault = Vault(
            backup_vault_name=backup_vault_name,
            encryption_key_arn=encryption_key_arn,
            creator_request_id=creator_request_id,
            backend=self,
        )
        if backup_vault_tags:
            self.tag_resource(vault.backup_vault_arn, backup_vault_tags)
        self.vaults[backup_vault_name] = vault
        return vault

    def list_backup_vaults(self) -> List[Vault]:
        """
        Pagination is not yet implemented
        """
        return list(self.vaults.values())

    def list_tags(self, resource_arn: str) -> Dict[str, str]:
        """
        Pagination is not yet implemented
        """
        return self.tagger.get_tag_dict_for_resource(resource_arn)

    def tag_resource(self, resource_arn: str, tags: Dict[str, str]) -> None:
        tags_input = TaggingService.convert_dict_to_tags_input(tags or {})
        self.tagger.tag_resource(resource_arn, tags_input)

    def untag_resource(self, resource_arn: str, tag_key_list: List[str]) -> None:
        self.tagger.untag_resource_using_names(resource_arn, tag_key_list)


backup_backends = BackendDict(BackupBackend, "backup")
