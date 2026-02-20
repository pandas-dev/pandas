"""Handles incoming backup requests, invokes methods, returns responses."""

import json
from urllib.parse import unquote

from moto.core.responses import ActionResult, BaseResponse, EmptyResult

from .models import BackupBackend, backup_backends


class BackupResponse(BaseResponse):
    """Handler for Backup requests and responses."""

    def __init__(self) -> None:
        super().__init__(service_name="backup")

    @property
    def backup_backend(self) -> BackupBackend:
        """Return backend instance specific for this region."""
        return backup_backends[self.current_account][self.region]

    def create_backup_plan(self) -> str:
        params = json.loads(self.body)
        backup_plan = params.get("BackupPlan")
        backup_plan_tags = params.get("BackupPlanTags")
        creator_request_id = params.get("CreatorRequestId")
        plan = self.backup_backend.create_backup_plan(
            backup_plan=backup_plan,
            backup_plan_tags=backup_plan_tags,
            creator_request_id=creator_request_id,
        )
        return json.dumps(dict(plan.to_dict()))

    def get_backup_plan(self) -> str:
        params = self._get_params()
        backup_plan_id = self.path.split("/plans/")[-1]
        backup_plan_id = backup_plan_id.replace("/", "")  # replace any trailing slash
        version_id = params.get("versionId")
        plan = self.backup_backend.get_backup_plan(
            backup_plan_id=backup_plan_id, version_id=version_id
        )
        return json.dumps(dict(plan.to_get_dict()))

    def delete_backup_plan(self) -> str:
        backup_plan_id = self.path.split("/")[-1]
        (
            backup_plan_id,
            backup_plan_arn,
            deletion_date,
            version_id,
        ) = self.backup_backend.delete_backup_plan(
            backup_plan_id=backup_plan_id,
        )
        return json.dumps(
            {
                "BackupPlanId": backup_plan_id,
                "BackupPlanArn": backup_plan_arn,
                "DeletionDate": deletion_date,
                "VersionId": version_id,
            }
        )

    def list_backup_plans(self) -> str:
        params = self._get_params()
        include_deleted = params.get("includeDeleted")
        backup_plans_list = self.backup_backend.list_backup_plans(
            include_deleted=include_deleted
        )
        return json.dumps(
            {"BackupPlansList": [p.to_list_dict() for p in backup_plans_list]}
        )

    def create_backup_vault(self) -> str:
        params = json.loads(self.body)
        backup_vault_name = self.path.split("/")[-1]
        backup_vault_tags = params.get("BackupVaultTags")
        encryption_key_arn = params.get("EncryptionKeyArn")
        creator_request_id = params.get("CreatorRequestId")
        backup_vault = self.backup_backend.create_backup_vault(
            backup_vault_name=backup_vault_name,
            backup_vault_tags=backup_vault_tags,
            encryption_key_arn=encryption_key_arn,
            creator_request_id=creator_request_id,
        )
        return json.dumps(dict(backup_vault.to_dict()))

    def delete_backup_vault(self) -> EmptyResult:
        backup_vault_name = self.path.split("/")[-1]
        self.backup_backend.delete_backup_vault(backup_vault_name)
        return EmptyResult()

    def describe_backup_vault(self) -> ActionResult:
        backup_vault_name = self.path.split("/")[-1]
        vault = self.backup_backend.describe_backup_vault(backup_vault_name)
        return ActionResult(result=vault)

    def list_backup_vaults(self) -> str:
        backup_vault_list = self.backup_backend.list_backup_vaults()
        return json.dumps(
            {"BackupVaultList": [v.to_list_dict() for v in backup_vault_list]}
        )

    def list_tags(self) -> str:
        resource_arn = unquote(self.path.split("/")[-2])
        tags = self.backup_backend.list_tags(
            resource_arn=resource_arn,
        )
        return json.dumps({"Tags": tags})

    def tag_resource(self) -> str:
        params = json.loads(self.body)
        resource_arn = unquote(self.path.split("/")[-1])
        tags = params.get("Tags")
        self.backup_backend.tag_resource(
            resource_arn=resource_arn,
            tags=tags,
        )
        return "{}"

    def untag_resource(self) -> str:
        params = json.loads(self.body)
        resource_arn = unquote(self.path.split("/")[-1])
        tag_key_list = params.get("TagKeyList")
        self.backup_backend.untag_resource(
            resource_arn=resource_arn,
            tag_key_list=tag_key_list,
        )
        return "{}"

    def put_backup_vault_lock_configuration(self) -> str:
        backup_vault_name = self.path.split("/")[-2]
        params = json.loads(self.body) if self.body else {}
        min_retention_days = params.get("MinRetentionDays")
        max_retention_days = params.get("MaxRetentionDays")
        changeable_for_days = params.get("ChangeableForDays")

        self.backup_backend.put_backup_vault_lock_configuration(
            backup_vault_name=backup_vault_name,
            min_retention_days=min_retention_days,
            max_retention_days=max_retention_days,
            changeable_for_days=changeable_for_days,
        )

        return "{}"

    def delete_backup_vault_lock_configuration(self) -> str:
        backup_vault_name = self.path.split("/")[-2]

        self.backup_backend.delete_backup_vault_lock_configuration(
            backup_vault_name=backup_vault_name,
        )

        return "{}"

    def list_report_plans(self) -> ActionResult:
        report_plans = self.backup_backend.list_report_plans()
        return ActionResult(result={"ReportPlans": report_plans})

    def create_report_plan(self) -> ActionResult:
        report_plan_name = self._get_param("ReportPlanName")
        report_plan_description = self._get_param("ReportPlanDescription")
        report_delivery_channel = self._get_param("ReportDeliveryChannel")
        report_setting = self._get_param("ReportSetting")
        report_plan = self.backup_backend.create_report_plan(
            report_plan_name=report_plan_name,
            report_plan_description=report_plan_description,
            report_delivery_channel=report_delivery_channel,
            report_setting=report_setting,
        )
        return ActionResult(result=report_plan)

    def describe_report_plan(self) -> ActionResult:
        report_plan_name = self._get_param("reportPlanName")
        report_plan = self.backup_backend.describe_report_plan(
            report_plan_name=report_plan_name
        )
        return ActionResult(result={"ReportPlan": report_plan})

    def delete_report_plan(self) -> EmptyResult:
        report_plan_name = self.path.split("/report-plans/")[-1]
        plan_name = report_plan_name.replace("/", "")  # replace any trailing slash
        self.backup_backend.delete_report_plan(report_plan_name=plan_name)
        return EmptyResult()
