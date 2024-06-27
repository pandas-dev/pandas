"""Handles incoming backup requests, invokes methods, returns responses."""

import json
from urllib.parse import unquote

from moto.core.responses import BaseResponse

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
        backup_plan_id = self.path.split("/")[-2]
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
            dict(
                BackupPlanId=backup_plan_id,
                BackupPlanArn=backup_plan_arn,
                DeletionDate=deletion_date,
                VersionId=version_id,
            )
        )

    def list_backup_plans(self) -> str:
        params = self._get_params()
        include_deleted = params.get("includeDeleted")
        backup_plans_list = self.backup_backend.list_backup_plans(
            include_deleted=include_deleted
        )
        return json.dumps(
            dict(BackupPlansList=[p.to_list_dict() for p in backup_plans_list])
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

    def list_backup_vaults(self) -> str:
        backup_vault_list = self.backup_backend.list_backup_vaults()
        return json.dumps(
            dict(BackupVaultList=[v.to_list_dict() for v in backup_vault_list])
        )

    def list_tags(self) -> str:
        resource_arn = unquote(self.path.split("/")[-2])
        tags = self.backup_backend.list_tags(
            resource_arn=resource_arn,
        )
        return json.dumps(dict(Tags=tags))

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
