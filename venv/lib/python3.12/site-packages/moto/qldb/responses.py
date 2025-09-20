"""Handles incoming qldb requests, invokes methods, returns responses."""

import json
from datetime import datetime
from typing import Dict, Optional, Union

from moto.core.responses import BaseResponse

from .models import QLDBBackend, qldb_backends


class QLDBResponse(BaseResponse):
    """Handler for QLDB requests and responses."""

    def __init__(self) -> None:
        super().__init__(service_name="qldb")

    @property
    def qldb_backend(self) -> QLDBBackend:
        """Return backend instance specific for this region."""
        return qldb_backends[self.current_account][self.region]

    def _format_encryption_description(
        self, encryption_description: Optional[Dict[str, Union[str, datetime]]]
    ) -> Dict[str, Optional[Union[str, datetime]]]:
        return {
            "KmsKeyArn": (encryption_description or {}).get("kms_key_arn"),
            "EncryptionStatus": (encryption_description or {}).get("encryption_status"),
            "InaccessibleKmsKeyDateTime": (encryption_description or {}).get(
                "inaccessible_kms_key_date_time"
            ),
        }

    def describe_ledger(self) -> str:
        name = self._get_param("name")
        (
            name,
            arn,
            state,
            creation_date_time,
            permissions_mode,
            deletion_protection,
            encryption_description,
        ) = self.qldb_backend.describe_ledger(
            name=name,
        )
        creation_date_time_string = (
            creation_date_time.strftime("%d/%m/%Y, %H:%M:%S")
            if creation_date_time is not None
            else ""
        )
        return json.dumps(
            dict(
                Name=name,
                Arn=arn,
                State=state,
                CreationDateTime=creation_date_time_string,
                PermissionsMode=permissions_mode,
                DeletionProtection=deletion_protection,
                EncryptionDescription=self._format_encryption_description(
                    encryption_description
                ),
            )
        )

    def create_ledger(self) -> str:
        params = json.loads(self.body)
        name = params.get("Name")
        tags = params.get("Tags")
        permissions_mode = params.get("PermissionsMode")
        deletion_protection = params.get("DeletionProtection")
        kms_key = params.get("KmsKey")
        (
            name,
            arn,
            state,
            creation_date_time,
            permissions_mode,
            deletion_protection,
            kms_key_arn,
        ) = self.qldb_backend.create_ledger(
            name=name,
            tags=tags,
            permissions_mode=permissions_mode,
            deletion_protection=deletion_protection,
            kms_key=kms_key,
        )
        return json.dumps(
            dict(
                Name=name,
                Arn=arn,
                State=state,
                CreationDateTime=creation_date_time.strftime("%d/%m/%Y, %H:%M:%S"),
                PermissionsMode=permissions_mode,
                DeletionProtection=deletion_protection,
                KmsKeyArn=kms_key_arn,
            )
        )

    def delete_ledger(self) -> str:
        name = self._get_param("name")
        self.qldb_backend.delete_ledger(
            name=name,
        )
        return json.dumps(dict())

    def update_ledger(self) -> str:
        name = self._get_param("name")
        params = json.loads(self.body)
        deletion_protection = params.get("DeletionProtection")
        kms_key = params.get("KmsKey")
        (
            name,
            arn,
            state,
            creation_date_time,
            deletion_protection,
            encryption_description,
        ) = self.qldb_backend.update_ledger(
            name=name,
            deletion_protection=deletion_protection,
            kms_key=kms_key,
        )
        creation_date_time_string = (
            creation_date_time.strftime("%d/%m/%Y, %H:%M:%S")
            if creation_date_time is not None
            else ""
        )
        return json.dumps(
            dict(
                name=name,
                arn=arn,
                state=state,
                creationDateTime=creation_date_time_string,
                deletionProtection=deletion_protection,
                encryptionDescription=self._format_encryption_description(
                    encryption_description
                ),
            )
        )

    def tag_resource(self) -> str:
        params = json.loads(self.body)
        resource_arn = params.get("ResourceArn")
        tags = params.get("Tags")
        self.qldb_backend.tag_resource(
            resource_arn=resource_arn,
            new_tags=tags,
        )
        return json.dumps(dict())

    def list_tags_for_resource(self) -> str:
        resource_arn = self._get_param("ResourceArn")
        tags = self.qldb_backend.list_tags_for_resource(
            resource_arn=resource_arn,
        )
        return json.dumps(dict(Tags=tags))
