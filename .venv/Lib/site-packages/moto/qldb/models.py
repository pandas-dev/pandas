"""QLDBBackend class with methods for supported APIs."""

from datetime import datetime
from typing import Any, Dict, Optional, Tuple, Union

from moto.core.base_backend import BackendDict, BaseBackend
from moto.qldb.exceptions import (
    LedgerNameAlreadyTakenException,
    LedgerNotFoundException,
    ResourceNotFoundException,
)


class QLDBBackend(BaseBackend):
    """Implementation of QLDB APIs."""

    def __init__(self, region_name: str, account_id: str) -> None:
        super().__init__(region_name, account_id)
        self.ledgers: Dict[str, Any] = dict()  # type ignore[misc]

    def describe_ledger(
        self, name: str
    ) -> Tuple[
        str,
        Optional[str],
        Optional[str],
        Optional[datetime],
        Optional[str],
        Optional[bool],
        Optional[Dict[str, Union[str, datetime]]],
    ]:
        if name not in self.ledgers:
            raise LedgerNotFoundException(name)
        ledger = self.ledgers[name]
        return (
            name,
            ledger.get("arn"),
            ledger.get("state"),
            ledger.get("creation_date_time"),
            ledger.get("permissions_mode"),
            ledger.get("deletion_protection"),
            ledger.get("encryption_description"),
        )

    def _get_kms_key_arn(self, key: str) -> str:
        return f"arn:aws:kms:us-east-1:123456789012:key/{key}"

    def create_ledger(
        self,
        name: str,
        tags: Optional[Dict[str, str]],
        permissions_mode: str,
        deletion_protection: Optional[bool],
        kms_key: Optional[str],
    ) -> Tuple[
        str,
        str,
        str,
        datetime,
        str,
        Optional[bool],
        str,
    ]:
        if name in self.ledgers:
            raise LedgerNameAlreadyTakenException(name)
        kms_key_arn = self._get_kms_key_arn(kms_key or name)
        encryption_description = {
            "kms_key_arn": kms_key_arn,
            "encryption_status": "ENABLED",
            "inaccessible_kms_key_date_time": datetime.now().strftime(
                "%d/%m/%Y, %H:%M:%S"
            ),
        }
        arn = f"arn:aws:qldb:us-east-1:123456789012:ledger/{name}"
        creation_date_time = datetime.now()
        state = "ACTIVE"
        ledger = {
            "name": name,
            "arn": arn,
            "state": state,
            "creation_date_time": creation_date_time,
            "permissions_mode": permissions_mode,
            "deletion_protection": deletion_protection,
            "encryption_description": encryption_description,
            "tags": tags,
        }
        self.ledgers[name] = ledger
        return (
            name,
            arn,
            state,
            creation_date_time,
            permissions_mode,
            deletion_protection,
            kms_key_arn,
        )

    def delete_ledger(self, name: str) -> None:
        if name not in self.ledgers:
            raise LedgerNotFoundException(name)
        del self.ledgers[name]
        return

    def update_ledger(
        self, name: str, deletion_protection: Optional[bool], kms_key: Optional[str]
    ) -> Tuple[
        str,
        Optional[str],
        Optional[str],
        Optional[datetime],
        Optional[bool],
        Optional[Dict[str, Union[str, datetime]]],
    ]:
        if name not in self.ledgers:
            raise LedgerNotFoundException(name)
        if deletion_protection is not None:
            self.ledgers[name]["deletion_protection"] = deletion_protection
        if kms_key and "encryption_description" in self.ledgers[name]:
            self.ledgers[name]["encryption_description"]["kms_key_arn"] = (
                self._get_kms_key_arn(kms_key)
            )

        ledger = self.ledgers[name]

        return (
            name,
            ledger.get("arn"),
            ledger.get("state"),
            ledger.get("creation_date_time"),
            ledger.get("deletion_protection"),
            ledger.get("encryption_description"),
        )

    def _get_ledger_by_resource_arn(self, resource_arn: str) -> Dict[str, Any]:
        for ledger in self.ledgers.values():
            if ledger.get("resource_arn") == resource_arn:
                return ledger
        raise ResourceNotFoundException(resource_arn)

    def tag_resource(self, resource_arn: str, new_tags: Dict[str, str]) -> None:
        ledger = self._get_ledger_by_resource_arn(resource_arn)
        tags = ledger.get("tags") or {}
        tags.update(new_tags)
        ledger["tags"] = tags
        return

    def list_tags_for_resource(self, resource_arn: str) -> Dict[str, str]:
        ledger = self._get_ledger_by_resource_arn(resource_arn)
        return ledger.get("tags") or {}


qldb_backends = BackendDict(QLDBBackend, "qldb")
