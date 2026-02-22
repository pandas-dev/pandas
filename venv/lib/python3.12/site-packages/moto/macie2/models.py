from typing import Any, Optional

from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.common_models import BaseModel
from moto.core.utils import utcnow
from moto.moto_api._internal import mock_random

from .exceptions import ResourceNotFoundException


class Invitation(BaseModel):
    def __init__(self, account_id: str, region_name: str, admin_account_id: str):
        self.account_id = account_id
        self.invitation_id = mock_random.get_random_hex()
        self.invited_at = utcnow()
        self.relationship_status = "Invited"
        self.arn = f"arn:aws:macie2:{region_name}:{admin_account_id}:invitation/{self.invitation_id}"

    def to_json(self) -> dict[str, Any]:
        return {
            "accountId": self.account_id,
            "invitationId": self.invitation_id,
            "invitedAt": self.invited_at.isoformat(),
            "relationshipStatus": self.relationship_status,
        }


class Member(BaseModel):
    def __init__(
        self,
        account_id: str,
        region_name: str,
        admin_account_id: str,
        invitation: Invitation,
    ):
        self.account_id = account_id
        self.relationship_status = "Enabled"
        self.updated_at = utcnow()
        self.arn = (
            f"arn:aws:macie2:{region_name}:{admin_account_id}:member/{self.account_id}"
        )
        self.administrator_account_id = admin_account_id
        self.invited_at = invitation.invited_at

    def to_json(self) -> dict[str, Any]:
        return {
            "accountId": self.account_id,
            "administratorAccountId": self.administrator_account_id,
            "arn": self.arn,
            "email": f"{self.account_id}@example.com",  # Moto does not have access to real emails
            "invitedAt": self.invited_at.isoformat(),
            "masterAccountId": self.administrator_account_id,  # masterAccountId is deprecated but still present
            "relationshipStatus": self.relationship_status,
            "updatedAt": self.updated_at.isoformat(),
        }


class MacieBackend(BaseBackend):
    def __init__(self, region_name: str, account_id: str):
        super().__init__(region_name, account_id)
        self.invitations: dict[str, Invitation] = {}
        self.members: dict[str, Member] = {}
        self.administrator_account: Optional[Member] = None
        self.organization_admin_account_id: Optional[str] = None
        self.macie_session: Optional[dict[str, Any]] = {
            "createdAt": utcnow(),
            "findingPublishingFrequency": "FIFTEEN_MINUTES",
            "serviceRole": f"arn:aws:iam::{account_id}:role/aws-service-role/macie.amazonaws.com/AWSServiceRoleForAmazonMacie",
            "status": "ENABLED",
            "updatedAt": utcnow(),
        }

    def create_invitations(self, account_ids: list[str]) -> None:
        for account_id in account_ids:
            invitation = Invitation(
                account_id=account_id,
                region_name=self.region_name,
                admin_account_id=self.account_id,
            )
            self.invitations[account_id] = invitation

    def list_invitations(self) -> list[Invitation]:
        return list(self.invitations.values())

    def decline_invitations(self, account_ids: list[str]) -> None:
        for account_id in account_ids:
            for backend_dict in macie2_backends.values():
                backend = backend_dict.get(self.region_name)
                if backend and account_id in backend.invitations:
                    backend.invitations.pop(account_id)

    def accept_invitation(
        self, administrator_account_id: str, invitation_id: str
    ) -> None:
        # Loop through every account backend to find the matching invitation
        # Invitations can only be sent to the same region
        # If the region does not exist, then the user hasn't used Macie in this account yet, and we can safely skip it
        for backend_dict in macie2_backends.values():
            backend = backend_dict.get(self.region_name)
            if not backend:
                continue

            invitation = backend.invitations.get(self.account_id)
            if invitation and invitation.invitation_id == invitation_id:
                # Create member relationship
                member = Member(
                    account_id=self.account_id,
                    region_name=self.region_name,
                    admin_account_id=administrator_account_id,
                    invitation=invitation,
                )
                backend.members[self.account_id] = member
                self.administrator_account = member
                backend.invitations.pop(self.account_id)
                return

    def list_members(self) -> list[Member]:
        return list(self.members.values())

    def get_administrator_account(self) -> Optional[Member]:
        return self.administrator_account

    def delete_member(self, member_account_id: str) -> None:
        if member_account_id not in self.members:
            raise ResourceNotFoundException(
                "The request failed because the resource doesn't exist."
            )

        self.members.pop(member_account_id)

        for backend_dict in macie2_backends.values():
            backend = backend_dict.get(self.region_name)
            if backend and backend.account_id == member_account_id:
                backend.administrator_account = None
                break

    def disassociate_member(self, member_account_id: str) -> None:
        if member_account_id not in self.members:
            raise ResourceNotFoundException(
                "The request failed because the resource doesn't exist."
            )

        self.members.pop(member_account_id)

        for backend_dict in macie2_backends.values():
            backend = backend_dict.get(self.region_name)
            if backend and backend.account_id == member_account_id:
                backend.administrator_account = None
                break

    def get_macie_session(self) -> dict[str, Any]:
        if not self.macie_session:
            raise ResourceNotFoundException(
                "The request failed because the specified resource doesn't exist."
            )
        return {
            "createdAt": self.macie_session["createdAt"].isoformat(),
            "findingPublishingFrequency": self.macie_session[
                "findingPublishingFrequency"
            ],
            "serviceRole": self.macie_session["serviceRole"],
            "status": self.macie_session["status"],
            "updatedAt": self.macie_session["updatedAt"].isoformat(),
        }

    def enable_macie(
        self,
        finding_publishing_frequency: str = "FIFTEEN_MINUTES",
        status: str = "ENABLED",
    ) -> None:
        now = utcnow()
        self.macie_session = {
            "createdAt": now,
            "findingPublishingFrequency": finding_publishing_frequency,
            "serviceRole": f"arn:aws:iam::{self.account_id}:role/aws-service-role/macie.amazonaws.com/AWSServiceRoleForAmazonMacie",
            "status": status,
            "updatedAt": now,
        }

    def enable_organization_admin_account(self, admin_account_id: str) -> None:
        self.organization_admin_account_id = admin_account_id

    def list_organization_admin_accounts(self) -> list[dict[str, str]]:
        if self.organization_admin_account_id:
            return [
                {
                    "accountId": self.organization_admin_account_id,
                    "status": "ENABLED",
                }
            ]
        return []

    def disable_macie(self) -> None:
        self.invitations.clear()
        self.members.clear()
        self.administrator_account = None
        self.macie_session = None


macie2_backends = BackendDict(MacieBackend, "macie2")
