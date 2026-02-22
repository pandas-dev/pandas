import json

from moto.core.responses import BaseResponse

from .models import MacieBackend, macie2_backends


class MacieResponse(BaseResponse):
    """Handles API requests for the Macie service."""

    def __init__(self) -> None:
        super().__init__(service_name="macie2")

    @property
    def macie_backend(self) -> MacieBackend:
        return macie2_backends[self.current_account][self.region]

    def create_invitations(self) -> str:
        account_ids = self._get_param("accountIds", [])
        self.macie_backend.create_invitations(account_ids=account_ids)
        return json.dumps({"unprocessedAccounts": []})

    def list_invitations(self) -> str:
        invitations = self.macie_backend.list_invitations()
        return json.dumps({"invitations": [inv.to_json() for inv in invitations]})

    def decline_invitations(self) -> str:
        account_ids = self._get_param("accountIds", [])
        self.macie_backend.decline_invitations(account_ids=account_ids)
        return json.dumps({"unprocessedAccounts": []})

    def accept_invitation(self) -> str:
        administrator_account_id = self._get_param("administratorAccountId")
        invitation_id = self._get_param("invitationId")
        self.macie_backend.accept_invitation(administrator_account_id, invitation_id)
        return "{}"

    def list_members(self) -> str:
        members = self.macie_backend.list_members()
        return json.dumps({"members": [m.to_json() for m in members]})

    def get_administrator_account(self) -> str:
        admin = self.macie_backend.get_administrator_account()
        if admin:
            response = {
                "administrator": {
                    "accountId": admin.administrator_account_id,
                    "relationshipStatus": "Enabled",
                }
            }
        else:
            # When there's no administrator, return an empty dict
            response = {}
        return json.dumps(response)

    def delete_member(self) -> str:
        member_account_id = self.path.split("/")[-1]
        self.macie_backend.delete_member(member_account_id)
        return "{}"

    def disassociate_member(self) -> str:
        member_account_id = self.path.split("/")[-1]
        self.macie_backend.disassociate_member(member_account_id)
        return "{}"

    def get_macie_session(self) -> str:
        session = self.macie_backend.get_macie_session()
        return json.dumps(session)

    def enable_macie(self) -> str:
        finding_publishing_frequency = self._get_param(
            "findingPublishingFrequency", "FIFTEEN_MINUTES"
        )
        status = self._get_param("status", "ENABLED")
        self.macie_backend.enable_macie(finding_publishing_frequency, status)
        return "{}"

    def disable_macie(self) -> str:
        self.macie_backend.disable_macie()
        return json.dumps({})

    def enable_organization_admin_account(self) -> str:
        admin_account_id = self._get_param("adminAccountId")
        self.macie_backend.enable_organization_admin_account(
            admin_account_id=admin_account_id
        )
        return json.dumps({})

    def list_organization_admin_accounts(self) -> str:
        admin_accounts = self.macie_backend.list_organization_admin_accounts()
        return json.dumps({"adminAccounts": admin_accounts})
