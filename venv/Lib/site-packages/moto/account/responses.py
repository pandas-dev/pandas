"""Handles incoming account requests, invokes methods, returns responses."""

from moto.core.responses import ActionResult, BaseResponse, EmptyResult

from .exceptions import UnknownContactType
from .models import ALLOWED_CONTACT_TYPES, account_backends


class AccountResponse(BaseResponse):
    def __init__(self) -> None:
        super().__init__(service_name="account")

    def put_alternate_contact(self) -> ActionResult:
        account_id = self._get_account_id()
        alternate_contact_type = self._get_contact_type(account_id)
        email_address = self._get_param("EmailAddress")
        name = self._get_param("Name")
        phone_number = self._get_param("PhoneNumber")
        title = self._get_param("Title")

        backend = account_backends[account_id or self.current_account][self.partition]

        backend.put_alternate_contact(
            alternate_contact_type=alternate_contact_type,
            email_address=email_address,
            name=name,
            phone_number=phone_number,
            title=title,
        )
        return EmptyResult()

    def get_alternate_contact(self) -> ActionResult:
        account_id = self._get_account_id()
        alternate_contact_type = self._get_contact_type(account_id)

        backend = account_backends[account_id][self.partition]

        contact = backend.get_alternate_contact(
            alternate_contact_type=alternate_contact_type
        )
        return ActionResult(result={"AlternateContact": contact})

    def delete_alternate_contact(self) -> EmptyResult:
        account_id = self._get_account_id()
        alternate_contact_type = self._get_contact_type(account_id)

        backend = account_backends[account_id][self.partition]

        backend.delete_alternate_contact(alternate_contact_type=alternate_contact_type)
        return EmptyResult()

    def _get_account_id(self) -> str:
        return self._get_param("AccountId") or self.current_account

    def _get_contact_type(self, account_id: str) -> str:
        alternate_contact_type = self._get_param("AlternateContactType")
        if alternate_contact_type not in ALLOWED_CONTACT_TYPES:
            from moto.sts.models import STSBackend, sts_backends

            access_key_id = self.get_access_key()
            sts_backend: STSBackend = sts_backends[account_id][self.partition]
            _, user_arn, _ = sts_backend.get_caller_identity(access_key_id, self.region)
            raise UnknownContactType(user_arn, action="GetAlternateContact")
        return alternate_contact_type
