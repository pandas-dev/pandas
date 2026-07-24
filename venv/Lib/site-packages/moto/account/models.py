"""AccountBackend class with methods for supported APIs."""

from dataclasses import dataclass

from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.common_models import BaseModel
from moto.utilities.utils import PARTITION_NAMES

from .exceptions import UnspecifiedContactType

ALLOWED_CONTACT_TYPES = ["SECURITY", "OPERATIONS", "BILLING"]


@dataclass
class AlternateContact(BaseModel):
    alternate_contact_type: str
    title: str
    name: str
    email_address: str
    phone_number: str


class AccountBackend(BaseBackend):
    def __init__(self, region_name: str, account_id: str):
        super().__init__(region_name, account_id)
        self._alternate_contacts: dict[str, AlternateContact] = {}

    def put_alternate_contact(
        self,
        alternate_contact_type: str,
        email_address: str,
        name: str,
        phone_number: str,
        title: str,
    ) -> None:
        self._alternate_contacts[alternate_contact_type] = AlternateContact(
            alternate_contact_type=alternate_contact_type,
            name=name,
            title=title,
            email_address=email_address,
            phone_number=phone_number,
        )

    def get_alternate_contact(self, alternate_contact_type: str) -> AlternateContact:
        if alternate_contact_type not in self._alternate_contacts:
            raise UnspecifiedContactType
        return self._alternate_contacts[alternate_contact_type]

    def delete_alternate_contact(self, alternate_contact_type: str) -> None:
        self._alternate_contacts.pop(alternate_contact_type, None)


account_backends = BackendDict(
    AccountBackend,
    "account",
    use_boto3_regions=False,
    additional_regions=PARTITION_NAMES,
)
