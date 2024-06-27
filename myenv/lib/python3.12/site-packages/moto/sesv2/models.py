"""SESV2Backend class with methods for supported APIs."""

from typing import Any, Dict, List

from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.common_models import BaseModel
from moto.core.utils import iso_8601_datetime_with_milliseconds

from ..ses.models import Message, RawMessage, ses_backends
from .exceptions import NotFoundException


class Contact(BaseModel):
    def __init__(
        self,
        contact_list_name: str,
        email_address: str,
        topic_preferences: List[Dict[str, str]],
        unsubscribe_all: bool,
    ) -> None:
        self.contact_list_name = contact_list_name
        self.email_address = email_address
        self.topic_default_preferences: List[Dict[str, str]] = []
        self.topic_preferences = topic_preferences
        self.unsubscribe_all = unsubscribe_all
        self.created_timestamp = iso_8601_datetime_with_milliseconds()
        self.last_updated_timestamp = iso_8601_datetime_with_milliseconds()

    @property
    def response_object(self) -> Dict[str, Any]:  # type: ignore[misc]
        return {
            "ContactListName": self.contact_list_name,
            "EmailAddress": self.email_address,
            "TopicDefaultPreferences": self.topic_default_preferences,
            "TopicPreferences": self.topic_preferences,
            "UnsubscribeAll": self.unsubscribe_all,
            "CreatedTimestamp": self.created_timestamp,
            "LastUpdatedTimestamp": self.last_updated_timestamp,
        }


class ContactList(BaseModel):
    def __init__(
        self,
        contact_list_name: str,
        description: str,
        topics: List[Dict[str, str]],
    ) -> None:
        self.contact_list_name = contact_list_name
        self.description = description
        self.topics = topics
        self.created_timestamp = iso_8601_datetime_with_milliseconds()
        self.last_updated_timestamp = iso_8601_datetime_with_milliseconds()
        self.contacts: Dict[str, Contact] = {}

    def create_contact(self, contact_list_name: str, params: Dict[str, Any]) -> None:
        email_address = params["EmailAddress"]
        topic_preferences = (
            [] if "TopicPreferences" not in params else params["TopicPreferences"]
        )
        unsubscribe_all = (
            False if "UnsubscribeAll" not in params else params["UnsubscribeAll"]
        )
        new_contact = Contact(
            contact_list_name, email_address, topic_preferences, unsubscribe_all
        )
        self.contacts[email_address] = new_contact

    def list_contacts(self) -> List[Contact]:
        return self.contacts.values()  # type: ignore[return-value]

    def get_contact(self, email: str) -> Contact:
        if email in self.contacts:
            return self.contacts[email]
        else:
            raise NotFoundException(f"{email} doesn't exist in List.")

    def delete_contact(self, email: str) -> None:
        # delete if contact exists, otherwise get_contact will throw appropriate exception
        if self.get_contact(email):
            del self.contacts[email]

    @property
    def response_object(self) -> Dict[str, Any]:  # type: ignore[misc]
        return {
            "ContactListName": self.contact_list_name,
            "Description": self.description,
            "Topics": self.topics,
            "CreatedTimestamp": self.created_timestamp,
            "LastUpdatedTimestamp": self.last_updated_timestamp,
        }


class SESV2Backend(BaseBackend):
    """Implementation of SESV2 APIs, piggy back on v1 SES"""

    def __init__(self, region_name: str, account_id: str):
        super().__init__(region_name, account_id)
        self.contacts: Dict[str, Contact] = {}
        self.contacts_lists: Dict[str, ContactList] = {}

    def create_contact_list(self, params: Dict[str, Any]) -> None:
        name = params["ContactListName"]
        description = params.get("Description")
        topics = [] if "Topics" not in params else params["Topics"]
        new_list = ContactList(name, str(description), topics)
        self.contacts_lists[name] = new_list

    def get_contact_list(self, contact_list_name: str) -> ContactList:
        if contact_list_name in self.contacts_lists:
            return self.contacts_lists[contact_list_name]
        else:
            raise NotFoundException(
                f"List with name: {contact_list_name} doesn't exist."
            )

    def list_contact_lists(self) -> List[ContactList]:
        return self.contacts_lists.values()  # type: ignore[return-value]

    def delete_contact_list(self, name: str) -> None:
        if name in self.contacts_lists:
            del self.contacts_lists[name]
        else:
            raise NotFoundException(f"List with name: {name} doesn't exist")

    def create_contact(self, contact_list_name: str, params: Dict[str, Any]) -> None:
        contact_list = self.get_contact_list(contact_list_name)
        contact_list.create_contact(contact_list_name, params)
        return

    def get_contact(self, email: str, contact_list_name: str) -> Contact:
        contact_list = self.get_contact_list(contact_list_name)
        contact = contact_list.get_contact(email)
        return contact

    def list_contacts(self, contact_list_name: str) -> List[Contact]:
        contact_list = self.get_contact_list(contact_list_name)
        contacts = contact_list.list_contacts()
        return contacts

    def delete_contact(self, email: str, contact_list_name: str) -> None:
        contact_list = self.get_contact_list(contact_list_name)
        contact_list.delete_contact(email)
        return

    def send_email(
        self, source: str, destinations: Dict[str, List[str]], subject: str, body: str
    ) -> Message:
        v1_backend = ses_backends[self.account_id][self.region_name]
        message = v1_backend.send_email(
            source=source,
            destinations=destinations,
            subject=subject,
            body=body,
        )
        return message

    def send_raw_email(
        self, source: str, destinations: List[str], raw_data: str
    ) -> RawMessage:
        v1_backend = ses_backends[self.account_id][self.region_name]
        message = v1_backend.send_raw_email(
            source=source, destinations=destinations, raw_data=raw_data
        )
        return message


sesv2_backends = BackendDict(SESV2Backend, "sesv2")
