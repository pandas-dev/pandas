"""SESV2Backend class with methods for supported APIs."""

from typing import Any, Dict, List, Optional

from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.common_models import BaseModel
from moto.core.utils import iso_8601_datetime_with_milliseconds
from moto.utilities.paginator import paginate

from ..ses.models import ConfigurationSet, Message, RawMessage, ses_backends
from .exceptions import NotFoundException

PAGINATION_MODEL = {
    "list_dedicated_ip_pools": {
        "input_token": "next_token",
        "limit_key": "page_size",
        "limit_default": 100,
        "unique_attribute": ["pool_name"],
    },
    "list_email_identities": {
        "input_token": "next_token",
        "limit_key": "page_size",
        "limit_default": 100,
        "unique_attribute": "IdentityName",
    },
    "list_configuration_sets": {
        "input_token": "next_token",
        "limit_key": "page_size",
        "limit_default": 100,
        "unique_attribute": "configuration_set_name",
    },
}


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


class EmailIdentity(BaseModel):
    def __init__(
        self,
        email_identity: str,
        tags: Optional[Dict[str, str]],
        dkim_signing_attributes: Optional[object],
        configuration_set_name: Optional[str],
    ) -> None:
        self.email_identity = email_identity
        self.tags = tags
        self.dkim_signing_attributes = dkim_signing_attributes
        self.configuration_set_name = configuration_set_name
        self.identity_type = "EMAIL_ADDRESS" if "@" in email_identity else "DOMAIN"
        self.verified_for_sending_status = False
        self.feedback_forwarding_status = False
        self.verification_status = "SUCCESS"
        self.sending_enabled = True
        self.dkim_attributes: Dict[str, Any] = {}
        if not self.dkim_signing_attributes or not isinstance(
            self.dkim_signing_attributes, dict
        ):
            self.dkim_attributes["SigningEnabled"] = False
            self.dkim_attributes["Status"] = "NOT_STARTED"
        else:
            self.dkim_attributes["SigningEnabled"] = True
            self.dkim_attributes["Status"] = "SUCCESS"
            self.dkim_attributes["NextSigningKeyLength"] = (
                self.dkim_signing_attributes.get("NextSigningKeyLength", "RSA_1024_BIT")
            )
            self.dkim_attributes["SigningAttributesOrigin"] = (
                self.dkim_signing_attributes.get("SigningAttributesOrigin", "AWS_SES")
            )
            self.dkim_attributes["LastKeyGenerationTimestamp"] = (
                iso_8601_datetime_with_milliseconds()
            )
        self.policies: Dict[str, Any] = {}

    @property
    def get_response_object(self) -> Dict[str, Any]:  # type: ignore[misc]
        return {
            "IdentityType": self.identity_type,
            "FeedbackForwardingStatus": self.feedback_forwarding_status,
            "VerifiedForSendingStatus": self.verified_for_sending_status,
            "DkimAttributes": self.dkim_attributes,
            "Policies": self.policies,
            "Tags": self.tags,
            "ConfigurationSetName": self.configuration_set_name,
            "VerificationStatus": self.verification_status,
        }

    @property
    def list_response_object(self) -> Dict[str, Any]:  # type: ignore[misc]
        return {
            "IdentityType": self.identity_type,
            "IdentityName": self.email_identity,
            "SendingEnabled": self.sending_enabled,
            "VerificationStatus": self.verification_status,
        }


class DedicatedIpPool(BaseModel):
    def __init__(
        self, pool_name: str, scaling_mode: str, tags: List[Dict[str, str]]
    ) -> None:
        self.pool_name = pool_name
        self.scaling_mode = scaling_mode
        self.tags = tags

    def to_dict(self) -> Dict[str, Any]:
        return {
            "PoolName": self.pool_name,
            "Tags": self.tags,
            "ScalingMode": self.scaling_mode,
        }


class SESV2Backend(BaseBackend):
    """Implementation of SESV2 APIs, piggy back on v1 SES"""

    def __init__(self, region_name: str, account_id: str):
        super().__init__(region_name, account_id)
        self.contacts: Dict[str, Contact] = {}
        self.contacts_lists: Dict[str, ContactList] = {}
        self.email_identities: Dict[str, EmailIdentity] = {}
        self.v1_backend = ses_backends[self.account_id][self.region_name]
        self.dedicated_ip_pools: Dict[str, DedicatedIpPool] = {}

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
        message = self.v1_backend.send_email(
            source=source,
            destinations=destinations,
            subject=subject,
            body=body,
        )
        return message

    def send_raw_email(
        self, source: str, destinations: List[str], raw_data: str
    ) -> RawMessage:
        message = self.v1_backend.send_raw_email(
            source=source, destinations=destinations, raw_data=raw_data
        )
        return message

    def create_email_identity(
        self,
        email_identity: str,
        tags: Optional[Dict[str, str]],
        dkim_signing_attributes: Optional[object],
        configuration_set_name: Optional[str],
    ) -> EmailIdentity:
        identity = EmailIdentity(
            email_identity=email_identity,
            tags=tags,
            dkim_signing_attributes=dkim_signing_attributes,
            configuration_set_name=configuration_set_name,
        )
        self.email_identities[email_identity] = identity
        return identity

    @paginate(pagination_model=PAGINATION_MODEL)
    def list_email_identities(self) -> List[EmailIdentity]:
        identities = [identity for identity in self.email_identities.values()]
        return identities

    def create_configuration_set(
        self,
        configuration_set_name: str,
        tracking_options: Dict[str, str],
        delivery_options: Dict[str, Any],
        reputation_options: Dict[str, Any],
        sending_options: Dict[str, bool],
        tags: List[Dict[str, str]],
        suppression_options: Dict[str, List[str]],
        vdm_options: Dict[str, Dict[str, str]],
    ) -> None:
        self.v1_backend.create_configuration_set_v2(
            configuration_set_name=configuration_set_name,
            tracking_options=tracking_options,
            delivery_options=delivery_options,
            reputation_options=reputation_options,
            sending_options=sending_options,
            tags=tags,
            suppression_options=suppression_options,
            vdm_options=vdm_options,
        )

    def delete_configuration_set(self, configuration_set_name: str) -> None:
        self.v1_backend.delete_configuration_set(configuration_set_name)

    def get_configuration_set(self, configuration_set_name: str) -> ConfigurationSet:
        config_set = self.v1_backend.describe_configuration_set(
            configuration_set_name=configuration_set_name
        )
        return config_set

    @paginate(pagination_model=PAGINATION_MODEL)
    def list_configuration_sets(self) -> List[ConfigurationSet]:
        return self.v1_backend._list_all_configuration_sets()

    def create_dedicated_ip_pool(
        self, pool_name: str, tags: List[Dict[str, str]], scaling_mode: str
    ) -> None:
        if pool_name not in self.dedicated_ip_pools:
            new_pool = DedicatedIpPool(
                pool_name=pool_name, tags=tags, scaling_mode=scaling_mode
            )
            self.dedicated_ip_pools[pool_name] = new_pool

    def delete_dedicated_ip_pool(self, pool_name: str) -> None:
        self.dedicated_ip_pools.pop(pool_name)

    @paginate(pagination_model=PAGINATION_MODEL)
    def list_dedicated_ip_pools(self) -> List[str]:
        return list(self.dedicated_ip_pools.keys())

    def get_dedicated_ip_pool(self, pool_name: str) -> DedicatedIpPool:
        if not self.dedicated_ip_pools.get(pool_name, None):
            raise NotFoundException(pool_name)
        return self.dedicated_ip_pools[pool_name]

    def get_email_identity(self, email_identity: str) -> EmailIdentity:
        if email_identity not in self.email_identities:
            raise NotFoundException(email_identity)
        return self.email_identities[email_identity]

    def create_email_identity_policy(
        self, email_identity: str, policy_name: str, policy: str
    ) -> None:
        email_id = self.get_email_identity(email_identity)

        email_id.policies[policy_name] = policy

        return

    def delete_email_identity_policy(
        self, email_identity: str, policy_name: str
    ) -> None:
        if email_identity not in self.email_identities:
            raise NotFoundException(email_identity)

        email_id = self.email_identities[email_identity]

        if policy_name in email_id.policies:
            del email_id.policies[policy_name]

        return

    def update_email_identity_policy(
        self, email_identity: str, policy_name: str, policy: str
    ) -> None:
        if email_identity not in self.email_identities:
            raise NotFoundException(email_identity)

        email_id = self.email_identities[email_identity]

        email_id.policies[policy_name] = policy

        return

    def get_email_identity_policies(self, email_identity: str) -> Dict[str, Any]:
        email_id = self.get_email_identity(email_identity)

        return email_id.policies


sesv2_backends = BackendDict(SESV2Backend, "sesv2")
