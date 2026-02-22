import copy
import datetime
import email
import json
import re
from email.encoders import encode_7or8bit
from email.mime.base import MIMEBase
from email.mime.multipart import MIMEMultipart
from email.utils import formataddr, getaddresses, parseaddr
from typing import Any, Literal, Optional

from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.common_models import BaseModel
from moto.core.utils import iso_8601_datetime_with_milliseconds, utcnow
from moto.sns.models import sns_backends
from moto.utilities.tagging_service import TaggingService
from moto.utilities.utils import get_partition

from .exceptions import (
    AlreadyExists,
    CannotDelete,
    ConfigurationSetAlreadyExists,
    ConfigurationSetDoesNotExist,
    EventDestinationAlreadyExists,
    InvalidLambdaFunctionException,
    InvalidParameterValue,
    InvalidRenderingParameterException,
    InvalidS3ConfigurationException,
    InvalidSnsTopicException,
    MessageRejectedError,
    NotFoundException,
    RuleDoesNotExist,
    RuleSetDoesNotExist,
    TemplateDoesNotExist,
    TemplateNameAlreadyExists,
    ValidationError,
)
from .feedback import BOUNCE, COMMON_MAIL, COMPLAINT, DELIVERY
from .template import parse_template
from .utils import get_arn, get_random_message_id, is_valid_address

RECIPIENT_LIMIT = 50


class SESFeedback(BaseModel):
    BOUNCE = "Bounce"
    COMPLAINT = "Complaint"
    DELIVERY = "Delivery"

    SUCCESS_ADDR = "success"
    BOUNCE_ADDR = "bounce"
    COMPLAINT_ADDR = "complaint"

    FEEDBACK_SUCCESS_MSG = {"test": "success"}
    FEEDBACK_BOUNCE_MSG = {"test": "bounce"}
    FEEDBACK_COMPLAINT_MSG = {"test": "complaint"}

    FORWARDING_ENABLED = "feedback_forwarding_enabled"

    @staticmethod
    def generate_message(account_id: str, msg_type: str) -> dict[str, Any]:  # type: ignore[misc]
        msg: dict[str, Any] = dict(COMMON_MAIL)
        msg["mail"]["sendingAccountId"] = account_id
        if msg_type == SESFeedback.BOUNCE:
            msg["bounce"] = BOUNCE
        elif msg_type == SESFeedback.COMPLAINT:
            msg["complaint"] = COMPLAINT
        elif msg_type == SESFeedback.DELIVERY:
            msg["delivery"] = DELIVERY

        return msg


class Message(BaseModel):
    def __init__(
        self,
        message_id: str,
        source: str,
        subject: str,
        body: str,
        destinations: dict[str, list[str]],
    ):
        self.id = message_id
        self.source = source
        self.subject = subject
        self.body = body
        self.destinations = destinations


class TemplateMessage(BaseModel):
    def __init__(
        self,
        message_id: str,
        source: str,
        template: str,
        template_data: str,
        destinations: Any,
    ):
        self.id = message_id
        self.source = source
        self.template = template
        self.template_data = template_data
        self.destinations = destinations


class BulkTemplateMessage(BaseModel):
    def __init__(
        self,
        message_ids: list[str],
        source: str,
        template: str,
        template_data: str,
        destinations: Any,
    ):
        self.ids = message_ids
        self.source = source
        self.template = template
        self.template_data = template_data
        self.destinations = destinations


class RawMessage(BaseModel):
    def __init__(
        self, message_id: str, source: str, destinations: list[str], raw_data: str
    ):
        self.id = message_id
        self.source = source
        self.destinations = destinations
        self.raw_data = raw_data


class SESQuota(BaseModel):
    def __init__(self, sent: int):
        self.sent_last24_hours = sent
        self.max24_hour_send = 200
        self.max_send_rate = 1


class ReceiptRuleSet(BaseModel):
    def __init__(self, name: str):
        self.name = name
        self.created_timestamp = datetime.datetime.now(
            datetime.timezone.utc
        ).isoformat()
        self.is_active = False  # By default, during creation
        self.rules: list[dict[str, Any]] = []

    @property
    def metadata(self) -> dict[str, str]:
        return {
            "CreatedTimestamp": self.created_timestamp,
            "Name": self.name,
        }


class ConfigurationSet(BaseModel):
    def __init__(
        self,
        configuration_set_name: str,
        tracking_options: Optional[dict[str, str]] = None,
        delivery_options: Optional[dict[str, Any]] = None,
        reputation_options: Optional[dict[str, Any]] = None,
        sending_options: Optional[dict[str, bool]] = None,
        tags: Optional[list[dict[str, str]]] = None,
        suppression_options: Optional[dict[str, list[str]]] = None,
        vdm_options: Optional[dict[str, dict[str, str]]] = None,
    ) -> None:
        # Shared between SES and SESv2
        self.configuration_set_name = configuration_set_name
        self.tracking_options = tracking_options or {}
        self.delivery_options = delivery_options or {}
        self.reputation_options = reputation_options or {}
        self.enabled = sending_options or {}  # Enabled in v1, SendingOptions in v2
        # SESv2 specific fields
        self.tags = tags or []
        self.suppression_options = suppression_options or {}
        self.vdm_options = vdm_options or {}

    def to_dict_v2(self) -> dict[str, Any]:
        return {
            "ConfigurationSetName": self.configuration_set_name,
            "TrackingOptions": self.tracking_options,
            "DeliveryOptions": self.delivery_options,
            "ReputationOptions": self.reputation_options,
            "SendingOptions": {"SendingEnabled": self.enabled},
            "Tags": self.tags,
            "SuppressionOptions": self.suppression_options,
            "VdmOptions": self.vdm_options,
        }


class Contact(BaseModel):
    def __init__(
        self,
        contact_list_name: str,
        email_address: str,
        topic_preferences: list[dict[str, str]],
        unsubscribe_all: bool,
    ) -> None:
        self.contact_list_name = contact_list_name
        self.email_address = email_address
        self.topic_default_preferences: list[dict[str, str]] = []
        self.topic_preferences = topic_preferences
        self.unsubscribe_all = unsubscribe_all
        self.created_timestamp = iso_8601_datetime_with_milliseconds()
        self.last_updated_timestamp = iso_8601_datetime_with_milliseconds()

    @property
    def response_object(self) -> dict[str, Any]:  # type: ignore[misc]
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
        topics: list[dict[str, str]],
    ) -> None:
        self.contact_list_name = contact_list_name
        self.description = description
        self.topics = topics
        self.created_timestamp = iso_8601_datetime_with_milliseconds()
        self.last_updated_timestamp = iso_8601_datetime_with_milliseconds()
        self.contacts: dict[str, Contact] = {}

    def create_contact(self, contact_list_name: str, params: dict[str, Any]) -> None:
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

    def list_contacts(self) -> list[Contact]:
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
    def response_object(self) -> dict[str, Any]:  # type: ignore[misc]
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
        tags: Optional[list[dict[str, str]]],
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
        self.dkim_attributes: dict[str, Any] = {}
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
        self.policies: dict[str, Any] = {}

    @property
    def get_response_object(self) -> dict[str, Any]:  # type: ignore[misc]
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
    def list_response_object(self) -> dict[str, Any]:  # type: ignore[misc]
        return {
            "IdentityType": self.identity_type,
            "IdentityName": self.email_identity,
            "SendingEnabled": self.sending_enabled,
            "VerificationStatus": self.verification_status,
        }


class DedicatedIpPool(BaseModel):
    def __init__(
        self, pool_name: str, scaling_mode: str, tags: list[dict[str, str]]
    ) -> None:
        self.pool_name = pool_name
        self.scaling_mode = scaling_mode
        self.tags = tags

    def to_dict(self) -> dict[str, Any]:
        return {
            "PoolName": self.pool_name,
            "Tags": self.tags,
            "ScalingMode": self.scaling_mode,
        }


class SESBackend(BaseBackend):
    """
    Responsible for mocking calls to SES.

    Sent messages are persisted in the backend. If you need to verify that a message was sent successfully, you can use the internal API to check:

    .. sourcecode:: python

        from moto.core import DEFAULT_ACCOUNT_ID
        from moto.ses import ses_backends
        ses_backend = ses_backends[DEFAULT_ACCOUNT_ID][region]
        messages = ses_backend.sent_messages # sent_messages is a List of Message objects

    Note that, as this is an internal API, the exact format may differ per versions.
    """

    __RULE_NAME_REGEX = r"^[a-zA-Z0-9_.-]+$"
    __RULE_SET_PARAM = "ruleSetName"
    __RULE_PARAM = "rule.name"

    def __init__(self, region_name: str, account_id: str):
        super().__init__(region_name, account_id)
        self.sent_messages: list[Any] = []
        self.sent_message_count = 0
        self.rejected_messages_count = 0
        self.sns_topics: dict[str, dict[str, Any]] = {}
        self.config_sets: dict[str, ConfigurationSet] = {}
        self.config_set_event_destination: dict[str, dict[str, Any]] = {}
        self.event_destinations: dict[str, int] = {}
        self.identity_mail_from_domains: dict[str, dict[str, Any]] = {}
        self.templates: dict[str, dict[str, str]] = {}
        self.receipt_rule_set: dict[str, ReceiptRuleSet] = {}
        self.dkim_tokens: dict[str, list[str]] = {}
        self.contacts: dict[str, Contact] = {}
        self.contacts_lists: dict[str, ContactList] = {}
        self.email_identities: dict[str, EmailIdentity] = {}
        self.dedicated_ip_pools: dict[str, DedicatedIpPool] = {}
        self.tagger = TaggingService()

    def _is_verified_address(self, source: str) -> bool:
        _, address = parseaddr(source)
        if address in self.email_identities:
            return True
        host = address.split("@", 1)[-1]
        return host in self.email_identities

    def verify_email_identity(self, address: str) -> None:
        _, address = parseaddr(address)
        if address not in self.email_identities:
            self.create_email_identity_v2(address, None, None, None)

    def verify_email_address(self, address: str) -> None:
        _, address = parseaddr(address)
        if address not in self.email_identities:
            self.create_email_identity_v2(address, None, None, None)

    def verify_domain(self, domain: str) -> None:
        if domain.lower() not in self.email_identities:
            self.create_email_identity_v2(domain.lower(), None, None, None)

    def create_email_identity_v2(
        self,
        email_identity: str,
        tags: Optional[list[dict[str, str]]],
        dkim_signing_attributes: Optional[object],
        configuration_set_name: Optional[str],
    ) -> EmailIdentity:
        identity = EmailIdentity(
            email_identity=email_identity,
            tags=tags,
            dkim_signing_attributes=dkim_signing_attributes,
            configuration_set_name=configuration_set_name,
        )

        if tags:
            self.tagger.tag_resource(
                arn=get_arn(self, "ses", f"email-identity/{email_identity}"),
                tags=tags,
            )
        self.email_identities[email_identity] = identity
        return identity

    def list_identities(
        self, identity_type: Optional[Literal["EmailAddress", "Domain"]]
    ) -> list[str]:
        if identity_type is not None:
            return [
                (value.email_identity)
                for _, value in self.email_identities.items()
                if value.identity_type
                == ("EMAIL_ADDRESS" if identity_type == "EmailAddress" else "DOMAIN")
            ]
        return [(value.email_identity) for _, value in self.email_identities.items()]

    def list_verified_email_addresses(self) -> list[str]:
        x = [
            (value)
            for _, value in self.email_identities.items()
            if value.identity_type == "EMAIL_ADDRESS"
        ]
        y = [(identity.email_identity) for identity in x]
        return y

    def delete_identity(self, identity: str) -> None:
        if self._is_verified_address(identity):
            del self.email_identities[identity]

    def send_email(
        self, source: str, subject: str, body: str, destinations: dict[str, list[str]]
    ) -> Message:
        recipient_count = sum(map(len, destinations.values()))
        if recipient_count > RECIPIENT_LIMIT:
            raise MessageRejectedError("Too many recipients.")
        if not self._is_verified_address(source):
            self.rejected_messages_count += 1
            raise MessageRejectedError(f"Email address not verified {source}")
        destination_addresses = [
            address for addresses in destinations.values() for address in addresses
        ]
        for address in [source, *destination_addresses]:
            msg = is_valid_address(address)
            if msg is not None:
                raise InvalidParameterValue(msg)

        self.__process_sns_feedback__(source, destinations)

        message_id = get_random_message_id()
        message = Message(message_id, source, subject, body, destinations)
        self.sent_messages.append(message)
        self.sent_message_count += recipient_count
        return message

    def send_bulk_templated_email(
        self,
        source: str,
        template: str,
        template_data: str,
        destinations: list[dict[str, dict[str, list[str]]]],
    ) -> BulkTemplateMessage:
        recipient_count = len(destinations)
        if recipient_count > RECIPIENT_LIMIT:
            raise MessageRejectedError("Too many destinations.")

        total_recipient_count = sum(
            sum(map(len, d["Destination"].values())) for d in destinations
        )
        if total_recipient_count > RECIPIENT_LIMIT:
            raise MessageRejectedError("Too many destinations.")

        if not self._is_verified_address(source):
            self.rejected_messages_count += 1
            raise MessageRejectedError(f"Email address not verified {source}")

        if not self.templates.get(template):
            raise TemplateDoesNotExist(f"Template ({template}) does not exist")

        self.__process_sns_feedback__(source, destinations)

        message_id = get_random_message_id()
        message = TemplateMessage(
            message_id, source, template, template_data, destinations
        )
        self.sent_messages.append(message)
        self.sent_message_count += total_recipient_count

        ids = [get_random_message_id() for x in range(len(destinations))]
        return BulkTemplateMessage(ids, source, template, template_data, destinations)

    def send_templated_email(
        self,
        source: str,
        template: str,
        template_data: str,
        destinations: dict[str, list[str]],
    ) -> TemplateMessage:
        recipient_count = sum(map(len, destinations.values()))
        if recipient_count > RECIPIENT_LIMIT:
            raise MessageRejectedError("Too many recipients.")
        if not self._is_verified_address(source):
            self.rejected_messages_count += 1
            raise MessageRejectedError(f"Email address not verified {source}")
        destination_addresses = [
            address for addresses in destinations.values() for address in addresses
        ]
        for address in [source, *destination_addresses]:
            msg = is_valid_address(address)
            if msg is not None:
                raise InvalidParameterValue(msg)

        if not self.templates.get(template):
            raise TemplateDoesNotExist(f"Template ({template}) does not exist")

        self.__process_sns_feedback__(source, destinations)

        message_id = get_random_message_id()
        message = TemplateMessage(
            message_id, source, template, template_data, destinations
        )
        self.sent_messages.append(message)
        self.sent_message_count += recipient_count
        return message

    def __type_of_message__(self, destinations: Any) -> Optional[str]:
        """Checks the destination for any special address that could indicate delivery,
        complaint or bounce like in SES simulator"""
        if isinstance(destinations, list):
            alladdress = destinations
        else:
            alladdress = (
                destinations.get("ToAddresses", [])
                + destinations.get("CcAddresses", [])
                + destinations.get("BccAddresses", [])
            )

        for addr in alladdress:
            if SESFeedback.SUCCESS_ADDR in addr:
                return SESFeedback.DELIVERY
            elif SESFeedback.COMPLAINT_ADDR in addr:
                return SESFeedback.COMPLAINT
            elif SESFeedback.BOUNCE_ADDR in addr:
                return SESFeedback.BOUNCE

        return None

    def __generate_feedback__(self, msg_type: str) -> dict[str, Any]:
        """Generates the SNS message for the feedback"""
        return SESFeedback.generate_message(self.account_id, msg_type)

    def __process_sns_feedback__(self, source: str, destinations: Any) -> None:
        domain = str(source)
        if "@" in domain:
            domain = domain.split("@")[1]
        if domain in self.sns_topics:
            msg_type = self.__type_of_message__(destinations)
            if msg_type is not None:
                sns_topic = self.sns_topics[domain].get(msg_type, None)
                if sns_topic is not None:
                    message = self.__generate_feedback__(msg_type)
                    if message:
                        sns_backends[self.account_id][self.region_name].publish(
                            message,  # type: ignore[arg-type]
                            arn=sns_topic,
                        )

    def send_raw_email(
        self, source: str, destinations: list[str], raw_data: str
    ) -> RawMessage:
        if source is not None:
            _, source_email_address = parseaddr(source)
            if not self._is_verified_address(source_email_address):
                raise MessageRejectedError(
                    f"Did not have authority to send from email {source_email_address}"
                )

        message = email.message_from_string(raw_data)
        if source is None:
            if message["from"] is None:
                raise MessageRejectedError("Source not specified")

            _, source = parseaddr(message["from"])
            if not self._is_verified_address(source):
                raise MessageRejectedError(
                    f"Did not have authority to send from email {source}"
                )

        fieldvalues = [
            message[header] for header in ["TO", "CC", "BCC"] if header in message
        ]
        destinations += [
            formataddr((realname, email_address))
            for realname, email_address in getaddresses(fieldvalues)
            if email_address
        ]
        if len(destinations) > RECIPIENT_LIMIT:
            raise MessageRejectedError("Too many recipients.")
        for address in [addr for addr in [source, *destinations] if addr is not None]:
            msg = is_valid_address(address)
            if msg is not None:
                raise InvalidParameterValue(msg)

        self.__process_sns_feedback__(source, destinations)

        self.sent_message_count += len(destinations)
        message_id = get_random_message_id()
        raw_message = RawMessage(message_id, source, destinations, raw_data)
        self.sent_messages.append(raw_message)
        return raw_message

    def get_send_quota(self) -> SESQuota:
        return SESQuota(self.sent_message_count)

    def get_identity_notification_attributes(
        self, identities: list[str]
    ) -> dict[str, dict[str, Any]]:
        response: dict[str, dict[str, Any]] = {}
        for identity in identities:
            config = self.sns_topics.get(identity, {})
            response[identity] = {
                "ForwardingEnabled": config.get("feedback_forwarding_enabled", True),
                "HeadersInBounceNotificationsEnabled": False,
                "HeadersInComplaintNotificationsEnabled": False,
                "HeadersInDeliveryNotificationsEnabled": False,
            }
        return response

    def set_identity_feedback_forwarding_enabled(
        self, identity: str, enabled: bool
    ) -> None:
        identity_sns_topics = self.sns_topics.get(identity, {})
        identity_sns_topics[SESFeedback.FORWARDING_ENABLED] = enabled
        self.sns_topics[identity] = identity_sns_topics

    def set_identity_notification_topic(
        self, identity: str, notification_type: str, sns_topic: Optional[str]
    ) -> None:
        identity_sns_topics = self.sns_topics.get(identity, {})
        if sns_topic is None:
            del identity_sns_topics[notification_type]
        else:
            identity_sns_topics[notification_type] = sns_topic

        self.sns_topics[identity] = identity_sns_topics

    def create_configuration_set(self, configuration_set_name: str) -> None:
        if configuration_set_name in self.config_sets:
            raise ConfigurationSetAlreadyExists(
                f"Configuration set <{configuration_set_name}> already exists"
            )
        config_set = ConfigurationSet(configuration_set_name=configuration_set_name)

        self.config_sets[configuration_set_name] = config_set

    def create_configuration_set_v2(
        self,
        configuration_set_name: str,
        tracking_options: dict[str, str],
        delivery_options: dict[str, Any],
        reputation_options: dict[str, Any],
        sending_options: dict[str, bool],
        tags: list[dict[str, str]],
        suppression_options: dict[str, list[str]],
        vdm_options: dict[str, dict[str, str]],
    ) -> None:
        if configuration_set_name in self.config_sets:
            raise ConfigurationSetAlreadyExists(
                f"Configuration set <{configuration_set_name}> already exists"
            )
        new_config_set = ConfigurationSet(
            configuration_set_name=configuration_set_name,
            tracking_options=tracking_options,
            delivery_options=delivery_options,
            reputation_options=reputation_options,
            sending_options=sending_options,
            tags=tags,
            suppression_options=suppression_options,
            vdm_options=vdm_options,
        )

        if tags:
            self.tagger.tag_resource(
                arn=get_arn(self, "ses", f"configuration-set/{configuration_set_name}"),
                tags=tags,
            )

        self.config_sets[configuration_set_name] = new_config_set

    def describe_configuration_set(
        self, configuration_set_name: str
    ) -> ConfigurationSet:
        if configuration_set_name not in self.config_sets:
            raise ConfigurationSetDoesNotExist(
                f"Configuration set <{configuration_set_name}> does not exist."
            )
        return self.config_sets[configuration_set_name]

    def delete_configuration_set(self, configuration_set_name: str) -> None:
        self.config_sets.pop(configuration_set_name)

    def list_configuration_sets(self) -> list[ConfigurationSet]:
        return list(self.config_sets.values())

    def _list_all_configuration_sets(self) -> list[ConfigurationSet]:
        return list(self.config_sets.values())

    def create_configuration_set_event_destination(
        self, configuration_set_name: str, event_destination: dict[str, Any]
    ) -> None:
        if self.config_sets.get(configuration_set_name) is None:
            raise ConfigurationSetDoesNotExist("Invalid Configuration Set Name.")

        if self.event_destinations.get(event_destination["Name"]):
            raise EventDestinationAlreadyExists("Duplicate Event destination Name.")

        self.config_set_event_destination[configuration_set_name] = event_destination
        self.event_destinations[event_destination["Name"]] = 1

    def get_send_statistics(self) -> dict[str, Any]:
        return {
            "DeliveryAttempts": self.sent_message_count,
            "Rejects": self.rejected_messages_count,
            "Complaints": 0,
            "Bounces": 0,
            "Timestamp": utcnow(),
        }

    def add_template(self, template_info: dict[str, str]) -> None:
        template_name = template_info["template_name"]
        if not template_name:
            raise ValidationError(
                "1 validation error detected: "
                "Value null at 'template.templateName'"
                "failed to satisfy constraint: Member must not be null"
            )

        if self.templates.get(template_name, None):
            raise TemplateNameAlreadyExists("Duplicate Template Name.")

        template_subject = template_info["subject_part"]
        if not template_subject:
            raise InvalidParameterValue("The subject must be specified.")
        self.templates[template_name] = template_info

    def update_template(self, template_info: dict[str, str]) -> None:
        template_name = template_info["template_name"]
        if not template_name:
            raise ValidationError(
                "1 validation error detected: "
                "Value null at 'template.templateName'"
                "failed to satisfy constraint: Member must not be null"
            )

        if not self.templates.get(template_name, None):
            raise TemplateDoesNotExist("Invalid Template Name.")

        template_subject = template_info["subject_part"]
        if not template_subject:
            raise InvalidParameterValue("The subject must be specified.")
        self.templates[template_name] = template_info

    def get_template(self, template_name: str) -> dict[str, str]:
        if not self.templates.get(template_name, None):
            raise TemplateDoesNotExist("Invalid Template Name.")
        return self.templates[template_name]

    def list_templates(self) -> list[dict[str, str]]:
        return list(self.templates.values())

    def render_template(self, render_data: dict[str, Any]) -> str:
        template_name = render_data.get("name", "")
        template = self.templates.get(template_name, None)
        if not template:
            raise TemplateDoesNotExist("Invalid Template Name.")

        template_data = render_data.get("data")
        try:
            template_data = json.loads(template_data)  # type: ignore
        except ValueError:
            raise InvalidRenderingParameterException(
                "Template rendering data is invalid"
            )

        subject_part = template["subject_part"]
        text_part = template["text_part"]
        html_part = template["html_part"]

        subject_part = parse_template(str(subject_part), template_data)
        text_part = parse_template(str(text_part), template_data)
        html_part = parse_template(str(html_part), template_data)

        email_obj = MIMEMultipart("alternative")

        mime_text = MIMEBase("text", "plain;charset=UTF-8")
        mime_text.set_payload(text_part.encode("utf-8"))
        encode_7or8bit(mime_text)
        email_obj.attach(mime_text)

        mime_html = MIMEBase("text", "html;charset=UTF-8")
        mime_html.set_payload(html_part.encode("utf-8"))
        encode_7or8bit(mime_html)
        email_obj.attach(mime_html)

        now = datetime.datetime.now().isoformat()

        rendered_template = (
            f"Date: {now}\r\nSubject: {subject_part}\r\n{email_obj.as_string()}"
        )
        return rendered_template

    def delete_template(self, name: str) -> None:
        self.templates.pop(name)

    def create_receipt_rule_set(self, rule_set_name: str) -> None:
        """
        We have to validate rule_set_name against the following conditions:
            1. Must match a given regex pattern
            2. Must start and end with a number or a character
            3. Contain 64 characters or lesser
        """
        # Boto3 throws an error with the same message for both failures, even though we could have very well combined it into one regex
        self._validate_name_param(self.__RULE_SET_PARAM, rule_set_name)
        if self.receipt_rule_set.get(rule_set_name) is not None:
            raise AlreadyExists(f"Rule set already exists: {rule_set_name}")
        self.receipt_rule_set[rule_set_name] = ReceiptRuleSet(name=rule_set_name)

    def create_receipt_rule(
        self, rule_set_name: str, rule: dict[str, Any], after: Optional[str]
    ) -> None:
        # Validate ruleSetName
        self._validate_name_param(self.__RULE_SET_PARAM, rule_set_name)
        # Validate the name of the rule
        self._validate_name_param(self.__RULE_PARAM, rule["Name"])
        # Start operating on the rule set, if it exists
        rule_set = self.receipt_rule_set.get(rule_set_name)
        if rule_set is None:
            raise RuleSetDoesNotExist(f"Rule set does not exist: {rule_set_name}")
        # If "after" is specified, then verify if a rule exists by that name under the rule set
        if after and not any(r["Name"] == after for r in rule_set.rules):
            raise RuleDoesNotExist(f"Rule does not exist: {after}")
        if rule in rule_set.rules:
            raise AlreadyExists(f"Rule already exists: {rule['Name']}")
        # Validate the actions as part of the receipt_rule
        self._validate_receipt_rule_actions(rule)
        if not after:
            # If after is not provided, the new rule is inserted at the beginning of the rule list
            rule_set.rules.insert(0, rule)
        else:
            # Find the position of the rule specified by after parameter
            for rule_idx, _ in enumerate(rule_set.rules):
                if rule_set.rules[rule_idx]["Name"] == after:
                    rule_set.rules.insert(rule_idx + 1, rule)
                    break
        self.receipt_rule_set[rule_set_name] = rule_set

    def clone_receipt_rule_set(
        self, original_rule_set_name: str, rule_set_name: str
    ) -> None:
        # Boto3 validates both the original and new rule set names
        self._validate_name_param(self.__RULE_SET_PARAM, original_rule_set_name)
        self._validate_name_param(self.__RULE_SET_PARAM, rule_set_name)

        # Check if original_rule_set_name exists
        if self.receipt_rule_set.get(original_rule_set_name) is None:
            raise RuleSetDoesNotExist(
                f"Rule set does not exist: {original_rule_set_name}"
            )

        # Check if rule_set_name already exists
        if self.receipt_rule_set.get(rule_set_name) is not None:
            raise AlreadyExists(f"Rule set already exists: {rule_set_name}")

        # Clone the original rule set
        original_rule_set = self.receipt_rule_set[original_rule_set_name]
        self.receipt_rule_set[rule_set_name] = ReceiptRuleSet(name=rule_set_name)
        self.receipt_rule_set[rule_set_name].rules = copy.deepcopy(
            original_rule_set.rules
        )

    def set_active_receipt_rule_set(
        self, rule_set_name: Optional[str]
    ) -> Optional[ReceiptRuleSet]:
        if not rule_set_name:
            # A null rule_set_name parameter (i.e., not passed in the request at all) means that all receipt rule sets should be marked inactive
            for rs in self.receipt_rule_set.values():
                rs.is_active = False
            return None
        self._validate_name_param(self.__RULE_SET_PARAM, rule_set_name)
        # Verify that the rule set exists
        if rule_set_name not in self.receipt_rule_set:
            raise RuleSetDoesNotExist(f"Rule set does not exist: {rule_set_name}")
        # Only one active rule set is allowed at a time
        for rs in self.receipt_rule_set.values():
            rs.is_active = False
        self.receipt_rule_set[rule_set_name].is_active = True
        return self.receipt_rule_set[rule_set_name]

    def describe_active_receipt_rule_set(self) -> Optional[ReceiptRuleSet]:
        for rs in self.receipt_rule_set.values():
            if rs.is_active:
                return rs
        return None

    def delete_receipt_rule_set(self, rule_set_name: str) -> None:
        self._validate_name_param(self.__RULE_SET_PARAM, rule_set_name)
        # If the rule set does not exist, boto3 silently returns with success response
        if rule_set_name not in self.receipt_rule_set:
            return
        # If the rule set is active, then raise a CannotDeleteException
        if self.receipt_rule_set[rule_set_name].is_active:
            raise CannotDelete(f"Cannot delete active rule set: {rule_set_name}")
        del self.receipt_rule_set[rule_set_name]

    def list_receipt_rule_sets(self) -> list[ReceiptRuleSet]:
        # The receipt rule sets are ordered by name
        return sorted(self.receipt_rule_set.values(), key=lambda rs: rs.name)

    def _validate_name_param(self, param_name: str, name: str) -> None:
        # Boto3 throws an error with the same message for both failures, even though we could have very well combined it into one regex
        if (
            not re.match(SESBackend.__RULE_NAME_REGEX, name)
            or not name[0].isalnum()
            or not name[-1].isalnum()
        ):
            raise ValidationError(
                f"Value at '{param_name}' failed to satisfy constraint: Member must satisfy regular expression pattern: ^[a-zA-Z0-9_.-]+$"
            )
        # A different message thrown for length of the rule_name_set string
        if len(name) > 64:
            raise ValidationError(f"Not a valid {param_name}: {name}")
        return

    def _validate_receipt_rule_actions(self, rule: dict[str, Any]) -> None:
        # Allowed to be empty
        actions = rule.get("Actions", [])
        for action in actions:
            if "S3Action" in action:
                self._validate_s3_action(action["S3Action"])
            if "BounceAction" in action:
                self._validate_sns_topic(action["BounceAction"].get("TopicArn"))
            if "WorkmailAction" in action:
                self._validate_sns_topic(action["WorkmailAction"].get("TopicArn"))
            if "LambdaAction" in action:
                self._validate_lambda_action(action["LambdaAction"])

    def _validate_s3_action(self, s3_action: dict[str, str]) -> None:
        from moto.s3.models import s3_backends

        # Raise an exception if the bucket does not exist
        try:
            partition = get_partition(self.region_name)
            s3_backends[self.account_id][partition].get_bucket(s3_action["BucketName"])
        except Exception:
            raise InvalidS3ConfigurationException(
                f"Could not write to bucket: {s3_action['BucketName']}"
            )

        if s3_action.get("KmsKeyArn"):
            self._validate_kms_key(s3_action["KmsKeyArn"])

        if s3_action.get("TopicArn"):
            self._validate_sns_topic(s3_action["TopicArn"])

        if s3_action.get("IamRoleArn"):
            self._validate_iam_role(s3_action["IamRoleArn"])

    def _validate_lambda_action(self, lambda_action: dict[str, str]) -> None:
        from moto.awslambda.models import lambda_backends

        # Raise an exception if the Lambda function does not exist
        try:
            _ = lambda_backends[self.account_id][self.region_name].get_function(
                lambda_action["FunctionArn"]
            )
        except Exception:
            raise InvalidLambdaFunctionException(
                f"Invalid Lambda function: {lambda_action['FunctionArn']}"
            )

        self._validate_sns_topic(lambda_action.get("TopicArn"))

    def _validate_kms_key(self, kms_key_arn: str) -> None:
        from moto.kms.models import kms_backends

        # Raise an exception if the KMS key does not exist
        try:
            region_kms_backend = kms_backends[self.account_id][self.region_name]
            _ = region_kms_backend.describe_key(kms_key_arn)
        except Exception:
            raise InvalidS3ConfigurationException(
                f"Unable to use AWS KMS key: {kms_key_arn}"
            )

    def _validate_sns_topic(self, topic_arn: Optional[str]) -> None:
        # Nothing to validate
        if not topic_arn:
            return

        from moto.sns.models import sns_backends

        # Raise an exception if the SNS topic does not exist
        if topic_arn not in sns_backends[self.account_id][self.region_name].topics:
            raise InvalidSnsTopicException(f"Invalid SNS topic: {topic_arn}")

    def _validate_iam_role(self, role_arn: str) -> None:
        from moto.iam.models import iam_backends

        # Raise an exception if the IAM role does not exist
        if (
            role_arn
            not in iam_backends[self.account_id][get_partition(self.region_name)].roles
        ):
            raise InvalidParameterValue("Could not assume the provided IAM role")

    def describe_receipt_rule_set(self, rule_set_name: str) -> ReceiptRuleSet:
        rule_set = self.receipt_rule_set.get(rule_set_name)

        if rule_set is None:
            raise RuleSetDoesNotExist(f"Rule set does not exist: {rule_set_name}")

        return rule_set

    def describe_receipt_rule(
        self, rule_set_name: str, rule_name: str
    ) -> dict[str, Any]:
        rule_set = self.receipt_rule_set.get(rule_set_name)

        if rule_set is None:
            raise RuleSetDoesNotExist(f"Rule set does not exist: {rule_set_name}")

        for receipt_rule in rule_set.rules:
            if receipt_rule["Name"] == rule_name:
                return receipt_rule

        raise RuleDoesNotExist(f"Rule does not exist: {rule_name}")

    def update_receipt_rule(self, rule_set_name: str, rule: dict[str, Any]) -> None:
        rule_set = self.receipt_rule_set.get(rule_set_name)

        if rule_set is None:
            raise RuleSetDoesNotExist(f"Rule set does not exist: {rule_set_name}")

        for i, receipt_rule in enumerate(rule_set.rules):
            if receipt_rule["Name"] == rule["Name"]:
                rule_set.rules[i] = rule
                break
        else:
            raise RuleDoesNotExist(f"Rule does not exist: {rule['Name']}")

    def set_identity_mail_from_domain(
        self,
        identity: str,
        mail_from_domain: Optional[str] = None,
        behavior_on_mx_failure: Optional[str] = None,
    ) -> None:
        if not self._is_verified_address(identity):
            raise InvalidParameterValue(f"Identity '{identity}' does not exist.")

        if mail_from_domain is None:
            self.identity_mail_from_domains.pop(identity)
            return

        if not mail_from_domain.endswith(identity.split("@")[-1]):
            raise InvalidParameterValue(
                f"Provided MAIL-FROM domain '{mail_from_domain}' is not subdomain of "
                f"the domain of the identity '{identity.split('@')[-1]}'."
            )

        if behavior_on_mx_failure not in (None, "RejectMessage", "UseDefaultValue"):
            raise ValidationError(
                "1 validation error detected: "
                f"Value '{behavior_on_mx_failure}' at 'behaviorOnMXFailure'"
                "failed to satisfy constraint: Member must satisfy enum value set: "
                "[RejectMessage, UseDefaultValue]"
            )

        self.identity_mail_from_domains[identity] = {
            "mail_from_domain": mail_from_domain,
            "behavior_on_mx_failure": behavior_on_mx_failure,
        }

    def get_identity_mail_from_domain_attributes(
        self, identities: Optional[list[str]] = None
    ) -> dict[str, dict[str, str]]:
        if identities is None:
            actual_identities = []
        else:
            actual_identities = [
                ident
                for ident in self.email_identities.values()
                if ident.email_identity in identities
            ]
        attributes_by_identity = {}
        for identity in actual_identities:
            value = self.identity_mail_from_domains.get(identity.email_identity, {})
            mail_from_domain = value.get("mail_from_domain")
            attributes_by_identity[identity.email_identity] = {
                "MailFromDomain": mail_from_domain,
                "MailFromDomainStatus": "Success" if mail_from_domain else None,
                "BehaviorOnMXFailure": value.get(
                    "behavior_on_mx_failure", "UseDefaultValue"
                ),
            }
        return attributes_by_identity

    def get_identity_verification_attributes(
        self, identities: Optional[list[str]] = None
    ) -> dict[str, dict[str, str]]:
        if identities is None:
            actual_identities = []
        else:
            actual_identities = [
                ident
                for ident in self.email_identities.values()
                if ident.email_identity in identities
            ]
        attributes_by_identity = {}
        for identity in actual_identities:
            attributes_by_identity[identity.email_identity] = {
                "VerificationStatus": "Success",
                "VerificationToken": "ILQMESfEW0p6i6gIJcEWvO65TP5hg6B99hGFZ2lxrIs=",
            }
        return attributes_by_identity

    def update_configuration_set_reputation_metrics_enabled(
        self, configuration_set_name: str, enabled: bool
    ) -> None:
        if configuration_set_name not in self.config_sets:
            raise ConfigurationSetDoesNotExist(
                f"Configuration set <{configuration_set_name}> does not exist"
            )
        config_set = self.config_sets[configuration_set_name]
        if config_set.reputation_options is None:
            config_set.reputation_options = {}
        config_set.reputation_options["ReputationMetricsEnabled"] = enabled

    def get_identity_dkim_attributes(
        self, identities: list[str]
    ) -> dict[str, dict[str, Any]]:
        result = {}
        actual_identities: list[EmailIdentity] = []
        for ident in identities:
            if ident in self.email_identities.keys():
                actual_identities.append(self.email_identities[ident])
            else:
                actual_identities.append(
                    EmailIdentity(
                        email_identity=ident,
                        tags=None,
                        dkim_signing_attributes=None,
                        configuration_set_name=None,
                    )
                )

        for identity in actual_identities:
            is_domain = identity.identity_type == "DOMAIN"
            dkim_enabled = True
            verification_status = (
                "Success"
                if identity.email_identity in self.email_identities.keys()
                else "NotStarted"
            )

            dkim_data = {
                "DkimEnabled": dkim_enabled,
                "DkimVerificationStatus": verification_status,
            }

            # Only include tokens for domain identities
            if is_domain and verification_status == "Success":
                if identity.email_identity not in self.dkim_tokens:
                    # Generate new DKIM tokens for the domain
                    self.dkim_tokens[identity.email_identity] = [
                        "vvjuipp74whm76gqoni7qmwwn4w4qusjiainivf6sf",
                        "3frqe7jn4obpuxjpwpolz6ipb3k5nvt2nhjpik2oy",
                        "wrqplteh7oodxnad7hsl4mixg2uavzneazxv5sxi2",
                    ]
                dkim_data["DkimTokens"] = self.dkim_tokens[identity.email_identity]

            result[identity.email_identity] = dkim_data

        return result


ses_backends = BackendDict(SESBackend, "ses")
