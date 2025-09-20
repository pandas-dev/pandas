import copy
import datetime
import email
import json
import re
from email.encoders import encode_7or8bit
from email.mime.base import MIMEBase
from email.mime.multipart import MIMEMultipart
from email.utils import formataddr, getaddresses, parseaddr
from typing import Any, Dict, List, Optional

from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.common_models import BaseModel
from moto.core.utils import utcnow
from moto.sns.models import sns_backends
from moto.utilities.paginator import paginate
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
    RuleDoesNotExist,
    RuleSetDoesNotExist,
    TemplateDoesNotExist,
    TemplateNameAlreadyExists,
    ValidationError,
)
from .feedback import BOUNCE, COMMON_MAIL, COMPLAINT, DELIVERY
from .template import parse_template
from .utils import get_random_message_id, is_valid_address

RECIPIENT_LIMIT = 50

PAGINATION_MODEL = {
    "list_configuration_sets": {
        "input_token": "next_token",
        "limit_key": "max_items",
        "limit_default": 100,
        "unique_attribute": "configuration_set_name",
    },
}


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
    def generate_message(account_id: str, msg_type: str) -> Dict[str, Any]:  # type: ignore[misc]
        msg: Dict[str, Any] = dict(COMMON_MAIL)
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
        destinations: Dict[str, List[str]],
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
        message_ids: List[str],
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
        self, message_id: str, source: str, destinations: List[str], raw_data: str
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
        self.rules: List[Dict[str, Any]] = []

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
        tracking_options: Optional[Dict[str, str]] = {},
        delivery_options: Optional[Dict[str, Any]] = {},
        reputation_options: Optional[Dict[str, Any]] = {},
        sending_options: Optional[Dict[str, bool]] = {},
        tags: Optional[List[Dict[str, str]]] = [],
        suppression_options: Optional[Dict[str, List[str]]] = {},
        vdm_options: Optional[Dict[str, Dict[str, str]]] = {},
    ) -> None:
        # Shared between SES and SESv2
        self.configuration_set_name = configuration_set_name
        self.tracking_options = tracking_options
        self.delivery_options = delivery_options
        self.reputation_options = reputation_options
        self.enabled = sending_options  # Enabled in v1, SendingOptions in v2
        # SESv2 specific fields
        self.tags = tags
        self.suppression_options = suppression_options
        self.vdm_options = vdm_options

    def to_dict_v2(self) -> Dict[str, Any]:
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
        self.addresses: List[str] = []
        self.email_addresses: List[str] = []
        self.domains: List[str] = []
        self.sent_messages: List[Any] = []
        self.sent_message_count = 0
        self.rejected_messages_count = 0
        self.sns_topics: Dict[str, Dict[str, Any]] = {}
        self.config_sets: Dict[str, ConfigurationSet] = {}
        self.config_set_event_destination: Dict[str, Dict[str, Any]] = {}
        self.event_destinations: Dict[str, int] = {}
        self.identity_mail_from_domains: Dict[str, Dict[str, Any]] = {}
        self.templates: Dict[str, Dict[str, str]] = {}
        self.receipt_rule_set: Dict[str, ReceiptRuleSet] = {}
        self.dkim_tokens: Dict[str, List[str]] = {}

    def _is_verified_address(self, source: str) -> bool:
        _, address = parseaddr(source)
        if address in self.addresses:
            return True
        if address in self.email_addresses:
            return True
        host = address.split("@", 1)[-1]
        return host in self.domains

    def verify_email_identity(self, address: str) -> None:
        _, address = parseaddr(address)
        if address not in self.addresses:
            self.addresses.append(address)

    def verify_email_address(self, address: str) -> None:
        _, address = parseaddr(address)
        self.email_addresses.append(address)

    def verify_domain(self, domain: str) -> None:
        if domain.lower() not in self.domains:
            self.domains.append(domain.lower())

    def list_identities(self, identity_type: str) -> List[str]:
        if identity_type == "Domain":
            return self.domains
        if identity_type == "EmailAddress":
            return self.addresses
        return self.domains + self.addresses

    def list_verified_email_addresses(self) -> List[str]:
        return self.email_addresses

    def delete_identity(self, identity: str) -> None:
        if "@" in identity:
            self.addresses.remove(identity)
        else:
            self.domains.remove(identity)

    def send_email(
        self, source: str, subject: str, body: str, destinations: Dict[str, List[str]]
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
        destinations: List[Dict[str, Dict[str, List[str]]]],
    ) -> BulkTemplateMessage:
        recipient_count = len(destinations)
        if recipient_count > RECIPIENT_LIMIT:
            raise MessageRejectedError("Too many destinations.")

        total_recipient_count = sum(
            map(lambda d: sum(map(len, d["Destination"].values())), destinations)
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

        ids = list(map(lambda x: get_random_message_id(), range(len(destinations))))
        return BulkTemplateMessage(ids, source, template, template_data, destinations)

    def send_templated_email(
        self,
        source: str,
        template: str,
        template_data: str,
        destinations: Dict[str, List[str]],
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

    def __generate_feedback__(self, msg_type: str) -> Dict[str, Any]:
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
        self, source: str, destinations: List[str], raw_data: str
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
        self, identities: List[str]
    ) -> Dict[str, Dict[str, Any]]:
        response: Dict[str, Dict[str, Any]] = {}
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
        tracking_options: Dict[str, str],
        delivery_options: Dict[str, Any],
        reputation_options: Dict[str, Any],
        sending_options: Dict[str, bool],
        tags: List[Dict[str, str]],
        suppression_options: Dict[str, List[str]],
        vdm_options: Dict[str, Dict[str, str]],
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
        self.config_sets[configuration_set_name] = new_config_set

    def describe_configuration_set(
        self, configuration_set_name: str
    ) -> ConfigurationSet:
        if configuration_set_name not in self.config_sets:
            raise ConfigurationSetDoesNotExist(
                f"Configuration set <{configuration_set_name}> does not exist"
            )
        return self.config_sets[configuration_set_name]

    def delete_configuration_set(self, configuration_set_name: str) -> None:
        self.config_sets.pop(configuration_set_name)

    @paginate(pagination_model=PAGINATION_MODEL)
    def list_configuration_sets(self) -> List[ConfigurationSet]:
        return list(self.config_sets.values())

    def _list_all_configuration_sets(self) -> List[ConfigurationSet]:
        return list(self.config_sets.values())

    def create_configuration_set_event_destination(
        self, configuration_set_name: str, event_destination: Dict[str, Any]
    ) -> None:
        if self.config_sets.get(configuration_set_name) is None:
            raise ConfigurationSetDoesNotExist("Invalid Configuration Set Name.")

        if self.event_destinations.get(event_destination["Name"]):
            raise EventDestinationAlreadyExists("Duplicate Event destination Name.")

        self.config_set_event_destination[configuration_set_name] = event_destination
        self.event_destinations[event_destination["Name"]] = 1

    def get_send_statistics(self) -> Dict[str, Any]:
        return {
            "DeliveryAttempts": self.sent_message_count,
            "Rejects": self.rejected_messages_count,
            "Complaints": 0,
            "Bounces": 0,
            "Timestamp": utcnow(),
        }

    def add_template(self, template_info: Dict[str, str]) -> None:
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

    def update_template(self, template_info: Dict[str, str]) -> None:
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

    def get_template(self, template_name: str) -> Dict[str, str]:
        if not self.templates.get(template_name, None):
            raise TemplateDoesNotExist("Invalid Template Name.")
        return self.templates[template_name]

    def list_templates(self) -> List[Dict[str, str]]:
        return list(self.templates.values())

    def render_template(self, render_data: Dict[str, Any]) -> str:
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
        self, rule_set_name: str, rule: Dict[str, Any], after: Optional[str]
    ) -> None:
        # Validate ruleSetName
        self._validate_name_param(self.__RULE_SET_PARAM, rule_set_name)
        # Validate the name of the rule
        self._validate_name_param(self.__RULE_PARAM, rule["name"])
        # Start operating on the rule set, if it exists
        rule_set = self.receipt_rule_set.get(rule_set_name)
        if rule_set is None:
            raise RuleSetDoesNotExist(f"Rule set does not exist: {rule_set_name}")
        # If "after" is specified, then verify if a rule exists by that name under the rule set
        if after and not any(r["name"] == after for r in rule_set.rules):
            raise RuleDoesNotExist(f"Rule does not exist: {after}")
        if rule in rule_set.rules:
            raise AlreadyExists(f"Rule already exists: {rule['name']}")
        # Validate the actions as part of the receipt_rule
        self._validate_receipt_rule_actions(rule)
        if not after:
            # If after is not provided, the new rule is inserted at the beginning of the rule list
            rule_set.rules.insert(0, rule)
        else:
            # Find the position of the rule specified by after parameter
            for rule_idx, _ in enumerate(rule_set.rules):
                if rule_set.rules[rule_idx]["name"] == after:
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

    def list_receipt_rule_sets(self) -> List[ReceiptRuleSet]:
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

    def _validate_receipt_rule_actions(self, rule: Dict[str, Any]) -> None:
        # Allowed to be empty
        actions = rule.get("actions", [])
        for action in actions:
            if "S3Action" in action:
                self._validate_s3_action(action["S3Action"])
            if "BounceAction" in action:
                self._validate_sns_topic(action["BounceAction"].get("topic_arn"))
            if "WorkmailAction" in action:
                self._validate_sns_topic(action["WorkmailAction"].get("topic_arn"))
            if "LambdaAction" in action:
                self._validate_lambda_action(action["LambdaAction"])

    def _validate_s3_action(self, s3_action: Dict[str, str]) -> None:
        from moto.s3.models import s3_backends

        # Raise an exception if the bucket does not exist
        try:
            partition = get_partition(self.region_name)
            s3_backends[self.account_id][partition].get_bucket(s3_action["bucket_name"])
        except Exception:
            raise InvalidS3ConfigurationException(
                f"Could not write to bucket: {s3_action['bucket_name']}"
            )

        if s3_action.get("kms_key_arn"):
            self._validate_kms_key(s3_action["kms_key_arn"])

        if s3_action.get("topic_arn"):
            self._validate_sns_topic(s3_action["topic_arn"])

        if s3_action.get("iam_role_arn"):
            self._validate_iam_role(s3_action["iam_role_arn"])

    def _validate_lambda_action(self, lambda_action: Dict[str, str]) -> None:
        from moto.awslambda.models import lambda_backends

        # Raise an exception if the Lambda function does not exist
        try:
            _ = lambda_backends[self.account_id][self.region_name].get_function(
                lambda_action["function_arn"]
            )
        except Exception:
            raise InvalidLambdaFunctionException(
                f"Invalid Lambda function: {lambda_action['function_arn']}"
            )

        self._validate_sns_topic(lambda_action.get("topic_arn"))

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
    ) -> Dict[str, Any]:
        rule_set = self.receipt_rule_set.get(rule_set_name)

        if rule_set is None:
            raise RuleSetDoesNotExist(f"Rule set does not exist: {rule_set_name}")

        for receipt_rule in rule_set.rules:
            if receipt_rule["name"] == rule_name:
                return receipt_rule

        raise RuleDoesNotExist(f"Rule does not exist: {rule_name}")

    def update_receipt_rule(self, rule_set_name: str, rule: Dict[str, Any]) -> None:
        rule_set = self.receipt_rule_set.get(rule_set_name)

        if rule_set is None:
            raise RuleSetDoesNotExist(f"Rule set does not exist: {rule_set_name}")

        for i, receipt_rule in enumerate(rule_set.rules):
            if receipt_rule["name"] == rule["name"]:
                rule_set.rules[i] = rule
                break
        else:
            raise RuleDoesNotExist(f"Rule does not exist: {rule['name']}")

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
        self, identities: Optional[List[str]] = None
    ) -> Dict[str, Dict[str, str]]:
        if identities is None:
            identities = []
        attributes_by_identity = {}
        for identity in identities:
            if identity in (self.domains + self.addresses):
                value = self.identity_mail_from_domains.get(identity, {})
                mail_from_domain = value.get("mail_from_domain")
                attributes_by_identity[identity] = {
                    "MailFromDomain": mail_from_domain,
                    "MailFromDomainStatus": "Success" if mail_from_domain else None,
                    "BehaviorOnMXFailure": value.get(
                        "behavior_on_mx_failure", "UseDefaultValue"
                    ),
                }
        return attributes_by_identity

    def get_identity_verification_attributes(
        self, identities: Optional[List[str]] = None
    ) -> Dict[str, Dict[str, str]]:
        if identities is None:
            identities = []
        attributes_by_identity = {}
        for identity in identities:
            if identity in (self.domains + self.addresses):
                attributes_by_identity[identity] = {
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
        self, identities: List[str]
    ) -> Dict[str, Dict[str, Any]]:
        result = {}
        for identity in identities:
            is_domain = "@" not in identity
            dkim_enabled = True
            verification_status = (
                "Success"
                if identity in (self.domains + self.addresses)
                else "NotStarted"
            )

            dkim_data = {
                "DkimEnabled": dkim_enabled,
                "DkimVerificationStatus": verification_status,
            }

            # Only include tokens for domain identities
            if is_domain and verification_status == "Success":
                if identity not in self.dkim_tokens:
                    # Generate new DKIM tokens for the domain
                    self.dkim_tokens[identity] = [
                        "vvjuipp74whm76gqoni7qmwwn4w4qusjiainivf6sf",
                        "3frqe7jn4obpuxjpwpolz6ipb3k5nvt2nhjpik2oy",
                        "wrqplteh7oodxnad7hsl4mixg2uavzneazxv5sxi2",
                    ]
                dkim_data["DkimTokens"] = self.dkim_tokens[identity]

            result[identity] = dkim_data

        return result


ses_backends = BackendDict(SESBackend, "ses")
