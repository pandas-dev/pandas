import contextlib
import json
import re
from collections import OrderedDict
from typing import Any, Dict, Iterable, List, Optional, Set, Tuple

import requests

from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.common_models import BaseModel, CloudFormationModel
from moto.core.utils import (
    camelcase_to_underscores,
    iso_8601_datetime_with_milliseconds,
)
from moto.moto_api._internal import mock_random
from moto.sqs import sqs_backends
from moto.sqs.exceptions import MissingParameter
from moto.utilities.arns import parse_arn
from moto.utilities.utils import get_partition

from .exceptions import (
    BatchEntryIdsNotDistinct,
    DuplicateSnsEndpointError,
    InternalError,
    InvalidParameterValue,
    ResourceNotFoundError,
    SnsEndpointDisabled,
    SNSInvalidParameter,
    SNSNotFoundError,
    TagLimitExceededError,
    TooManyEntriesInBatchRequest,
    TopicNotFound,
)
from .utils import (
    FilterPolicyMatcher,
    is_e164,
    make_arn_for_subscription,
    make_arn_for_topic,
)

DEFAULT_PAGE_SIZE = 100
MAXIMUM_MESSAGE_LENGTH = 262144  # 256 KiB
MAXIMUM_SMS_MESSAGE_BYTES = 1600  # Amazon limit for a single publish SMS action


class Topic(CloudFormationModel):
    def __init__(self, name: str, sns_backend: "SNSBackend"):
        self.name = name
        self.sns_backend = sns_backend
        self.account_id = sns_backend.account_id
        self.display_name = ""
        self.delivery_policy = ""
        self.kms_master_key_id = ""
        self.effective_delivery_policy = json.dumps(DEFAULT_EFFECTIVE_DELIVERY_POLICY)
        self.arn = make_arn_for_topic(self.account_id, name, sns_backend.region_name)

        self.subscriptions_pending = 0
        self.subscriptions_confimed = 0
        self.subscriptions_deleted = 0
        self.sent_notifications: List[
            Tuple[str, str, Optional[str], Optional[Dict[str, Any]], Optional[str]]
        ] = []

        self._policy_json = self._create_default_topic_policy(
            sns_backend.region_name, self.account_id, name
        )
        self._tags: Dict[str, str] = {}
        self.fifo_topic = "false"
        self.content_based_deduplication = "false"

    def publish(
        self,
        message: str,
        subject: Optional[str] = None,
        message_attributes: Optional[Dict[str, Any]] = None,
        group_id: Optional[str] = None,
        deduplication_id: Optional[str] = None,
    ) -> str:
        message_id = str(mock_random.uuid4())
        subscriptions, _ = self.sns_backend.list_subscriptions_by_topic(
            topic_arn=self.arn
        )
        for subscription in subscriptions:
            subscription.publish(
                message,
                message_id,
                subject=subject,
                message_attributes=message_attributes,
                group_id=group_id,
                deduplication_id=deduplication_id,
            )
        self.sent_notifications.append(
            (message_id, message, subject, message_attributes, group_id)
        )
        return message_id

    @classmethod
    def has_cfn_attr(cls, attr: str) -> bool:
        return attr in ["TopicName"]

    def get_cfn_attribute(self, attribute_name: str) -> str:
        from moto.cloudformation.exceptions import UnformattedGetAttTemplateException

        if attribute_name == "TopicName":
            return self.name
        raise UnformattedGetAttTemplateException()

    @property
    def physical_resource_id(self) -> str:
        return self.arn

    @property
    def policy(self) -> str:
        return json.dumps(self._policy_json, separators=(",", ":"))

    @policy.setter
    def policy(self, policy: Any) -> None:
        self._policy_json = json.loads(policy)

    @staticmethod
    def cloudformation_name_type() -> str:
        return "TopicName"

    @staticmethod
    def cloudformation_type() -> str:
        # https://docs.aws.amazon.com/AWSCloudFormation/latest/UserGuide/aws-resource-sns-topic.html
        return "AWS::SNS::Topic"

    @classmethod
    def create_from_cloudformation_json(  # type: ignore[misc]
        cls,
        resource_name: str,
        cloudformation_json: Any,
        account_id: str,
        region_name: str,
        **kwargs: Any,
    ) -> "Topic":
        sns_backend = sns_backends[account_id][region_name]
        properties = cloudformation_json["Properties"]

        topic = sns_backend.create_topic(resource_name)
        for subscription in properties.get("Subscription", []):
            sns_backend.subscribe(
                topic.arn, subscription["Endpoint"], subscription["Protocol"]
            )
        return topic

    @classmethod
    def update_from_cloudformation_json(  # type: ignore[misc]
        cls,
        original_resource: "Topic",
        new_resource_name: str,
        cloudformation_json: Any,
        account_id: str,
        region_name: str,
    ) -> "Topic":
        original_resource.delete(account_id, region_name)
        return cls.create_from_cloudformation_json(
            new_resource_name, cloudformation_json, account_id, region_name
        )

    def delete(self, account_id: str, region_name: str) -> None:
        sns_backend: SNSBackend = sns_backends[account_id][region_name]
        sns_backend.delete_topic(self.arn)

    def _create_default_topic_policy(
        self, region_name: str, account_id: str, name: str
    ) -> Dict[str, Any]:
        return {
            "Version": "2008-10-17",
            "Id": "__default_policy_ID",
            "Statement": [
                {
                    "Effect": "Allow",
                    "Sid": "__default_statement_ID",
                    "Principal": {"AWS": "*"},
                    "Action": [
                        "SNS:GetTopicAttributes",
                        "SNS:SetTopicAttributes",
                        "SNS:AddPermission",
                        "SNS:RemovePermission",
                        "SNS:DeleteTopic",
                        "SNS:Subscribe",
                        "SNS:ListSubscriptionsByTopic",
                        "SNS:Publish",
                        "SNS:Receive",
                    ],
                    "Resource": make_arn_for_topic(self.account_id, name, region_name),
                    "Condition": {"StringEquals": {"AWS:SourceOwner": str(account_id)}},
                }
            ],
        }


class Subscription(BaseModel):
    def __init__(self, account_id: str, topic: Topic, endpoint: str, protocol: str):
        self.account_id = account_id
        self.topic = topic
        self.endpoint = endpoint
        self.protocol = protocol
        self.arn = make_arn_for_subscription(self.topic.arn)
        self.attributes: Dict[str, Any] = {}
        self._filter_policy = None  # filter policy as a dict, not json.
        self._filter_policy_matcher = None
        self.confirmed = False

    def publish(
        self,
        message: str,
        message_id: str,
        subject: Optional[str] = None,
        message_attributes: Optional[Dict[str, Any]] = None,
        group_id: Optional[str] = None,
        deduplication_id: Optional[str] = None,
    ) -> None:
        if self._filter_policy_matcher is not None:
            if not self._filter_policy_matcher.matches(message_attributes, message):
                return

        if self.protocol == "sqs":
            queue_name = self.endpoint.split(":")[-1]
            region = self.endpoint.split(":")[3]
            if self.attributes.get("RawMessageDelivery") != "true":
                sqs_backends[self.account_id][region].send_message(
                    queue_name,
                    json.dumps(
                        self.get_post_data(
                            message,
                            message_id,
                            subject,
                            message_attributes=message_attributes,
                        ),
                        sort_keys=True,
                        indent=2,
                        separators=(",", ": "),
                    ),
                    deduplication_id=deduplication_id,
                    group_id=group_id,
                )
            else:
                raw_message_attributes = {}
                for key, value in message_attributes.items():  # type: ignore
                    attr_type = "string_value"
                    type_value = value["Value"]
                    if value["Type"].startswith("Binary"):
                        attr_type = "binary_value"
                    elif value["Type"].startswith("Number"):
                        type_value = str(value["Value"])

                    raw_message_attributes[key] = {
                        "data_type": value["Type"],
                        attr_type: type_value,
                    }

                sqs_backends[self.account_id][region].send_message(
                    queue_name,
                    message,
                    message_attributes=raw_message_attributes,
                    deduplication_id=deduplication_id,
                    group_id=group_id,
                )
        elif self.protocol in ["http", "https"]:
            post_data = self.get_post_data(message, message_id, subject)
            requests.post(
                self.endpoint,
                json=post_data,
                headers={"Content-Type": "text/plain; charset=UTF-8"},
            )
        elif self.protocol == "lambda":
            # TODO: support bad function name
            # http://docs.aws.amazon.com/general/latest/gr/aws-arns-and-namespaces.html
            arr = self.endpoint.split(":")
            region = arr[3]
            qualifier = None
            if len(arr) == 7:
                assert arr[5] == "function"
                function_name = arr[-1]
            elif len(arr) == 8:
                assert arr[5] == "function"
                qualifier = arr[-1]
                function_name = arr[-2]
            else:
                assert False

            from moto.awslambda.utils import get_backend

            get_backend(self.account_id, region).send_sns_message(
                function_name, message, subject=subject, qualifier=qualifier
            )

    def get_post_data(
        self,
        message: str,
        message_id: str,
        subject: Optional[str],
        message_attributes: Optional[Dict[str, Any]] = None,
    ) -> Dict[str, Any]:
        post_data: Dict[str, Any] = {
            "Type": "Notification",
            "MessageId": message_id,
            "TopicArn": self.topic.arn,
            "Message": message,
            "Timestamp": iso_8601_datetime_with_milliseconds(),
            "SignatureVersion": "1",
            "Signature": "EXAMPLElDMXvB8r9R83tGoNn0ecwd5UjllzsvSvbItzfaMpN2nk5HVSw7XnOn/49IkxDKz8YrlH2qJXj2iZB0Zo2O71c4qQk1fMUDi3LGpij7RCW7AW9vYYsSqIKRnFS94ilu7NFhUzLiieYr4BKHpdTmdD6c0esKEYBpabxDSc=",
            "SigningCertURL": "https://sns.us-east-1.amazonaws.com/SimpleNotificationService-f3ecfb7224c7233fe7bb5f59f96de52f.pem",
            "UnsubscribeURL": f"https://sns.us-east-1.amazonaws.com/?Action=Unsubscribe&SubscriptionArn=arn:aws:sns:us-east-1:{self.account_id}:some-topic:2bcfbf39-05c3-41de-beaa-fcfcc21c8f55",
        }
        if subject:
            post_data["Subject"] = subject
        if message_attributes:
            post_data["MessageAttributes"] = message_attributes
        return post_data


class PlatformApplication(BaseModel):
    def __init__(
        self,
        account_id: str,
        region: str,
        name: str,
        platform: str,
        attributes: Dict[str, str],
    ):
        self.region = region
        self.name = name
        self.platform = platform
        self.attributes = attributes
        self.arn = f"arn:{get_partition(region)}:sns:{region}:{account_id}:app/{platform}/{name}"


class PlatformEndpoint(BaseModel):
    def __init__(
        self,
        account_id: str,
        region: str,
        application: PlatformApplication,
        custom_user_data: str,
        token: str,
        attributes: Dict[str, str],
    ):
        self.region = region
        self.application = application
        self.custom_user_data = custom_user_data
        self.token = token
        self.attributes = attributes
        self.id = mock_random.uuid4()
        self.arn = f"arn:{get_partition(region)}:sns:{region}:{account_id}:endpoint/{self.application.platform}/{self.application.name}/{self.id}"
        self.messages: Dict[str, str] = OrderedDict()
        self.__fixup_attributes()

    def __fixup_attributes(self) -> None:
        # When AWS returns the attributes dict, it always contains these two elements, so we need to
        # automatically ensure they exist as well.
        if "Token" not in self.attributes:
            self.attributes["Token"] = self.token
        if "Enabled" in self.attributes:
            enabled = self.attributes["Enabled"]
            self.attributes["Enabled"] = enabled.lower()
        else:
            self.attributes["Enabled"] = "true"

    @property
    def enabled(self) -> bool:
        return json.loads(self.attributes.get("Enabled", "true").lower())

    def publish(self, message: str) -> str:
        if not self.enabled:
            raise SnsEndpointDisabled("Endpoint is disabled")

        # This is where we would actually send a message
        message_id = str(mock_random.uuid4())
        self.messages[message_id] = message
        return message_id


class SNSBackend(BaseBackend):
    """
    Responsible for mocking calls to SNS. Integration with SQS/HTTP/etc is supported.

    Messages published to a topic are persisted in the backend. If you need to verify that a message was published successfully, you can use the internal API to check the message was published successfully:

    .. sourcecode:: python

        from moto.core import DEFAULT_ACCOUNT_ID
        from moto.sns import sns_backends
        sns_backend = sns_backends[DEFAULT_ACCOUNT_ID]["us-east-1"]  # Use the appropriate account/region
        all_send_notifications = sns_backend.topics[topic_arn].sent_notifications

    Note that, as this is an internal API, the exact format may differ per versions.
    """

    def __init__(self, region_name: str, account_id: str):
        super().__init__(region_name, account_id)
        self.topics: Dict[str, Topic] = OrderedDict()
        self.subscriptions: OrderedDict[str, Subscription] = OrderedDict()
        self.applications: Dict[str, PlatformApplication] = {}
        self.platform_endpoints: Dict[str, PlatformEndpoint] = {}
        self.region_name = region_name
        self.sms_attributes: Dict[str, str] = {}
        self.sms_messages: Dict[str, Tuple[str, str]] = OrderedDict()
        self.opt_out_numbers = [
            "+447420500600",
            "+447420505401",
            "+447632960543",
            "+447632960028",
            "+447700900149",
            "+447700900550",
            "+447700900545",
            "+447700900907",
        ]

    def get_sms_attributes(self, filter_list: Set[str]) -> Dict[str, str]:
        if len(filter_list) > 0:
            return {k: v for k, v in self.sms_attributes.items() if k in filter_list}
        else:
            return self.sms_attributes

    def set_sms_attributes(self, attrs: Dict[str, str]) -> None:
        self.sms_attributes.update(attrs)

    def create_topic(
        self,
        name: str,
        attributes: Optional[Dict[str, str]] = None,
        tags: Optional[Dict[str, str]] = None,
    ) -> Topic:
        if attributes is None:
            attributes = {}
        if attributes.get("FifoTopic") and attributes["FifoTopic"].lower() == "true":
            fails_constraints = not re.match(r"^[a-zA-Z0-9_-]{1,256}\.fifo$", name)
            msg = "Fifo Topic names must end with .fifo and must be made up of only uppercase and lowercase ASCII letters, numbers, underscores, and hyphens, and must be between 1 and 256 characters long."

        else:
            fails_constraints = not re.match(r"^[a-zA-Z0-9_-]{1,256}$", name)
            msg = "Topic names must be made up of only uppercase and lowercase ASCII letters, numbers, underscores, and hyphens, and must be between 1 and 256 characters long."

        if fails_constraints:
            raise InvalidParameterValue(msg)

        candidate_topic = Topic(name, self)
        if attributes:
            for attribute in attributes:
                setattr(
                    candidate_topic,
                    camelcase_to_underscores(attribute),
                    attributes[attribute],
                )
        if tags:
            candidate_topic._tags = tags
        if candidate_topic.arn in self.topics:
            return self.topics[candidate_topic.arn]
        else:
            self.topics[candidate_topic.arn] = candidate_topic
            return candidate_topic

    def _get_values_nexttoken(
        self, values_map: Dict[str, Any], next_token: Optional[str] = None
    ) -> Tuple[List[Any], Optional[int]]:
        i_next_token = int(next_token or "0")
        values = list(values_map.values())[
            i_next_token : i_next_token + DEFAULT_PAGE_SIZE
        ]
        if len(values) == DEFAULT_PAGE_SIZE:
            i_next_token = i_next_token + DEFAULT_PAGE_SIZE
        else:
            i_next_token = None  # type: ignore
        return values, i_next_token

    def _get_topic_subscriptions(self, topic: Topic) -> List[Subscription]:
        return [sub for sub in self.subscriptions.values() if sub.topic == topic]

    def list_topics(
        self, next_token: Optional[str] = None
    ) -> Tuple[List[Topic], Optional[int]]:
        return self._get_values_nexttoken(self.topics, next_token)

    def delete_topic_subscriptions(self, topic: Topic) -> None:
        for key, value in dict(self.subscriptions).items():
            if value.topic == topic:
                self.subscriptions.pop(key)

    def delete_topic(self, arn: str) -> None:
        with contextlib.suppress(TopicNotFound):
            topic = self.get_topic(arn)
            self.delete_topic_subscriptions(topic)
            parsed_arn = parse_arn(arn)
            sns_backends[parsed_arn.account][parsed_arn.region].topics.pop(arn, None)

    def get_topic(self, arn: str) -> Topic:
        parsed_arn = parse_arn(arn)
        try:
            return sns_backends[parsed_arn.account][self.region_name].topics[arn]
        except KeyError:
            raise TopicNotFound

    def set_topic_attribute(
        self, topic_arn: str, attribute_name: str, attribute_value: str
    ) -> None:
        topic = self.get_topic(topic_arn)
        setattr(topic, attribute_name, attribute_value)

    def subscribe(self, topic_arn: str, endpoint: str, protocol: str) -> Subscription:
        if protocol == "sms":
            if re.search(r"[./-]{2,}", endpoint) or re.search(
                r"(^[./-]|[./-]$)", endpoint
            ):
                raise SNSInvalidParameter(f"Invalid SMS endpoint: {endpoint}")

            reduced_endpoint = re.sub(r"[./-]", "", endpoint)

            if not is_e164(reduced_endpoint):
                raise SNSInvalidParameter(f"Invalid SMS endpoint: {endpoint}")

        # AWS doesn't create duplicates
        old_subscription = self._find_subscription(topic_arn, endpoint, protocol)
        if old_subscription:
            return old_subscription

        if protocol == "application":
            try:
                # Validate the endpoint exists, assuming it's a new subscription
                # Existing subscriptions are still found if an endpoint is deleted, at least for a short period
                self.get_endpoint(endpoint)
            except SNSNotFoundError:
                # Space between `arn{endpoint}` is lacking in AWS as well
                raise SNSInvalidParameter(
                    f"Invalid parameter: Endpoint Reason: Endpoint does not exist for endpoint arn{endpoint}"
                )

        topic = self.get_topic(topic_arn)
        subscription = Subscription(self.account_id, topic, endpoint, protocol)
        attributes = {
            "PendingConfirmation": "false",
            "ConfirmationWasAuthenticated": "true",
            "Endpoint": endpoint,
            "TopicArn": topic_arn,
            "Protocol": protocol,
            "SubscriptionArn": subscription.arn,
            "Owner": self.account_id,
            "RawMessageDelivery": "false",
        }

        if protocol in ["http", "https"]:
            attributes["EffectiveDeliveryPolicy"] = topic.effective_delivery_policy

        subscription.attributes = attributes
        self.subscriptions[subscription.arn] = subscription
        return subscription

    def _find_subscription(
        self, topic_arn: str, endpoint: str, protocol: str
    ) -> Optional[Subscription]:
        for subscription in self.subscriptions.values():
            if (
                subscription.topic.arn == topic_arn
                and subscription.endpoint == endpoint
                and subscription.protocol == protocol
            ):
                return subscription
        return None

    def unsubscribe(self, subscription_arn: str) -> None:
        self.subscriptions.pop(subscription_arn, None)

    def list_subscriptions(
        self, next_token: Optional[str] = None
    ) -> Tuple[List[Subscription], Optional[int]]:
        return self._get_values_nexttoken(self.subscriptions, next_token)

    def list_subscriptions_by_topic(
        self, topic_arn: str, next_token: Optional[str] = None
    ) -> Tuple[List[Subscription], Optional[int]]:
        topic = self.get_topic(topic_arn)
        filtered = OrderedDict(
            [(sub.arn, sub) for sub in self._get_topic_subscriptions(topic)]
        )
        return self._get_values_nexttoken(filtered, next_token)

    def publish(
        self,
        message: str,
        arn: Optional[str],
        phone_number: Optional[str] = None,
        subject: Optional[str] = None,
        message_attributes: Optional[Dict[str, Any]] = None,
        group_id: Optional[str] = None,
        deduplication_id: Optional[str] = None,
    ) -> str:
        if subject is not None and len(subject) > 100:
            # Note that the AWS docs around length are wrong: https://github.com/getmoto/moto/issues/1503
            raise ValueError("Subject must be less than 100 characters")

        if phone_number:
            # This is only an approximation. In fact, we should try to use GSM-7 or UCS-2 encoding to count used bytes
            if len(message) > MAXIMUM_SMS_MESSAGE_BYTES:
                raise ValueError("SMS message must be less than 1600 bytes")

            message_id = str(mock_random.uuid4())
            self.sms_messages[message_id] = (phone_number, message)
            return message_id

        if len(message) > MAXIMUM_MESSAGE_LENGTH:
            raise InvalidParameterValue(
                "An error occurred (InvalidParameter) when calling the Publish operation: Invalid parameter: Message too long"
            )

        try:
            topic = self.get_topic(arn)  # type: ignore

            fifo_topic = topic.fifo_topic == "true"
            if fifo_topic:
                if not group_id:
                    # MessageGroupId is a mandatory parameter for all
                    # messages in a fifo queue
                    raise MissingParameter("MessageGroupId")
                deduplication_id_required = topic.content_based_deduplication == "false"
                if not deduplication_id and deduplication_id_required:
                    raise InvalidParameterValue(
                        "The topic should either have ContentBasedDeduplication enabled or MessageDeduplicationId provided explicitly"
                    )
            elif group_id or deduplication_id:
                parameter = "MessageGroupId" if group_id else "MessageDeduplicationId"
                raise InvalidParameterValue(
                    f"Invalid parameter: {parameter} "
                    f"Reason: The request includes {parameter} parameter that is not valid for this topic type"
                )

            message_id = topic.publish(
                message,
                subject=subject,
                message_attributes=message_attributes,
                group_id=group_id,
                deduplication_id=deduplication_id,
            )
        except SNSNotFoundError:
            endpoint = self.get_endpoint(arn)  # type: ignore
            message_id = endpoint.publish(message)
        return message_id

    def create_platform_application(
        self, name: str, platform: str, attributes: Dict[str, str]
    ) -> PlatformApplication:
        application = PlatformApplication(
            self.account_id, self.region_name, name, platform, attributes
        )
        self.applications[application.arn] = application
        return application

    def get_application(self, arn: str) -> PlatformApplication:
        try:
            return self.applications[arn]
        except KeyError:
            raise SNSNotFoundError("PlatformApplication does not exist")

    def set_platform_application_attributes(
        self, arn: str, attributes: Dict[str, Any]
    ) -> PlatformApplication:
        application = self.get_application(arn)
        application.attributes.update(attributes)
        return application

    def list_platform_applications(self) -> Iterable[PlatformApplication]:
        return self.applications.values()

    def delete_platform_application(self, platform_arn: str) -> None:
        self.applications.pop(platform_arn)
        endpoints = self.list_endpoints_by_platform_application(platform_arn)
        for endpoint in endpoints:
            self.platform_endpoints.pop(endpoint.arn)

    def create_platform_endpoint(
        self,
        application: PlatformApplication,
        custom_user_data: str,
        token: str,
        attributes: Dict[str, str],
    ) -> PlatformEndpoint:
        for endpoint in self.platform_endpoints.values():
            if token == endpoint.token:
                same_user_data = custom_user_data == endpoint.custom_user_data
                same_attrs = (
                    attributes.get("Enabled", "true").lower()
                    == endpoint.attributes["Enabled"]
                )

                if same_user_data and same_attrs:
                    return endpoint
                raise DuplicateSnsEndpointError(
                    f"Invalid parameter: Token Reason: Endpoint {endpoint.arn} already exists with the same Token, but different attributes."
                )
        platform_endpoint = PlatformEndpoint(
            self.account_id,
            self.region_name,
            application,
            custom_user_data,
            token,
            attributes,
        )
        self.platform_endpoints[platform_endpoint.arn] = platform_endpoint
        return platform_endpoint

    def list_endpoints_by_platform_application(
        self, application_arn: str
    ) -> List[PlatformEndpoint]:
        return [
            endpoint
            for endpoint in self.platform_endpoints.values()
            if endpoint.application.arn == application_arn
        ]

    def get_endpoint(self, arn: str) -> PlatformEndpoint:
        try:
            return self.platform_endpoints[arn]
        except KeyError:
            raise SNSNotFoundError("Endpoint does not exist")

    def set_endpoint_attributes(
        self, arn: str, attributes: Dict[str, Any]
    ) -> PlatformEndpoint:
        endpoint = self.get_endpoint(arn)
        if "Enabled" in attributes:
            attributes["Enabled"] = attributes["Enabled"].lower()
        endpoint.attributes.update(attributes)
        return endpoint

    def delete_endpoint(self, arn: str) -> None:
        try:
            del self.platform_endpoints[arn]
        except KeyError:
            pass  # idempotent operation

    def get_subscription_attributes(self, arn: str) -> Dict[str, Any]:
        subscription = self.subscriptions.get(arn)

        if not subscription:
            raise SNSNotFoundError(
                "Subscription does not exist", template="wrapped_single_error"
            )
        # AWS does not return the FilterPolicy scope if the FilterPolicy is not set
        # if the FilterPolicy is set and not the FilterPolicyScope, it returns the default value
        attributes = {**subscription.attributes}
        if "FilterPolicyScope" in attributes and not attributes.get("FilterPolicy"):
            attributes.pop("FilterPolicyScope", None)
            attributes.pop("FilterPolicy", None)

        elif "FilterPolicy" in attributes and "FilterPolicyScope" not in attributes:
            attributes["FilterPolicyScope"] = "MessageAttributes"

        return attributes

    def set_subscription_attributes(self, arn: str, name: str, value: Any) -> None:
        if name not in [
            "RawMessageDelivery",
            "DeliveryPolicy",
            "FilterPolicy",
            "FilterPolicyScope",
            "RedrivePolicy",
            "SubscriptionRoleArn",
        ]:
            raise SNSInvalidParameter("AttributeName")

        # TODO: should do validation
        _subscription = [_ for _ in self.subscriptions.values() if _.arn == arn]
        if not _subscription:
            raise SNSNotFoundError(f"Subscription with arn {arn} not found")
        subscription = _subscription[0]

        if name == "FilterPolicy":
            if value:
                filter_policy = json.loads(value)
                # we validate the filter policy differently depending on the scope
                # we need to always set the scope first
                filter_policy_scope = subscription.attributes.get("FilterPolicyScope")
                self._validate_filter_policy(filter_policy, scope=filter_policy_scope)
                subscription._filter_policy = filter_policy
                subscription._filter_policy_matcher = FilterPolicyMatcher(
                    filter_policy, filter_policy_scope
                )
            else:
                subscription._filter_policy = None
                subscription._filter_policy_matcher = None

        subscription.attributes[name] = value

    def _validate_filter_policy(self, value: Any, scope: str) -> None:
        combinations = 1

        def aggregate_rules(
            filter_policy: Dict[str, Any], depth: int = 1
        ) -> List[List[Any]]:
            """
            This method evaluate the filter policy recursively, and returns only a list of lists of rules.
            It also calculates the combinations of rules, calculated depending on the nesting of the rules.
            Example:
            nested_filter_policy = {
                "key_a": {
                    "key_b": {
                        "key_c": ["value_one", "value_two", "value_three", "value_four"]
                    }
                },
                "key_d": {
                    "key_e": ["value_one", "value_two", "value_three"]
                }
            }
            This function then iterates on the values of the top level keys of the filter policy: ("key_a", "key_d")
            If the iterated value is not a list, it means it is a nested property. If the scope is `MessageBody`, it is
            allowed, we call this method on the value, adding a level to the depth to keep track on how deep the key is.
            If the value is a list, it means it contains rules: we will append this list of rules in _rules, and
            calculate the combinations it adds.
            For the example filter policy containing nested properties, we calculate it this way
            The first array has four values in a three-level nested key, and the second has three values in a two-level
            nested key. 3 x 4 x 2 x 3 = 72
            The return value would be:
            [["value_one", "value_two", "value_three", "value_four"], ["value_one", "value_two", "value_three"]]
            It allows us to later iterate of the list of rules in an easy way, to verify its conditions.

            :param filter_policy: a dict, starting at the FilterPolicy
            :param depth: the depth/level of the rules we are evaluating
            :return: a list of lists of rules
            """
            nonlocal combinations
            _rules = []
            for key, _value in filter_policy.items():
                if isinstance(_value, dict):
                    if scope == "MessageBody":
                        # From AWS docs: "unlike attribute-based policies, payload-based policies support property nesting."
                        _rules.extend(aggregate_rules(_value, depth=depth + 1))
                    else:
                        raise SNSInvalidParameter(
                            "Invalid parameter: Filter policy scope MessageAttributes does not support nested filter policy"
                        )
                elif isinstance(_value, list):
                    if key == "$or":
                        for val in _value:
                            _rules.extend(aggregate_rules(val, depth=depth + 1))
                    else:
                        _rules.append(_value)
                    combinations = combinations * len(_value) * depth
                else:
                    raise SNSInvalidParameter(
                        f'Invalid parameter: FilterPolicy: "{key}" must be an object or an array'
                    )
            return _rules

        # A filter policy can have a maximum of five attribute names. For a nested policy, only parent keys are counted.
        if len(value.values()) > 5:
            raise SNSInvalidParameter(
                "Invalid parameter: FilterPolicy: Filter policy can not have more than 5 keys"
            )

        aggregated_rules = aggregate_rules(value)
        # For the complexity of the filter policy, the total combination of values must not exceed 150.
        # https://docs.aws.amazon.com/sns/latest/dg/subscription-filter-policy-constraints.html
        if combinations > 150:
            raise SNSInvalidParameter(
                "Invalid parameter: FilterPolicy: Filter policy is too complex"
            )

        for rules in aggregated_rules:
            for rule in rules:
                if rule is None:
                    continue
                if isinstance(rule, str):
                    continue
                if isinstance(rule, bool):
                    continue
                if isinstance(rule, (int, float)):
                    if rule <= -1000000000 or rule >= 1000000000:
                        raise InternalError("Unknown")
                    continue
                if isinstance(rule, dict):
                    keyword = list(rule.keys())[0]
                    attributes = list(rule.values())[0]
                    if keyword in ["anything-but", "equals-ignore-case"]:
                        continue
                    elif keyword == "exists":
                        if not isinstance(attributes, bool):
                            raise SNSInvalidParameter(
                                "Invalid parameter: FilterPolicy: exists match pattern must be either true or false."
                            )
                        continue
                    elif keyword == "numeric":
                        # TODO: All of the exceptions listed below contain column pointing where the error is (in AWS response)
                        # Example: 'Value of < must be numeric\n at [Source: (String)"{"price":[{"numeric":["<","100"]}]}"; line: 1, column: 28]'
                        # While it probably can be implemented, it doesn't feel as important as the general parameter checking

                        attributes_copy = attributes[:]
                        if not attributes_copy:
                            raise SNSInvalidParameter(
                                "Invalid parameter: Attributes Reason: FilterPolicy: Invalid member in numeric match: ]\n at ..."
                            )

                        operator = attributes_copy.pop(0)

                        if not isinstance(operator, str):
                            raise SNSInvalidParameter(
                                f"Invalid parameter: Attributes Reason: FilterPolicy: Invalid member in numeric match: {(str(operator))}\n at ..."
                            )

                        if operator not in ("<", "<=", "=", ">", ">="):
                            raise SNSInvalidParameter(
                                f"Invalid parameter: Attributes Reason: FilterPolicy: Unrecognized numeric range operator: {(str(operator))}\n at ..."
                            )

                        try:
                            value = attributes_copy.pop(0)
                        except IndexError:
                            value = None

                        if value is None or not isinstance(value, (int, float)):
                            raise SNSInvalidParameter(
                                f"Invalid parameter: Attributes Reason: FilterPolicy: Value of {(str(operator))} must be numeric\n at ..."
                            )

                        if not attributes_copy:
                            continue

                        if operator not in (">", ">="):
                            raise SNSInvalidParameter(
                                "Invalid parameter: Attributes Reason: FilterPolicy: Too many elements in numeric expression\n at ..."
                            )

                        second_operator = attributes_copy.pop(0)

                        if second_operator not in ("<", "<="):
                            raise SNSInvalidParameter(
                                f"Invalid parameter: Attributes Reason: FilterPolicy: Bad numeric range operator: {(str(second_operator))}\n at ..."
                            )

                        try:
                            second_value = attributes_copy.pop(0)
                        except IndexError:
                            second_value = None

                        if second_value is None or not isinstance(
                            second_value, (int, float)
                        ):
                            raise SNSInvalidParameter(
                                f"Invalid parameter: Attributes Reason: FilterPolicy: Value of {(str(second_operator))} must be numeric\n at ..."
                            )

                        if second_value <= value:
                            raise SNSInvalidParameter(
                                "Invalid parameter: Attributes Reason: FilterPolicy: Bottom must be less than top\n at ..."
                            )

                        continue
                    elif keyword in ["prefix", "suffix"]:
                        continue
                    else:
                        raise SNSInvalidParameter(
                            f"Invalid parameter: FilterPolicy: Unrecognized match type {keyword}"
                        )

                raise SNSInvalidParameter(
                    "Invalid parameter: FilterPolicy: Match value must be String, number, true, false, or null"
                )

    def add_permission(
        self,
        region_name: str,
        topic_arn: str,
        label: str,
        aws_account_ids: List[str],
        action_names: List[str],
    ) -> None:
        topic = self.get_topic(topic_arn)
        policy = topic._policy_json
        statement = next(
            (
                statement
                for statement in policy["Statement"]
                if statement["Sid"] == label
            ),
            None,
        )

        if statement:
            raise SNSInvalidParameter("Statement already exists")

        if any(action_name not in VALID_POLICY_ACTIONS for action_name in action_names):
            raise SNSInvalidParameter("Policy statement action out of service scope!")

        principals = [
            f"arn:{get_partition(region_name)}:iam::{account_id}:root"
            for account_id in aws_account_ids
        ]
        actions = [f"SNS:{action_name}" for action_name in action_names]

        statement = {
            "Sid": label,
            "Effect": "Allow",
            "Principal": {"AWS": principals[0] if len(principals) == 1 else principals},
            "Action": actions[0] if len(actions) == 1 else actions,
            "Resource": topic_arn,
        }

        topic._policy_json["Statement"].append(statement)

    def remove_permission(self, topic_arn: str, label: str) -> None:
        topic = self.get_topic(topic_arn)
        statements = topic._policy_json["Statement"]
        statements = [
            statement for statement in statements if statement["Sid"] != label
        ]

        topic._policy_json["Statement"] = statements

    def list_tags_for_resource(self, resource_arn: str) -> Dict[str, str]:
        if resource_arn not in self.topics:
            raise ResourceNotFoundError

        return self.topics[resource_arn]._tags

    def tag_resource(self, resource_arn: str, tags: Dict[str, str]) -> None:
        if resource_arn not in self.topics:
            raise ResourceNotFoundError

        updated_tags = self.topics[resource_arn]._tags.copy()
        updated_tags.update(tags)

        if len(updated_tags) > 50:
            raise TagLimitExceededError

        self.topics[resource_arn]._tags = updated_tags

    def untag_resource(self, resource_arn: str, tag_keys: List[str]) -> None:
        if resource_arn not in self.topics:
            raise ResourceNotFoundError

        for key in tag_keys:
            self.topics[resource_arn]._tags.pop(key, None)

    def publish_batch(
        self, topic_arn: str, publish_batch_request_entries: List[Dict[str, Any]]
    ) -> Tuple[List[Dict[str, str]], List[Dict[str, Any]]]:
        """
        The MessageStructure and MessageDeduplicationId-parameters have not yet been implemented.
        """
        topic = self.get_topic(topic_arn)

        if len(publish_batch_request_entries) > 10:
            raise TooManyEntriesInBatchRequest

        ids = [m["Id"] for m in publish_batch_request_entries]
        if len(set(ids)) != len(ids):
            raise BatchEntryIdsNotDistinct

        fifo_topic = topic.fifo_topic == "true"
        if fifo_topic:
            if not all(
                ["MessageGroupId" in entry for entry in publish_batch_request_entries]
            ):
                raise SNSInvalidParameter(
                    "Invalid parameter: The MessageGroupId parameter is required for FIFO topics"
                )

        successful: List[Dict[str, str]] = []
        failed: List[Dict[str, Any]] = []

        for entry in publish_batch_request_entries:
            try:
                message_id = self.publish(
                    message=entry["Message"],
                    arn=topic_arn,
                    subject=entry.get("Subject"),
                    message_attributes=entry.get("MessageAttributes", {}),
                    group_id=entry.get("MessageGroupId"),
                    deduplication_id=entry.get("MessageDeduplicationId"),
                )
                successful.append({"MessageId": message_id, "Id": entry["Id"]})
            except Exception as e:
                if isinstance(e, InvalidParameterValue):
                    failed.append(
                        {
                            "Id": entry["Id"],
                            "Code": "InvalidParameter",
                            "Message": e.message,
                            "SenderFault": True,
                        }
                    )
        return successful, failed

    def check_if_phone_number_is_opted_out(self, number: str) -> bool:
        """
        Current implementation returns True for all numbers ending in '99'
        """
        return number.endswith("99")

    def list_phone_numbers_opted_out(self) -> List[str]:
        return self.opt_out_numbers

    def opt_in_phone_number(self, number: str) -> None:
        try:
            self.opt_out_numbers.remove(number)
        except ValueError:
            pass

    def confirm_subscription(self) -> None:
        pass

    def get_endpoint_attributes(self, arn: str) -> Dict[str, str]:
        endpoint = self.get_endpoint(arn)
        return endpoint.attributes

    def get_platform_application_attributes(self, arn: str) -> Dict[str, str]:
        application = self.get_application(arn)
        return application.attributes

    def get_topic_attributes(self) -> None:
        pass


sns_backends = BackendDict(SNSBackend, "sns")


DEFAULT_EFFECTIVE_DELIVERY_POLICY = {
    "defaultHealthyRetryPolicy": {
        "numNoDelayRetries": 0,
        "numMinDelayRetries": 0,
        "minDelayTarget": 20,
        "maxDelayTarget": 20,
        "numMaxDelayRetries": 0,
        "numRetries": 3,
        "backoffFunction": "linear",
    },
    "sicklyRetryPolicy": None,
    "throttlePolicy": None,
    "guaranteed": False,
}


VALID_POLICY_ACTIONS = [
    "GetTopicAttributes",
    "SetTopicAttributes",
    "AddPermission",
    "RemovePermission",
    "DeleteTopic",
    "Subscribe",
    "ListSubscriptionsByTopic",
    "Publish",
    "Receive",
]
