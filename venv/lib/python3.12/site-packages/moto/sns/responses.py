import base64
import re
from collections import defaultdict
from typing import Any, Dict

from moto.core.responses import TYPE_RESPONSE, ActionResult, BaseResponse, EmptyResult
from moto.core.utils import camelcase_to_underscores

from ..core.parsers import XFormedDict
from .exceptions import InvalidParameterValue, SNSInvalidParameter, SNSNotFoundError
from .models import SNSBackend, sns_backends
from .utils import is_e164


def transform_tags(tags: dict[str, str]) -> list[dict[str, str]]:
    return [{"Key": key, "Value": value} for key, value in tags.items()]


class SNSResponse(BaseResponse):
    SMS_ATTR_REGEX = re.compile(
        r"^attributes\.entry\.(?P<index>\d+)\.(?P<type>key|value)$"
    )
    OPT_OUT_PHONE_NUMBER_REGEX = re.compile(r"^\+?\d+$")
    RESPONSE_KEY_PATH_TO_TRANSFORMER = {
        "ListTagsForResourceResponse.Tags": transform_tags,
    }

    def __init__(self) -> None:
        super().__init__(service_name="sns")
        self.automated_parameter_parsing = True

    @property
    def backend(self) -> SNSBackend:
        return sns_backends[self.current_account][self.region]

    def _get_attributes(self) -> Dict[str, Any]:
        attributes = self._get_param("Attributes", XFormedDict())
        return attributes.original_dict()

    def _get_tags(self) -> Dict[str, str]:
        tags = self._get_param("Tags", [])
        return {tag["key"]: tag["value"] for tag in tags}

    def _parse_message_attributes(self) -> Dict[str, Any]:
        message_attributes = self._get_param("MessageAttributes", XFormedDict())
        message_attributes = message_attributes.original_dict()
        return self._transform_message_attributes(message_attributes)

    def _transform_message_attributes(
        self, message_attributes: Dict[str, Any]
    ) -> Dict[str, Any]:
        # SNS converts some key names before forwarding messages
        # DataType -> Type, StringValue -> Value, BinaryValue -> Value
        transformed_message_attributes = {}
        for name, value in message_attributes.items():
            # validation
            data_type = value["DataType"]
            if not data_type:
                raise InvalidParameterValue(
                    f"The message attribute '{name}' must contain non-empty message attribute value."
                )

            data_type_parts = data_type.split(".")
            if len(data_type_parts) > 2 or data_type_parts[0] not in [
                "String",
                "Binary",
                "Number",
            ]:
                raise InvalidParameterValue(
                    f"The message attribute '{name}' has an invalid message "
                    "attribute type, the set of supported type prefixes is "
                    "Binary, Number, and String."
                )

            transform_value = None
            if "StringValue" in value:
                transform_value = value["StringValue"]
                if data_type == "Number":
                    try:
                        int(transform_value)
                    except ValueError:
                        try:
                            float(transform_value)
                        except ValueError:
                            raise InvalidParameterValue(
                                "An error occurred (ParameterValueInvalid) "
                                "when calling the Publish operation: "
                                f"Could not cast message attribute '{name}' value to number."
                            )
            elif "BinaryValue" in value:
                transform_value = base64.b64encode(value["BinaryValue"]).decode()
            if transform_value == "":
                raise InvalidParameterValue(
                    f"The message attribute '{name}' must contain non-empty "
                    "message attribute value for message attribute "
                    f"type '{data_type[0]}'."
                )

            # transformation
            transformed_message_attributes[name] = {
                "Type": data_type,
                "Value": transform_value,
            }

        return transformed_message_attributes

    def create_topic(self) -> ActionResult:
        name = self._get_param("Name")
        attributes = self._get_attributes()
        tags = self._get_tags()
        topic = self.backend.create_topic(name, attributes, tags)
        result = {"TopicArn": topic.arn}
        return ActionResult(result)

    def list_topics(self) -> ActionResult:
        next_token = self._get_param("NextToken")
        topics, next_token = self.backend.list_topics(next_token=next_token)
        result = {"Topics": topics, "NextToken": next_token}
        return ActionResult(result)

    def delete_topic(self) -> ActionResult:
        topic_arn = self._get_param("TopicArn")
        self.backend.delete_topic(topic_arn)
        return EmptyResult()

    def get_topic_attributes(self) -> ActionResult:
        topic_arn = self._get_param("TopicArn")
        topic = self.backend.get_topic(topic_arn)
        attributes = {
            "Owner": topic.account_id,
            "Policy": topic.policy,
            "TopicArn": topic.arn,
            "DisplayName": topic.display_name,
            "SubscriptionsPending": topic.subscriptions_pending,
            "SubscriptionsConfirmed": topic.subscriptions_confimed,
            "SubscriptionsDeleted": topic.subscriptions_deleted,
            "DeliveryPolicy": topic.delivery_policy,
            "EffectiveDeliveryPolicy": topic.effective_delivery_policy,
        }
        if topic.kms_master_key_id:
            attributes["KmsMasterKeyId"] = topic.kms_master_key_id
        if topic.fifo_topic == "true":
            attributes["FifoTopic"] = topic.fifo_topic
            attributes["ContentBasedDeduplication"] = topic.content_based_deduplication
        result = {"Attributes": attributes}
        return ActionResult(result)

    def set_topic_attributes(self) -> ActionResult:
        topic_arn = self._get_param("TopicArn")
        attribute_name = self._get_param("AttributeName")
        attribute_name = camelcase_to_underscores(attribute_name)
        attribute_value = self._get_param("AttributeValue")
        self.backend.set_topic_attribute(topic_arn, attribute_name, attribute_value)
        return EmptyResult()

    def subscribe(self) -> ActionResult:
        topic_arn = self._get_param("TopicArn")
        endpoint = self._get_param("Endpoint")
        protocol = self._get_param("Protocol")
        attributes = self._get_attributes()

        subscription = self.backend.subscribe(topic_arn, endpoint, protocol)

        if attributes is not None:
            # We need to set the FilterPolicyScope first, as the validation of the FilterPolicy will depend on it
            if "FilterPolicyScope" in attributes:
                filter_policy_scope = attributes.pop("FilterPolicyScope")
                self.backend.set_subscription_attributes(
                    subscription.arn, "FilterPolicyScope", filter_policy_scope
                )

            for attr_name, attr_value in attributes.items():
                self.backend.set_subscription_attributes(
                    subscription.arn, attr_name, attr_value
                )

        result = {"SubscriptionArn": subscription.arn}
        return ActionResult(result)

    def unsubscribe(self) -> ActionResult:
        subscription_arn = self._get_param("SubscriptionArn")
        self.backend.unsubscribe(subscription_arn)
        return EmptyResult()

    def list_subscriptions(self) -> ActionResult:
        next_token = self._get_param("NextToken")
        subscriptions, next_token = self.backend.list_subscriptions(
            next_token=next_token
        )
        result = {"Subscriptions": subscriptions, "NextToken": next_token}
        return ActionResult(result)

    def list_subscriptions_by_topic(self) -> ActionResult:
        topic_arn = self._get_param("TopicArn")
        next_token = self._get_param("NextToken")
        subscriptions, next_token = self.backend.list_subscriptions_by_topic(
            topic_arn, next_token=next_token
        )
        result = {"Subscriptions": subscriptions, "NextToken": next_token}
        return ActionResult(result)

    def publish(self) -> ActionResult:
        target_arn = self._get_param("TargetArn")
        topic_arn = self._get_param("TopicArn")
        phone_number = self._get_param("PhoneNumber")
        subject = self._get_param("Subject")
        message_group_id = self._get_param("MessageGroupId")
        message_deduplication_id = self._get_param("MessageDeduplicationId")
        message_structure = self._get_param("MessageStructure")

        message_attributes = self._parse_message_attributes()

        arn = None
        if phone_number is not None:
            # Check phone is correct syntax (e164)
            if not is_e164(phone_number):
                raise SNSInvalidParameter("Phone number does not meet the E164 format")
        elif target_arn is not None:
            arn = target_arn
        else:
            arn = topic_arn

        message = self._get_param("Message")

        try:
            message_id = self.backend.publish(
                message,
                arn=arn,
                phone_number=phone_number,
                subject=subject,
                message_attributes=message_attributes,
                group_id=message_group_id,
                deduplication_id=message_deduplication_id,
                message_structure=message_structure,
            )
        except ValueError as err:
            raise SNSInvalidParameter(str(err))

        result = {"MessageId": message_id}
        return ActionResult(result)

    def publish_batch(self) -> ActionResult:
        topic_arn = self._get_param("TopicArn")
        publish_batch_request_entries = self._get_param(
            "PublishBatchRequestEntries", []
        )
        for entry in publish_batch_request_entries:
            if "MessageAttributes" in entry:
                # Use the same validation/processing as the regular publish-method
                entry["MessageAttributes"] = self._transform_message_attributes(
                    entry["MessageAttributes"]
                )
        successful, failed = self.backend.publish_batch(
            topic_arn=topic_arn,
            publish_batch_request_entries=publish_batch_request_entries,
        )
        result = {"Successful": successful, "Failed": failed}
        return ActionResult(result)

    def create_platform_application(self) -> ActionResult:
        name = self._get_param("Name")
        platform = self._get_param("Platform")
        attributes = self._get_attributes()
        platform_application = self.backend.create_platform_application(
            name, platform, attributes
        )
        result = {"PlatformApplicationArn": platform_application.arn}
        return ActionResult(result)

    def get_platform_application_attributes(self) -> ActionResult:
        arn = self._get_param("PlatformApplicationArn")
        attributes = self.backend.get_platform_application_attributes(arn)
        result = {"Attributes": attributes}
        return ActionResult(result)

    def set_platform_application_attributes(self) -> ActionResult:
        arn = self._get_param("PlatformApplicationArn")
        attributes = self._get_attributes()

        self.backend.set_platform_application_attributes(arn, attributes)
        return EmptyResult()

    def list_platform_applications(self) -> ActionResult:
        applications = self.backend.list_platform_applications()
        result = {"PlatformApplications": applications}
        return ActionResult(result)

    def delete_platform_application(self) -> ActionResult:
        platform_arn = self._get_param("PlatformApplicationArn")
        self.backend.delete_platform_application(platform_arn)
        return EmptyResult()

    def create_platform_endpoint(self) -> ActionResult:
        application_arn = self._get_param("PlatformApplicationArn")
        application = self.backend.get_application(application_arn)
        custom_user_data = self._get_param("CustomUserData")
        token = self._get_param("Token")
        attributes = self._get_attributes()
        platform_endpoint = self.backend.create_platform_endpoint(
            application, custom_user_data, token, attributes
        )
        result = {"EndpointArn": platform_endpoint.arn}
        return ActionResult(result)

    def list_endpoints_by_platform_application(self) -> ActionResult:
        application_arn = self._get_param("PlatformApplicationArn")
        self.backend.get_application(application_arn)
        endpoints = self.backend.list_endpoints_by_platform_application(application_arn)
        result = {"Endpoints": endpoints}
        return ActionResult(result)

    def get_endpoint_attributes(self) -> ActionResult:
        arn = self._get_param("EndpointArn")
        attributes = self.backend.get_endpoint_attributes(arn)
        result = {"Attributes": attributes}
        return ActionResult(result)

    def set_endpoint_attributes(self) -> ActionResult:
        arn = self._get_param("EndpointArn")
        attributes = self._get_attributes()
        self.backend.set_endpoint_attributes(arn, attributes)
        return EmptyResult()

    def delete_endpoint(self) -> ActionResult:
        arn = self._get_param("EndpointArn")
        self.backend.delete_endpoint(arn)
        return EmptyResult()

    def get_subscription_attributes(self) -> ActionResult:
        arn = self._get_param("SubscriptionArn")
        attributes = self.backend.get_subscription_attributes(arn)
        result = {"Attributes": attributes}
        return ActionResult(result)

    def set_subscription_attributes(self) -> ActionResult:
        arn = self._get_param("SubscriptionArn")
        attr_name = self._get_param("AttributeName")
        attr_value = self._get_param("AttributeValue")
        self.backend.set_subscription_attributes(arn, attr_name, attr_value)
        return EmptyResult()

    def set_sms_attributes(self) -> ActionResult:
        # attributes.entry.1.key
        # attributes.entry.1.value
        # to
        # 1: {key:X, value:Y}
        temp_dict: Dict[str, Any] = defaultdict(dict)
        for key, value in self.querystring.items():
            match = self.SMS_ATTR_REGEX.match(key)
            if match is not None:
                temp_dict[match.group("index")][match.group("type")] = value[0]

        # 1: {key:X, value:Y}
        # to
        # X: Y
        # All of this, just to take into account when people provide invalid stuff.
        result = {}
        for item in temp_dict.values():
            if "key" in item and "value" in item:
                result[item["key"]] = item["value"]

        self.backend.set_sms_attributes(result)

        return EmptyResult()

    def get_sms_attributes(self) -> ActionResult:
        filter_list = set()
        for key, value in self.querystring.items():
            if key.startswith("attributes.member.1"):
                filter_list.add(value[0])
        attributes = self.backend.get_sms_attributes(filter_list)
        result = {"attributes": attributes}
        return ActionResult(result)

    def check_if_phone_number_is_opted_out(
        self,
    ) -> ActionResult:
        number = self._get_param("phoneNumber")
        if self.OPT_OUT_PHONE_NUMBER_REGEX.match(number) is None:
            raise SNSInvalidParameter(
                "Invalid parameter: PhoneNumber Reason: input incorrectly formatted"
            )

        x = self.backend.check_if_phone_number_is_opted_out(number)
        result = {"isOptedOut": x}
        return ActionResult(result)

    def list_phone_numbers_opted_out(self) -> ActionResult:
        numbers = self.backend.list_phone_numbers_opted_out()
        result = {"phoneNumbers": numbers}
        return ActionResult(result)

    def opt_in_phone_number(self) -> ActionResult:
        number = self._get_param("phoneNumber")
        self.backend.opt_in_phone_number(number)
        return EmptyResult()

    def add_permission(self) -> ActionResult:
        topic_arn = self._get_param("TopicArn")
        label = self._get_param("Label")
        aws_account_ids = self._get_param("AWSAccountId", [])
        action_names = self._get_param("ActionName", [])
        self.backend.add_permission(
            region_name=self.region,
            topic_arn=topic_arn,
            label=label,
            aws_account_ids=aws_account_ids,
            action_names=action_names,
        )
        return EmptyResult()

    def remove_permission(self) -> ActionResult:
        topic_arn = self._get_param("TopicArn")
        label = self._get_param("Label")
        self.backend.remove_permission(topic_arn, label)
        return EmptyResult()

    def confirm_subscription(self) -> ActionResult:
        arn = self._get_param("TopicArn")

        if arn not in self.backend.topics:
            raise SNSNotFoundError("Topic does not exist")

        # Once Tokens are stored by the `subscribe` endpoint and distributed
        # to the client somehow, then we can check validity of tokens
        # presented to this method. The following code works, all thats
        # needed is to perform a token check and assign that value to the
        # `already_subscribed` variable.
        #
        # token = self._get_param('Token')
        # auth = self._get_param('AuthenticateOnUnsubscribe')
        # if already_subscribed:
        #     error_response = self._error(
        #         code='AuthorizationError',
        #         message='Subscription already confirmed'
        #     )
        #     return error_response, dict(status=400)

        result = {"SubscriptionArn": f"{arn}:68762e72-e9b1-410a-8b3b-903da69ee1d5"}
        return ActionResult(result)

    def list_tags_for_resource(self) -> ActionResult:
        arn = self._get_param("ResourceArn")
        tags = self.backend.list_tags_for_resource(arn)
        result = {"Tags": tags}
        return ActionResult(result)

    def tag_resource(self) -> ActionResult:
        arn = self._get_param("ResourceArn")
        tags = self._get_tags()
        self.backend.tag_resource(arn, tags)
        return EmptyResult()

    def untag_resource(self) -> ActionResult:
        arn = self._get_param("ResourceArn")
        tag_keys = self._get_param("TagKeys", [])
        self.backend.untag_resource(arn, tag_keys)
        return EmptyResult()

    @staticmethod
    def serve_pem(request: Any, full_url: str, headers: Any) -> TYPE_RESPONSE:
        sns = SNSResponse()
        sns.automated_parameter_parsing = False
        sns.setup_class(request, full_url, headers)
        key_name = full_url.split("/")[-1]
        key = sns.backend._message_public_keys[key_name]
        return 200, {}, key
