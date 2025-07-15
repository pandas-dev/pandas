import json
import re
from collections import defaultdict
from typing import Any, Dict, Tuple, Union

from moto.core.responses import TYPE_RESPONSE, BaseResponse
from moto.core.utils import camelcase_to_underscores

from .exceptions import InvalidParameterValue, SNSNotFoundError
from .models import SNSBackend, sns_backends
from .utils import is_e164


class SNSResponse(BaseResponse):
    SMS_ATTR_REGEX = re.compile(
        r"^attributes\.entry\.(?P<index>\d+)\.(?P<type>key|value)$"
    )
    OPT_OUT_PHONE_NUMBER_REGEX = re.compile(r"^\+?\d+$")

    def __init__(self) -> None:
        super().__init__(service_name="sns")

    @property
    def backend(self) -> SNSBackend:
        return sns_backends[self.current_account][self.region]

    def _error(self, code: str, message: str, sender: str = "Sender") -> str:
        template = self.response_template(ERROR_RESPONSE)
        return template.render(code=code, message=message, sender=sender)

    def _get_attributes(self) -> Dict[str, str]:
        attributes = self._get_list_prefix("Attributes.entry")
        return dict((attribute["key"], attribute["value"]) for attribute in attributes)

    def _get_tags(self) -> Dict[str, str]:
        tags = self._get_list_prefix("Tags.member")
        return {tag["key"]: tag["value"] for tag in tags}

    def _parse_message_attributes(self) -> Dict[str, Any]:
        message_attributes = self._get_object_map(
            "MessageAttributes.entry", name="Name", value="Value"
        )
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
                transform_value = value["BinaryValue"]
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

    def create_topic(self) -> str:
        name = self._get_param("Name")
        attributes = self._get_attributes()
        tags = self._get_tags()
        topic = self.backend.create_topic(name, attributes, tags)

        if self.request_json:
            return json.dumps(
                {
                    "CreateTopicResponse": {
                        "CreateTopicResult": {"TopicArn": topic.arn},
                        "ResponseMetadata": {
                            "RequestId": "a8dec8b3-33a4-11df-8963-01868b7c937a"
                        },
                    }
                }
            )

        template = self.response_template(CREATE_TOPIC_TEMPLATE)
        return template.render(topic=topic)

    def list_topics(self) -> str:
        next_token = self._get_param("NextToken")
        topics, next_token = self.backend.list_topics(next_token=next_token)

        if self.request_json:
            return json.dumps(
                {
                    "ListTopicsResponse": {
                        "ListTopicsResult": {
                            "Topics": [{"TopicArn": topic.arn} for topic in topics],
                            "NextToken": next_token,
                        }
                    },
                    "ResponseMetadata": {
                        "RequestId": "a8dec8b3-33a4-11df-8963-01868b7c937a"
                    },
                }
            )

        template = self.response_template(LIST_TOPICS_TEMPLATE)
        return template.render(topics=topics, next_token=next_token)

    def delete_topic(self) -> str:
        topic_arn = self._get_param("TopicArn")
        self.backend.delete_topic(topic_arn)

        if self.request_json:
            return json.dumps(
                {
                    "DeleteTopicResponse": {
                        "ResponseMetadata": {
                            "RequestId": "a8dec8b3-33a4-11df-8963-01868b7c937a"
                        }
                    }
                }
            )

        template = self.response_template(DELETE_TOPIC_TEMPLATE)
        return template.render()

    def get_topic_attributes(self) -> str:
        topic_arn = self._get_param("TopicArn")
        topic = self.backend.get_topic(topic_arn)

        if self.request_json:
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
                attributes["ContentBasedDeduplication"] = (
                    topic.content_based_deduplication
                )
            response = {
                "GetTopicAttributesResponse": {
                    "GetTopicAttributesResult": {"Attributes": attributes},
                    "ResponseMetadata": {
                        "RequestId": "057f074c-33a7-11df-9540-99d0768312d3"
                    },
                }
            }
            return json.dumps(response)

        template = self.response_template(GET_TOPIC_ATTRIBUTES_TEMPLATE)
        return template.render(topic=topic)

    def set_topic_attributes(self) -> str:
        topic_arn = self._get_param("TopicArn")
        attribute_name = self._get_param("AttributeName")
        attribute_name = camelcase_to_underscores(attribute_name)
        attribute_value = self._get_param("AttributeValue")
        self.backend.set_topic_attribute(topic_arn, attribute_name, attribute_value)

        if self.request_json:
            return json.dumps(
                {
                    "SetTopicAttributesResponse": {
                        "ResponseMetadata": {
                            "RequestId": "a8763b99-33a7-11df-a9b7-05d48da6f042"
                        }
                    }
                }
            )

        template = self.response_template(SET_TOPIC_ATTRIBUTES_TEMPLATE)
        return template.render()

    def subscribe(self) -> str:
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

        if self.request_json:
            return json.dumps(
                {
                    "SubscribeResponse": {
                        "SubscribeResult": {"SubscriptionArn": subscription.arn},
                        "ResponseMetadata": {
                            "RequestId": "a8763b99-33a7-11df-a9b7-05d48da6f042"
                        },
                    }
                }
            )

        template = self.response_template(SUBSCRIBE_TEMPLATE)
        return template.render(subscription=subscription)

    def unsubscribe(self) -> str:
        subscription_arn = self._get_param("SubscriptionArn")
        self.backend.unsubscribe(subscription_arn)

        if self.request_json:
            return json.dumps(
                {
                    "UnsubscribeResponse": {
                        "ResponseMetadata": {
                            "RequestId": "a8763b99-33a7-11df-a9b7-05d48da6f042"
                        }
                    }
                }
            )

        template = self.response_template(UNSUBSCRIBE_TEMPLATE)
        return template.render()

    def list_subscriptions(self) -> str:
        next_token = self._get_param("NextToken")
        subscriptions, next_token = self.backend.list_subscriptions(
            next_token=next_token
        )

        if self.request_json:
            return json.dumps(
                {
                    "ListSubscriptionsResponse": {
                        "ListSubscriptionsResult": {
                            "Subscriptions": [
                                {
                                    "TopicArn": subscription.topic.arn,
                                    "Protocol": subscription.protocol,
                                    "SubscriptionArn": subscription.arn,
                                    "Owner": subscription.topic.account_id,
                                    "Endpoint": subscription.endpoint,
                                }
                                for subscription in subscriptions
                            ],
                            "NextToken": next_token,
                        },
                        "ResponseMetadata": {
                            "RequestId": "384ac68d-3775-11df-8963-01868b7c937a"
                        },
                    }
                }
            )

        template = self.response_template(LIST_SUBSCRIPTIONS_TEMPLATE)
        return template.render(subscriptions=subscriptions, next_token=next_token)

    def list_subscriptions_by_topic(self) -> str:
        topic_arn = self._get_param("TopicArn")
        next_token = self._get_param("NextToken")
        subscriptions, next_token = self.backend.list_subscriptions_by_topic(
            topic_arn, next_token=next_token
        )

        if self.request_json:
            return json.dumps(
                {
                    "ListSubscriptionsByTopicResponse": {
                        "ListSubscriptionsByTopicResult": {
                            "Subscriptions": [
                                {
                                    "TopicArn": subscription.topic.arn,
                                    "Protocol": subscription.protocol,
                                    "SubscriptionArn": subscription.arn,
                                    "Owner": subscription.topic.account_id,
                                    "Endpoint": subscription.endpoint,
                                }
                                for subscription in subscriptions
                            ],
                            "NextToken": next_token,
                        },
                        "ResponseMetadata": {
                            "RequestId": "384ac68d-3775-11df-8963-01868b7c937a"
                        },
                    }
                }
            )

        template = self.response_template(LIST_SUBSCRIPTIONS_BY_TOPIC_TEMPLATE)
        return template.render(subscriptions=subscriptions, next_token=next_token)

    def publish(self) -> Union[str, Tuple[str, Dict[str, int]]]:
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
                return (
                    self._error(
                        "InvalidParameter", "Phone number does not meet the E164 format"
                    ),
                    dict(status=400),
                )
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
            error_response = self._error("InvalidParameter", str(err))
            return error_response, dict(status=400)

        if self.request_json:
            return json.dumps(
                {
                    "PublishResponse": {
                        "PublishResult": {"MessageId": message_id},
                        "ResponseMetadata": {
                            "RequestId": "384ac68d-3775-11df-8963-01868b7c937a"
                        },
                    }
                }
            )

        template = self.response_template(PUBLISH_TEMPLATE)
        return template.render(message_id=message_id)

    def publish_batch(self) -> str:
        topic_arn = self._get_param("TopicArn")
        publish_batch_request_entries = self._get_multi_param(
            "PublishBatchRequestEntries.member"
        )
        for entry in publish_batch_request_entries:
            if "MessageAttributes" in entry:
                # Convert into the same format as the regular publish-method
                # FROM: [{'Name': 'a', 'Value': {'DataType': 'String', 'StringValue': 'v'}}]
                # TO  : {'name': {'DataType': 'Number', 'StringValue': '123'}}
                msg_attrs = {y["Name"]: y["Value"] for y in entry["MessageAttributes"]}
                # Use the same validation/processing as the regular publish-method
                entry["MessageAttributes"] = self._transform_message_attributes(
                    msg_attrs
                )
        successful, failed = self.backend.publish_batch(
            topic_arn=topic_arn,
            publish_batch_request_entries=publish_batch_request_entries,
        )
        template = self.response_template(PUBLISH_BATCH_TEMPLATE)
        return template.render(successful=successful, failed=failed)

    def create_platform_application(self) -> str:
        name = self._get_param("Name")
        platform = self._get_param("Platform")
        attributes = self._get_attributes()
        platform_application = self.backend.create_platform_application(
            name, platform, attributes
        )

        if self.request_json:
            return json.dumps(
                {
                    "CreatePlatformApplicationResponse": {
                        "CreatePlatformApplicationResult": {
                            "PlatformApplicationArn": platform_application.arn
                        },
                        "ResponseMetadata": {
                            "RequestId": "384ac68d-3775-11df-8963-01868b7c937b"
                        },
                    }
                }
            )

        template = self.response_template(CREATE_PLATFORM_APPLICATION_TEMPLATE)
        return template.render(platform_application=platform_application)

    def get_platform_application_attributes(self) -> str:
        arn = self._get_param("PlatformApplicationArn")
        attributes = self.backend.get_platform_application_attributes(arn)

        if self.request_json:
            return json.dumps(
                {
                    "GetPlatformApplicationAttributesResponse": {
                        "GetPlatformApplicationAttributesResult": {
                            "Attributes": attributes
                        },
                        "ResponseMetadata": {
                            "RequestId": "384ac68d-3775-11df-8963-01868b7c937f"
                        },
                    }
                }
            )

        template = self.response_template(GET_PLATFORM_APPLICATION_ATTRIBUTES_TEMPLATE)
        return template.render(attributes=attributes)

    def set_platform_application_attributes(self) -> str:
        arn = self._get_param("PlatformApplicationArn")
        attributes = self._get_attributes()

        self.backend.set_platform_application_attributes(arn, attributes)

        if self.request_json:
            return json.dumps(
                {
                    "SetPlatformApplicationAttributesResponse": {
                        "ResponseMetadata": {
                            "RequestId": "384ac68d-3775-12df-8963-01868b7c937f"
                        }
                    }
                }
            )

        template = self.response_template(SET_PLATFORM_APPLICATION_ATTRIBUTES_TEMPLATE)
        return template.render()

    def list_platform_applications(self) -> str:
        applications = self.backend.list_platform_applications()

        if self.request_json:
            return json.dumps(
                {
                    "ListPlatformApplicationsResponse": {
                        "ListPlatformApplicationsResult": {
                            "PlatformApplications": [
                                {
                                    "PlatformApplicationArn": application.arn,
                                    "attributes": application.attributes,
                                }
                                for application in applications
                            ],
                            "NextToken": None,
                        },
                        "ResponseMetadata": {
                            "RequestId": "384ac68d-3775-11df-8963-01868b7c937c"
                        },
                    }
                }
            )

        template = self.response_template(LIST_PLATFORM_APPLICATIONS_TEMPLATE)
        return template.render(applications=applications)

    def delete_platform_application(self) -> str:
        platform_arn = self._get_param("PlatformApplicationArn")
        self.backend.delete_platform_application(platform_arn)

        if self.request_json:
            return json.dumps(
                {
                    "DeletePlatformApplicationResponse": {
                        "ResponseMetadata": {
                            "RequestId": "384ac68d-3775-11df-8963-01868b7c937e"
                        }
                    }
                }
            )

        template = self.response_template(DELETE_PLATFORM_APPLICATION_TEMPLATE)
        return template.render()

    def create_platform_endpoint(self) -> str:
        application_arn = self._get_param("PlatformApplicationArn")
        application = self.backend.get_application(application_arn)

        custom_user_data = self._get_param("CustomUserData")
        token = self._get_param("Token")
        attributes = self._get_attributes()

        platform_endpoint = self.backend.create_platform_endpoint(
            application, custom_user_data, token, attributes
        )

        if self.request_json:
            return json.dumps(
                {
                    "CreatePlatformEndpointResponse": {
                        "CreatePlatformEndpointResult": {
                            "EndpointArn": platform_endpoint.arn
                        },
                        "ResponseMetadata": {
                            "RequestId": "384ac68d-3779-11df-8963-01868b7c937b"
                        },
                    }
                }
            )

        template = self.response_template(CREATE_PLATFORM_ENDPOINT_TEMPLATE)
        return template.render(platform_endpoint=platform_endpoint)

    def list_endpoints_by_platform_application(self) -> str:
        application_arn = self._get_param("PlatformApplicationArn")
        self.backend.get_application(application_arn)
        endpoints = self.backend.list_endpoints_by_platform_application(application_arn)

        if self.request_json:
            return json.dumps(
                {
                    "ListEndpointsByPlatformApplicationResponse": {
                        "ListEndpointsByPlatformApplicationResult": {
                            "Endpoints": [
                                {
                                    "Attributes": endpoint.attributes,
                                    "EndpointArn": endpoint.arn,
                                }
                                for endpoint in endpoints
                            ],
                            "NextToken": None,
                        },
                        "ResponseMetadata": {
                            "RequestId": "384ac68d-3775-11df-8963-01868b7c937a"
                        },
                    }
                }
            )

        template = self.response_template(
            LIST_ENDPOINTS_BY_PLATFORM_APPLICATION_TEMPLATE
        )
        return template.render(endpoints=endpoints)

    def get_endpoint_attributes(self) -> Union[str, Tuple[str, Dict[str, int]]]:
        arn = self._get_param("EndpointArn")
        try:
            attributes = self.backend.get_endpoint_attributes(arn)

            if self.request_json:
                return json.dumps(
                    {
                        "GetEndpointAttributesResponse": {
                            "GetEndpointAttributesResult": {"Attributes": attributes},
                            "ResponseMetadata": {
                                "RequestId": "384ac68d-3775-11df-8963-01868b7c937f"
                            },
                        }
                    }
                )

            template = self.response_template(GET_ENDPOINT_ATTRIBUTES_TEMPLATE)
            return template.render(attributes=attributes)
        except SNSNotFoundError:
            error_response = self._error("NotFound", "Endpoint does not exist")
            return error_response, dict(status=404)

    def set_endpoint_attributes(self) -> Union[str, Tuple[str, Dict[str, int]]]:
        arn = self._get_param("EndpointArn")
        attributes = self._get_attributes()

        self.backend.set_endpoint_attributes(arn, attributes)

        if self.request_json:
            return json.dumps(
                {
                    "SetEndpointAttributesResponse": {
                        "ResponseMetadata": {
                            "RequestId": "384bc68d-3775-12df-8963-01868b7c937f"
                        }
                    }
                }
            )

        template = self.response_template(SET_ENDPOINT_ATTRIBUTES_TEMPLATE)
        return template.render()

    def delete_endpoint(self) -> str:
        arn = self._get_param("EndpointArn")
        self.backend.delete_endpoint(arn)

        if self.request_json:
            return json.dumps(
                {
                    "DeleteEndpointResponse": {
                        "ResponseMetadata": {
                            "RequestId": "384bc68d-3775-12df-8963-01868b7c937f"
                        }
                    }
                }
            )

        template = self.response_template(DELETE_ENDPOINT_TEMPLATE)
        return template.render()

    def get_subscription_attributes(self) -> str:
        arn = self._get_param("SubscriptionArn")
        attributes = self.backend.get_subscription_attributes(arn)
        template = self.response_template(GET_SUBSCRIPTION_ATTRIBUTES_TEMPLATE)
        return template.render(attributes=attributes)

    def set_subscription_attributes(self) -> str:
        arn = self._get_param("SubscriptionArn")
        attr_name = self._get_param("AttributeName")
        attr_value = self._get_param("AttributeValue")
        self.backend.set_subscription_attributes(arn, attr_name, attr_value)
        template = self.response_template(SET_SUBSCRIPTION_ATTRIBUTES_TEMPLATE)
        return template.render()

    def set_sms_attributes(self) -> str:
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

        template = self.response_template(SET_SMS_ATTRIBUTES_TEMPLATE)
        return template.render()

    def get_sms_attributes(self) -> str:
        filter_list = set()
        for key, value in self.querystring.items():
            if key.startswith("attributes.member.1"):
                filter_list.add(value[0])

        result = self.backend.get_sms_attributes(filter_list)

        template = self.response_template(GET_SMS_ATTRIBUTES_TEMPLATE)
        return template.render(attributes=result)

    def check_if_phone_number_is_opted_out(
        self,
    ) -> Union[str, Tuple[str, Dict[str, int]]]:
        number = self._get_param("phoneNumber")
        if self.OPT_OUT_PHONE_NUMBER_REGEX.match(number) is None:
            error_response = self._error(
                code="InvalidParameter",
                message="Invalid parameter: PhoneNumber Reason: input incorrectly formatted",
            )
            return error_response, dict(status=400)

        x = self.backend.check_if_phone_number_is_opted_out(number)
        template = self.response_template(CHECK_IF_OPTED_OUT_TEMPLATE)
        return template.render(opt_out=str(x).lower())

    def list_phone_numbers_opted_out(self) -> str:
        numbers = self.backend.list_phone_numbers_opted_out()
        template = self.response_template(LIST_OPTOUT_TEMPLATE)
        return template.render(opt_outs=numbers)

    def opt_in_phone_number(self) -> str:
        number = self._get_param("phoneNumber")

        self.backend.opt_in_phone_number(number)

        template = self.response_template(OPT_IN_NUMBER_TEMPLATE)
        return template.render()

    def add_permission(self) -> str:
        topic_arn = self._get_param("TopicArn")
        label = self._get_param("Label")
        aws_account_ids = self._get_multi_param("AWSAccountId.member.")
        action_names = self._get_multi_param("ActionName.member.")

        self.backend.add_permission(
            region_name=self.region,
            topic_arn=topic_arn,
            label=label,
            aws_account_ids=aws_account_ids,
            action_names=action_names,
        )

        template = self.response_template(ADD_PERMISSION_TEMPLATE)
        return template.render()

    def remove_permission(self) -> str:
        topic_arn = self._get_param("TopicArn")
        label = self._get_param("Label")

        self.backend.remove_permission(topic_arn, label)

        template = self.response_template(DEL_PERMISSION_TEMPLATE)
        return template.render()

    def confirm_subscription(self) -> Union[str, Tuple[str, Dict[str, int]]]:
        arn = self._get_param("TopicArn")

        if arn not in self.backend.topics:
            error_response = self._error("NotFound", "Topic does not exist")
            return error_response, dict(status=404)

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

        template = self.response_template(CONFIRM_SUBSCRIPTION_TEMPLATE)
        return template.render(sub_arn=f"{arn}:68762e72-e9b1-410a-8b3b-903da69ee1d5")

    def list_tags_for_resource(self) -> str:
        arn = self._get_param("ResourceArn")

        result = self.backend.list_tags_for_resource(arn)

        template = self.response_template(LIST_TAGS_FOR_RESOURCE_TEMPLATE)
        return template.render(tags=result)

    def tag_resource(self) -> str:
        arn = self._get_param("ResourceArn")
        tags = self._get_tags()

        self.backend.tag_resource(arn, tags)

        return self.response_template(TAG_RESOURCE_TEMPLATE).render()

    def untag_resource(self) -> str:
        arn = self._get_param("ResourceArn")
        tag_keys = self._get_multi_param("TagKeys.member")

        self.backend.untag_resource(arn, tag_keys)

        return self.response_template(UNTAG_RESOURCE_TEMPLATE).render()

    @staticmethod
    def serve_pem(request: Any, full_url: str, headers: Any) -> TYPE_RESPONSE:  # type: ignore[misc]
        sns = SNSResponse()
        sns.setup_class(request, full_url, headers)
        key_name = full_url.split("/")[-1]
        key = sns.backend._message_public_keys[key_name]
        return 200, {}, key


CREATE_TOPIC_TEMPLATE = """<CreateTopicResponse xmlns="http://sns.amazonaws.com/doc/2010-03-31/">
     <CreateTopicResult>
       <TopicArn>{{ topic.arn }}</TopicArn>
     </CreateTopicResult>
     <ResponseMetadata>
       <RequestId>a8dec8b3-33a4-11df-8963-01868b7c937a</RequestId>
     </ResponseMetadata>
   </CreateTopicResponse>"""

LIST_TOPICS_TEMPLATE = """<ListTopicsResponse xmlns="http://sns.amazonaws.com/doc/2010-03-31/">
  <ListTopicsResult>
    <Topics>
    {% for topic in topics %}
      <member>
        <TopicArn>{{ topic.arn }}</TopicArn>
      </member>
    {% endfor %}
    </Topics>
    {% if next_token  %}
    <NextToken>{{ next_token }}</NextToken>
    {% endif %}
  </ListTopicsResult>
  <ResponseMetadata>
    <RequestId>3f1478c7-33a9-11df-9540-99d0768312d3</RequestId>
  </ResponseMetadata>
</ListTopicsResponse>"""

DELETE_TOPIC_TEMPLATE = """<DeleteTopicResponse xmlns="http://sns.amazonaws.com/doc/2010-03-31/">
  <ResponseMetadata>
    <RequestId>f3aa9ac9-3c3d-11df-8235-9dab105e9c32</RequestId>
  </ResponseMetadata>
</DeleteTopicResponse>"""

GET_TOPIC_ATTRIBUTES_TEMPLATE = """<GetTopicAttributesResponse xmlns="http://sns.amazonaws.com/doc/2010-03-31/">
  <GetTopicAttributesResult>
    <Attributes>
      <entry>
        <key>Owner</key>
        <value>{{ topic.account_id }}</value>
      </entry>
      <entry>
        <key>Policy</key>
        <value>{{ topic.policy }}</value>
      </entry>
      <entry>
        <key>TopicArn</key>
        <value>{{ topic.arn }}</value>
      </entry>
      <entry>
        <key>DisplayName</key>
        <value>{{ topic.display_name }}</value>
      </entry>
      <entry>
        <key>SubscriptionsPending</key>
        <value>{{ topic.subscriptions_pending }}</value>
      </entry>
      <entry>
        <key>SubscriptionsConfirmed</key>
        <value>{{ topic.subscriptions_confimed }}</value>
      </entry>
      <entry>
        <key>SubscriptionsDeleted</key>
        <value>{{ topic.subscriptions_deleted }}</value>
      </entry>
      <entry>
        <key>DeliveryPolicy</key>
        <value>{{ topic.delivery_policy }}</value>
      </entry>
      <entry>
        <key>EffectiveDeliveryPolicy</key>
        <value>{{ topic.effective_delivery_policy }}</value>
      </entry>
      {% if topic.kms_master_key_id %}
      <entry>
        <key>KmsMasterKeyId</key>
        <value>{{ topic.kms_master_key_id }}</value>
      </entry>
      {% endif %}
      {% if topic.fifo_topic == 'true' %}
      <entry>
        <key>FifoTopic</key>
        <value>{{ topic.fifo_topic }}</value>
      </entry>
      <entry>
        <key>ContentBasedDeduplication</key>
        <value>{{ topic.content_based_deduplication }}</value>
      </entry>
      {% endif %}
    </Attributes>
  </GetTopicAttributesResult>
  <ResponseMetadata>
    <RequestId>057f074c-33a7-11df-9540-99d0768312d3</RequestId>
  </ResponseMetadata>
</GetTopicAttributesResponse>"""

SET_TOPIC_ATTRIBUTES_TEMPLATE = """<SetTopicAttributesResponse xmlns="http://sns.amazonaws.com/doc/2010-03-31/">
  <ResponseMetadata>
    <RequestId>a8763b99-33a7-11df-a9b7-05d48da6f042</RequestId>
  </ResponseMetadata>
</SetTopicAttributesResponse>"""

CREATE_PLATFORM_APPLICATION_TEMPLATE = """<CreatePlatformApplicationResponse xmlns="http://sns.amazonaws.com/doc/2010-03-31/">
  <CreatePlatformApplicationResult>
    <PlatformApplicationArn>{{ platform_application.arn }}</PlatformApplicationArn>
  </CreatePlatformApplicationResult>
  <ResponseMetadata>
    <RequestId>b6f0e78b-e9d4-5a0e-b973-adc04e8a4ff9</RequestId>
  </ResponseMetadata>
</CreatePlatformApplicationResponse>"""

CREATE_PLATFORM_ENDPOINT_TEMPLATE = """<CreatePlatformEndpointResponse xmlns="http://sns.amazonaws.com/doc/2010-03-31/">
  <CreatePlatformEndpointResult>
    <EndpointArn>{{ platform_endpoint.arn }}</EndpointArn>
  </CreatePlatformEndpointResult>
  <ResponseMetadata>
    <RequestId>6613341d-3e15-53f7-bf3c-7e56994ba278</RequestId>
  </ResponseMetadata>
</CreatePlatformEndpointResponse>"""

LIST_PLATFORM_APPLICATIONS_TEMPLATE = """<ListPlatformApplicationsResponse xmlns="http://sns.amazonaws.com/doc/2010-03-31/">
  <ListPlatformApplicationsResult>
    <PlatformApplications>
    {% for application in applications %}
      <member>
      <PlatformApplicationArn>{{ application.arn }}</PlatformApplicationArn>
        <Attributes>
        {% for attribute in application.attributes %}
          <entry>
            <key>{{ attribute }}</key>
            <value>{{ application.attributes[attribute] }}</value>
          </entry>
        {% endfor %}
        </Attributes>
      </member>
    {% endfor %}
    </PlatformApplications>
  </ListPlatformApplicationsResult>
  <ResponseMetadata>
    <RequestId>315a335e-85d8-52df-9349-791283cbb529</RequestId>
  </ResponseMetadata>
</ListPlatformApplicationsResponse>"""

DELETE_PLATFORM_APPLICATION_TEMPLATE = """<DeletePlatformApplicationResponse xmlns="http://sns.amazonaws.com/doc/2010-03-31/">
  <ResponseMetadata>
    <RequestId>097dac18-7a77-5823-a8dd-e65476dcb037</RequestId>
  </ResponseMetadata>
</DeletePlatformApplicationResponse>"""

GET_ENDPOINT_ATTRIBUTES_TEMPLATE = """<GetEndpointAttributesResponse xmlns="http://sns.amazonaws.com/doc/2010-03-31/">
  <GetEndpointAttributesResult>
    <Attributes>
    {% for attribute in attributes %}
      <entry>
        <key>{{ attribute }}</key>
        <value>{{ attributes[attribute] }}</value>
      </entry>
    {% endfor %}
    </Attributes>
  </GetEndpointAttributesResult>
  <ResponseMetadata>
    <RequestId>6c725a19-a142-5b77-94f9-1055a9ea04e7</RequestId>
  </ResponseMetadata>
</GetEndpointAttributesResponse>"""

LIST_ENDPOINTS_BY_PLATFORM_APPLICATION_TEMPLATE = """<ListEndpointsByPlatformApplicationResponse xmlns="http://sns.amazonaws.com/doc/2010-03-31/">
  <ListEndpointsByPlatformApplicationResult>
    <Endpoints>
      {% for endpoint in endpoints %}
      <member>
        <EndpointArn>{{ endpoint.arn }}</EndpointArn>
        <Attributes>
          {% for attribute in endpoint.attributes %}
          <entry>
            <key>{{ attribute }}</key>
            <value>{{ endpoint.attributes[attribute] }}</value>
          </entry>
          {% endfor %}
        </Attributes>
      </member>
      {% endfor %}
    </Endpoints>
  </ListEndpointsByPlatformApplicationResult>
  <ResponseMetadata>
    <RequestId>9a48768c-dac8-5a60-aec0-3cc27ea08d96</RequestId>
  </ResponseMetadata>
</ListEndpointsByPlatformApplicationResponse>"""

GET_PLATFORM_APPLICATION_ATTRIBUTES_TEMPLATE = """<GetPlatformApplicationAttributesResponse xmlns="http://sns.amazonaws.com/doc/2010-03-31/">
  <GetPlatformApplicationAttributesResult>
    <Attributes>
    {% for attribute in attributes %}
      <entry>
        <key>{{ attribute }}</key>
        <value>{{ attributes[attribute] }}</value>
      </entry>
    {% endfor %}
    </Attributes>
  </GetPlatformApplicationAttributesResult>
  <ResponseMetadata>
    <RequestId>74848df2-87f6-55ed-890c-c7be80442462</RequestId>
  </ResponseMetadata>
</GetPlatformApplicationAttributesResponse>"""

PUBLISH_TEMPLATE = """<PublishResponse xmlns="http://sns.amazonaws.com/doc/2010-03-31/">
  <PublishResult>
    <MessageId>{{ message_id }}</MessageId>
  </PublishResult>
  <ResponseMetadata>
    <RequestId>f187a3c1-376f-11df-8963-01868b7c937a</RequestId>
  </ResponseMetadata>
</PublishResponse>"""

SET_ENDPOINT_ATTRIBUTES_TEMPLATE = """<SetEndpointAttributesResponse xmlns="http://sns.amazonaws.com/doc/2010-03-31/">
  <ResponseMetadata>
    <RequestId>2fe0bfc7-3e85-5ee5-a9e2-f58b35e85f6a</RequestId>
  </ResponseMetadata>
</SetEndpointAttributesResponse>"""

DELETE_ENDPOINT_TEMPLATE = """<DeleteEndpointResponse xmlns="http://sns.amazonaws.com/doc/2010-03-31/">
 <ResponseMetadata>
 <RequestId>c1d2b191-353c-5a5f-8969-fbdd3900afa8</RequestId>
 </ResponseMetadata>
</DeleteEndpointResponse>"""


SET_PLATFORM_APPLICATION_ATTRIBUTES_TEMPLATE = """<SetPlatformApplicationAttributesResponse xmlns="http://sns.amazonaws.com/doc/2010-03-31/">
  <ResponseMetadata>
    <RequestId>cf577bcc-b3dc-5463-88f1-3180b9412395</RequestId>
  </ResponseMetadata>
</SetPlatformApplicationAttributesResponse>"""

SUBSCRIBE_TEMPLATE = """<SubscribeResponse xmlns="http://sns.amazonaws.com/doc/2010-03-31/">
  <SubscribeResult>
    <SubscriptionArn>{{ subscription.arn }}</SubscriptionArn>
  </SubscribeResult>
  <ResponseMetadata>
    <RequestId>c4407779-24a4-56fa-982c-3d927f93a775</RequestId>
  </ResponseMetadata>
</SubscribeResponse>"""

UNSUBSCRIBE_TEMPLATE = """<UnsubscribeResponse xmlns="http://sns.amazonaws.com/doc/2010-03-31/">
  <ResponseMetadata>
    <RequestId>18e0ac39-3776-11df-84c0-b93cc1666b84</RequestId>
  </ResponseMetadata>
</UnsubscribeResponse>"""

LIST_SUBSCRIPTIONS_TEMPLATE = """<ListSubscriptionsResponse xmlns="http://sns.amazonaws.com/doc/2010-03-31/">
  <ListSubscriptionsResult>
    <Subscriptions>
    {% for subscription in subscriptions %}
      <member>
        <TopicArn>{{ subscription.topic.arn }}</TopicArn>
        <Protocol>{{ subscription.protocol }}</Protocol>
        <SubscriptionArn>{{ subscription.arn }}</SubscriptionArn>
        <Owner>{{ subscription.account_id }}</Owner>
        <Endpoint>{{ subscription.endpoint }}</Endpoint>
      </member>
    {% endfor %}
    </Subscriptions>
    {% if next_token  %}
    <NextToken>{{ next_token }}</NextToken>
    {% endif %}
  </ListSubscriptionsResult>
  <ResponseMetadata>
    <RequestId>384ac68d-3775-11df-8963-01868b7c937a</RequestId>
  </ResponseMetadata>
</ListSubscriptionsResponse>"""

LIST_SUBSCRIPTIONS_BY_TOPIC_TEMPLATE = """<ListSubscriptionsByTopicResponse xmlns="http://sns.amazonaws.com/doc/2010-03-31/">
  <ListSubscriptionsByTopicResult>
    <Subscriptions>
    {% for subscription in subscriptions %}
      <member>
        <TopicArn>{{ subscription.topic.arn }}</TopicArn>
        <Protocol>{{ subscription.protocol }}</Protocol>
        <SubscriptionArn>{{ subscription.arn }}</SubscriptionArn>
        <Owner>{{ subscription.account_id }}</Owner>
        <Endpoint>{{ subscription.endpoint }}</Endpoint>
      </member>
    {% endfor %}
    </Subscriptions>
    {% if next_token  %}
    <NextToken>{{ next_token }}</NextToken>
    {% endif %}
  </ListSubscriptionsByTopicResult>
  <ResponseMetadata>
    <RequestId>384ac68d-3775-11df-8963-01868b7c937a</RequestId>
  </ResponseMetadata>
</ListSubscriptionsByTopicResponse>"""


# Not responding aws system attribetus like 'Owner' and 'SubscriptionArn'
GET_SUBSCRIPTION_ATTRIBUTES_TEMPLATE = """<GetSubscriptionAttributesResponse xmlns="http://sns.amazonaws.com/doc/2010-03-31/">
  <GetSubscriptionAttributesResult>
    <Attributes>
      {% for name, value in attributes.items() %}
      <entry>
        <key>{{ name }}</key>
        <value>{{ value }}</value>
      </entry>
      {% endfor %}
    </Attributes>
  </GetSubscriptionAttributesResult>
  <ResponseMetadata>
    <RequestId>057f074c-33a7-11df-9540-99d0768312d3</RequestId>
  </ResponseMetadata>
</GetSubscriptionAttributesResponse>"""


SET_SUBSCRIPTION_ATTRIBUTES_TEMPLATE = """<SetSubscriptionAttributesResponse xmlns="http://sns.amazonaws.com/doc/2010-03-31/">
  <ResponseMetadata>
    <RequestId>a8763b99-33a7-11df-a9b7-05d48da6f042</RequestId>
  </ResponseMetadata>
</SetSubscriptionAttributesResponse>"""

SET_SMS_ATTRIBUTES_TEMPLATE = """<SetSMSAttributesResponse xmlns="http://sns.amazonaws.com/doc/2010-03-31/">
  <SetSMSAttributesResult/>
  <ResponseMetadata>
    <RequestId>26332069-c04a-5428-b829-72524b56a364</RequestId>
  </ResponseMetadata>
</SetSMSAttributesResponse>"""

GET_SMS_ATTRIBUTES_TEMPLATE = """<GetSMSAttributesResponse xmlns="http://sns.amazonaws.com/doc/2010-03-31/">
  <GetSMSAttributesResult>
    <attributes>
      {% for name, value in attributes.items() %}
      {% if value %}
      <entry>
        <key>{{ name }}</key>
        <value>{{ value }}</value>
      </entry>
      {% endif %}
      {% endfor %}
    </attributes>
  </GetSMSAttributesResult>
  <ResponseMetadata>
    <RequestId>287f9554-8db3-5e66-8abc-c76f0186db7e</RequestId>
  </ResponseMetadata>
</GetSMSAttributesResponse>"""

CHECK_IF_OPTED_OUT_TEMPLATE = """<CheckIfPhoneNumberIsOptedOutResponse xmlns="http://sns.amazonaws.com/doc/2010-03-31/">
  <CheckIfPhoneNumberIsOptedOutResult>
    <isOptedOut>{{ opt_out }}</isOptedOut>
  </CheckIfPhoneNumberIsOptedOutResult>
  <ResponseMetadata>
    <RequestId>287f9554-8db3-5e66-8abc-c76f0186db7e</RequestId>
  </ResponseMetadata>
</CheckIfPhoneNumberIsOptedOutResponse>"""

ERROR_RESPONSE = """<ErrorResponse xmlns="http://sns.amazonaws.com/doc/2010-03-31/">
  <Error>
    <Type>{{ sender }}</Type>
    <Code>{{ code }}</Code>
    <Message>{{ message }}</Message>
  </Error>
  <RequestId>9dd01905-5012-5f99-8663-4b3ecd0dfaef</RequestId>
</ErrorResponse>"""

LIST_OPTOUT_TEMPLATE = """<ListPhoneNumbersOptedOutResponse xmlns="http://sns.amazonaws.com/doc/2010-03-31/">
  <ListPhoneNumbersOptedOutResult>
    <phoneNumbers>
      {% for item in opt_outs %}
      <member>{{ item }}</member>
      {% endfor %}
    </phoneNumbers>
  </ListPhoneNumbersOptedOutResult>
  <ResponseMetadata>
    <RequestId>985e196d-a237-51b6-b33a-4b5601276b38</RequestId>
  </ResponseMetadata>
</ListPhoneNumbersOptedOutResponse>"""

OPT_IN_NUMBER_TEMPLATE = """<OptInPhoneNumberResponse xmlns="http://sns.amazonaws.com/doc/2010-03-31/">
  <OptInPhoneNumberResult/>
  <ResponseMetadata>
    <RequestId>4c61842c-0796-50ef-95ac-d610c0bc8cf8</RequestId>
  </ResponseMetadata>
</OptInPhoneNumberResponse>"""

ADD_PERMISSION_TEMPLATE = """<AddPermissionResponse xmlns="http://sns.amazonaws.com/doc/2010-03-31/">
  <ResponseMetadata>
    <RequestId>c046e713-c5ff-5888-a7bc-b52f0e4f1299</RequestId>
  </ResponseMetadata>
</AddPermissionResponse>"""

DEL_PERMISSION_TEMPLATE = """<RemovePermissionResponse xmlns="http://sns.amazonaws.com/doc/2010-03-31/">
  <ResponseMetadata>
    <RequestId>e767cc9f-314b-5e1b-b283-9ea3fd4e38a3</RequestId>
  </ResponseMetadata>
</RemovePermissionResponse>"""

CONFIRM_SUBSCRIPTION_TEMPLATE = """<ConfirmSubscriptionResponse xmlns="http://sns.amazonaws.com/doc/2010-03-31/">
  <ConfirmSubscriptionResult>
    <SubscriptionArn>{{ sub_arn }}</SubscriptionArn>
  </ConfirmSubscriptionResult>
  <ResponseMetadata>
    <RequestId>16eb4dde-7b3c-5b3e-a22a-1fe2a92d3293</RequestId>
  </ResponseMetadata>
</ConfirmSubscriptionResponse>"""

LIST_TAGS_FOR_RESOURCE_TEMPLATE = """<ListTagsForResourceResponse xmlns="http://sns.amazonaws.com/doc/2010-03-31/">
  <ListTagsForResourceResult>
    <Tags>
      {% for name, value in tags.items() %}
      <member>
        <Key>{{ name }}</Key>
        <Value>{{ value }}</Value>
      </member>
      {% endfor %}
    </Tags>
  </ListTagsForResourceResult>
  <ResponseMetadata>
    <RequestId>97fa763f-861b-5223-a946-20251f2a42e2</RequestId>
  </ResponseMetadata>
</ListTagsForResourceResponse>"""

TAG_RESOURCE_TEMPLATE = """<TagResourceResponse xmlns="http://sns.amazonaws.com/doc/2010-03-31/">
    <TagResourceResult/>
    <ResponseMetadata>
        <RequestId>fd4ab1da-692f-50a7-95ad-e7c665877d98</RequestId>
    </ResponseMetadata>
</TagResourceResponse>"""

UNTAG_RESOURCE_TEMPLATE = """<UntagResourceResponse xmlns="http://sns.amazonaws.com/doc/2010-03-31/">
    <UntagResourceResult/>
    <ResponseMetadata>
        <RequestId>14eb7b1a-4cbd-5a56-80db-2d06412df769</RequestId>
    </ResponseMetadata>
</UntagResourceResponse>"""

PUBLISH_BATCH_TEMPLATE = """<PublishBatchResponse xmlns="http://sns.amazonaws.com/doc/2010-03-31/">
  <ResponseMetadata>
    <RequestId>1549581b-12b7-11e3-895e-1334aEXAMPLE</RequestId>
  </ResponseMetadata>
  <PublishBatchResult>
    <Successful>
{% for successful in successful %}
      <member>
        <Id>{{ successful["Id"] }}</Id>
        <MessageId>{{ successful["MessageId"] }}</MessageId>
      </member>
{% endfor %}
    </Successful>
    <Failed>
{% for failed in failed %}
      <member>
        <Id>{{ failed["Id"] }}</Id>
        <Code>{{ failed["Code"] }}</Code>
        <Message>{{ failed["Message"] }}</Message>
        <SenderFault>{{'true' if failed["SenderFault"] else 'false'}}</SenderFault>
      </member>
{% endfor %}
    </Failed>
  </PublishBatchResult>
</PublishBatchResponse>"""
