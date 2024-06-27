import base64
from typing import Any, Dict, List

from moto.core.responses import BaseResponse
from moto.core.utils import utcnow

from .exceptions import ValidationError
from .models import SESBackend, ses_backends


class EmailResponse(BaseResponse):
    def __init__(self) -> None:
        super().__init__(service_name="ses")

    @property
    def backend(self) -> SESBackend:
        return ses_backends[self.current_account][self.region]

    def verify_email_identity(self) -> str:
        address = self.querystring.get("EmailAddress")[0]  # type: ignore
        self.backend.verify_email_identity(address)
        template = self.response_template(VERIFY_EMAIL_IDENTITY)
        return template.render()

    def verify_email_address(self) -> str:
        address = self.querystring.get("EmailAddress")[0]  # type: ignore
        self.backend.verify_email_address(address)
        template = self.response_template(VERIFY_EMAIL_ADDRESS)
        return template.render()

    def list_identities(self) -> str:
        identity_type = self._get_param("IdentityType")
        if identity_type not in [None, "EmailAddress", "Domain"]:
            raise ValidationError(
                f"Value '{identity_type}' at 'identityType' failed to satisfy constraint: Member must satisfy enum value set: [Domain, EmailAddress]"
            )
        identities = self.backend.list_identities(identity_type)
        template = self.response_template(LIST_IDENTITIES_RESPONSE)
        return template.render(identities=identities)

    def list_verified_email_addresses(self) -> str:
        email_addresses = self.backend.list_verified_email_addresses()
        template = self.response_template(LIST_VERIFIED_EMAIL_RESPONSE)
        return template.render(email_addresses=email_addresses)

    def verify_domain_dkim(self) -> str:
        domain = self.querystring.get("Domain")[0]  # type: ignore
        self.backend.verify_domain(domain)
        template = self.response_template(VERIFY_DOMAIN_DKIM_RESPONSE)
        return template.render()

    def verify_domain_identity(self) -> str:
        domain = self.querystring.get("Domain")[0]  # type: ignore
        self.backend.verify_domain(domain)
        template = self.response_template(VERIFY_DOMAIN_IDENTITY_RESPONSE)
        return template.render()

    def delete_identity(self) -> str:
        domain = self.querystring.get("Identity")[0]  # type: ignore
        self.backend.delete_identity(domain)
        template = self.response_template(DELETE_IDENTITY_RESPONSE)
        return template.render()

    def send_email(self) -> str:
        bodydatakey = "Message.Body.Text.Data"
        if "Message.Body.Html.Data" in self.querystring:
            bodydatakey = "Message.Body.Html.Data"
        body = self.querystring.get(bodydatakey)[0]  # type: ignore
        source = self.querystring.get("Source")[0]  # type: ignore
        subject = self.querystring.get("Message.Subject.Data")[0]  # type: ignore
        destinations: Dict[str, List[str]] = {
            "ToAddresses": [],
            "CcAddresses": [],
            "BccAddresses": [],
        }
        for dest_type in destinations:
            # consume up to 51 to allow exception
            for i in range(1, 52):
                field = f"Destination.{dest_type}.member.{i}"
                address = self.querystring.get(field)
                if address is None:
                    break
                destinations[dest_type].append(address[0])

        message = self.backend.send_email(source, subject, body, destinations)
        template = self.response_template(SEND_EMAIL_RESPONSE)
        return template.render(message=message)

    def send_templated_email(self) -> str:
        source = self.querystring.get("Source")[0]  # type: ignore
        template: List[str] = self.querystring.get("Template")  # type: ignore
        template_data: List[str] = self.querystring.get("TemplateData")  # type: ignore

        destinations: Dict[str, List[str]] = {
            "ToAddresses": [],
            "CcAddresses": [],
            "BccAddresses": [],
        }
        for dest_type in destinations:
            # consume up to 51 to allow exception
            for i in range(1, 52):
                field = f"Destination.{dest_type}.member.{i}"
                address = self.querystring.get(field)
                if address is None:
                    break
                destinations[dest_type].append(address[0])

        message = self.backend.send_templated_email(
            source, template, template_data, destinations
        )
        return self.response_template(SEND_TEMPLATED_EMAIL_RESPONSE).render(
            message=message
        )

    def send_bulk_templated_email(self) -> str:
        source = self.querystring.get("Source")[0]  # type: ignore
        template = self.querystring.get("Template")
        template_data = self.querystring.get("DefaultTemplateData")

        destinations = []
        for i in range(1, 52):
            destination_field = (
                f"Destinations.member.{i}.Destination.ToAddresses.member.1"
            )
            if self.querystring.get(destination_field) is None:
                break
            destination: Dict[str, List[str]] = {
                "ToAddresses": [],
                "CcAddresses": [],
                "BccAddresses": [],
            }
            for dest_type in destination:
                # consume up to 51 to allow exception
                for j in range(1, 52):
                    field = (
                        f"Destinations.member.{i}.Destination.{dest_type}.member.{j}"
                    )
                    address = self.querystring.get(field)
                    if address is None:
                        break
                    destination[dest_type].append(address[0])
            destinations.append({"Destination": destination})

        message = self.backend.send_bulk_templated_email(
            source,
            template,  # type: ignore
            template_data,  # type: ignore
            destinations,
        )
        template = self.response_template(SEND_BULK_TEMPLATED_EMAIL_RESPONSE)
        result = template.render(message=message)
        return result

    def send_raw_email(self) -> str:
        source = self.querystring.get("Source")
        if source is not None:
            (source,) = source

        raw_data = self.querystring.get("RawMessage.Data")[0]  # type: ignore
        raw_data = base64.b64decode(raw_data)
        raw_data = raw_data.decode("utf-8")
        destinations = []
        # consume up to 51 to allow exception
        for i in range(1, 52):
            field = f"Destinations.member.{i}"
            address = self.querystring.get(field)
            if address is None:
                break
            destinations.append(address[0])

        message = self.backend.send_raw_email(source, destinations, raw_data)
        template = self.response_template(SEND_RAW_EMAIL_RESPONSE)
        return template.render(message=message)

    def get_send_quota(self) -> str:
        quota = self.backend.get_send_quota()
        template = self.response_template(GET_SEND_QUOTA_RESPONSE)
        return template.render(quota=quota)

    def get_identity_notification_attributes(self) -> str:
        identities = self._get_params()["Identities"]
        identities = self.backend.get_identity_notification_attributes(identities)
        template = self.response_template(GET_IDENTITY_NOTIFICATION_ATTRIBUTES)
        return template.render(identities=identities)

    def set_identity_feedback_forwarding_enabled(self) -> str:
        identity = self._get_param("Identity")
        enabled = self._get_bool_param("ForwardingEnabled")
        self.backend.set_identity_feedback_forwarding_enabled(identity, enabled)
        template = self.response_template(SET_IDENTITY_FORWARDING_ENABLED_RESPONSE)
        return template.render()

    def set_identity_notification_topic(self) -> str:
        identity = self.querystring.get("Identity")[0]  # type: ignore
        not_type = self.querystring.get("NotificationType")[0]  # type: ignore
        sns_topic = self.querystring.get("SnsTopic")
        if sns_topic:
            sns_topic = sns_topic[0]

        self.backend.set_identity_notification_topic(identity, not_type, sns_topic)
        template = self.response_template(SET_IDENTITY_NOTIFICATION_TOPIC_RESPONSE)
        return template.render()

    def get_send_statistics(self) -> str:
        statistics = self.backend.get_send_statistics()
        template = self.response_template(GET_SEND_STATISTICS)
        return template.render(all_statistics=[statistics])

    def create_configuration_set(self) -> str:
        configuration_set_name = self.querystring.get("ConfigurationSet.Name")[0]  # type: ignore
        self.backend.create_configuration_set(
            configuration_set_name=configuration_set_name
        )
        template = self.response_template(CREATE_CONFIGURATION_SET)
        return template.render()

    def describe_configuration_set(self) -> str:
        configuration_set_name = self.querystring.get("ConfigurationSetName")[0]  # type: ignore
        self.backend.describe_configuration_set(configuration_set_name)
        template = self.response_template(DESCRIBE_CONFIGURATION_SET)
        return template.render(name=configuration_set_name)

    def create_configuration_set_event_destination(self) -> str:
        configuration_set_name = self._get_param("ConfigurationSetName")
        is_configuration_event_enabled = self.querystring.get(
            "EventDestination.Enabled"
        )[0]  # type: ignore
        configuration_event_name = self.querystring.get("EventDestination.Name")[0]  # type: ignore
        event_topic_arn = self.querystring.get(  # type: ignore
            "EventDestination.SNSDestination.TopicARN"
        )[0]
        event_matching_types = self._get_multi_param(
            "EventDestination.MatchingEventTypes.member"
        )

        event_destination = {
            "Name": configuration_event_name,
            "Enabled": is_configuration_event_enabled,
            "EventMatchingTypes": event_matching_types,
            "SNSDestination": event_topic_arn,
        }

        self.backend.create_configuration_set_event_destination(
            configuration_set_name=configuration_set_name,
            event_destination=event_destination,
        )

        template = self.response_template(CREATE_CONFIGURATION_SET_EVENT_DESTINATION)
        return template.render()

    def create_template(self) -> str:
        template_data = self._get_dict_param("Template")
        template_info = {}
        template_info["text_part"] = template_data.get("._text_part", "")
        template_info["html_part"] = template_data.get("._html_part", "")
        template_info["template_name"] = template_data.get("._name", "")
        template_info["subject_part"] = template_data.get("._subject_part", "")
        template_info["Timestamp"] = utcnow()
        self.backend.add_template(template_info=template_info)
        template = self.response_template(CREATE_TEMPLATE)
        return template.render()

    def update_template(self) -> str:
        template_data = self._get_dict_param("Template")
        template_info = {}
        template_info["text_part"] = template_data.get("._text_part", "")
        template_info["html_part"] = template_data.get("._html_part", "")
        template_info["template_name"] = template_data.get("._name", "")
        template_info["subject_part"] = template_data.get("._subject_part", "")
        template_info["Timestamp"] = utcnow()
        self.backend.update_template(template_info=template_info)
        template = self.response_template(UPDATE_TEMPLATE)
        return template.render()

    def get_template(self) -> str:
        template_name = self._get_param("TemplateName")
        template_data = self.backend.get_template(template_name)
        template = self.response_template(GET_TEMPLATE)
        return template.render(template_data=template_data)

    def list_templates(self) -> str:
        email_templates = self.backend.list_templates()
        template = self.response_template(LIST_TEMPLATES)
        return template.render(templates=email_templates)

    def test_render_template(self) -> str:
        render_info = self._get_dict_param("Template")
        rendered_template = self.backend.render_template(render_info)
        template = self.response_template(RENDER_TEMPLATE)
        return template.render(template=rendered_template)

    def delete_template(self) -> str:
        name = self._get_param("TemplateName")
        self.backend.delete_template(name)
        return self.response_template(DELETE_TEMPLATE).render()

    def create_receipt_rule_set(self) -> str:
        rule_set_name = self._get_param("RuleSetName")
        self.backend.create_receipt_rule_set(rule_set_name)
        template = self.response_template(CREATE_RECEIPT_RULE_SET)
        return template.render()

    def create_receipt_rule(self) -> str:
        rule_set_name = self._get_param("RuleSetName")
        rule = self._get_dict_param("Rule.")
        self.backend.create_receipt_rule(rule_set_name, rule)
        template = self.response_template(CREATE_RECEIPT_RULE)
        return template.render()

    def describe_receipt_rule_set(self) -> str:
        rule_set_name = self._get_param("RuleSetName")

        rule_set = self.backend.describe_receipt_rule_set(rule_set_name)

        for i, rule in enumerate(rule_set):
            formatted_rule: Dict[str, Any] = {}

            for k, v in rule.items():
                self._parse_param(k, v, formatted_rule)

            rule_set[i] = formatted_rule

        template = self.response_template(DESCRIBE_RECEIPT_RULE_SET)

        return template.render(rule_set=rule_set, rule_set_name=rule_set_name)

    def describe_receipt_rule(self) -> str:
        rule_set_name = self._get_param("RuleSetName")
        rule_name = self._get_param("RuleName")

        receipt_rule = self.backend.describe_receipt_rule(rule_set_name, rule_name)

        rule: Dict[str, Any] = {}

        for k, v in receipt_rule.items():
            self._parse_param(k, v, rule)

        template = self.response_template(DESCRIBE_RECEIPT_RULE)
        return template.render(rule=rule)

    def update_receipt_rule(self) -> str:
        rule_set_name = self._get_param("RuleSetName")
        rule = self._get_dict_param("Rule.")

        self.backend.update_receipt_rule(rule_set_name, rule)

        template = self.response_template(UPDATE_RECEIPT_RULE)
        return template.render()

    def set_identity_mail_from_domain(self) -> str:
        identity = self._get_param("Identity")
        mail_from_domain = self._get_param("MailFromDomain")
        behavior_on_mx_failure = self._get_param("BehaviorOnMXFailure")

        self.backend.set_identity_mail_from_domain(
            identity, mail_from_domain, behavior_on_mx_failure
        )

        template = self.response_template(SET_IDENTITY_MAIL_FROM_DOMAIN)
        return template.render()

    def get_identity_mail_from_domain_attributes(self) -> str:
        identities = self._get_multi_param("Identities.member.")
        attributes_by_identity = self.backend.get_identity_mail_from_domain_attributes(
            identities
        )
        template = self.response_template(GET_IDENTITY_MAIL_FROM_DOMAIN_ATTRIBUTES)

        return template.render(identities=attributes_by_identity)

    def get_identity_verification_attributes(self) -> str:
        params = self._get_params()
        identities = params.get("Identities")
        verification_attributes = self.backend.get_identity_verification_attributes(
            identities=identities,
        )

        template = self.response_template(GET_IDENTITY_VERIFICATION_ATTRIBUTES_TEMPLATE)
        return template.render(verification_attributes=verification_attributes)


VERIFY_EMAIL_IDENTITY = """<VerifyEmailIdentityResponse xmlns="http://ses.amazonaws.com/doc/2010-12-01/">
  <VerifyEmailIdentityResult/>
  <ResponseMetadata>
    <RequestId>47e0ef1a-9bf2-11e1-9279-0100e8cf109a</RequestId>
  </ResponseMetadata>
</VerifyEmailIdentityResponse>"""

VERIFY_EMAIL_ADDRESS = """<VerifyEmailAddressResponse xmlns="http://ses.amazonaws.com/doc/2010-12-01/">
  <VerifyEmailAddressResult/>
  <ResponseMetadata>
    <RequestId>47e0ef1a-9bf2-11e1-9279-0100e8cf109a</RequestId>
  </ResponseMetadata>
</VerifyEmailAddressResponse>"""

LIST_IDENTITIES_RESPONSE = """<ListIdentitiesResponse xmlns="http://ses.amazonaws.com/doc/2010-12-01/">
  <ListIdentitiesResult>
    <Identities>
        {% for identity in identities %}
          <member>{{ identity }}</member>
        {% endfor %}
    </Identities>
  </ListIdentitiesResult>
  <ResponseMetadata>
    <RequestId>cacecf23-9bf1-11e1-9279-0100e8cf109a</RequestId>
  </ResponseMetadata>
</ListIdentitiesResponse>"""

LIST_VERIFIED_EMAIL_RESPONSE = """<ListVerifiedEmailAddressesResponse xmlns="http://ses.amazonaws.com/doc/2010-12-01/">
  <ListVerifiedEmailAddressesResult>
    <VerifiedEmailAddresses>
        {% for email in email_addresses %}
          <member>{{ email }}</member>
        {% endfor %}
    </VerifiedEmailAddresses>
  </ListVerifiedEmailAddressesResult>
  <ResponseMetadata>
    <RequestId>cacecf23-9bf1-11e1-9279-0100e8cf109a</RequestId>
  </ResponseMetadata>
</ListVerifiedEmailAddressesResponse>"""

VERIFY_DOMAIN_DKIM_RESPONSE = """<VerifyDomainDkimResponse xmlns="http://ses.amazonaws.com/doc/2010-12-01/">
  <VerifyDomainDkimResult>
    <DkimTokens>
      <member>vvjuipp74whm76gqoni7qmwwn4w4qusjiainivf6sf</member>
      <member>3frqe7jn4obpuxjpwpolz6ipb3k5nvt2nhjpik2oy</member>
      <member>wrqplteh7oodxnad7hsl4mixg2uavzneazxv5sxi2</member>
    </DkimTokens>
    </VerifyDomainDkimResult>
    <ResponseMetadata>
      <RequestId>9662c15b-c469-11e1-99d1-797d6ecd6414</RequestId>
    </ResponseMetadata>
</VerifyDomainDkimResponse>"""

VERIFY_DOMAIN_IDENTITY_RESPONSE = """\
<VerifyDomainIdentityResponse xmlns="http://ses.amazonaws.com/doc/2010-12-01/">
  <VerifyDomainIdentityResult>
    <VerificationToken>QTKknzFg2J4ygwa+XvHAxUl1hyHoY0gVfZdfjIedHZ0=</VerificationToken>
  </VerifyDomainIdentityResult>
  <ResponseMetadata>
    <RequestId>94f6368e-9bf2-11e1-8ee7-c98a0037a2b6</RequestId>
  </ResponseMetadata>
</VerifyDomainIdentityResponse>"""

DELETE_IDENTITY_RESPONSE = """<DeleteIdentityResponse xmlns="http://ses.amazonaws.com/doc/2010-12-01/">
  <DeleteIdentityResult/>
  <ResponseMetadata>
    <RequestId>d96bd874-9bf2-11e1-8ee7-c98a0037a2b6</RequestId>
  </ResponseMetadata>
</DeleteIdentityResponse>"""

SEND_EMAIL_RESPONSE = """<SendEmailResponse xmlns="http://ses.amazonaws.com/doc/2010-12-01/">
  <SendEmailResult>
    <MessageId>{{ message.id }}</MessageId>
  </SendEmailResult>
  <ResponseMetadata>
    <RequestId>d5964849-c866-11e0-9beb-01a62d68c57f</RequestId>
  </ResponseMetadata>
</SendEmailResponse>"""

SEND_TEMPLATED_EMAIL_RESPONSE = """<SendTemplatedEmailResponse xmlns="http://ses.amazonaws.com/doc/2010-12-01/">
  <SendTemplatedEmailResult>
    <MessageId>{{ message.id }}</MessageId>
  </SendTemplatedEmailResult>
  <ResponseMetadata>
    <RequestId>d5964849-c866-11e0-9beb-01a62d68c57f</RequestId>
  </ResponseMetadata>
</SendTemplatedEmailResponse>"""

SEND_BULK_TEMPLATED_EMAIL_RESPONSE = """<SendBulkTemplatedEmailResponse xmlns="http://ses.amazonaws.com/doc/2010-12-01/">
  <SendBulkTemplatedEmailResult>
    {% for id in message.ids %}
        <BulkEmailDestinationStatus>
            <MessageId>{{ id }}</MessageId>
        </BulkEmailDestinationStatus>
    {% endfor %}
  </SendBulkTemplatedEmailResult>
  <ResponseMetadata>
    <RequestId>d5964849-c866-11e0-9beb-01a62d68c57f</RequestId>
  </ResponseMetadata>
</SendBulkTemplatedEmailResponse>"""

SEND_RAW_EMAIL_RESPONSE = """<SendRawEmailResponse xmlns="http://ses.amazonaws.com/doc/2010-12-01/">
  <SendRawEmailResult>
    <MessageId>{{ message.id }}</MessageId>
  </SendRawEmailResult>
  <ResponseMetadata>
    <RequestId>e0abcdfa-c866-11e0-b6d0-273d09173b49</RequestId>
  </ResponseMetadata>
</SendRawEmailResponse>"""

GET_SEND_QUOTA_RESPONSE = """<GetSendQuotaResponse xmlns="http://ses.amazonaws.com/doc/2010-12-01/">
  <GetSendQuotaResult>
    <SentLast24Hours>{{ quota.sent_past_24 }}</SentLast24Hours>
    <Max24HourSend>200.0</Max24HourSend>
    <MaxSendRate>1.0</MaxSendRate>
  </GetSendQuotaResult>
  <ResponseMetadata>
    <RequestId>273021c6-c866-11e0-b926-699e21c3af9e</RequestId>
  </ResponseMetadata>
</GetSendQuotaResponse>"""

GET_IDENTITY_NOTIFICATION_ATTRIBUTES = """<GetIdentityNotificationAttributesResponse xmlns="http://ses.amazonaws.com/doc/2010-12-01/">
  <GetIdentityNotificationAttributesResult>
    <NotificationAttributes>
    {% for identity, config in identities.items() %}
      <entry>
        <key>{{ identity }}</key>
        <value>
          <HeadersInBounceNotificationsEnabled>false</HeadersInBounceNotificationsEnabled>
          <HeadersInDeliveryNotificationsEnabled>false</HeadersInDeliveryNotificationsEnabled>
          <HeadersInComplaintNotificationsEnabled>false</HeadersInComplaintNotificationsEnabled>
          {% if config.get("feedback_forwarding_enabled", True) == False %}
          <ForwardingEnabled>false</ForwardingEnabled>
          {% else %}
          <ForwardingEnabled>true</ForwardingEnabled>
          {% endif %}
        </value>
      </entry>
      {% endfor %}
    </NotificationAttributes>
  </GetIdentityNotificationAttributesResult>
  <ResponseMetadata>
    <RequestId>46c90cfc-9055-4b84-96e3-4d6a309a8b9b</RequestId>
  </ResponseMetadata>
</GetIdentityNotificationAttributesResponse>"""

SET_IDENTITY_FORWARDING_ENABLED_RESPONSE = """<SetIdentityFeedbackForwardingEnabledResponse xmlns="http://ses.amazonaws.com/doc/2010-12-01/">
  <SetIdentityFeedbackForwardingEnabledResult/>
  <ResponseMetadata>
    <RequestId>47e0ef1a-9bf2-11e1-9279-0100e8cf109a</RequestId>
  </ResponseMetadata>
</SetIdentityFeedbackForwardingEnabledResponse>"""

SET_IDENTITY_NOTIFICATION_TOPIC_RESPONSE = """<SetIdentityNotificationTopicResponse xmlns="http://ses.amazonaws.com/doc/2010-12-01/">
  <SetIdentityNotificationTopicResult/>
  <ResponseMetadata>
    <RequestId>47e0ef1a-9bf2-11e1-9279-0100e8cf109a</RequestId>
  </ResponseMetadata>
</SetIdentityNotificationTopicResponse>"""

GET_SEND_STATISTICS = """<GetSendStatisticsResponse xmlns="http://ses.amazonaws.com/doc/2010-12-01/">
  <GetSendStatisticsResult>
      <SendDataPoints>
        {% for statistics in all_statistics %}
            <member>
                <DeliveryAttempts>{{ statistics["DeliveryAttempts"] }}</DeliveryAttempts>
                <Rejects>{{ statistics["Rejects"] }}</Rejects>
                <Bounces>{{ statistics["Bounces"] }}</Bounces>
                <Complaints>{{ statistics["Complaints"] }}</Complaints>
                <Timestamp>{{ statistics["Timestamp"].strftime('%Y-%m-%dT%H:%M:%S.%fZ') }}</Timestamp>
            </member>
        {% endfor %}
      </SendDataPoints>
  </GetSendStatisticsResult>
  <ResponseMetadata>
    <RequestId>e0abcdfa-c866-11e0-b6d0-273d09173z49</RequestId>
  </ResponseMetadata>
</GetSendStatisticsResponse>"""

CREATE_CONFIGURATION_SET = """<CreateConfigurationSetResponse xmlns="http://ses.amazonaws.com/doc/2010-12-01/">
  <CreateConfigurationSetResult/>
  <ResponseMetadata>
    <RequestId>47e0ef1a-9bf2-11e1-9279-0100e8cf109a</RequestId>
  </ResponseMetadata>
</CreateConfigurationSetResponse>"""

DESCRIBE_CONFIGURATION_SET = """<DescribeConfigurationSetResponse xmlns="http://ses.amazonaws.com/doc/2010-12-01/">
  <DescribeConfigurationSetResult>
    <ConfigurationSet>
      <Name>{{ name }}</Name>
    </ConfigurationSet>
  </DescribeConfigurationSetResult>
  <ResponseMetadata>
    <RequestId>8e410745-c1bd-4450-82e0-f968cf2105f2</RequestId>
  </ResponseMetadata>
</DescribeConfigurationSetResponse>"""

CREATE_CONFIGURATION_SET_EVENT_DESTINATION = """<CreateConfigurationSetEventDestinationResponse xmlns="http://ses.amazonaws.com/doc/2010-12-01/">
  <CreateConfigurationSetEventDestinationResult/>
  <ResponseMetadata>
    <RequestId>67e0ef1a-9bf2-11e1-9279-0100e8cf109a</RequestId>
  </ResponseMetadata>
</CreateConfigurationSetEventDestinationResponse>"""

CREATE_TEMPLATE = """<CreateTemplateResponse xmlns="http://ses.amazonaws.com/doc/2010-12-01/">
  <CreateTemplateResult/>
  <ResponseMetadata>
    <RequestId>47e0ef1a-9bf2-11e1-9279-0100e8cf12ba</RequestId>
  </ResponseMetadata>
</CreateTemplateResponse>"""

UPDATE_TEMPLATE = """<UpdateTemplateResponse xmlns="http://ses.amazonaws.com/doc/2010-12-01/">
  <UpdateTemplateResult/>
  <ResponseMetadata>
    <RequestId>47e0ef1a-9bf2-11e1-9279-0100e8cf12ba</RequestId>
  </ResponseMetadata>
</UpdateTemplateResponse>"""

GET_TEMPLATE = """<GetTemplateResponse xmlns="http://ses.amazonaws.com/doc/2010-12-01/">
    <GetTemplateResult>
        <Template>
            <TemplateName>{{ template_data["template_name"] }}</TemplateName>
            <SubjectPart>{{ template_data["subject_part"] }}</SubjectPart>
            <HtmlPart><![CDATA[{{ template_data["html_part"] }}]]></HtmlPart>
            <TextPart>{{ template_data["text_part"] }}</TextPart>
        </Template>
    </GetTemplateResult>
    <ResponseMetadata>
        <RequestId>47e0ef1a-9bf2-11e1-9279-0100e8cf12ba</RequestId>
    </ResponseMetadata>
</GetTemplateResponse>"""

LIST_TEMPLATES = """<ListTemplatesResponse xmlns="http://ses.amazonaws.com/doc/2010-12-01/">
    <ListTemplatesResult>
        <TemplatesMetadata>
            {% for template in templates %}
                <Item>
                    <Name>{{ template["template_name"] }}</Name>
                    <CreatedTimestamp>{{ template["Timestamp"] }}</CreatedTimestamp>
                </Item>
            {% endfor %}
        </TemplatesMetadata>
    </ListTemplatesResult>
    <ResponseMetadata>
        <RequestId>47e0ef1a-9bf2-11e1-9279-0100e8cf12ba</RequestId>
    </ResponseMetadata>
</ListTemplatesResponse>"""

RENDER_TEMPLATE = """
<TestRenderTemplateResponse xmlns="http://ses.amazonaws.com/doc/2010-12-01/">
    <TestRenderTemplateResult>
      <RenderedTemplate>
      {{template | e}}
      </RenderedTemplate>
    </TestRenderTemplateResult>
    <ResponseMetadata>
        <RequestId>47e0ef1a-9bf2-11e1-9279-0100e8cf12ba</RequestId>
    </ResponseMetadata>
</TestRenderTemplateResponse>
"""

DELETE_TEMPLATE = """<DeleteTemplateResponse xmlns="http://ses.amazonaws.com/doc/2010-12-01/">
    <DeleteTemplateResult>
    </DeleteTemplateResult>
    <ResponseMetadata>
        <RequestId>47e0ef1a-9bf2-11e1-9279-0100e8cf12ba</RequestId>
    </ResponseMetadata>
</DeleteTemplateResponse>"""

CREATE_RECEIPT_RULE_SET = """<CreateReceiptRuleSetResponse xmlns="http://ses.amazonaws.com/doc/2010-12-01/">
  <CreateReceiptRuleSetResult/>
  <ResponseMetadata>
    <RequestId>47e0ef1a-9bf2-11e1-9279-01ab88cf109a</RequestId>
  </ResponseMetadata>
</CreateReceiptRuleSetResponse>"""

CREATE_RECEIPT_RULE = """<CreateReceiptRuleResponse xmlns="http://ses.amazonaws.com/doc/2010-12-01/">
  <CreateReceiptRuleResult/>
  <ResponseMetadata>
    <RequestId>15e0ef1a-9bf2-11e1-9279-01ab88cf109a</RequestId>
  </ResponseMetadata>
</CreateReceiptRuleResponse>"""

DESCRIBE_RECEIPT_RULE_SET = """<DescribeReceiptRuleSetResponse xmlns="http://ses.amazonaws.com/doc/2010-12-01/">
  <DescribeReceiptRuleSetResult>
    <Rules>
      {% for rule in rule_set %}
      <member>
        <Recipients>
          {% for recipient in rule["recipients"] %}
          <member>{{recipient}}</member>
          {% endfor %}
        </Recipients>
        <Name>{{rule["name"]}}</Name>
        <Actions>
          {% for action in rule["actions"] %}
          <member>
            {% if action["_s3_action"] %}
            <S3Action>
              <BucketName>{{action["_s3_action"]["_bucket_name"]}}</BucketName>
              <KmsKeyArn>{{action["_s3_action"]["_kms_key_arn"]}}</KmsKeyArn>
              <ObjectKeyPrefix>{{action["_s3_action"]["_object_key_prefix"]}}</ObjectKeyPrefix>
              <TopicArn>{{action["_s3_action"]["_topic_arn"]}}</TopicArn>
            </S3Action>
            {% endif %}
            {% if action["_bounce_action"] %}
            <BounceAction>
              <TopicArn>{{action["_bounce_action"]["_topic_arn"]}}</TopicArn>
              <SmtpReplyCode>{{action["_bounce_action"]["_smtp_reply_code"]}}</SmtpReplyCode>
              <StatusCode>{{action["_bounce_action"]["_status_code"]}}</StatusCode>
              <Message>{{action["_bounce_action"]["_message"]}}</Message>
              <Sender>{{action["_bounce_action"]["_sender"]}}</Sender>
            </BounceAction>
            {% endif %}
          </member>
          {% endfor %}
        </Actions>
        <TlsPolicy>{{rule["tls_policy"]}}</TlsPolicy>
        <ScanEnabled>{{rule["scan_enabled"]}}</ScanEnabled>
        <Enabled>{{rule["enabled"]}}</Enabled>
      </member>
      {% endfor %}
    </Rules>
    <Metadata>
      <Name>{{rule_set_name}}</Name>
      <CreatedTimestamp>2021-10-31</CreatedTimestamp>
    </Metadata>
  </DescribeReceiptRuleSetResult>
  <ResponseMetadata>
    <RequestId>15e0ef1a-9bf2-11e1-9279-01ab88cf109a</RequestId>
  </ResponseMetadata>
</DescribeReceiptRuleSetResponse>
"""

DESCRIBE_RECEIPT_RULE = """<DescribeReceiptRuleResponse xmlns="http://ses.amazonaws.com/doc/2010-12-01/">
  <DescribeReceiptRuleResult>
    <Rule>
      <Recipients>
        {% for recipient in rule["recipients"] %}
        <member>{{recipient}}</member>
        {% endfor %}
      </Recipients>
      <Name>{{rule["name"]}}</Name>
      <Actions>
        {% for action in rule["actions"] %}
        <member>
          {% if action["_s3_action"] %}
          <S3Action>
            <BucketName>{{action["_s3_action"]["_bucket_name"]}}</BucketName>
            <KmsKeyArn>{{action["_s3_action"]["_kms_key_arn"]}}</KmsKeyArn>
            <ObjectKeyPrefix>{{action["_s3_action"]["_object_key_prefix"]}}</ObjectKeyPrefix>
            <TopicArn>{{action["_s3_action"]["_topic_arn"]}}</TopicArn>
          </S3Action>
          {% endif %}
          {% if action["_bounce_action"] %}
          <BounceAction>
            <TopicArn>{{action["_bounce_action"]["_topic_arn"]}}</TopicArn>
            <SmtpReplyCode>{{action["_bounce_action"]["_smtp_reply_code"]}}</SmtpReplyCode>
            <StatusCode>{{action["_bounce_action"]["_status_code"]}}</StatusCode>
            <Message>{{action["_bounce_action"]["_message"]}}</Message>
            <Sender>{{action["_bounce_action"]["_sender"]}}</Sender>
          </BounceAction>
          {% endif %}
        </member>
        {% endfor %}
      </Actions>
      <TlsPolicy>{{rule["tls_policy"]}}</TlsPolicy>
      <ScanEnabled>{{rule["scan_enabled"]}}</ScanEnabled>
      <Enabled>{{rule["enabled"]}}</Enabled>
    </Rule>
  </DescribeReceiptRuleResult>
  <ResponseMetadata>
    <RequestId>15e0ef1a-9bf2-11e1-9279-01ab88cf109a</RequestId>
  </ResponseMetadata>
</DescribeReceiptRuleResponse>
"""

UPDATE_RECEIPT_RULE = """<UpdateReceiptRuleResponse xmlns="http://ses.amazonaws.com/doc/2010-12-01/">
  <UpdateReceiptRuleResult/>
  <ResponseMetadata>
    <RequestId>15e0ef1a-9bf2-11e1-9279-01ab88cf109a</RequestId>
  </ResponseMetadata>
</UpdateReceiptRuleResponse>"""

SET_IDENTITY_MAIL_FROM_DOMAIN = """<SetIdentityMailFromDomainResponse xmlns="http://ses.amazonaws.com/doc/2010-12-01/">
  <SetIdentityMailFromDomainResult/>
  <ResponseMetadata>
    <RequestId>47e0ef1a-9bf2-11e1-9279-0100e8cf109a</RequestId>
  </ResponseMetadata>
</SetIdentityMailFromDomainResponse>"""

GET_IDENTITY_MAIL_FROM_DOMAIN_ATTRIBUTES = """<GetIdentityMailFromDomainAttributesResponse xmlns="http://ses.amazonaws.com/doc/2010-12-01/">
  <GetIdentityMailFromDomainAttributesResult>
    {% if identities.items()|length > 0 %}
    <MailFromDomainAttributes>
    {% for name, value in identities.items() %}
      <entry>
        <key>{{ name }}</key>
        <value>
          {% if 'mail_from_domain' in value %}
          <MailFromDomain>{{ value.get("mail_from_domain") }}</MailFromDomain>
          <MailFromDomainStatus>Success</MailFromDomainStatus>
          {% endif %}
          <BehaviorOnMXFailure>{{ value.get("behavior_on_mx_failure") }}</BehaviorOnMXFailure>
        </value>
      </entry>
    {% endfor %}
    </MailFromDomainAttributes>
    {% else %}
    <MailFromDomainAttributes/>
    {% endif %}
  </GetIdentityMailFromDomainAttributesResult>
  <ResponseMetadata>
    <RequestId>47e0ef1a-9bf2-11e1-9279-0100e8cf109a</RequestId>
  </ResponseMetadata>
</GetIdentityMailFromDomainAttributesResponse>"""

GET_IDENTITY_VERIFICATION_ATTRIBUTES_TEMPLATE = """<GetIdentityVerificationAttributesResponse xmlns="http://ses.amazonaws.com/doc/2010-12-01/">
  <GetIdentityVerificationAttributesResult>
    <VerificationAttributes>
      {% for name, value in verification_attributes.items() %}
      <entry>
        <key>{{ name }}</key>
        <value>
          <VerificationStatus>{{ value }}</VerificationStatus>
          <VerificationToken>ILQMESfEW0p6i6gIJcEWvO65TP5hg6B99hGFZ2lxrIs=</VerificationToken>
        </value>
      </entry>
      {% endfor %}
    </VerificationAttributes>
  </GetIdentityVerificationAttributesResult>
  <ResponseMetadata>
    <RequestId>d435c1b8-a225-4b89-acff-81fcf7ef9236</RequestId>
  </ResponseMetadata>
</GetIdentityVerificationAttributesResponse>"""
