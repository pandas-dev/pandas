from typing import Any, Dict
from urllib.parse import unquote

import xmltodict

from moto.core.responses import TYPE_RESPONSE, BaseResponse

from .models import CloudFrontBackend, cloudfront_backends

XMLNS = "http://cloudfront.amazonaws.com/doc/2020-05-31/"


class CloudFrontResponse(BaseResponse):
    def __init__(self) -> None:
        super().__init__(service_name="cloudfront")

    def _get_xml_body(self) -> Dict[str, Any]:
        return xmltodict.parse(self.body, dict_constructor=dict, force_list="Path")

    @property
    def backend(self) -> CloudFrontBackend:
        return cloudfront_backends[self.current_account][self.partition]

    def create_distribution(self) -> TYPE_RESPONSE:
        params = self._get_xml_body()
        if "DistributionConfigWithTags" in params:
            config = params.get("DistributionConfigWithTags")
            tags = (config.get("Tags", {}).get("Items") or {}).get("Tag", [])  # type: ignore[union-attr]
            if not isinstance(tags, list):
                tags = [tags]
        else:
            config = params
            tags = []
        distribution_config = config.get("DistributionConfig")  # type: ignore[union-attr]
        distribution, location, e_tag = self.backend.create_distribution(
            distribution_config=distribution_config,
            tags=tags,
        )
        template = self.response_template(CREATE_DISTRIBUTION_TEMPLATE)
        response = template.render(distribution=distribution, xmlns=XMLNS)
        headers = {"ETag": e_tag, "Location": location}
        return 200, headers, response

    def list_distributions(self) -> TYPE_RESPONSE:
        distributions = self.backend.list_distributions()
        template = self.response_template(LIST_TEMPLATE)
        response = template.render(distributions=distributions)
        return 200, {}, response

    def delete_distribution(self) -> TYPE_RESPONSE:
        distribution_id = self.path.split("/")[-1]
        if_match = self._get_param("If-Match")
        self.backend.delete_distribution(distribution_id, if_match)
        return 204, {"status": 204}, ""

    def get_distribution(self) -> TYPE_RESPONSE:
        distribution_id = self.path.split("/")[-1]
        dist, etag = self.backend.get_distribution(distribution_id)
        template = self.response_template(GET_DISTRIBUTION_TEMPLATE)
        response = template.render(distribution=dist, xmlns=XMLNS)
        return 200, {"ETag": etag}, response

    def get_distribution_config(self) -> TYPE_RESPONSE:
        dist_id = self.path.split("/")[-2]
        distribution_config, etag = self.backend.get_distribution_config(dist_id)
        template = self.response_template(GET_DISTRIBUTION_CONFIG_TEMPLATE)
        response = template.render(distribution=distribution_config, xmlns=XMLNS)
        return 200, {"ETag": etag}, response

    def update_distribution(self) -> TYPE_RESPONSE:
        dist_id = self.path.split("/")[-2]
        params = self._get_xml_body()
        dist_config = params.get("DistributionConfig")
        if_match = self.headers["If-Match"]

        dist, location, e_tag = self.backend.update_distribution(
            dist_config=dist_config,  # type: ignore[arg-type]
            _id=dist_id,
            if_match=if_match,
        )
        template = self.response_template(UPDATE_DISTRIBUTION_TEMPLATE)
        response = template.render(distribution=dist, xmlns=XMLNS)
        headers = {"ETag": e_tag, "Location": location}
        return 200, headers, response

    def create_invalidation(self) -> TYPE_RESPONSE:
        dist_id = self.path.split("/")[-2]
        params = self._get_xml_body()["InvalidationBatch"]
        paths = ((params.get("Paths") or {}).get("Items") or {}).get("Path") or []
        caller_ref = params.get("CallerReference")

        invalidation = self.backend.create_invalidation(dist_id, paths, caller_ref)  # type: ignore[arg-type]
        template = self.response_template(CREATE_INVALIDATION_TEMPLATE)
        response = template.render(invalidation=invalidation, xmlns=XMLNS)

        return 200, {"Location": invalidation.location}, response

    def list_invalidations(self) -> TYPE_RESPONSE:
        dist_id = self.path.split("/")[-2]
        invalidations = self.backend.list_invalidations(dist_id)
        template = self.response_template(INVALIDATIONS_TEMPLATE)
        response = template.render(invalidations=invalidations, xmlns=XMLNS)

        return 200, {}, response

    def list_tags_for_resource(self) -> TYPE_RESPONSE:
        resource = unquote(self._get_param("Resource"))
        tags = self.backend.list_tags_for_resource(resource=resource)["Tags"]
        template = self.response_template(TAGS_TEMPLATE)
        response = template.render(tags=tags, xmlns=XMLNS)
        return 200, {}, response

    def create_origin_access_control(self) -> TYPE_RESPONSE:
        config = self._get_xml_body().get("OriginAccessControlConfig", {})
        config.pop("@xmlns", None)
        control = self.backend.create_origin_access_control(config)
        template = self.response_template(ORIGIN_ACCESS_CONTROl)
        return 200, {}, template.render(control=control)

    def get_origin_access_control(self) -> TYPE_RESPONSE:
        control_id = self.path.split("/")[-1]
        control = self.backend.get_origin_access_control(control_id)
        template = self.response_template(ORIGIN_ACCESS_CONTROl)
        return 200, {"ETag": control.etag}, template.render(control=control)

    def list_origin_access_controls(self) -> TYPE_RESPONSE:
        controls = self.backend.list_origin_access_controls()
        template = self.response_template(LIST_ORIGIN_ACCESS_CONTROl)
        return 200, {}, template.render(controls=controls)

    def update_origin_access_control(self) -> TYPE_RESPONSE:
        control_id = self.path.split("/")[-2]
        config = self._get_xml_body().get("OriginAccessControlConfig", {})
        config.pop("@xmlns", None)
        control = self.backend.update_origin_access_control(control_id, config)
        template = self.response_template(ORIGIN_ACCESS_CONTROl)
        return 200, {"ETag": control.etag}, template.render(control=control)

    def delete_origin_access_control(self) -> TYPE_RESPONSE:
        control_id = self.path.split("/")[-1]
        self.backend.delete_origin_access_control(control_id)
        return 200, {}, "{}"


DIST_META_TEMPLATE = """
    <Id>{{ distribution.distribution_id }}</Id>
    <ARN>{{ distribution.arn }}</ARN>
    <Status>{{ distribution.status }}</Status>
    <LastModifiedTime>{{ distribution.last_modified_time }}</LastModifiedTime>
    <InProgressInvalidationBatches>{{ distribution.in_progress_invalidation_batches }}</InProgressInvalidationBatches>
    <DomainName>{{ distribution.domain_name }}</DomainName>
"""


DIST_CONFIG_TEMPLATE = """
      <CallerReference>{{ distribution.distribution_config.caller_reference }}</CallerReference>
      <Aliases>
        <Quantity>{{ distribution.distribution_config.aliases|length }}</Quantity>
        <Items>
          {% for alias  in distribution.distribution_config.aliases %}
            <CNAME>{{ alias }}</CNAME>
          {% endfor %}
        </Items>
      </Aliases>
      <DefaultRootObject>{{ distribution.distribution_config.default_root_object }}</DefaultRootObject>
      <Origins>
        <Quantity>{{ distribution.distribution_config.origins|length }}</Quantity>
        <Items>
          {% for origin  in distribution.distribution_config.origins %}
          <Origin>
            <Id>{{ origin.id }}</Id>
            <DomainName>{{ origin.domain_name }}</DomainName>
            <OriginPath>{{ origin.origin_path }}</OriginPath>
            <CustomHeaders>
              <Quantity>{{ origin.custom_headers|length }}</Quantity>
              <Items>
                {% for header  in origin.custom_headers %}
                  <OriginCustomHeader>
                  <HeaderName>{{ header['HeaderName'] }}</HeaderName>
                  <HeaderValue>{{ header['HeaderValue'] }}</HeaderValue>
                  </OriginCustomHeader>
                {% endfor %}
              </Items>
            </CustomHeaders>
            {% if origin.s3_access_identity %}
            <S3OriginConfig>
              <OriginAccessIdentity>{{ origin.s3_access_identity }}</OriginAccessIdentity>
            </S3OriginConfig>
            {% endif %}
            {% if origin.custom_origin %}
            <CustomOriginConfig>
              <HTTPPort>{{ origin.custom_origin.http_port }}</HTTPPort>
              <HTTPSPort>{{ origin.custom_origin.https_port }}</HTTPSPort>
              <OriginProtocolPolicy>{{ origin.custom_origin.protocol_policy }}</OriginProtocolPolicy>
              <OriginSslProtocols>
                <Quantity>{{ origin.custom_origin.ssl_protocols|length }}</Quantity>
                <Items>
                  {% for protocol  in origin.custom_origin.ssl_protocols %}
                  <SslProtocol>{{ protocol }}</SslProtocol>
                  {% endfor %}
                </Items>
              </OriginSslProtocols>
              <OriginReadTimeout>{{ origin.custom_origin.read_timeout }}</OriginReadTimeout>
              <OriginKeepaliveTimeout>{{ origin.custom_origin.keep_alive }}</OriginKeepaliveTimeout>
            </CustomOriginConfig>
            {% endif %}
            <ConnectionAttempts>{{ origin.connection_attempts }}</ConnectionAttempts>
            <ConnectionTimeout>{{ origin.connection_timeout }}</ConnectionTimeout>
            {% if origin.origin_shield %}
            <OriginShield>
              <Enabled>{{ origin.origin_shield.get("Enabled") }}</Enabled>
              <OriginShieldRegion>{{ origin.origin_shield.get("OriginShieldRegion") }}</OriginShieldRegion>
            </OriginShield>
            {% else %}
            <OriginShield>
              <Enabled>false</Enabled>
            </OriginShield>
            {% endif %}
            </Origin>
          {% endfor %}
        </Items>
      </Origins>
      <OriginGroups>
        <Quantity>{{ distribution.distribution_config.origin_groups|length }}</Quantity>
        {% if distribution.distribution_config.origin_groups %}
        <Items>
          {% for origin_group  in distribution.distribution_config.origin_groups %}
            <Id>{{ origin_group.id }}</Id>
            <FailoverCriteria>
              <StatusCodes>
                <Quantity>{{ origin_group.failover_criteria.status_codes.quantity }}</Quantity>
                <Items>
                  {% for status_code_list  in origin_group_list.failover_criteria.status_codes.StatusCodeList %}
                    <StatusCode>{{ status_code_list.status_code }}</StatusCode>
                  {% endfor %}
                </Items>
              </StatusCodes>
            </FailoverCriteria>
            <Members>
              <Quantity>{{ origin_group.members.quantity }}</Quantity>
              <Items>
                {% for origin_group_member_list  in origin_group.members.OriginGroupMemberList %}
                  <OriginId>{{ origin_group_member_list.origin_id }}</OriginId>
                {% endfor %}
              </Items>
            </Members>
          {% endfor %}
        </Items>
        {% endif %}
      </OriginGroups>
      <DefaultCacheBehavior>
        <TargetOriginId>{{ distribution.distribution_config.default_cache_behavior.target_origin_id }}</TargetOriginId>
        <TrustedSigners>
          <Enabled>{{ distribution.distribution_config.default_cache_behavior.trusted_signers_enabled }}</Enabled>
          <Quantity>{{ distribution.distribution_config.default_cache_behavior.trusted_signers|length }}</Quantity>
          <Items>
            {% for aws_account_number  in distribution.distribution_config.default_cache_behavior.trusted_signers %}
              <AwsAccountNumber>{{ aws_account_number }}</AwsAccountNumber>
            {% endfor %}
          </Items>
        </TrustedSigners>
        <TrustedKeyGroups>
          <Enabled>{{ distribution.distribution_config.default_cache_behavior.trusted_key_groups_enabled }}</Enabled>
          <Quantity>{{ distribution.distribution_config.default_cache_behavior.trusted_key_groups|length }}</Quantity>
          <Items>
            {% for key_group  in distribution.distribution_config.default_cache_behavior.trusted_key_groups %}
              <KeyGroup>{{ key_group }}</KeyGroup>
            {% endfor %}
          </Items>
        </TrustedKeyGroups>
        <ViewerProtocolPolicy>{{ distribution.distribution_config.default_cache_behavior.viewer_protocol_policy }}</ViewerProtocolPolicy>
        <AllowedMethods>
          <Quantity>{{ distribution.distribution_config.default_cache_behavior.allowed_methods|length }}</Quantity>
          <Items>
            {% for method in distribution.distribution_config.default_cache_behavior.allowed_methods %}
            <Method>{{ method }}</Method>
            {% endfor %}
          </Items>
          <CachedMethods>
            <Quantity>{{ distribution.distribution_config.default_cache_behavior.cached_methods|length }}</Quantity>
            <Items>
              {% for method in distribution.distribution_config.default_cache_behavior.cached_methods %}
              <Method>{{ method }}</Method>
              {% endfor %}
            </Items>
          </CachedMethods>
        </AllowedMethods>
        <SmoothStreaming>{{ distribution.distribution_config.default_cache_behavior.smooth_streaming }}</SmoothStreaming>
        <Compress>{{ 'true' if distribution.distribution_config.default_cache_behavior.compress else 'false' }}</Compress>
        <LambdaFunctionAssociations>
          <Quantity>{{ distribution.distribution_config.default_cache_behavior.lambda_function_associations|length }}</Quantity>
          {% if distribution.distribution_config.default_cache_behavior.lambda_function_associations %}
          <Items>
            {% for func in distribution.distribution_config.default_cache_behavior.lambda_function_associations %}
              <LambdaFunctionARN>{{ func.arn }}</LambdaFunctionARN>
              <EventType>{{ func.event_type }}</EventType>
              <IncludeBody>{{ func.include_body }}</IncludeBody>
            {% endfor %}
          </Items>
          {% endif %}
        </LambdaFunctionAssociations>
        <FunctionAssociations>
          <Quantity>{{ distribution.distribution_config.default_cache_behavior.function_associations|length }}</Quantity>
          {% if distribution.distribution_config.default_cache_behavior.function_associations %}
          <Items>
            {% for func in distribution.distribution_config.default_cache_behavior.function_associations %}
              <FunctionARN>{{ func.arn }}</FunctionARN>
              <EventType>{{ func.event_type }}</EventType>
            {% endfor %}
          </Items>
          {% endif %}
        </FunctionAssociations>
        <FieldLevelEncryptionId>{{ distribution.distribution_config.default_cache_behavior.field_level_encryption_id }}</FieldLevelEncryptionId>
        <RealtimeLogConfigArn>{{ distribution.distribution_config.default_cache_behavior.realtime_log_config_arn }}</RealtimeLogConfigArn>
        <CachePolicyId>{{ distribution.distribution_config.default_cache_behavior.cache_policy_id }}</CachePolicyId>
        <OriginRequestPolicyId>{{ distribution.distribution_config.default_cache_behavior.origin_request_policy_id }}</OriginRequestPolicyId>
        <ResponseHeadersPolicyId>{{ distribution.distribution_config.default_cache_behavior.response_headers_policy_id }}</ResponseHeadersPolicyId>
        <ForwardedValues>
          <QueryString>{{ distribution.distribution_config.default_cache_behavior.forwarded_values.query_string }}</QueryString>
          <Cookies>
            <Forward>{{ distribution.distribution_config.default_cache_behavior.forwarded_values.cookie_forward }}</Forward>
            <WhitelistedNames>
              <Quantity>{{ distribution.distribution_config.default_cache_behavior.forwarded_values.whitelisted_names|length }}</Quantity>
              <Items>
                {% for name  in distribution.distribution_config.default_cache_behavior.forwarded_values.whitelisted_names %}
                  <Name>{{ name }}</Name>
                {% endfor %}
              </Items>
            </WhitelistedNames>
          </Cookies>
          <Headers>
            <Quantity>{{ distribution.distribution_config.default_cache_behavior.forwarded_values.headers|length }}</Quantity>
            <Items>
              {% for h  in distribution.distribution_config.default_cache_behavior.forwarded_values.headers %}
                <Name>{{ h }}</Name>
              {% endfor %}
            </Items>
          </Headers>
          <QueryStringCacheKeys>
            <Quantity>{{ distribution.distribution_config.default_cache_behavior.forwarded_values.query_string_cache_keys|length }}</Quantity>
            <Items>
              {% for key  in distribution.distribution_config.default_cache_behavior.forwarded_values.query_string_cache_keys %}
                <Name>{{ key }}</Name>
              {% endfor %}
            </Items>
          </QueryStringCacheKeys>
        </ForwardedValues>
        <MinTTL>{{ distribution.distribution_config.default_cache_behavior.min_ttl }}</MinTTL>
        <DefaultTTL>{{ distribution.distribution_config.default_cache_behavior.default_ttl }}</DefaultTTL>
        <MaxTTL>{{ distribution.distribution_config.default_cache_behavior.max_ttl }}</MaxTTL>
      </DefaultCacheBehavior>
      <CacheBehaviors>
        <Quantity>{{ distribution.distribution_config.cache_behaviors|length }}</Quantity>
        {% if distribution.distribution_config.cache_behaviors %}
        <Items>
          {% for behaviour in distribution.distribution_config.cache_behaviors %}
            <CacheBehavior>
                <PathPattern>{{ behaviour.path_pattern }}</PathPattern>
                <TargetOriginId>{{ behaviour.target_origin_id }}</TargetOriginId>
                <TrustedSigners>
                  <Enabled>{{ behaviour.trusted_signers.enabled }}</Enabled>
                  <Quantity>{{ behaviour.trusted_signers | length }}</Quantity>
                  <Items>
                    {% for account_nr  in behaviour.trusted_signers %}
                      <AwsAccountNumber>{{ account_nr }}</AwsAccountNumber>
                    {% endfor %}
                  </Items>
                </TrustedSigners>
                <TrustedKeyGroups>
                  <Enabled>{{ behaviour.trusted_key_groups.enabled }}</Enabled>
                  <Quantity>{{ behaviour.trusted_key_groups | length }}</Quantity>
                  <Items>
                    {% for trusted_key_group_id_list  in behaviour.trusted_key_groups.TrustedKeyGroupIdList %}
                      <KeyGroup>{{ trusted_key_group_id_list.key_group }}</KeyGroup>
                    {% endfor %}
                  </Items>
                </TrustedKeyGroups>
                <ViewerProtocolPolicy>{{ behaviour.viewer_protocol_policy }}</ViewerProtocolPolicy>
                <AllowedMethods>
                  <Quantity>{{ behaviour.allowed_methods | length }}</Quantity>
                  <Items>
                    {% for method in behaviour.allowed_methods %}<Method>{{ method }}</Method>{% endfor %}
                  </Items>
                  <CachedMethods>
                    <Quantity>{{ behaviour.cached_methods|length }}</Quantity>
                    <Items>
                      {% for method in behaviour.cached_methods %}<Method>{{ method }}</Method>{% endfor %}
                    </Items>
                  </CachedMethods>
                </AllowedMethods>
                <SmoothStreaming>{{ behaviour.smooth_streaming }}</SmoothStreaming>
                <Compress>{{ behaviour.compress }}</Compress>
                <LambdaFunctionAssociations>
                  <Quantity>{{ behaviour.lambda_function_associations | length }}</Quantity>
                  <Items>
                    {% for lambda_function_association_list in behaviour.lambda_function_associations.LambdaFunctionAssociationList %}
                      <LambdaFunctionARN>{{ LambdaFunctionARN }}</LambdaFunctionARN>
                      <EventType>{{ EventType }}</EventType>
                      <IncludeBody>{{ lambda_function_association_list.include_body }}</IncludeBody>
                    {% endfor %}
                  </Items>
                </LambdaFunctionAssociations>
                <FunctionAssociations>
                  <Quantity>{{ behaviour.function_associations | length }}</Quantity>
                  <Items>
                    {% for function_association_list  in behaviour.function_associations.FunctionAssociationList %}
                      <FunctionARN>{{ FunctionARN }}</FunctionARN>
                      <EventType>{{ EventType }}</EventType>
                    {% endfor %}
                  </Items>
                </FunctionAssociations>
                <FieldLevelEncryptionId>{{ behaviour.field_level_encryption_id }}</FieldLevelEncryptionId>
                <RealtimeLogConfigArn>{{ behaviour.realtime_log_config_arn }}</RealtimeLogConfigArn>
                <CachePolicyId>{{ behaviour.cache_policy_id }}</CachePolicyId>
                <OriginRequestPolicyId>{{ behaviour.origin_request_policy_id }}</OriginRequestPolicyId>
                <ResponseHeadersPolicyId>{{ behaviour.response_headers_policy_id }}</ResponseHeadersPolicyId>
                <ForwardedValues>
                  <QueryString>{{ behaviour.forwarded_values.query_string }}</QueryString>
                  <Cookies>
                    <Forward>{{ ItemSelection }}</Forward>
                    <WhitelistedNames>
                      <Quantity>{{ behaviour.forwarded_values.cookies.whitelisted_names| length }}</Quantity>
                      <Items>
                        {% for wl_name in behaviour.forwarded_values.cookies.whitelisted_names %}
                          <Name>{{ wl_name }}</Name>
                        {% endfor %}
                      </Items>
                    </WhitelistedNames>
                  </Cookies>
                  <Headers>
                    <Quantity>{{ behaviour.forwarded_values.headers | length }}</Quantity>
                    <Items>
                      {% for header_list in behaviour.forwarded_values.headers.HeaderList %}
                        <Name>{{ header_list.name }}</Name>
                      {% endfor %}
                    </Items>
                  </Headers>
                  <QueryStringCacheKeys>
                    <Quantity>{{ behaviour.forwarded_values.query_string_cache_keys | length }}</Quantity>
                    <Items>
                      {% for query_string_cache_keys_list in behaviour.forwarded_values.query_string_cache_keys.QueryStringCacheKeysList %}
                        <Name>{{ query_string_cache_keys_list.name }}</Name>
                      {% endfor %}
                    </Items>
                  </QueryStringCacheKeys>
                </ForwardedValues>
                <MinTTL>{{ behaviour.min_ttl }}</MinTTL>
                <DefaultTTL>{{ behaviour.default_ttl }}</DefaultTTL>
                <MaxTTL>{{ behaviour.max_ttl }}</MaxTTL>
            </CacheBehavior>
          {% endfor %}
        </Items>
        {% endif %}
      </CacheBehaviors>
      <CustomErrorResponses>
        <Quantity>{{ distribution.distribution_config.custom_error_responses|length }}</Quantity>
        {% if distribution.distribution_config.custom_error_responses %}
        <Items>
          {% for response  in distribution.distribution_config.custom_error_responses %}
            <ErrorCode>{{ response.error_code }}</ErrorCode>
            <ResponsePagePath>{{ response.response_page_path }}</ResponsePagePath>
            <ResponseCode>{{ response.response_code }}</ResponseCode>
            <ErrorCachingMinTTL>{{ response.error_caching_min_ttl }}</ErrorCachingMinTTL>
          {% endfor %}
        </Items>
        {% endif %}
      </CustomErrorResponses>
      <Comment>{{ distribution.distribution_config.comment }}</Comment>
      <Logging>
        <Enabled>{{ distribution.distribution_config.logging.enabled }}</Enabled>
        <IncludeCookies>{{ distribution.distribution_config.logging.include_cookies }}</IncludeCookies>
        <Bucket>{{ distribution.distribution_config.logging.bucket }}</Bucket>
        <Prefix>{{ distribution.distribution_config.logging.prefix }}</Prefix>
      </Logging>
      <PriceClass>{{ distribution.distribution_config.price_class }}</PriceClass>
      <Enabled>{{ distribution.distribution_config.enabled }}</Enabled>
      <ViewerCertificate>
        <CloudFrontDefaultCertificate>{{ 'true' if distribution.distribution_config.viewer_certificate.cloud_front_default_certificate else 'false' }}</CloudFrontDefaultCertificate>
        <IAMCertificateId>{{ distribution.distribution_config.viewer_certificate.iam_certificate_id }}</IAMCertificateId>
        <ACMCertificateArn>{{ distribution.distribution_config.viewer_certificate.acm_certificate_arn }}</ACMCertificateArn>
        <SSLSupportMethod>{{ SSLSupportMethod }}</SSLSupportMethod>
        <MinimumProtocolVersion>{{ distribution.distribution_config.viewer_certificate.min_protocol_version }}</MinimumProtocolVersion>
        <Certificate>{{ distribution.distribution_config.viewer_certificate.certificate }}</Certificate>
        <CertificateSource>{{ distribution.distribution_config.viewer_certificate.certificate_source }}</CertificateSource>
      </ViewerCertificate>
      <Restrictions>
        <GeoRestriction>
          <RestrictionType>{{ distribution.distribution_config.geo_restriction._type }}</RestrictionType>
          <Quantity>{{ distribution.distribution_config.geo_restriction.restrictions|length }}</Quantity>
          {% if distribution.distribution_config.geo_restriction.restrictions %}
          <Items>
            {% for location  in distribution.distribution_config.geo_restriction.restrictions %}
              <Location>{{ location }}</Location>
            {% endfor %}
          </Items>
          {% endif %}
        </GeoRestriction>
      </Restrictions>
      <WebACLId>{{ distribution.distribution_config.web_acl_id }}</WebACLId>
      <HttpVersion>{{ distribution.distribution_config.http_version }}</HttpVersion>
      <IsIPV6Enabled>{{ 'true' if distribution.distribution_config.is_ipv6_enabled else 'false' }}</IsIPV6Enabled>
"""


DISTRIBUTION_TEMPLATE = (
    DIST_META_TEMPLATE
    + """
    <ActiveTrustedSigners>
      <Enabled>{{ distribution.active_trusted_signers.enabled }}</Enabled>
      <Quantity>{{ distribution.active_trusted_signers.quantity }}</Quantity>
      <Items>
        {% for signer  in distribution.active_trusted_signers.signers %}
          <AwsAccountNumber>{{ signer.aws_account_number }}</AwsAccountNumber>
          <KeyPairIds>
            <Quantity>{{ signer.key_pair_ids.quantity }}</Quantity>
            <Items>
              {% for key_pair_id_list  in signer.key_pair_ids.KeyPairIdList %}
                <KeyPairId>{{ key_pair_id_list.key_pair_id }}</KeyPairId>
              {% endfor %}
            </Items>
          </KeyPairIds>
        {% endfor %}
      </Items>
    </ActiveTrustedSigners>
    <ActiveTrustedKeyGroups>
      <Enabled>{{ distribution.active_trusted_key_groups.enabled }}</Enabled>
      <Quantity>{{ distribution.active_trusted_key_groups.quantity }}</Quantity>
      <Items>
        {% for kg_key_pair_id  in distribution.active_trusted_key_groups.kg_key_pair_ids %}
          <KeyGroupId>{{ kg_key_pair_id.key_group_id }}</KeyGroupId>
          <KeyPairIds>
            <Quantity>{{ kg_key_pair_id.key_pair_ids.quantity }}</Quantity>
            <Items>
              {% for key_pair_id_list  in kg_key_pair_ids_list.key_pair_ids.KeyPairIdList %}
                <KeyPairId>{{ key_pair_id_list.key_pair_id }}</KeyPairId>
              {% endfor %}
            </Items>
          </KeyPairIds>
        {% endfor %}
      </Items>
    </ActiveTrustedKeyGroups>
    <DistributionConfig>
      """
    + DIST_CONFIG_TEMPLATE
    + """
    </DistributionConfig>
    <AliasICPRecordals>
      {% for a  in distribution.alias_icp_recordals %}
        <CNAME>{{ a.cname }}</CNAME>
        <ICPRecordalStatus>{{ a.status }}</ICPRecordalStatus>
      {% endfor %}
    </AliasICPRecordals>"""
)

CREATE_DISTRIBUTION_TEMPLATE = (
    """<?xml version="1.0"?>
  <CreateDistributionResult xmlns="{{ xmlns }}">
"""
    + DISTRIBUTION_TEMPLATE
    + """
  </CreateDistributionResult>
"""
)

GET_DISTRIBUTION_TEMPLATE = (
    """<?xml version="1.0"?>
  <Distribution xmlns="{{ xmlns }}">
"""
    + DISTRIBUTION_TEMPLATE
    + """
  </Distribution>
"""
)

GET_DISTRIBUTION_CONFIG_TEMPLATE = (
    """<?xml version="1.0"?>
  <DistributionConfig>
"""
    + DIST_CONFIG_TEMPLATE
    + """
  </DistributionConfig>
"""
)


LIST_TEMPLATE = (
    """<?xml version="1.0"?>
<DistributionList xmlns="http://cloudfront.amazonaws.com/doc/2020-05-31/">
  <Marker></Marker>
  <MaxItems>100</MaxItems>
  <IsTruncated>false</IsTruncated>
  <Quantity>{{ distributions|length }}</Quantity>
  {% if distributions %}
  <Items>
      {% for distribution in distributions %}
      <DistributionSummary>
      """
    + DIST_META_TEMPLATE
    + """
      """
    + DIST_CONFIG_TEMPLATE
    + """
      </DistributionSummary>
      {% endfor %}
  </Items>
  {% endif %}
</DistributionList>"""
)

UPDATE_DISTRIBUTION_TEMPLATE = (
    """<?xml version="1.0"?>
  <Distribution xmlns="{{ xmlns }}">
"""
    + DISTRIBUTION_TEMPLATE
    + """
  </Distribution>
"""
)

CREATE_INVALIDATION_TEMPLATE = """<?xml version="1.0"?>
<Invalidation>
  <Id>{{ invalidation.invalidation_id }}</Id>
  <Status>{{ invalidation.status }}</Status>
  <CreateTime>{{ invalidation.create_time }}</CreateTime>
  <InvalidationBatch>
    <CallerReference>{{ invalidation.caller_ref }}</CallerReference>
    <Paths>
      <Quantity>{{ invalidation.paths|length }}</Quantity>
      <Items>
        {% for path in invalidation.paths %}<Path>{{ path }}</Path>{% endfor %}
      </Items>
    </Paths>
  </InvalidationBatch>
</Invalidation>
"""

INVALIDATIONS_TEMPLATE = """<?xml version="1.0" encoding="UTF-8"?>
<InvalidationList>
   <IsTruncated>false</IsTruncated>
   {% if invalidations %}
   <Items>
      {% for invalidation in invalidations %}
      <InvalidationSummary>
         <CreateTime>{{ invalidation.create_time }}</CreateTime>
         <Id>{{ invalidation.invalidation_id }}</Id>
         <Status>{{ invalidation.status }}</Status>
      </InvalidationSummary>
      {% endfor %}
   </Items>
   {% endif %}
   <Marker></Marker>
   <MaxItems>100</MaxItems>
   <Quantity>{{ invalidations|length }}</Quantity>
</InvalidationList>
"""

TAGS_TEMPLATE = """<?xml version="1.0"?>
<Tags>
  <Items>
    {% for tag in tags %}
      <Tag>
      <Key>{{ tag["Key"] }}</Key>
      <Value>{{ tag["Value"] }}</Value>
      </Tag>
    {% endfor %}
  </Items>
</Tags>
"""


ORIGIN_ACCESS_CONTROl = """<?xml version="1.0"?>
<OriginAccessControl>
  <Id>{{ control.id }}</Id>
  <OriginAccessControlConfig>
    <Name>{{ control.name }}</Name>
    {% if control.description %}
    <Description>{{ control.description }}</Description>
    {% endif %}
    <SigningProtocol>{{ control.signing_protocol }}</SigningProtocol>
    <SigningBehavior>{{ control.signing_behaviour }}</SigningBehavior>
    <OriginAccessControlOriginType>{{ control.origin_type }}</OriginAccessControlOriginType>
  </OriginAccessControlConfig>
</OriginAccessControl>
"""


LIST_ORIGIN_ACCESS_CONTROl = """<?xml version="1.0"?>
<OriginAccessControlList>
  <Items>
  {% for control in controls %}
    <OriginAccessControlSummary>
      <Id>{{ control.id }}</Id>
      <Name>{{ control.name }}</Name>
      {% if control.description %}
      <Description>{{ control.description }}</Description>
      {% endif %}
      <SigningProtocol>{{ control.signing_protocol }}</SigningProtocol>
      <SigningBehavior>{{ control.signing_behaviour }}</SigningBehavior>
      <OriginAccessControlOriginType>{{ control.origin_type }}</OriginAccessControlOriginType>
    </OriginAccessControlSummary>
  {% endfor %}
  </Items>
</OriginAccessControlList>
"""
