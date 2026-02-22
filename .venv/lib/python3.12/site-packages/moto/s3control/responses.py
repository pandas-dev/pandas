import json
import re
from typing import Any
from urllib.parse import unquote

import xmltodict

from moto.core.common_types import TYPE_RESPONSE
from moto.core.responses import (
    ActionResult,
    BaseResponse,
    EmptyResult,
    _get_method_urls,
)
from moto.s3.responses import S3_PUBLIC_ACCESS_BLOCK_CONFIGURATION

from .models import S3ControlBackend, s3control_backends


class S3ControlResponse(BaseResponse):
    def __init__(self) -> None:
        super().__init__(service_name="s3control")

    def _get_action_from_method_and_request_uri(
        self, method: str, request_uri: str
    ) -> str:
        """
        Override to sort patterns by length (descending) so more specific
        patterns match before more general ones.

        This fixes issues where patterns like:
        - /mrap/instances/{name+} (GetMultiRegionAccessPoint)
        - /mrap/instances/{name+}/policy (GetMultiRegionAccessPointPolicy)

        Both match the same URL, but we want the more specific one to win.
        """
        methods_url = _get_method_urls(self.service_name, self.region)
        regexp_and_names = methods_url[method]
        sorted_patterns = sorted(
            regexp_and_names.items(), key=lambda x: len(x[0]), reverse=True
        )
        for regexp, name in sorted_patterns:
            match = re.match(regexp, request_uri)
            self.uri_match = match
            if match:
                return name
        return None  # type: ignore[return-value]

    @property
    def backend(self) -> S3ControlBackend:
        return s3control_backends[self.current_account][self.partition]

    def get_public_access_block(self) -> str:
        account_id = self.headers.get("x-amz-account-id")
        public_block_config = self.backend.get_public_access_block(
            account_id=account_id
        )
        template = self.response_template(S3_PUBLIC_ACCESS_BLOCK_CONFIGURATION)
        return template.render(public_block_config=public_block_config)

    def put_public_access_block(self) -> TYPE_RESPONSE:
        account_id = self.headers.get("x-amz-account-id")
        pab_config = self._parse_pab_config(self.body)
        self.backend.put_public_access_block(
            account_id, pab_config["PublicAccessBlockConfiguration"]
        )
        return 201, {"status": 201}, json.dumps({})

    def delete_public_access_block(self) -> TYPE_RESPONSE:
        account_id = self.headers.get("x-amz-account-id")
        self.backend.delete_public_access_block(account_id=account_id)
        return 204, {"status": 204}, json.dumps({})

    def _parse_pab_config(self, body: str) -> dict[str, Any]:
        parsed_xml = xmltodict.parse(body)
        parsed_xml["PublicAccessBlockConfiguration"].pop("@xmlns", None)

        return parsed_xml

    def create_access_point(self) -> str:
        account_id, name = self._get_accountid_and_name_from_accesspoint(self.uri)
        params = xmltodict.parse(self.body)["CreateAccessPointRequest"]
        bucket = params["Bucket"]
        vpc_configuration = params.get("VpcConfiguration")
        public_access_block_configuration = params.get("PublicAccessBlockConfiguration")
        access_point = self.backend.create_access_point(
            account_id=account_id,
            name=name,
            bucket=bucket,
            vpc_configuration=vpc_configuration,
            public_access_block_configuration=public_access_block_configuration,
        )
        template = self.response_template(CREATE_ACCESS_POINT_TEMPLATE)
        return template.render(access_point=access_point)

    def get_access_point(self) -> str:
        account_id, name = self._get_accountid_and_name_from_accesspoint(self.uri)

        access_point = self.backend.get_access_point(account_id=account_id, name=name)
        template = self.response_template(GET_ACCESS_POINT_TEMPLATE)
        return template.render(access_point=access_point)

    def delete_access_point(self) -> TYPE_RESPONSE:
        account_id, name = self._get_accountid_and_name_from_accesspoint(self.uri)
        self.backend.delete_access_point(account_id=account_id, name=name)
        return 204, {"status": 204}, ""

    def put_access_point_policy(self) -> str:
        account_id, name = self._get_accountid_and_name_from_policy(self.uri)
        params = xmltodict.parse(self.body)
        policy = params["PutAccessPointPolicyRequest"]["Policy"]
        self.backend.put_access_point_policy(account_id, name, policy)
        return ""

    def get_access_point_policy(self) -> str:
        account_id, name = self._get_accountid_and_name_from_policy(self.uri)
        policy = self.backend.get_access_point_policy(account_id, name)
        template = self.response_template(GET_ACCESS_POINT_POLICY_TEMPLATE)
        return template.render(policy=policy)

    def delete_access_point_policy(self) -> TYPE_RESPONSE:
        account_id, name = self._get_accountid_and_name_from_policy(self.uri)
        self.backend.delete_access_point_policy(account_id=account_id, name=name)
        return 204, {"status": 204}, ""

    def get_access_point_policy_status(self) -> str:
        account_id, name = self._get_accountid_and_name_from_policy(self.uri)
        self.backend.get_access_point_policy_status(account_id, name)
        template = self.response_template(GET_ACCESS_POINT_POLICY_STATUS_TEMPLATE)
        return template.render()

    def _get_accountid_and_name_from_accesspoint(
        self, full_url: str
    ) -> tuple[str, str]:
        url = full_url
        if full_url.startswith("http"):
            url = full_url.split("://")[1]
        account_id = url.split(".")[0]
        name = url.split("v20180820/accesspoint/")[-1]
        return account_id, name

    def _get_accountid_and_name_from_policy(self, full_url: str) -> tuple[str, str]:
        url = full_url
        if full_url.startswith("http"):
            url = full_url.split("://")[1]
        account_id = url.split(".")[0]
        name = self.path.split("/")[-2]
        return account_id, name

    def put_storage_lens_configuration(self) -> str:
        account_id = self.headers.get("x-amz-account-id")
        config_id = self.path.split("/")[-1]
        request = xmltodict.parse(self.body)["PutStorageLensConfigurationRequest"]
        storage_lens_configuration = request.get("StorageLensConfiguration")
        tags = request.get("Tags")
        self.backend.put_storage_lens_configuration(
            config_id=config_id,
            account_id=account_id,
            storage_lens_configuration=storage_lens_configuration,
            tags=tags,
        )
        return ""

    def get_storage_lens_configuration(self) -> str:
        account_id = self.headers.get("x-amz-account-id")
        config_id = self.path.split("/")[-1]
        storage_lens_configuration = self.backend.get_storage_lens_configuration(
            config_id=config_id,
            account_id=account_id,
        )
        template = self.response_template(GET_STORAGE_LENS_CONFIGURATION_TEMPLATE)
        return template.render(config=storage_lens_configuration.config)

    def delete_storage_lens_configuration(self) -> TYPE_RESPONSE:
        account_id = self.headers.get("x-amz-account-id")
        config_id = self.path.split("/")[-1]
        self.backend.delete_storage_lens_configuration(
            config_id=config_id,
            account_id=account_id,
        )
        return 204, {"status": 204}, ""

    def list_storage_lens_configurations(self) -> str:
        account_id = self.headers.get("x-amz-account-id")
        params = self._get_params()
        next_token = params.get("nextToken")
        storage_lens_configuration_list, next_token = (
            self.backend.list_storage_lens_configurations(
                account_id=account_id,
                next_token=next_token,
            )
        )
        template = self.response_template(LIST_STORAGE_LENS_CONFIGURATIONS_TEMPLATE)
        return template.render(
            next_token=next_token, configs=storage_lens_configuration_list
        )

    def put_storage_lens_configuration_tagging(self) -> str:
        account_id = self.headers.get("x-amz-account-id")
        config_id = self.path.split("/")[-2]
        request = xmltodict.parse(self.body)[
            "PutStorageLensConfigurationTaggingRequest"
        ]
        tags = request.get("Tags")
        self.backend.put_storage_lens_configuration_tagging(
            config_id=config_id,
            account_id=account_id,
            tags=tags,
        )
        return ""

    def get_storage_lens_configuration_tagging(self) -> str:
        account_id = self.headers.get("x-amz-account-id")
        config_id = self.path.split("/")[-2]
        storage_lens_tags = self.backend.get_storage_lens_configuration_tagging(
            config_id=config_id,
            account_id=account_id,
        )
        template = self.response_template(
            GET_STORAGE_LENS_CONFIGURATION_TAGGING_TEMPLATE
        )
        return template.render(tags=storage_lens_tags)

    def list_access_points(self) -> str:
        account_id = self.headers.get("x-amz-account-id")

        params = self._get_params()
        max_results = params.get("maxResults")
        if max_results:
            max_results = int(max_results)

        access_points, next_token = self.backend.list_access_points(
            account_id=account_id,
            bucket=params.get("bucket"),
            max_results=max_results,
            next_token=params.get("nextToken"),
        )

        template = self.response_template(LIST_ACCESS_POINTS_TEMPLATE)
        return template.render(access_points=access_points, next_token=next_token)

    def create_multi_region_access_point(self) -> str:
        account_id = self.headers.get("x-amz-account-id")
        params = xmltodict.parse(self.body, process_namespaces=False)[
            "CreateMultiRegionAccessPointRequest"
        ]

        details = params["Details"]
        name = details.get("Name")

        regions_data = details.get("Regions", {})
        if regions_data and "Region" in regions_data:
            regions_list = regions_data["Region"]
            if isinstance(regions_list, dict):
                regions_list = [regions_list]
        else:
            regions_list = []

        public_access_block = details.get("PublicAccessBlock", {})

        operation = self.backend.create_multi_region_access_point(
            account_id=account_id,
            name=name,
            public_access_block=public_access_block,
            regions=regions_list,
            region_name=self.region,
        )

        template = self.response_template(CREATE_MULTI_REGION_ACCESS_POINT_TEMPLATE)
        return template.render(request_token=operation.request_token_arn)

    def delete_multi_region_access_point(self) -> str:
        account_id = self.headers.get("x-amz-account-id")
        params = xmltodict.parse(self.body, process_namespaces=False)[
            "DeleteMultiRegionAccessPointRequest"
        ]

        details = params["Details"]
        name = details.get("Name")

        operation = self.backend.delete_multi_region_access_point(
            account_id=account_id,
            name=name,
            region_name=self.region,
        )

        template = self.response_template(DELETE_MULTI_REGION_ACCESS_POINT_TEMPLATE)
        return template.render(request_token=operation.request_token_arn)

    def describe_multi_region_access_point_operation(self) -> str:
        account_id = self.headers.get("x-amz-account-id")
        _prefix = "/async-requests/mrap/"
        if _prefix in self.path:
            request_token = unquote(self.path.partition(_prefix)[2])
        else:
            request_token = self.path.split("/")[-1]

        operation = self.backend.describe_multi_region_access_point_operation(
            account_id=account_id,
            request_token_arn=request_token,
        )

        template = self.response_template(
            DESCRIBE_MULTI_REGION_ACCESS_POINT_OPERATION_TEMPLATE
        )
        return template.render(operation=operation.to_dict())

    def get_multi_region_access_point(self) -> str:
        account_id = self.headers.get("x-amz-account-id")
        name = self.path.split("/")[-1]

        mrap = self.backend.get_multi_region_access_point(
            account_id=account_id,
            name=name,
        )

        template = self.response_template(GET_MULTI_REGION_ACCESS_POINT_TEMPLATE)
        return template.render(mrap=mrap.to_dict())

    def get_multi_region_access_point_policy(self) -> str:
        account_id = self.headers.get("x-amz-account-id")
        name = self.path.split("/")[-2]

        policy = self.backend.get_multi_region_access_point_policy(
            account_id=account_id,
            name=name,
        )

        template = self.response_template(GET_MULTI_REGION_ACCESS_POINT_POLICY_TEMPLATE)
        return template.render(policy=policy)

    def get_multi_region_access_point_policy_status(self) -> str:
        account_id = self.headers.get("x-amz-account-id")
        name = self.path.split("/")[-2]

        policy_status = self.backend.get_multi_region_access_point_policy_status(
            account_id=account_id,
            name=name,
        )

        template = self.response_template(
            GET_MULTI_REGION_ACCESS_POINT_POLICY_STATUS_TEMPLATE
        )
        return template.render(is_public=policy_status["IsPublic"])

    def list_multi_region_access_points(self) -> str:
        account_id = self.headers.get("x-amz-account-id")
        params = self._get_params()

        max_results = params.get("maxResults")
        if max_results:
            max_results = int(max_results)

        mraps, next_token = self.backend.list_multi_region_access_points(
            account_id=account_id,
            max_results=max_results,
            next_token=params.get("nextToken"),
        )

        template = self.response_template(LIST_MULTI_REGION_ACCESS_POINTS_TEMPLATE)
        return template.render(
            mraps=[mrap.to_dict() for mrap in mraps], next_token=next_token
        )

    def put_multi_region_access_point_policy(self) -> str:
        account_id = self.headers.get("x-amz-account-id")
        params = xmltodict.parse(self.body, process_namespaces=False)[
            "PutMultiRegionAccessPointPolicyRequest"
        ]

        details = params["Details"]
        name = details.get("Name")
        policy = details.get("Policy")

        operation = self.backend.put_multi_region_access_point_policy(
            account_id=account_id,
            name=name,
            policy=policy,
            region_name=self.region,
        )

        template = self.response_template(PUT_MULTI_REGION_ACCESS_POINT_POLICY_TEMPLATE)
        return template.render(request_token=operation.request_token_arn)

    def list_tags_for_resource(self) -> ActionResult:
        resource_arn = unquote(self.parsed_url.path.split("/tags/")[-1])
        tags = self.backend.list_tags_for_resource(resource_arn)
        return ActionResult(result={"Tags": tags})

    def tag_resource(self) -> EmptyResult:
        resource_arn = unquote(self.parsed_url.path.split("/tags/")[-1])
        tags = (
            xmltodict.parse(self.raw_body, force_list={"Tag": True})
            .get("TagResourceRequest", {})
            .get("Tags", {})["Tag"]
        )
        self.backend.tag_resource(resource_arn, tags=tags)
        return EmptyResult()

    def untag_resource(self) -> EmptyResult:
        resource_arn = unquote(self.parsed_url.path.split("/tags/")[-1])
        tag_keys = self.querystring.get("tagKeys", [])
        self.backend.untag_resource(resource_arn, tag_keys=tag_keys)
        return EmptyResult()


CREATE_ACCESS_POINT_TEMPLATE = """<CreateAccessPointResult>
  <ResponseMetadata>
    <RequestId>1549581b-12b7-11e3-895e-1334aEXAMPLE</RequestId>
  </ResponseMetadata>
  <Alias>{{ access_point.alias }}</Alias>
  <AccessPointArn>{{ access_point.arn }}</AccessPointArn>
</CreateAccessPointResult>
"""


GET_ACCESS_POINT_TEMPLATE = """<GetAccessPointResult>
  <ResponseMetadata>
    <RequestId>1549581b-12b7-11e3-895e-1334aEXAMPLE</RequestId>
  </ResponseMetadata>
  <Name>{{ access_point.name }}</Name>
  <Bucket>{{ access_point.bucket }}</Bucket>
  <NetworkOrigin>{{ access_point.network_origin }}</NetworkOrigin>
  {% if access_point.vpc_id %}
  <VpcConfiguration>
      <VpcId>{{ access_point.vpc_id }}</VpcId>
  </VpcConfiguration>
  {% endif %}
  <PublicAccessBlockConfiguration>
      <BlockPublicAcls>{{ access_point.pubc["BlockPublicAcls"] }}</BlockPublicAcls>
      <IgnorePublicAcls>{{ access_point.pubc["IgnorePublicAcls"] }}</IgnorePublicAcls>
      <BlockPublicPolicy>{{ access_point.pubc["BlockPublicPolicy"] }}</BlockPublicPolicy>
      <RestrictPublicBuckets>{{ access_point.pubc["RestrictPublicBuckets"] }}</RestrictPublicBuckets>
  </PublicAccessBlockConfiguration>
  <CreationDate>{{ access_point.created }}</CreationDate>
  <Alias>{{ access_point.alias }}</Alias>
  <AccessPointArn>{{ access_point.arn }}</AccessPointArn>
  <Endpoints>
      <entry>
          <key>ipv4</key>
          <value>s3-accesspoint.us-east-1.amazonaws.com</value>
      </entry>
      <entry>
          <key>fips</key>
          <value>s3-accesspoint-fips.us-east-1.amazonaws.com</value>
      </entry>
      <entry>
          <key>fips_dualstack</key>
          <value>s3-accesspoint-fips.dualstack.us-east-1.amazonaws.com</value>
      </entry>
      <entry>
          <key>dualstack</key>
          <value>s3-accesspoint.dualstack.us-east-1.amazonaws.com</value>
      </entry>
  </Endpoints>
</GetAccessPointResult>
"""


GET_ACCESS_POINT_POLICY_TEMPLATE = """<GetAccessPointPolicyResult>
  <ResponseMetadata>
    <RequestId>1549581b-12b7-11e3-895e-1334aEXAMPLE</RequestId>
  </ResponseMetadata>
  <Policy>{{ policy }}</Policy>
</GetAccessPointPolicyResult>
"""


GET_ACCESS_POINT_POLICY_STATUS_TEMPLATE = """<GetAccessPointPolicyResult>
  <ResponseMetadata>
    <RequestId>1549581b-12b7-11e3-895e-1334aEXAMPLE</RequestId>
  </ResponseMetadata>
  <PolicyStatus>
      <IsPublic>true</IsPublic>
  </PolicyStatus>
</GetAccessPointPolicyResult>
"""

XMLNS = 'xmlns="http://awss3control.amazonaws.com/doc/2018-08-20/"'

CREATE_MULTI_REGION_ACCESS_POINT_TEMPLATE = f"""<CreateMultiRegionAccessPointResult {XMLNS}>
  <RequestTokenARN>{{{{ request_token }}}}</RequestTokenARN>
</CreateMultiRegionAccessPointResult>
"""

DELETE_MULTI_REGION_ACCESS_POINT_TEMPLATE = f"""<DeleteMultiRegionAccessPointResult {XMLNS}>
  <RequestTokenARN>{{{{ request_token }}}}</RequestTokenARN>
</DeleteMultiRegionAccessPointResult>
"""

DESCRIBE_MULTI_REGION_ACCESS_POINT_OPERATION_TEMPLATE = f"""<DescribeMultiRegionAccessPointOperationResult {XMLNS}>
  <AsyncOperation>
    <CreationTime>{{{{ operation.CreationTime }}}}</CreationTime>
    <Operation>{{{{ operation.Operation }}}}</Operation>
    <RequestTokenARN>{{{{ operation.RequestTokenARN }}}}</RequestTokenARN>
    <RequestStatus>{{{{ operation.RequestStatus }}}}</RequestStatus>
    {{% if operation.RequestParameters %}}
    <RequestParameters>
      {{% if operation.RequestParameters.CreateMultiRegionAccessPointRequest %}}
      <CreateMultiRegionAccessPointRequest>
        <Name>{{{{ operation.RequestParameters.CreateMultiRegionAccessPointRequest.Name }}}}</Name>
        {{% if operation.RequestParameters.CreateMultiRegionAccessPointRequest.PublicAccessBlock %}}
        <PublicAccessBlock>
          <BlockPublicAcls>{{{{ operation.RequestParameters.CreateMultiRegionAccessPointRequest.PublicAccessBlock.BlockPublicAcls|default("true") }}}}</BlockPublicAcls>
          <IgnorePublicAcls>{{{{ operation.RequestParameters.CreateMultiRegionAccessPointRequest.PublicAccessBlock.IgnorePublicAcls|default("true") }}}}</IgnorePublicAcls>
          <BlockPublicPolicy>{{{{ operation.RequestParameters.CreateMultiRegionAccessPointRequest.PublicAccessBlock.BlockPublicPolicy|default("true") }}}}</BlockPublicPolicy>
          <RestrictPublicBuckets>{{{{ operation.RequestParameters.CreateMultiRegionAccessPointRequest.PublicAccessBlock.RestrictPublicBuckets|default("true") }}}}</RestrictPublicBuckets>
        </PublicAccessBlock>
        {{% endif %}}
        {{% if operation.RequestParameters.CreateMultiRegionAccessPointRequest.Regions %}}
        <Regions>
          {{% for region in operation.RequestParameters.CreateMultiRegionAccessPointRequest.Regions %}}
          <Region>
            <Bucket>{{{{ region.Bucket }}}}</Bucket>
          </Region>
          {{% endfor %}}
        </Regions>
        {{% endif %}}
      </CreateMultiRegionAccessPointRequest>
      {{% endif %}}
      {{% if operation.RequestParameters.DeleteMultiRegionAccessPointRequest %}}
      <DeleteMultiRegionAccessPointRequest>
        <Name>{{{{ operation.RequestParameters.DeleteMultiRegionAccessPointRequest.Name }}}}</Name>
      </DeleteMultiRegionAccessPointRequest>
      {{% endif %}}
      {{% if operation.RequestParameters.PutMultiRegionAccessPointPolicyRequest %}}
      <PutMultiRegionAccessPointPolicyRequest>
        <Name>{{{{ operation.RequestParameters.PutMultiRegionAccessPointPolicyRequest.Name }}}}</Name>
        <Policy>{{{{ operation.RequestParameters.PutMultiRegionAccessPointPolicyRequest.Policy }}}}</Policy>
      </PutMultiRegionAccessPointPolicyRequest>
      {{% endif %}}
    </RequestParameters>
    {{% endif %}}
    {{% if operation.ResponseDetails %}}
    <ResponseDetails>
      {{% if operation.ResponseDetails.MultiRegionAccessPointDetails %}}
      <MultiRegionAccessPointDetails>
        {{% if operation.ResponseDetails.MultiRegionAccessPointDetails.Regions %}}
        <Regions>
          {{% for region in operation.ResponseDetails.MultiRegionAccessPointDetails.Regions %}}
          <Region>
            <Name>{{{{ region.Name }}}}</Name>
            <RequestStatus>{{{{ region.RequestStatus }}}}</RequestStatus>
          </Region>
          {{% endfor %}}
        </Regions>
        {{% endif %}}
      </MultiRegionAccessPointDetails>
      {{% endif %}}
    </ResponseDetails>
    {{% endif %}}
  </AsyncOperation>
</DescribeMultiRegionAccessPointOperationResult>
"""

GET_MULTI_REGION_ACCESS_POINT_TEMPLATE = f"""<GetMultiRegionAccessPointResult {XMLNS}>
  <AccessPoint>
    <Name>{{{{ mrap.Name }}}}</Name>
    <Alias>{{{{ mrap.Alias }}}}</Alias>
    <CreatedAt>{{{{ mrap.CreatedAt }}}}</CreatedAt>
    <PublicAccessBlock>
      <BlockPublicAcls>{{{{ mrap.PublicAccessBlock.BlockPublicAcls|default("true") }}}}</BlockPublicAcls>
      <IgnorePublicAcls>{{{{ mrap.PublicAccessBlock.IgnorePublicAcls|default("true") }}}}</IgnorePublicAcls>
      <BlockPublicPolicy>{{{{ mrap.PublicAccessBlock.BlockPublicPolicy|default("true") }}}}</BlockPublicPolicy>
      <RestrictPublicBuckets>{{{{ mrap.PublicAccessBlock.RestrictPublicBuckets|default("true") }}}}</RestrictPublicBuckets>
    </PublicAccessBlock>
    <Status>{{{{ mrap.Status }}}}</Status>
    <Regions>
      {{% for region in mrap.Regions %}}
      <Region>
        <Bucket>{{{{ region.Bucket }}}}</Bucket>
        <Region>{{{{ region.Region }}}}</Region>
      </Region>
      {{% endfor %}}
    </Regions>
  </AccessPoint>
</GetMultiRegionAccessPointResult>
"""

GET_MULTI_REGION_ACCESS_POINT_POLICY_TEMPLATE = f"""<GetMultiRegionAccessPointPolicyResult {XMLNS}>
  <Policy>
    <Established>
      <Policy>{{{{ policy }}}}</Policy>
    </Established>
  </Policy>
</GetMultiRegionAccessPointPolicyResult>
"""

GET_MULTI_REGION_ACCESS_POINT_POLICY_STATUS_TEMPLATE = f"""<GetMultiRegionAccessPointPolicyStatusResult {XMLNS}>
  <Established>
    <IsPublic>{{{{ is_public|lower }}}}</IsPublic>
  </Established>
</GetMultiRegionAccessPointPolicyStatusResult>
"""

LIST_MULTI_REGION_ACCESS_POINTS_TEMPLATE = f"""<ListMultiRegionAccessPointsResult {XMLNS}>
  <AccessPoints>
    {{% for mrap in mraps %}}
    <AccessPoint>
      <Name>{{{{ mrap.Name }}}}</Name>
      <Alias>{{{{ mrap.Alias }}}}</Alias>
      <CreatedAt>{{{{ mrap.CreatedAt }}}}</CreatedAt>
      <PublicAccessBlock>
        <BlockPublicAcls>{{{{ mrap.PublicAccessBlock.BlockPublicAcls|default("true") }}}}</BlockPublicAcls>
        <IgnorePublicAcls>{{{{ mrap.PublicAccessBlock.IgnorePublicAcls|default("true") }}}}</IgnorePublicAcls>
        <BlockPublicPolicy>{{{{ mrap.PublicAccessBlock.BlockPublicPolicy|default("true") }}}}</BlockPublicPolicy>
        <RestrictPublicBuckets>{{{{ mrap.PublicAccessBlock.RestrictPublicBuckets|default("true") }}}}</RestrictPublicBuckets>
      </PublicAccessBlock>
      <Status>{{{{ mrap.Status }}}}</Status>
      <Regions>
        {{% for region in mrap.Regions %}}
        <Region>
          <Bucket>{{{{ region.Bucket }}}}</Bucket>
          <Region>{{{{ region.Region }}}}</Region>
        </Region>
        {{% endfor %}}
      </Regions>
    </AccessPoint>
    {{% endfor %}}
  </AccessPoints>
  {{% if next_token %}}
  <NextToken>{{{{ next_token }}}}</NextToken>
  {{% endif %}}
</ListMultiRegionAccessPointsResult>
"""

PUT_MULTI_REGION_ACCESS_POINT_POLICY_TEMPLATE = f"""<PutMultiRegionAccessPointPolicyResult {XMLNS}>
  <RequestTokenARN>{{{{ request_token }}}}</RequestTokenARN>
</PutMultiRegionAccessPointPolicyResult>
"""

GET_STORAGE_LENS_CONFIGURATION_TEMPLATE = """
<StorageLensConfiguration>
   <Id>{{config.get("Id")}}</Id>
   {% if config.get("DataExport") %}
   <DataExport>
      {% if config["DataExport"]["S3BucketDestination"] %}
      <S3BucketDestination>
         <AccountId>{{config["DataExport"]["S3BucketDestination"]["AccountId"]}}</AccountId>
         <Arn>{{config["DataExport"]["S3BucketDestination"]["Arn"]}}</Arn>
         {% if config["DataExport"]["S3BucketDestination"].get("Encryption") %}
         <Encryption>
            {% if config["DataExport"]["S3BucketDestination"]["Encryption"].get("SSEKMS") %}
            <SSE-KMS>
               <KeyId>config["DataExport"]["S3BucketDestination"]["Encryption"]["KeyId"]</KeyId>
            </SSE-KMS>
            {% endif %}
            {% if "SSE-S3" in config["DataExport"]["S3BucketDestination"]["Encryption"] %}
            <SSE-S3>
            </SSE-S3>
            {% endif %}
         </Encryption>
         {% endif %}
      </S3BucketDestination>
      {% endif %}
   </DataExport>
   {% endif %}
   <IsEnabled>{{config["IsEnabled"]}}</IsEnabled>
   <AccountLevel>
        <ActivityMetrics>
            <IsEnabled>{{config["AccountLevel"]["ActivityMetrics"]["IsEnabled"]}}</IsEnabled>
        </ActivityMetrics>
        <BucketLevel>
            <ActivityMetrics>
                <IsEnabled>{{config["AccountLevel"]["BucketLevel"]["ActivityMetrics"]["IsEnabled"]}}</IsEnabled>
            </ActivityMetrics>
            <PrefixLevel>
                <StorageMetrics>
                    <IsEnabled>{{config["AccountLevel"]["BucketLevel"]["PrefixLevel"]["StorageMetrics"]["IsEnabled"]}}</IsEnabled>
                    <SelectionCriteria>
                        <Delimiter>{{config["AccountLevel"]["BucketLevel"]["PrefixLevel"]["StorageMetrics"]["SelectionCriteria"]["Delimiter"]}}</Delimiter>
                        <MaxDepth>{{config["AccountLevel"]["BucketLevel"]["PrefixLevel"]["StorageMetrics"]["SelectionCriteria"]["MaxDepth"]}}</MaxDepth>
                        <MinStorageBytesPercentage>{{config["AccountLevel"]["BucketLevel"]["PrefixLevel"]["StorageMetrics"]["SelectionCriteria"]["MinStorageBytesPercentage"]}}</MinStorageBytesPercentage>
                    </SelectionCriteria>
                </StorageMetrics>
            </PrefixLevel>
            <DetailedStatusCodesMetrics>
                <IsEnabled>{{config["AccountLevel"]["BucketLevel"]["DetailedStatusCodesMetrics"]["IsEnabled"]}}</IsEnabled>
            </DetailedStatusCodesMetrics>
        </BucketLevel>
        <AdvancedDataProtectionMetrics>
            <IsEnabled>{{config["AccountLevel"]["AdvancedDataProtectionMetrics"]["IsEnabled"]}}</IsEnabled>
        </AdvancedDataProtectionMetrics>
        <DetailedStatusCodesMetrics>
            <IsEnabled>{{config["AccountLevel"]["DetailedStatusCodesMetrics"]["IsEnabled"]}}</IsEnabled>
        </DetailedStatusCodesMetrics>
   </AccountLevel>
   <AwsOrg>
        <Arn>{{config.get("AwsOrg", {}).get("Arn", "")}}</Arn>
    </AwsOrg>
    <StorageLensArn>{{config.get("StorageLensArn")}}</StorageLensArn>
</StorageLensConfiguration>
"""


LIST_STORAGE_LENS_CONFIGURATIONS_TEMPLATE = """
<ListStorageLensConfigurationsResult>
   {% if next_token %}
   <NextToken>{{ next_token }}</NextToken>
   {% endif %}
   {% for config in configs %}
   <StorageLensConfiguration>
      <HomeRegion></HomeRegion>
      <Id>{{ config.config.get("Id") }}</Id>
      <IsEnabled>{{ config.config.get("IsEnabled") }}</IsEnabled>
      <StorageLensArn>{{ config.arn }}</StorageLensArn>
    </StorageLensConfiguration>
    {% endfor %}
</ListStorageLensConfigurationsResult>
"""


GET_STORAGE_LENS_CONFIGURATION_TAGGING_TEMPLATE = """
<GetStorageLensConfigurationTaggingResult>
   <Tags>
      {% for tag in tags["Tag"] %}
      <Tag>
         <Key>{{ tag["Key"] }}</Key>
         <Value>{{ tag["Value"] }}</Value>
      </Tag>
      {% endfor %}
   </Tags>
</GetStorageLensConfigurationTaggingResult>

"""
LIST_ACCESS_POINTS_TEMPLATE = """<ListAccessPointsResult>
  <AccessPointList>
    {% for access_point in access_points %}
    <AccessPoint>
      <Name>{{ access_point.name }}</Name>
      <NetworkOrigin>{{ access_point.network_origin }}</NetworkOrigin>
      {% if access_point.vpc_id %}
      <VpcConfiguration>
        <VpcId>{{ access_point.vpc_id }}</VpcId>
      </VpcConfiguration>
      {% endif %}
      <Bucket>{{ access_point.bucket }}</Bucket>
      <AccessPointArn>{{ access_point.arn }}</AccessPointArn>
      <Alias>{{ access_point.alias }}</Alias>
    </AccessPoint>
    {% endfor %}
  </AccessPointList>
  {% if next_token %}
  <NextToken>{{ next_token }}</NextToken>
  {% endif %}
</ListAccessPointsResult>"""
