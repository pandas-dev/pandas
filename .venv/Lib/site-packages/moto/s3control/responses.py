import json
from typing import Any, Dict, Tuple

import xmltodict

from moto.core.common_types import TYPE_RESPONSE
from moto.core.responses import BaseResponse
from moto.s3.responses import S3_PUBLIC_ACCESS_BLOCK_CONFIGURATION

from .models import S3ControlBackend, s3control_backends


class S3ControlResponse(BaseResponse):
    def __init__(self) -> None:
        super().__init__(service_name="s3control")

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

    def _parse_pab_config(self, body: str) -> Dict[str, Any]:
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
    ) -> Tuple[str, str]:
        url = full_url
        if full_url.startswith("http"):
            url = full_url.split("://")[1]
        account_id = url.split(".")[0]
        name = url.split("v20180820/accesspoint/")[-1]
        return account_id, name

    def _get_accountid_and_name_from_policy(self, full_url: str) -> Tuple[str, str]:
        url = full_url
        if full_url.startswith("http"):
            url = full_url.split("://")[1]
        account_id = url.split(".")[0]
        name = self.path.split("/")[-2]
        return account_id, name


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
