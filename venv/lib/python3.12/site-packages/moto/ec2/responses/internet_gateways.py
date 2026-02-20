from moto.core.responses import ActionResult, EmptyResult

from ._base_response import EC2BaseResponse


class InternetGateways(EC2BaseResponse):
    def attach_internet_gateway(self) -> ActionResult:
        igw_id = self._get_param("InternetGatewayId")
        vpc_id = self._get_param("VpcId")

        self.error_on_dryrun()

        self.ec2_backend.attach_internet_gateway(igw_id, vpc_id)
        return EmptyResult()

    def create_internet_gateway(self) -> str:
        self.error_on_dryrun()

        tags = self._get_param("TagSpecifications", [])
        if tags:
            tags = tags[0].get("Tags") or []
        igw = self.ec2_backend.create_internet_gateway(tags=tags)
        return self.response_template(CREATE_INTERNET_GATEWAY_RESPONSE).render(
            internet_gateway=igw
        )

    def delete_internet_gateway(self) -> ActionResult:
        igw_id = self._get_param("InternetGatewayId")
        self.error_on_dryrun()

        self.ec2_backend.delete_internet_gateway(igw_id)
        return EmptyResult()

    def describe_internet_gateways(self) -> str:
        filter_dict = self._filters_from_querystring()
        igw_ids = self._get_param("InternetGatewayIds", None)
        igws = self.ec2_backend.describe_internet_gateways(igw_ids, filters=filter_dict)
        template = self.response_template(DESCRIBE_INTERNET_GATEWAYS_RESPONSE)
        return template.render(internet_gateways=igws)

    def detach_internet_gateway(self) -> ActionResult:
        # TODO validate no instances with EIPs in VPC before detaching
        # raise else DependencyViolationError()
        igw_id = self._get_param("InternetGatewayId")
        vpc_id = self._get_param("VpcId")
        self.error_on_dryrun()

        self.ec2_backend.detach_internet_gateway(igw_id, vpc_id)
        return EmptyResult()


CREATE_INTERNET_GATEWAY_RESPONSE = """<CreateInternetGatewayResponse xmlns="http://ec2.amazonaws.com/doc/2013-10-15/">
  <requestId>59dbff89-35bd-4eac-99ed-be587EXAMPLE</requestId>
  <internetGateway>
    <internetGatewayId>{{ internet_gateway.id }}</internetGatewayId>
    <attachmentSet/>
    <ownerId>{{ internet_gateway.owner_id }}</ownerId>
    <tagSet>
      {% for tag in internet_gateway.get_tags() %}
        <item>
          <key>{{ tag.key }}</key>
          <value>{{ tag.value }}</value>
        </item>
      {% endfor %}
    </tagSet>
  </internetGateway>
</CreateInternetGatewayResponse>"""


DESCRIBE_INTERNET_GATEWAYS_RESPONSE = """<DescribeInternetGatewaysResponse xmlns="http://ec2.amazonaws.com/doc/2013-10-
15/">
  <requestId>59dbff89-35bd-4eac-99ed-be587EXAMPLE</requestId>
  <internetGatewaySet>
    {% for igw in internet_gateways %}
    <item>
      <internetGatewayId>{{ igw.id }}</internetGatewayId>
      <ownerId>{{ igw.owner_id or none }}</ownerId>
      {% if igw.vpc  %}
        <attachmentSet>
          <item>
            <vpcId>{{ igw.vpc.id }}</vpcId>
            <state>available</state>
          </item>
        </attachmentSet>
      {% else %}
        <attachmentSet/>
      {% endif %}
      <tagSet>
        {% for tag in igw.get_tags() %}
          <item>
            <key>{{ tag.key }}</key>
            <value>{{ tag.value }}</value>
          </item>
        {% endfor %}
      </tagSet>
    </item>
    {% endfor %}
  </internetGatewaySet>
</DescribeInternetGatewaysResponse>"""
