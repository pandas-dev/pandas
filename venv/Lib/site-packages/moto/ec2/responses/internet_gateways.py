from ._base_response import EC2BaseResponse


class InternetGateways(EC2BaseResponse):
    def attach_internet_gateway(self) -> str:
        igw_id = self._get_param("InternetGatewayId")
        vpc_id = self._get_param("VpcId")

        self.error_on_dryrun()

        self.ec2_backend.attach_internet_gateway(igw_id, vpc_id)
        return self.response_template(ATTACH_INTERNET_GATEWAY_RESPONSE).render()

    def create_internet_gateway(self) -> str:
        self.error_on_dryrun()

        tags = self._get_multi_param("TagSpecification", skip_result_conversion=True)
        if tags:
            tags = tags[0].get("Tag") or []
        igw = self.ec2_backend.create_internet_gateway(tags=tags)
        return self.response_template(CREATE_INTERNET_GATEWAY_RESPONSE).render(
            internet_gateway=igw
        )

    def delete_internet_gateway(self) -> str:
        igw_id = self._get_param("InternetGatewayId")
        self.error_on_dryrun()

        self.ec2_backend.delete_internet_gateway(igw_id)
        return self.response_template(DELETE_INTERNET_GATEWAY_RESPONSE).render()

    def describe_internet_gateways(self) -> str:
        filter_dict = self._filters_from_querystring()
        if "InternetGatewayId.1" in self.querystring:
            igw_ids = self._get_multi_param("InternetGatewayId")
            igws = self.ec2_backend.describe_internet_gateways(
                igw_ids, filters=filter_dict
            )
        else:
            igws = self.ec2_backend.describe_internet_gateways(filters=filter_dict)

        template = self.response_template(DESCRIBE_INTERNET_GATEWAYS_RESPONSE)
        return template.render(internet_gateways=igws)

    def detach_internet_gateway(self) -> str:
        # TODO validate no instances with EIPs in VPC before detaching
        # raise else DependencyViolationError()
        igw_id = self._get_param("InternetGatewayId")
        vpc_id = self._get_param("VpcId")
        self.error_on_dryrun()

        self.ec2_backend.detach_internet_gateway(igw_id, vpc_id)
        return self.response_template(DETACH_INTERNET_GATEWAY_RESPONSE).render()


ATTACH_INTERNET_GATEWAY_RESPONSE = """<AttachInternetGatewayResponse xmlns="http://ec2.amazonaws.com/doc/2013-10-15/">
  <requestId>59dbff89-35bd-4eac-99ed-be587EXAMPLE</requestId>
  <return>true</return>
</AttachInternetGatewayResponse>"""

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

DELETE_INTERNET_GATEWAY_RESPONSE = """<DeleteInternetGatewayResponse xmlns="http://ec2.amazonaws.com/doc/2013-10-15/">
    <requestId>59dbff89-35bd-4eac-99ed-be587EXAMPLE</requestId>
    <return>true</return>
</DeleteInternetGatewayResponse>"""

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

DETACH_INTERNET_GATEWAY_RESPONSE = """<DetachInternetGatewayResponse xmlns="http://ec2.amazonaws.com/doc/2013-10-15/">
  <requestId>59dbff89-35bd-4eac-99ed-be587EXAMPLE</requestId>
  <return>true</return>
</DetachInternetGatewayResponse>"""
