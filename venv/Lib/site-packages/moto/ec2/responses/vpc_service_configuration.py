from ..exceptions import NoLoadBalancersProvided
from ._base_response import EC2BaseResponse


class VPCEndpointServiceConfiguration(EC2BaseResponse):
    def create_vpc_endpoint_service_configuration(self) -> str:
        gateway_lbs = self._get_multi_param("GatewayLoadBalancerArn")
        network_lbs = self._get_multi_param("NetworkLoadBalancerArn")
        if not gateway_lbs and not network_lbs:
            raise NoLoadBalancersProvided

        tags = self._get_multi_param("TagSpecification")
        if tags:
            tags = tags[0].get("Tag")
        acceptance_required = (
            str(self._get_param("AcceptanceRequired", "true")).lower() == "true"
        )
        private_dns_name = self._get_param("PrivateDnsName")

        config = self.ec2_backend.create_vpc_endpoint_service_configuration(
            gateway_lbs or network_lbs,
            acceptance_required=acceptance_required,
            private_dns_name=private_dns_name,
            tags=tags,
        )
        template = self.response_template(CREATE_VPC_ENDPOINT_SERVICE_CONFIGURATION)
        return template.render(config=config)

    def describe_vpc_endpoint_service_configurations(self) -> str:
        service_ids = self._get_multi_param("ServiceId")

        configs = self.ec2_backend.describe_vpc_endpoint_service_configurations(
            service_ids
        )

        template = self.response_template(DESCRIBE_VPC_ENDPOINT_SERVICE_CONFIGURATION)
        return template.render(configs=configs)

    def delete_vpc_endpoint_service_configurations(self) -> str:
        service_ids = self._get_multi_param("ServiceId")
        missing_configs = self.ec2_backend.delete_vpc_endpoint_service_configurations(
            service_ids
        )

        template = self.response_template(DELETE_VPC_ENDPOINT_SERVICE_CONFIGURATION)
        return template.render(missing=missing_configs)

    def describe_vpc_endpoint_service_permissions(self) -> str:
        service_id = self._get_param("ServiceId")

        principals = self.ec2_backend.describe_vpc_endpoint_service_permissions(
            service_id
        )

        template = self.response_template(DESCRIBE_VPC_ENDPOINT_SERVICE_PERMISSIONS)
        return template.render(principals=principals)

    def modify_vpc_endpoint_service_configuration(self) -> str:
        service_id = self._get_param("ServiceId")
        private_dns_name = self._get_param("PrivateDnsName")
        acceptance_required = self._get_param("AcceptanceRequired")
        add_network_lbs = self._get_multi_param("AddNetworkLoadBalancerArn")
        remove_network_lbs = self._get_multi_param("RemoveNetworkLoadBalancerArn")
        add_gateway_lbs = self._get_multi_param("AddGatewayLoadBalancerArn")
        remove_gateway_lbs = self._get_multi_param("RemoveGatewayLoadBalancerArn")

        self.ec2_backend.modify_vpc_endpoint_service_configuration(
            service_id,
            acceptance_required=acceptance_required,
            private_dns_name=private_dns_name,
            add_network_lbs=add_network_lbs,
            remove_network_lbs=remove_network_lbs,
            add_gateway_lbs=add_gateway_lbs,
            remove_gateway_lbs=remove_gateway_lbs,
        )

        return MODIFY_VPC_ENDPOINT_SERVICE_CONFIGURATION

    def modify_vpc_endpoint_service_permissions(self) -> str:
        service_id = self._get_param("ServiceId")
        add_principals = self._get_multi_param("AddAllowedPrincipals")
        remove_principals = self._get_multi_param("RemoveAllowedPrincipals")

        self.ec2_backend.modify_vpc_endpoint_service_permissions(
            service_id, add_principals, remove_principals
        )

        return MODIFY_VPC_ENDPOINT_SERVICE_PERMISSIONS


CREATE_VPC_ENDPOINT_SERVICE_CONFIGURATION = """
<CreateVpcEndpointServiceConfigurationResult xmlns="http://ec2.amazonaws.com/doc/2013-10-15/">
  <requestId>7a62c49f-347e-4fc4-9331-6e8eEXAMPLE</requestId>
  <serviceConfiguration>
      <serviceType>
        <item><serviceType>{{ config.service_type }}</serviceType></item>
      </serviceType>
      <serviceId>{{ config.id }}</serviceId>
      <serviceName>{{ config.service_name }}</serviceName>
      <serviceState>{{ config.service_state }}</serviceState>
      <availabilityZoneSet>
        {% for zone in config.availability_zones %}<item>{{ zone }}</item>{% endfor %}
      </availabilityZoneSet>
      <acceptanceRequired>{{ 'true' if config.acceptance_required else 'false' }}</acceptanceRequired>
      <managesVpcEndpoints>{{ 'true' if config.manages_vpc_endpoints else 'false' }}</managesVpcEndpoints>
      {%- if config.network_load_balancer_arns %}
      <networkLoadBalancerArnSet>
        {% for lb in config.network_load_balancer_arns %}<item>{{ lb }}</item>{% endfor %}
      </networkLoadBalancerArnSet>
      {% endif -%}
      {%- if config.gateway_load_balancer_arns %}
      <gatewayLoadBalancerArnSet>
        {% for lb in config.gateway_load_balancer_arns %}<item>{{ lb }}</item>{% endfor %}
      </gatewayLoadBalancerArnSet>
      {% endif -%}
      <baseEndpointDnsNameSet><item>{{ config.endpoint_dns_name }}</item></baseEndpointDnsNameSet>
      <privateDnsName>{{ config.private_dns_name }}</privateDnsName>
      <privateDnsNameConfiguration>
      {% if config.private_dns_name %}
        <state>verified</state>
        <type>TXT</type>
        <value>val</value>
        <name>n</name>
      {% endif %}
      </privateDnsNameConfiguration>
  </serviceConfiguration>
</CreateVpcEndpointServiceConfigurationResult>
"""


DESCRIBE_VPC_ENDPOINT_SERVICE_CONFIGURATION = """
<DescribeVpcEndpointServiceConfigurationsResult>
  <serviceConfigurationSet>
    {% for config in configs %}
      <item>
          <serviceType>
            <item><serviceType>{{ config.service_type }}</serviceType></item>
          </serviceType>
          <serviceId>{{ config.id }}</serviceId>
          <serviceName>{{ config.service_name }}</serviceName>
          <serviceState>{{ config.service_state }}</serviceState>
          <availabilityZoneSet>
            {% for zone in config.availability_zones %}<item>{{ zone }}</item>{% endfor %}
          </availabilityZoneSet>
          <acceptanceRequired>{{ 'true' if config.acceptance_required else 'false' }}</acceptanceRequired>
          <managesVpcEndpoints>{{ 'true' if config.manages_vpc_endpoints else 'false' }}</managesVpcEndpoints>
          {%- if config.network_load_balancer_arns %}
          <networkLoadBalancerArnSet>
            {% for lb in config.network_load_balancer_arns %}<item>{{ lb }}</item>{% endfor %}
          </networkLoadBalancerArnSet>
          {% endif -%}
          {%- if config.gateway_load_balancer_arns %}
          <gatewayLoadBalancerArnSet>
            {% for lb in config.gateway_load_balancer_arns %}<item>{{ lb }}</item>{% endfor %}
          </gatewayLoadBalancerArnSet>
          {% endif -%}
          <baseEndpointDnsNameSet><item>{{ config.endpoint_dns_name }}</item></baseEndpointDnsNameSet>
          <privateDnsName>{{ config.private_dns_name }}</privateDnsName>
          <privateDnsNameConfiguration>
          {% if config.private_dns_name %}
            <state>verified</state>
            <type>TXT</type>
            <value>val</value>
            <name>n</name>
          {% endif %}
          </privateDnsNameConfiguration>
          <tagSet>
                {% for tag in config.get_tags() %}
                    <item>
                        <key>{{ tag.key }}</key>
                        <value>{{ tag.value }}</value>
                    </item>
                {% endfor %}
            </tagSet>
      </item>
    {% endfor %}
  </serviceConfigurationSet>
</DescribeVpcEndpointServiceConfigurationsResult>
"""


DELETE_VPC_ENDPOINT_SERVICE_CONFIGURATION = """
<DeleteVpcEndpointServiceConfigurationsResult>
  <unsuccessful>
    {% for m in missing %}
    <item>
      <error>
        <code>InvalidVpcEndpointService.NotFound</code>
        <message>The VpcEndpointService Id '{{ m }}' does not exist</message>
      </error>
      <resourceId>{{ m }}</resourceId>
    </item>
    {% endfor %}
  </unsuccessful>
</DeleteVpcEndpointServiceConfigurationsResult>
"""


DESCRIBE_VPC_ENDPOINT_SERVICE_PERMISSIONS = """
<DescribeVpcEndpointServicePermissionsResult>
  <allowedPrincipals>
    {% for principal in principals %}
      <item>
        <principal>{{ principal }}</principal>
      </item>
    {% endfor %}
  </allowedPrincipals>
</DescribeVpcEndpointServicePermissionsResult>
"""

MODIFY_VPC_ENDPOINT_SERVICE_PERMISSIONS = """
<ModifyVpcEndpointServicePermissionsResult>
<return>true</return>
</ModifyVpcEndpointServicePermissionsResult>
"""


MODIFY_VPC_ENDPOINT_SERVICE_CONFIGURATION = """
<ModifyVpcEndpointServiceConfigurationResult>
<return>true</return>
</ModifyVpcEndpointServiceConfigurationResult>
"""
