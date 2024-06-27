from ._base_response import EC2BaseResponse


class DHCPOptions(EC2BaseResponse):
    def associate_dhcp_options(self) -> str:
        dhcp_opt_id = self._get_param("DhcpOptionsId")
        vpc_id = self._get_param("VpcId")

        vpc = self.ec2_backend.get_vpc(vpc_id)

        if dhcp_opt_id == "default":
            self.ec2_backend.disassociate_dhcp_options(vpc)
        else:
            dhcp_opt = self.ec2_backend.describe_dhcp_options([dhcp_opt_id])[0]
            self.ec2_backend.associate_dhcp_options(dhcp_opt, vpc)

        template = self.response_template(ASSOCIATE_DHCP_OPTIONS_RESPONSE)
        return template.render()

    def create_dhcp_options(self) -> str:
        provided_config = self._get_multi_param("DhcpConfiguration")
        flat_config = {f["Key"]: f["Value"] for f in provided_config}

        # TODO validate we only got the options we know about

        domain_name_servers = flat_config.get("domain-name-servers", None)
        domain_name = flat_config.get("domain-name", None)
        ntp_servers = flat_config.get("ntp-servers", None)
        netbios_name_servers = flat_config.get("netbios-name-servers", None)
        netbios_node_type = flat_config.get("netbios-node-type", None)

        dhcp_options_set = self.ec2_backend.create_dhcp_options(
            domain_name_servers=domain_name_servers,
            domain_name=domain_name,
            ntp_servers=ntp_servers,
            netbios_name_servers=netbios_name_servers,
            netbios_node_type=netbios_node_type,
        )

        template = self.response_template(CREATE_DHCP_OPTIONS_RESPONSE)
        return template.render(dhcp_options_set=dhcp_options_set)

    def delete_dhcp_options(self) -> str:
        dhcp_opt_id = self._get_param("DhcpOptionsId")
        self.ec2_backend.delete_dhcp_options_set(dhcp_opt_id)
        template = self.response_template(DELETE_DHCP_OPTIONS_RESPONSE)
        return template.render(delete_status="true")

    def describe_dhcp_options(self) -> str:
        dhcp_opt_ids = self._get_multi_param("DhcpOptionsId")
        filters = self._filters_from_querystring()
        dhcp_opts = self.ec2_backend.describe_dhcp_options(dhcp_opt_ids, filters)
        template = self.response_template(DESCRIBE_DHCP_OPTIONS_RESPONSE)
        return template.render(dhcp_options=dhcp_opts)


CREATE_DHCP_OPTIONS_RESPONSE = """
<CreateDhcpOptionsResponse xmlns="http://ec2.amazonaws.com/doc/2013-10-15/">
  <requestId>7a62c49f-347e-4fc4-9331-6e8eEXAMPLE</requestId>
  <dhcpOptions>
      <dhcpOptionsId>{{ dhcp_options_set.id }}</dhcpOptionsId>
      <dhcpConfigurationSet>
      {% for key, values in dhcp_options_set.options.items() %}
        {{ values }}
        {% if values %}
        <item>
          <key>{{key}}</key>
          <valueSet>
            {% for value in values %}
            <item>
              <value>{{ value }}</value>
            </item>
            {% endfor %}
          </valueSet>
        </item>
        {% endif %}
      {% endfor %}
      </dhcpConfigurationSet>
      <tagSet>
        {% for tag in dhcp_options_set.get_tags() %}
          <item>
            <resourceId>{{ tag.resource_id }}</resourceId>
            <resourceType>{{ tag.resource_type }}</resourceType>
            <key>{{ tag.key }}</key>
            <value>{{ tag.value }}</value>
          </item>
        {% endfor %}
      </tagSet>
  </dhcpOptions>
</CreateDhcpOptionsResponse>
"""

DELETE_DHCP_OPTIONS_RESPONSE = """
<DeleteDhcpOptionsResponse xmlns="http://ec2.amazonaws.com/doc/2013-10-15/">
  <requestId>7a62c49f-347e-4fc4-9331-6e8eEXAMPLE</requestId>
  <return>{{delete_status}}</return>
</DeleteDhcpOptionsResponse>
"""

DESCRIBE_DHCP_OPTIONS_RESPONSE = """
<DescribeDhcpOptionsResponse xmlns="http://ec2.amazonaws.com/doc/2013-10-15/">
  <requestId>7a62c49f-347e-4fc4-9331-6e8eEXAMPLE</requestId>
  <dhcpOptionsSet>
  {% for dhcp_options_set in dhcp_options %}
    <item>
      <dhcpOptionsId>{{ dhcp_options_set.id }}</dhcpOptionsId>
      <dhcpConfigurationSet>
      {% for key, values in dhcp_options_set.options.items() %}
        {{ values }}
        {% if values %}
        <item>
          <key>{{ key }}</key>
          <valueSet>
            {% for value in values %}
            <item>
              <value>{{ value }}</value>
            </item>
            {% endfor %}
          </valueSet>
        </item>
        {% endif %}
      {% endfor %}
      </dhcpConfigurationSet>
      <tagSet>
        {% for tag in dhcp_options_set.get_tags() %}
          <item>
            <resourceId>{{ tag.resource_id }}</resourceId>
            <resourceType>{{ tag.resource_type }}</resourceType>
            <key>{{ tag.key }}</key>
            <value>{{ tag.value }}</value>
          </item>
        {% endfor %}
      </tagSet>
    </item>
  {% endfor %}
  </dhcpOptionsSet>
</DescribeDhcpOptionsResponse>
"""

ASSOCIATE_DHCP_OPTIONS_RESPONSE = """
<AssociateDhcpOptionsResponse xmlns="http://ec2.amazonaws.com/doc/2013-10-15/">
<requestId>7a62c49f-347e-4fc4-9331-6e8eEXAMPLE</requestId>
<return>true</return>
</AssociateDhcpOptionsResponse>
"""
