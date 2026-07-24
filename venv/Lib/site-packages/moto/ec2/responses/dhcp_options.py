from moto.core.responses import ActionResult, EmptyResult

from ._base_response import EC2BaseResponse


class DHCPOptions(EC2BaseResponse):
    def associate_dhcp_options(self) -> ActionResult:
        dhcp_opt_id = self._get_param("DhcpOptionsId")
        vpc_id = self._get_param("VpcId")

        vpc = self.ec2_backend.get_vpc(vpc_id)

        if dhcp_opt_id == "default":
            self.ec2_backend.disassociate_dhcp_options(vpc)
        else:
            dhcp_opt = self.ec2_backend.describe_dhcp_options([dhcp_opt_id])[0]
            self.ec2_backend.associate_dhcp_options(dhcp_opt, vpc)

        return EmptyResult()

    def create_dhcp_options(self) -> ActionResult:
        provided_config = self._get_param("DhcpConfigurations", [])
        flat_config = {f["Key"]: f["Values"] for f in provided_config}

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

        result = {"DhcpOptions": dhcp_options_set}
        return ActionResult(result)

    def delete_dhcp_options(self) -> ActionResult:
        dhcp_opt_id = self._get_param("DhcpOptionsId")
        self.ec2_backend.delete_dhcp_options_set(dhcp_opt_id)
        return EmptyResult()

    def describe_dhcp_options(self) -> ActionResult:
        dhcp_opt_ids = self._get_param("DhcpOptionsIds", [])
        filters = self._filters_from_querystring()
        dhcp_opts = self.ec2_backend.describe_dhcp_options(dhcp_opt_ids, filters)
        result = {"DhcpOptions": dhcp_opts}
        return ActionResult(result)
