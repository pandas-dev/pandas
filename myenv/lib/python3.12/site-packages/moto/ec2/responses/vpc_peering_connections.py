from ._base_response import EC2BaseResponse


class VPCPeeringConnections(EC2BaseResponse):
    def create_vpc_peering_connection(self) -> str:
        tags = self._parse_tag_specification().get("vpc-peering-connection", {})

        account_id = self._get_param("PeerOwnerId") or self.current_account
        region_name = self._get_param("PeerRegion") or self.region

        vpc = self.ec2_backend.get_vpc(self._get_param("VpcId"))

        # Peer VPC could belong to another account or region
        from moto.ec2.models import ec2_backends

        peer_vpc = ec2_backends[account_id][region_name].get_vpc(
            self._get_param("PeerVpcId")
        )

        vpc_pcx = self.ec2_backend.create_vpc_peering_connection(vpc, peer_vpc, tags)
        template = self.response_template(CREATE_VPC_PEERING_CONNECTION_RESPONSE)
        return template.render(vpc_pcx=vpc_pcx)

    def delete_vpc_peering_connection(self) -> str:
        vpc_pcx_id = self._get_param("VpcPeeringConnectionId")
        vpc_pcx = self.ec2_backend.delete_vpc_peering_connection(vpc_pcx_id)
        template = self.response_template(DELETE_VPC_PEERING_CONNECTION_RESPONSE)
        return template.render(vpc_pcx=vpc_pcx)

    def describe_vpc_peering_connections(self) -> str:
        ids = self._get_multi_param("VpcPeeringConnectionId")
        vpc_pcxs = self.ec2_backend.describe_vpc_peering_connections(
            vpc_peering_ids=ids
        )
        template = self.response_template(DESCRIBE_VPC_PEERING_CONNECTIONS_RESPONSE)
        return template.render(vpc_pcxs=vpc_pcxs)

    def accept_vpc_peering_connection(self) -> str:
        vpc_pcx_id = self._get_param("VpcPeeringConnectionId")
        vpc_pcx = self.ec2_backend.accept_vpc_peering_connection(vpc_pcx_id)
        template = self.response_template(ACCEPT_VPC_PEERING_CONNECTION_RESPONSE)
        return template.render(vpc_pcx=vpc_pcx)

    def reject_vpc_peering_connection(self) -> str:
        vpc_pcx_id = self._get_param("VpcPeeringConnectionId")
        self.ec2_backend.reject_vpc_peering_connection(vpc_pcx_id)
        template = self.response_template(REJECT_VPC_PEERING_CONNECTION_RESPONSE)
        return template.render()

    def modify_vpc_peering_connection_options(self) -> str:
        vpc_pcx_id = self._get_param("VpcPeeringConnectionId")
        accepter_options = self._get_multi_param_dict(
            "AccepterPeeringConnectionOptions"
        )
        requester_options = self._get_multi_param_dict(
            "RequesterPeeringConnectionOptions"
        )
        self.ec2_backend.modify_vpc_peering_connection_options(
            vpc_pcx_id, accepter_options, requester_options
        )
        template = self.response_template(MODIFY_VPC_PEERING_CONNECTION_RESPONSE)
        return template.render(
            accepter_options=accepter_options, requester_options=requester_options
        )


# we are assuming that the owner id for accepter and requester vpc are same
# as we are checking for the vpc exsistance
CREATE_VPC_PEERING_CONNECTION_RESPONSE = """
<CreateVpcPeeringConnectionResponse xmlns="http://ec2.amazonaws.com/doc/2016-11-15/">
 <requestId>7a62c49f-347e-4fc4-9331-6e8eEXAMPLE</requestId>
 <vpcPeeringConnection>
   <vpcPeeringConnectionId>{{ vpc_pcx.id }}</vpcPeeringConnectionId>
   <requesterVpcInfo>
     <ownerId>{{ vpc_pcx.vpc.owner_id }}</ownerId>
     <region>{{ vpc_pcx.vpc.region }}</region>
     <vpcId>{{ vpc_pcx.vpc.id }}</vpcId>
     <cidrBlock>{{ vpc_pcx.vpc.cidr_block }}</cidrBlock>
     <cidrBlockSet></cidrBlockSet>
     <ipv6CidrBlockSet></ipv6CidrBlockSet>
     <peeringOptions>
       <allowEgressFromLocalClassicLinkToRemoteVpc>{{ vpc_pcx.requester_options.AllowEgressFromLocalClassicLinkToRemoteVpc or '' }}</allowEgressFromLocalClassicLinkToRemoteVpc>
       <allowEgressFromLocalVpcToRemoteClassicLink>{{ vpc_pcx.requester_options.AllowEgressFromLocalVpcToRemoteClassicLink or '' }}</allowEgressFromLocalVpcToRemoteClassicLink>
       <allowDnsResolutionFromRemoteVpc>{{ vpc_pcx.requester_options.AllowDnsResolutionFromRemoteVpc or '' }}</allowDnsResolutionFromRemoteVpc>
     </peeringOptions>
   </requesterVpcInfo>
   <accepterVpcInfo>
     <ownerId>{{ vpc_pcx.peer_vpc.owner_id }}</ownerId>
     <region>{{ vpc_pcx.peer_vpc.region }}</region>
     <vpcId>{{ vpc_pcx.peer_vpc.id }}</vpcId>
     <cidrBlock>{{ vpc_pcx.peer_vpc.cidr_block }}</cidrBlock>
     <cidrBlockSet></cidrBlockSet>
     <ipv6CidrBlockSet></ipv6CidrBlockSet>
     <peeringOptions>
       <allowEgressFromLocalClassicLinkToRemoteVpc>{{ vpc_pcx.accepter_options.AllowEgressFromLocalClassicLinkToRemoteVpc or '' }}</allowEgressFromLocalClassicLinkToRemoteVpc>
       <allowEgressFromLocalVpcToRemoteClassicLink>{{ vpc_pcx.accepter_options.AllowEgressFromLocalVpcToRemoteClassicLink or '' }}</allowEgressFromLocalVpcToRemoteClassicLink>
       <allowDnsResolutionFromRemoteVpc>{{ vpc_pcx.accepter_options.AllowDnsResolutionFromRemoteVpc or '' }}</allowDnsResolutionFromRemoteVpc>
     </peeringOptions>
   </accepterVpcInfo>
   <status>
     <code>initiating-request</code>
     <message>Initiating Request to {{ vpc_pcx.peer_vpc.owner_id }}</message>
   </status>
   <expirationTime>2014-02-18T14:37:25.000Z</expirationTime>
   <tagSet>
   {% for tag in vpc_pcx.get_tags() %}
     <item>
       <key>{{ tag.key }}</key>
       <value>{{ tag.value }}</value>
     </item>
   {% endfor %}
   </tagSet>
 </vpcPeeringConnection>
</CreateVpcPeeringConnectionResponse>
"""

DESCRIBE_VPC_PEERING_CONNECTIONS_RESPONSE = """
<DescribeVpcPeeringConnectionsResponse xmlns="http://ec2.amazonaws.com/doc/2016-11-15/">
<requestId>7a62c49f-347e-4fc4-9331-6e8eEXAMPLE</requestId>
 <vpcPeeringConnectionSet>
 {% for vpc_pcx in vpc_pcxs %}
 <item>
  <vpcPeeringConnectionId>{{ vpc_pcx.id }}</vpcPeeringConnectionId>
    <requesterVpcInfo>
     <ownerId>{{ vpc_pcx.vpc.owner_id }}</ownerId>
     <region>{{ vpc_pcx.vpc.region }}</region>
     <vpcId>{{ vpc_pcx.vpc.id }}</vpcId>
     <cidrBlock>{{ vpc_pcx.vpc.cidr_block }}</cidrBlock>
     <cidrBlockSet></cidrBlockSet>
     <ipv6CidrBlockSet></ipv6CidrBlockSet>
     <peeringOptions>
      <allowEgressFromLocalClassicLinkToRemoteVpc>{{ vpc_pcx.requester_options.AllowEgressFromLocalClassicLinkToRemoteVpc or '' }}</allowEgressFromLocalClassicLinkToRemoteVpc>
      <allowEgressFromLocalVpcToRemoteClassicLink>{{ vpc_pcx.requester_options.AllowEgressFromLocalVpcToRemoteClassicLink or '' }}</allowEgressFromLocalVpcToRemoteClassicLink>
      <allowDnsResolutionFromRemoteVpc>{{ vpc_pcx.requester_options.AllowDnsResolutionFromRemoteVpc or '' }}</allowDnsResolutionFromRemoteVpc>
     </peeringOptions>
    </requesterVpcInfo>
    <accepterVpcInfo>
     <ownerId>{{ vpc_pcx.peer_vpc.owner_id }}</ownerId>
     <region>{{ vpc_pcx.peer_vpc.region }}</region>
     <vpcId>{{ vpc_pcx.peer_vpc.id }}</vpcId>
     <cidrBlock>{{ vpc_pcx.peer_vpc.cidr_block }}</cidrBlock>
     <cidrBlockSet></cidrBlockSet>
     <ipv6CidrBlockSet></ipv6CidrBlockSet>
     <peeringOptions>
      <allowEgressFromLocalClassicLinkToRemoteVpc>{{ vpc_pcx.accepter_options.AllowEgressFromLocalClassicLinkToRemoteVpc or '' }}</allowEgressFromLocalClassicLinkToRemoteVpc>
      <allowEgressFromLocalVpcToRemoteClassicLink>{{ vpc_pcx.accepter_options.AllowEgressFromLocalVpcToRemoteClassicLink or '' }}</allowEgressFromLocalVpcToRemoteClassicLink>
      <allowDnsResolutionFromRemoteVpc>{{ vpc_pcx.accepter_options.AllowDnsResolutionFromRemoteVpc or '' }}</allowDnsResolutionFromRemoteVpc>
     </peeringOptions>
    </accepterVpcInfo>
     <status>
      <code>{{ vpc_pcx._status.code }}</code>
      <message>{{ vpc_pcx._status.message }}</message>
     </status>
     <tagSet>
     {% for tag in vpc_pcx.get_tags() %}
       <item>
         <key>{{ tag.key }}</key>
         <value>{{ tag.value }}</value>
       </item>
     {% endfor %}
     </tagSet>
 </item>
 {% endfor %}
 </vpcPeeringConnectionSet>
</DescribeVpcPeeringConnectionsResponse>
"""

DELETE_VPC_PEERING_CONNECTION_RESPONSE = """
<DeleteVpcPeeringConnectionResponse xmlns="http://ec2.amazonaws.com/doc/2013-10-15/">
  <requestId>7a62c49f-347e-4fc4-9331-6e8eEXAMPLE</requestId>
  <return>true</return>
</DeleteVpcPeeringConnectionResponse>
"""

ACCEPT_VPC_PEERING_CONNECTION_RESPONSE = """
<AcceptVpcPeeringConnectionResponse xmlns="http://ec2.amazonaws.com/doc/2016-11-15/">
  <requestId>7a62c49f-347e-4fc4-9331-6e8eEXAMPLE</requestId>
  <vpcPeeringConnection>
    <vpcPeeringConnectionId>{{ vpc_pcx.id }}</vpcPeeringConnectionId>
    <requesterVpcInfo>
     <ownerId>{{ vpc_pcx.vpc.owner_id }}</ownerId>
     <region>{{ vpc_pcx.vpc.region }}</region>
     <vpcId>{{ vpc_pcx.vpc.id }}</vpcId>
     <cidrBlock>{{ vpc_pcx.vpc.cidr_block }}</cidrBlock>
     <cidrBlockSet></cidrBlockSet>
     <ipv6CidrBlockSet></ipv6CidrBlockSet>
     <peeringOptions>
      <allowEgressFromLocalClassicLinkToRemoteVpc>{{ vpc_pcx.requester_options.AllowEgressFromLocalClassicLinkToRemoteVpc or '' }}</allowEgressFromLocalClassicLinkToRemoteVpc>
      <allowEgressFromLocalVpcToRemoteClassicLink>{{ vpc_pcx.requester_options.AllowEgressFromLocalVpcToRemoteClassicLink or '' }}</allowEgressFromLocalVpcToRemoteClassicLink>
      <allowDnsResolutionFromRemoteVpc>{{ vpc_pcx.requester_options.AllowDnsResolutionFromRemoteVpc or '' }}</allowDnsResolutionFromRemoteVpc>
     </peeringOptions>
    </requesterVpcInfo>
    <accepterVpcInfo>
      <ownerId>{{ vpc_pcx.peer_vpc.owner_id }}</ownerId>
      <region>{{ vpc_pcx.peer_vpc.region }}</region>
      <vpcId>{{ vpc_pcx.peer_vpc.id }}</vpcId>
      <cidrBlock>{{ vpc_pcx.peer_vpc.cidr_block }}</cidrBlock>
      <cidrBlockSet></cidrBlockSet>
      <ipv6CidrBlockSet></ipv6CidrBlockSet>
      <peeringOptions>
        <allowEgressFromLocalClassicLinkToRemoteVpc>{{ vpc_pcx.accepter_options.AllowEgressFromLocalClassicLinkToRemoteVpc or '' }}</allowEgressFromLocalClassicLinkToRemoteVpc>
        <allowEgressFromLocalVpcToRemoteClassicLink>{{ vpc_pcx.accepter_options.AllowEgressFromLocalVpcToRemoteClassicLink or '' }}</allowEgressFromLocalVpcToRemoteClassicLink>
        <allowDnsResolutionFromRemoteVpc>{{ vpc_pcx.accepter_options.AllowDnsResolutionFromRemoteVpc or '' }}</allowDnsResolutionFromRemoteVpc>
      </peeringOptions>
    </accepterVpcInfo>
    <status>
      <code>{{ vpc_pcx._status.code }}</code>
      <message>{{ vpc_pcx._status.message }}</message>
    </status>
    <tagSet>
    {% for tag in vpc_pcx.get_tags() %}
      <item>
        <key>{{ tag.key }}</key>
        <value>{{ tag.value }}</value>
      </item>
    {% endfor %}
    </tagSet>
  </vpcPeeringConnection>
</AcceptVpcPeeringConnectionResponse>
"""

REJECT_VPC_PEERING_CONNECTION_RESPONSE = """
<RejectVpcPeeringConnectionResponse xmlns="http://ec2.amazonaws.com/doc/2013-10-15/">
  <requestId>7a62c49f-347e-4fc4-9331-6e8eEXAMPLE</requestId>
  <return>true</return>
</RejectVpcPeeringConnectionResponse>
"""

MODIFY_VPC_PEERING_CONNECTION_RESPONSE = """
<ModifyVpcPeeringConnectionOptionsResponse xmlns="http://ec2.amazonaws.com/doc/2016-11-15/">
  <requestId>8d977c82-8aba-4cd1-81ca-example</requestId>
  {% if requester_options %}
  <requesterPeeringConnectionOptions>
    <allowEgressFromLocalClassicLinkToRemoteVpc>{{ requester_options.AllowEgressFromLocalClassicLinkToRemoteVpc or '' }}</allowEgressFromLocalClassicLinkToRemoteVpc>
    <allowEgressFromLocalVpcToRemoteClassicLink>{{ requester_options.AllowEgressFromLocalVpcToRemoteClassicLink or '' }}</allowEgressFromLocalVpcToRemoteClassicLink>
    <allowDnsResolutionFromRemoteVpc>{{ requester_options.AllowDnsResolutionFromRemoteVpc or '' }}</allowDnsResolutionFromRemoteVpc>
  </requesterPeeringConnectionOptions>
  {% endif %}
  {% if accepter_options %}
  <accepterPeeringConnectionOptions>
    <allowEgressFromLocalClassicLinkToRemoteVpc>{{ accepter_options.AllowEgressFromLocalClassicLinkToRemoteVpc or '' }}</allowEgressFromLocalClassicLinkToRemoteVpc>
    <allowEgressFromLocalVpcToRemoteClassicLink>{{ accepter_options.AllowEgressFromLocalVpcToRemoteClassicLink or '' }}</allowEgressFromLocalVpcToRemoteClassicLink>
    <allowDnsResolutionFromRemoteVpc>{{ accepter_options.AllowDnsResolutionFromRemoteVpc or '' }}</allowDnsResolutionFromRemoteVpc>
  </accepterPeeringConnectionOptions>
  {% endif %}
</ModifyVpcPeeringConnectionOptionsResponse>
"""
