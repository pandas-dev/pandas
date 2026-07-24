from moto.core.responses import ActionResult

from ._base_response import EC2BaseResponse


class VPCPeeringConnections(EC2BaseResponse):
    def create_vpc_peering_connection(self) -> ActionResult:
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
        # Use a response transformer to return initial connection status.
        self.RESPONSE_KEY_PATH_TO_TRANSFORMER[
            "CreateVpcPeeringConnectionResult.VpcPeeringConnection.Status"
        ] = lambda _: {
            "Code": "initiating-request",
            "Message": f"Initiating Request to {vpc_pcx.peer_vpc.owner_id}",
        }
        return ActionResult({"VpcPeeringConnection": vpc_pcx})

    def delete_vpc_peering_connection(self) -> ActionResult:
        vpc_pcx_id = self._get_param("VpcPeeringConnectionId")
        self.ec2_backend.delete_vpc_peering_connection(vpc_pcx_id)
        return ActionResult({"Return": True})

    def describe_vpc_peering_connections(self) -> ActionResult:
        ids = self._get_param("VpcPeeringConnectionIds", [])
        vpc_pcxs = self.ec2_backend.describe_vpc_peering_connections(
            vpc_peering_ids=ids
        )
        return ActionResult({"VpcPeeringConnections": vpc_pcxs})

    def accept_vpc_peering_connection(self) -> ActionResult:
        vpc_pcx_id = self._get_param("VpcPeeringConnectionId")
        vpc_pcx = self.ec2_backend.accept_vpc_peering_connection(vpc_pcx_id)
        return ActionResult({"VpcPeeringConnection": vpc_pcx})

    def reject_vpc_peering_connection(self) -> ActionResult:
        vpc_pcx_id = self._get_param("VpcPeeringConnectionId")
        self.ec2_backend.reject_vpc_peering_connection(vpc_pcx_id)
        return ActionResult({"Return": True})

    def modify_vpc_peering_connection_options(self) -> ActionResult:
        vpc_pcx_id = self._get_param("VpcPeeringConnectionId")
        accepter_options = self._get_param("AccepterPeeringConnectionOptions", {})
        requester_options = self._get_param("RequesterPeeringConnectionOptions", {})
        vpc_pcx = self.ec2_backend.modify_vpc_peering_connection_options(
            vpc_pcx_id, accepter_options, requester_options
        )
        return ActionResult(
            {
                "AccepterPeeringConnectionOptions": vpc_pcx.accepter_options,
                "RequesterPeeringConnectionOptions": vpc_pcx.requester_options,
            }
        )
