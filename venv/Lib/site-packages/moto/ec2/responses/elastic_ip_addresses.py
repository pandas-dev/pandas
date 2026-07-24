from moto.core.responses import ActionResult, EmptyResult
from moto.ec2.utils import add_tag_specification

from ._base_response import EC2BaseResponse


class ElasticIPAddresses(EC2BaseResponse):
    def allocate_address(self) -> ActionResult:
        domain = self._get_param("Domain", if_none=None)
        reallocate_address = self._get_param("Address", if_none=None)
        tag_param = self._get_param("TagSpecifications", [])
        tags = add_tag_specification(tag_param)

        self.error_on_dryrun()

        if reallocate_address:
            address = self.ec2_backend.allocate_address(
                domain, address=reallocate_address, tags=tags
            )
        else:
            address = self.ec2_backend.allocate_address(domain, tags=tags)
        return ActionResult(
            {
                "PublicIp": address.public_ip,
                "Domain": address.domain,
                "AllocationId": address.allocation_id,
            }
        )

    def associate_address(self) -> ActionResult:
        instance = eni = None

        if "InstanceId" in self._get_params():
            instance = self.ec2_backend.get_instance(self._get_param("InstanceId"))
        elif "NetworkInterfaceId" in self._get_params():
            eni = self.ec2_backend.get_network_interface(
                self._get_param("NetworkInterfaceId")
            )
        else:
            self.ec2_backend.raise_error(
                "MissingParameter",
                "Invalid request, expect InstanceId/NetworkId parameter.",
            )

        reassociate = self._get_param("AllowReassociation", False)

        self.error_on_dryrun()

        if instance or eni:
            if "PublicIp" in self._get_params():
                eip = self.ec2_backend.associate_address(
                    instance=instance,
                    eni=eni,
                    address=self._get_param("PublicIp"),
                    reassociate=reassociate,
                )
            elif "AllocationId" in self._get_params():
                eip = self.ec2_backend.associate_address(
                    instance=instance,
                    eni=eni,
                    allocation_id=self._get_param("AllocationId"),
                    reassociate=reassociate,
                )
            else:
                self.ec2_backend.raise_error(
                    "MissingParameter",
                    "Invalid request, expect PublicIp/AllocationId parameter.",
                )
        else:
            self.ec2_backend.raise_error(
                "MissingParameter",
                "Invalid request, expect either instance or ENI.",
            )

        return ActionResult({"AssociationId": eip.association_id})

    def describe_addresses(self) -> ActionResult:
        self.error_on_dryrun()
        allocation_ids = self._get_param("AllocationIds", [])
        public_ips = self._get_param("PublicIps", [])
        filters = self._filters_from_querystring()
        addresses = self.ec2_backend.describe_addresses(
            allocation_ids, public_ips, filters
        )
        return ActionResult({"Addresses": addresses})

    def describe_addresses_attribute(self) -> ActionResult:
        self.error_on_dryrun()
        allocation_ids = self._get_param("AllocationIds", [])
        addresses = self.ec2_backend.describe_addresses_attribute(allocation_ids)
        return ActionResult({"Addresses": addresses})

    def disassociate_address(self) -> ActionResult:
        if (
            "PublicIp" not in self._get_params()
            and "AssociationId" not in self._get_params()
        ):
            self.ec2_backend.raise_error(
                "MissingParameter",
                "Invalid request, expect PublicIp/AssociationId parameter.",
            )

        self.error_on_dryrun()

        if "PublicIp" in self._get_params():
            self.ec2_backend.disassociate_address(address=self._get_param("PublicIp"))
        elif "AssociationId" in self._get_params():
            self.ec2_backend.disassociate_address(
                association_id=self._get_param("AssociationId")
            )

        return EmptyResult()

    def release_address(self) -> ActionResult:
        if (
            "PublicIp" not in self._get_params()
            and "AllocationId" not in self._get_params()
        ):
            self.ec2_backend.raise_error(
                "MissingParameter",
                "Invalid request, expect PublicIp/AllocationId parameter.",
            )

        self.error_on_dryrun()

        if "PublicIp" in self._get_params():
            self.ec2_backend.release_address(address=self._get_param("PublicIp"))
        elif "AllocationId" in self._get_params():
            self.ec2_backend.release_address(
                allocation_id=self._get_param("AllocationId")
            )

        return EmptyResult()
