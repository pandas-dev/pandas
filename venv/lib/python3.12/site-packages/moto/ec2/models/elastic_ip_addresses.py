from typing import Any, Dict, List, Optional

from moto.core.common_models import CloudFormationModel

from ..exceptions import (
    FilterNotImplementedError,
    InvalidAddressError,
    InvalidAllocationIdError,
    InvalidAssociationIdError,
    ResourceAlreadyAssociatedError,
)
from ..utils import (
    generic_filter,
    random_eip_allocation_id,
    random_eip_association_id,
    random_ip,
)
from .core import TaggedEC2Resource


class ElasticAddress(TaggedEC2Resource, CloudFormationModel):
    def __init__(
        self,
        ec2_backend: Any,
        domain: str,
        address: Optional[str] = None,
        tags: Optional[Dict[str, str]] = None,
    ):
        self.ec2_backend = ec2_backend
        if address:
            self.public_ip = address
        else:
            self.public_ip = random_ip()
        self.allocation_id = random_eip_allocation_id() if domain == "vpc" else None
        self.id = self.allocation_id
        self.domain = domain
        self.instance = None
        self.eni = None
        self.association_id: Optional[str] = None
        self.add_tags(tags or {})

    @staticmethod
    def cloudformation_name_type() -> str:
        return ""

    @staticmethod
    def cloudformation_type() -> str:
        # https://docs.aws.amazon.com/AWSCloudFormation/latest/UserGuide/aws-properties-ec2-eip.html
        return "AWS::EC2::EIP"

    @classmethod
    def create_from_cloudformation_json(  # type: ignore[misc]
        cls,
        resource_name: str,
        cloudformation_json: Any,
        account_id: str,
        region_name: str,
        **kwargs: Any,
    ) -> "ElasticAddress":
        from ..models import ec2_backends

        ec2_backend = ec2_backends[account_id][region_name]

        properties = cloudformation_json.get("Properties")
        instance_id = None
        if properties:
            domain = properties.get("Domain")
            # TODO: support tags from cloudformation template
            eip = ec2_backend.allocate_address(domain=domain if domain else "standard")
            instance_id = properties.get("InstanceId")
        else:
            eip = ec2_backend.allocate_address(domain="standard")

        if instance_id:
            instance = ec2_backend.get_instance_by_id(instance_id)
            ec2_backend.associate_address(instance, address=eip.public_ip)

        return eip

    @property
    def physical_resource_id(self) -> str:
        return self.public_ip

    @classmethod
    def has_cfn_attr(cls, attr: str) -> bool:
        return attr in ["AllocationId"]

    def get_cfn_attribute(self, attribute_name: str) -> Any:
        from moto.cloudformation.exceptions import UnformattedGetAttTemplateException

        if attribute_name == "AllocationId":
            return self.allocation_id
        raise UnformattedGetAttTemplateException()

    def get_filter_value(
        self, filter_name: str, method_name: Optional[str] = None
    ) -> Any:
        if filter_name == "allocation-id":
            return self.allocation_id
        elif filter_name == "association-id":
            return self.association_id
        elif filter_name == "domain":
            return self.domain
        elif filter_name == "instance-id":
            if self.instance:
                return self.instance.id
            return None
        elif filter_name == "network-interface-id":
            if self.eni:
                return self.eni.id
            return None
        elif filter_name == "private-ip-address":
            if self.eni:
                return self.eni.private_ip_address
            return None
        elif filter_name == "public-ip":
            return self.public_ip
        elif filter_name == "network-interface-owner-id":
            # TODO: implement network-interface-owner-id
            raise FilterNotImplementedError(filter_name, "DescribeAddresses")
        else:
            return super().get_filter_value(filter_name, "DescribeAddresses")


class ElasticAddressBackend:
    def __init__(self) -> None:
        self.addresses: List[ElasticAddress] = []

    def allocate_address(
        self,
        domain: str,
        address: Optional[str] = None,
        tags: Optional[Dict[str, str]] = None,
    ) -> ElasticAddress:
        if domain not in ["standard", "vpc"]:
            domain = "vpc"
        if address:
            ea = ElasticAddress(self, domain=domain, address=address, tags=tags)
        else:
            ea = ElasticAddress(self, domain=domain, tags=tags)
        self.addresses.append(ea)
        return ea

    def address_by_ip(
        self, ips: List[str], fail_if_not_found: bool = True
    ) -> List[ElasticAddress]:
        eips = [
            address for address in self.addresses.copy() if address.public_ip in ips
        ]

        # TODO: Trim error message down to specific invalid address.
        if (not eips or len(ips) > len(eips)) and fail_if_not_found:
            raise InvalidAddressError(ips)

        return eips

    def address_by_allocation(self, allocation_ids: List[str]) -> List[ElasticAddress]:
        eips = [
            address
            for address in self.addresses
            if address.allocation_id in allocation_ids
        ]

        # TODO: Trim error message down to specific invalid id.
        if not eips or len(allocation_ids) > len(eips):
            raise InvalidAllocationIdError(allocation_ids)

        return eips

    def address_by_association(
        self, association_ids: List[str]
    ) -> List[ElasticAddress]:
        eips = [
            address
            for address in self.addresses
            if address.association_id in association_ids
        ]

        # TODO: Trim error message down to specific invalid id.
        if not eips or len(association_ids) > len(eips):
            raise InvalidAssociationIdError(association_ids)

        return eips

    def associate_address(
        self,
        instance: Any = None,
        eni: Any = None,
        address: Optional[str] = None,
        allocation_id: Optional[str] = None,
        reassociate: bool = False,
    ) -> ElasticAddress:
        eips = []
        if address:
            eips = self.address_by_ip([address])
        elif allocation_id:
            eips = self.address_by_allocation([allocation_id])
        eip = eips[0]

        new_instance_association = bool(
            instance and (not eip.instance or eip.instance.id == instance.id)
        )
        new_eni_association = bool(eni and (not eip.eni or eni.id == eip.eni.id))

        if new_instance_association or new_eni_association or reassociate:
            eip.instance = instance
            eip.eni = eni
            if not eip.eni and instance:
                # default to primary network interface
                eip.eni = instance.nics[0]
            if eip.eni:
                eip.eni.public_ip = eip.public_ip
            if eip.domain == "vpc":
                eip.association_id = random_eip_association_id()

            return eip

        raise ResourceAlreadyAssociatedError(eip.public_ip)

    def describe_addresses(
        self,
        allocation_ids: Optional[List[str]] = None,
        public_ips: Optional[List[str]] = None,
        filters: Any = None,
    ) -> List[ElasticAddress]:
        matches = self.addresses.copy()
        if allocation_ids:
            matches = [addr for addr in matches if addr.allocation_id in allocation_ids]
            if len(allocation_ids) > len(matches):
                unknown_ids = set(allocation_ids) - set(matches)  # type: ignore[arg-type]
                raise InvalidAllocationIdError(unknown_ids)
        if public_ips:
            matches = [addr for addr in matches if addr.public_ip in public_ips]
            if len(public_ips) > len(matches):
                unknown_ips = set(public_ips) - set(matches)  # type: ignore[arg-type]
                raise InvalidAddressError(unknown_ips)
        if filters:
            matches = generic_filter(filters, matches)

        return matches

    def describe_addresses_attribute(
        self, allocation_ids: Optional[List[str]] = None
    ) -> List[ElasticAddress]:
        return self.describe_addresses(allocation_ids)

    def disassociate_address(
        self, address: Optional[str] = None, association_id: Optional[str] = None
    ) -> None:
        eips = []
        if address:
            eips = self.address_by_ip([address])
        elif association_id:
            eips = self.address_by_association([association_id])
        eip = eips[0]

        if eip.eni:
            eip.eni.public_ip = None
            if eip.eni.instance and eip.eni.instance._state.name == "running":
                eip.eni.check_auto_public_ip()
            eip.eni = None

        eip.instance = None
        eip.association_id = None

    def release_address(
        self, address: Optional[str] = None, allocation_id: Optional[str] = None
    ) -> None:
        eips = []
        if address:
            eips = self.address_by_ip([address])
        elif allocation_id:
            eips = self.address_by_allocation([allocation_id])
        eip = eips[0]

        self.disassociate_address(address=eip.public_ip)
        eip.allocation_id = None
        self.addresses.remove(eip)
