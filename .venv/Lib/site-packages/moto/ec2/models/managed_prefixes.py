from typing import Any, Dict, List, Optional

from moto.utilities.utils import filter_resources, get_partition

from ..utils import describe_tag_filter, random_managed_prefix_list_id
from .core import TaggedEC2Resource


class ManagedPrefixList(TaggedEC2Resource):
    def __init__(
        self,
        backend: Any,
        region: str,
        address_family: Optional[str] = None,
        entry: Optional[List[Dict[str, str]]] = None,
        max_entries: Optional[str] = None,
        prefix_list_name: Optional[str] = None,
        tags: Optional[Dict[str, str]] = None,
        owner_id: Optional[str] = None,
    ):
        self.ec2_backend = backend
        self.address_family = address_family
        self.max_entries = max_entries
        self.id = random_managed_prefix_list_id()
        self.prefix_list_name = prefix_list_name
        self.state = "create-complete"
        self.state_message = "create complete"
        self.add_tags(tags or {})
        self.version: Optional[int] = 1
        self.entries = {self.version: entry} if entry else {}
        self.resource_owner_id = owner_id if owner_id else None
        self.prefix_list_arn = self.arn(region, self.owner_id)
        self.delete_counter = 1

    def arn(self, region: str, owner_id: str) -> str:
        return (
            f"arn:{get_partition(region)}:ec2:{region}:{owner_id}:prefix-list/{self.id}"
        )

    @property
    def owner_id(self) -> str:
        return self.resource_owner_id or self.ec2_backend.account_id


class ManagedPrefixListBackend:
    def __init__(self) -> None:
        self.managed_prefix_lists: Dict[str, ManagedPrefixList] = {}
        self.create_default_pls()

    def create_managed_prefix_list(
        self,
        address_family: Optional[str] = None,
        entry: Optional[List[Dict[str, str]]] = None,
        max_entries: Optional[str] = None,
        prefix_list_name: Optional[str] = None,
        tags: Optional[Dict[str, str]] = None,
        owner_id: Optional[str] = None,
    ) -> ManagedPrefixList:
        managed_prefix_list = ManagedPrefixList(
            self,
            address_family=address_family,
            entry=entry,
            max_entries=max_entries,
            prefix_list_name=prefix_list_name,
            region=self.region_name,  # type: ignore[attr-defined]
            tags=tags,
            owner_id=owner_id,
        )
        self.managed_prefix_lists[managed_prefix_list.id] = managed_prefix_list
        return managed_prefix_list

    def describe_managed_prefix_lists(
        self, prefix_list_ids: Optional[List[str]] = None, filters: Any = None
    ) -> List[ManagedPrefixList]:
        managed_prefix_lists = list(self.managed_prefix_lists.values())
        attr_pairs = (
            ("owner-id", "owner_id"),
            ("prefix-list-id", "id"),
            ("prefix-list-name", "prefix_list_name"),
        )

        if prefix_list_ids:
            managed_prefix_lists = [
                managed_prefix_list
                for managed_prefix_list in managed_prefix_lists
                if managed_prefix_list.id in prefix_list_ids
            ]

        result = managed_prefix_lists
        if filters:
            result = filter_resources(result, filters, attr_pairs)
            result = describe_tag_filter(filters, result)

        for item in result.copy():
            if not item.delete_counter:
                self.managed_prefix_lists.pop(item.id, None)
                result.remove(item)
            if item.state == "delete-complete":
                item.delete_counter -= 1
        return result

    def get_managed_prefix_list_entries(
        self, prefix_list_id: str
    ) -> Optional[ManagedPrefixList]:
        return self.managed_prefix_lists.get(prefix_list_id)

    def delete_managed_prefix_list(self, prefix_list_id: str) -> ManagedPrefixList:
        managed_prefix_list: ManagedPrefixList = self.managed_prefix_lists.get(
            prefix_list_id
        )  # type: ignore
        managed_prefix_list.state = "delete-complete"
        return managed_prefix_list

    def modify_managed_prefix_list(
        self,
        add_entry: List[Dict[str, str]],
        remove_entry: List[Dict[str, str]],
        prefix_list_id: Optional[str] = None,
        current_version: Optional[str] = None,
        prefix_list_name: Optional[str] = None,
    ) -> ManagedPrefixList:
        managed_pl: ManagedPrefixList = self.managed_prefix_lists.get(prefix_list_id)  # type: ignore
        managed_pl.prefix_list_name = prefix_list_name
        if remove_entry or add_entry:
            latest_version = managed_pl.entries.get(managed_pl.version)  # type: ignore[arg-type]
            entries = (
                managed_pl.entries.get(current_version, latest_version).copy()  # type: ignore
                if managed_pl.entries
                else []
            )
            for item in entries.copy():
                if item.get("Cidr", "") in remove_entry:
                    entries.remove(item)

            for item in add_entry:
                if item not in entries.copy():
                    entries.append(item)
            managed_pl.version += 1  # type: ignore[operator]
            managed_pl.entries[managed_pl.version] = entries
        managed_pl.state = "modify-complete"
        return managed_pl

    def _create_aws_managed_prefix_list(
        self, name: str, address_family: str, entries: List[Dict[str, str]]
    ) -> None:
        managed_prefix_list = self.create_managed_prefix_list(
            address_family=address_family,
            entry=entries,
            prefix_list_name=name,
            owner_id="aws",
        )
        managed_prefix_list.version = None
        managed_prefix_list.max_entries = None
        self.managed_prefix_lists[managed_prefix_list.id] = managed_prefix_list

    def create_default_pls(self) -> None:
        # See https://docs.aws.amazon.com/vpc/latest/userguide/working-with-aws-managed-prefix-lists.html

        # S3
        self._create_aws_managed_prefix_list(
            name=f"com.amazonaws.{self.region_name}.s3",  # type: ignore[attr-defined]
            address_family="IPv4",
            entries=[
                {"Cidr": "52.216.0.0/15", "Description": "default"},
                {"Cidr": "3.5.0.0/19", "Description": "default"},
                {"Cidr": "54.231.0.0/16", "Description": "default"},
            ],
        )

        # DynamoDB
        self._create_aws_managed_prefix_list(
            name=f"com.amazonaws.{self.region_name}.dynamodb",  # type: ignore[attr-defined]
            address_family="IPv4",
            entries=[
                {"Cidr": "3.218.182.0/24", "Description": "default"},
                {"Cidr": "3.218.180.0/23", "Description": "default"},
                {"Cidr": "52.94.0.0/22", "Description": "default"},
                {"Cidr": "52.119.224.0/20", "Description": "default"},
            ],
        )

        # CloudFront
        self._create_aws_managed_prefix_list(
            name="com.amazonaws.global.cloudfront.origin-facing",
            address_family="IPv4",
            entries=[
                {"Cidr": "13.124.199.0/24", "Description": "default"},
                {"Cidr": "130.176.0.0/18", "Description": "default"},
                {"Cidr": "15.158.0.0/16", "Description": "default"},
                {"Cidr": "18.68.0.0/16", "Description": "default"},
                {"Cidr": "204.246.166.0/24", "Description": "default"},
                {"Cidr": "205.251.218.0/24", "Description": "default"},
                {"Cidr": "3.172.0.0/18", "Description": "default"},
                {"Cidr": "54.239.208.0/21", "Description": "default"},
                {"Cidr": "64.252.64.0/18", "Description": "default"},
                {"Cidr": "70.132.0.0/18", "Description": "default"},
            ],
        )

        # Ground Station
        self._create_aws_managed_prefix_list(
            name="com.amazonaws.global.groundstation",
            address_family="IPv4",
            entries=[{"Cidr": "3.2.16.0/20", "Description": "default"}],
        )

        # VPC Lattice
        self._create_aws_managed_prefix_list(
            name=f"com.amazonaws.{self.region_name}.vpc-lattice",  # type: ignore[attr-defined]
            address_family="IPv4",
            entries=[{"Cidr": "169.254.171.0/24", "Description": "default"}],
        )

        # VPC Lattice ipv6
        self._create_aws_managed_prefix_list(
            name=f"com.amazonaws.{self.region_name}.ipv6.vpc-lattice",  # type: ignore[attr-defined]
            address_family="IPv6",
            entries=[{"Cidr": "fd00:ec2:80::/64", "Description": "default"}],
        )
