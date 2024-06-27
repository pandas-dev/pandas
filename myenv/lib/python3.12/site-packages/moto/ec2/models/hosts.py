from typing import Any, Dict, List, Optional

from moto.core.utils import unix_time

from ..utils import generic_filter, random_dedicated_host_id
from .core import TaggedEC2Resource


class Host(TaggedEC2Resource):
    def __init__(
        self,
        host_recovery: str,
        zone: str,
        instance_type: str,
        instance_family: str,
        auto_placement: str,
        backend: Any,
    ):
        self.id = random_dedicated_host_id()
        self.state = "available"
        self.host_recovery = host_recovery or "off"
        self.zone = zone
        self.instance_type: Optional[str] = instance_type
        self.instance_family: Optional[str] = instance_family
        self.auto_placement = auto_placement or "on"
        self.ec2_backend = backend
        self.allocation_time = unix_time()

    def release(self) -> None:
        self.state = "released"

    def get_filter_value(
        self, filter_name: str, method_name: Optional[str] = None
    ) -> Any:
        if filter_name == "availability-zone":
            return self.zone
        if filter_name == "state":
            return self.state
        if filter_name == "tag-key":
            return [t["key"] for t in self.get_tags()]
        if filter_name == "instance-type":
            return self.instance_type
        return None


class HostsBackend:
    def __init__(self) -> None:
        self.hosts: Dict[str, Host] = {}

    def allocate_hosts(
        self,
        quantity: int,
        host_recovery: str,
        zone: str,
        instance_type: str,
        instance_family: str,
        auto_placement: str,
        tags: Dict[str, str],
    ) -> List[str]:
        hosts = [
            Host(
                host_recovery,
                zone,
                instance_type,
                instance_family,
                auto_placement,
                self,
            )
            for _ in range(quantity)
        ]
        for host in hosts:
            self.hosts[host.id] = host
            if tags:
                host.add_tags(tags)
        return [host.id for host in hosts]

    def describe_hosts(
        self, host_ids: List[str], filters: Dict[str, Any]
    ) -> List[Host]:
        """
        Pagination is not yet implemented
        """
        results = list(self.hosts.values())
        if host_ids:
            results = [r for r in results if r.id in host_ids]
        if filters:
            results = generic_filter(filters, results)
        return results

    def modify_hosts(
        self,
        host_ids: List[str],
        auto_placement: str,
        host_recovery: str,
        instance_type: str,
        instance_family: str,
    ) -> None:
        for _id in host_ids:
            host = self.hosts[_id]
            if auto_placement is not None:
                host.auto_placement = auto_placement
            if host_recovery is not None:
                host.host_recovery = host_recovery
            if instance_type is not None:
                host.instance_type = instance_type
                host.instance_family = None
            if instance_family is not None:
                host.instance_family = instance_family
                host.instance_type = None

    def release_hosts(self, host_ids: List[str]) -> None:
        for host_id in host_ids:
            self.hosts[host_id].release()
