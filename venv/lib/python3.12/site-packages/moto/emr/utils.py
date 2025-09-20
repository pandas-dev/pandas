import copy
import re
import string
from typing import Any, Dict, Iterator, List, Tuple

from moto.core.utils import iso_8601_datetime_with_milliseconds
from moto.moto_api._internal import mock_random as random


def random_id(size: int = 13) -> str:
    chars = list(range(10)) + list(string.ascii_uppercase)
    return "".join(str(random.choice(chars)) for x in range(size))


def random_cluster_id() -> str:
    return f"j-{random_id()}"


def random_step_id() -> str:
    return f"s-{random_id()}"


def random_instance_group_id() -> str:
    return f"i-{random_id()}"


class ReleaseLabel:
    version_re = re.compile(r"^emr-(\d+)\.(\d+)\.(\d+)$")

    def __init__(self, release_label: str):
        major, minor, patch = self.parse(release_label)

        self.major = major
        self.minor = minor
        self.patch = patch

    @classmethod
    def parse(cls, release_label: str) -> Tuple[int, int, int]:
        if not release_label:
            raise ValueError(f"Invalid empty ReleaseLabel: {release_label}")

        match = cls.version_re.match(release_label)
        if not match:
            raise ValueError(f"Invalid ReleaseLabel: {release_label}")

        major, minor, patch = match.groups()

        return int(major), int(minor), int(patch)

    def __str__(self) -> str:
        version = f"emr-{self.major}.{self.minor}.{self.patch}"
        return version

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}({str(self)})"

    def __iter__(self) -> Iterator[Tuple[int, int, int]]:
        return iter((self.major, self.minor, self.patch))  # type: ignore

    def __eq__(self, other: Any) -> bool:
        if not isinstance(other, self.__class__):
            return NotImplemented
        return (
            self.major == other.major
            and self.minor == other.minor
            and self.patch == other.patch
        )

    def __ne__(self, other: Any) -> bool:
        if not isinstance(other, self.__class__):
            return NotImplemented
        return tuple(self) != tuple(other)

    def __lt__(self, other: Any) -> bool:
        if not isinstance(other, self.__class__):
            return NotImplemented
        return tuple(self) < tuple(other)

    def __le__(self, other: Any) -> bool:
        if not isinstance(other, self.__class__):
            return NotImplemented
        return tuple(self) <= tuple(other)

    def __gt__(self, other: Any) -> bool:
        if not isinstance(other, self.__class__):
            return NotImplemented
        return tuple(self) > tuple(other)

    def __ge__(self, other: Any) -> bool:
        if not isinstance(other, self.__class__):
            return NotImplemented
        return tuple(self) >= tuple(other)


class EmrManagedSecurityGroup:
    class Kind:
        MASTER = "Master"
        SLAVE = "Slave"
        SERVICE = "Service"

    kind = ""

    group_name = ""
    short_name = ""
    desc_fmt = "{short_name} for Elastic MapReduce created on {created}"

    @classmethod
    def description(cls) -> str:
        created = iso_8601_datetime_with_milliseconds()
        return cls.desc_fmt.format(short_name=cls.short_name, created=created)


class EmrManagedMasterSecurityGroup(EmrManagedSecurityGroup):
    kind = EmrManagedSecurityGroup.Kind.MASTER
    group_name = "ElasticMapReduce-Master-Private"
    short_name = "Master"


class EmrManagedSlaveSecurityGroup(EmrManagedSecurityGroup):
    kind = EmrManagedSecurityGroup.Kind.SLAVE
    group_name = "ElasticMapReduce-Slave-Private"
    short_name = "Slave"


class EmrManagedServiceAccessSecurityGroup(EmrManagedSecurityGroup):
    kind = EmrManagedSecurityGroup.Kind.SERVICE
    group_name = "ElasticMapReduce-ServiceAccess"
    short_name = "Service access"


class EmrSecurityGroupManager:
    MANAGED_RULES_EGRESS = [
        {
            "group_name_or_id": EmrManagedSecurityGroup.Kind.MASTER,
            "from_port": None,
            "ip_protocol": "-1",
            "ip_ranges": [{"CidrIp": "0.0.0.0/0"}],
            "to_port": None,
            "source_groups": [],
        },
        {
            "group_name_or_id": EmrManagedSecurityGroup.Kind.SLAVE,
            "from_port": None,
            "ip_protocol": "-1",
            "ip_ranges": [{"CidrIp": "0.0.0.0/0"}],
            "to_port": None,
            "source_groups": [],
        },
        {
            "group_name_or_id": EmrManagedSecurityGroup.Kind.SERVICE,
            "from_port": 8443,
            "ip_protocol": "tcp",
            "ip_ranges": [],
            "to_port": 8443,
            "source_groups": [
                {"GroupId": EmrManagedSecurityGroup.Kind.MASTER},
                {"GroupId": EmrManagedSecurityGroup.Kind.SLAVE},
            ],
        },
    ]

    MANAGED_RULES_INGRESS = [
        {
            "group_name_or_id": EmrManagedSecurityGroup.Kind.MASTER,
            "from_port": 0,
            "ip_protocol": "tcp",
            "ip_ranges": [],
            "to_port": 65535,
            "source_groups": [
                {"GroupId": EmrManagedSecurityGroup.Kind.MASTER},
                {"GroupId": EmrManagedSecurityGroup.Kind.SLAVE},
            ],
        },
        {
            "group_name_or_id": EmrManagedSecurityGroup.Kind.MASTER,
            "from_port": 8443,
            "ip_protocol": "tcp",
            "ip_ranges": [],
            "to_port": 8443,
            "source_groups": [{"GroupId": EmrManagedSecurityGroup.Kind.SERVICE}],
        },
        {
            "group_name_or_id": EmrManagedSecurityGroup.Kind.MASTER,
            "from_port": 0,
            "ip_protocol": "udp",
            "ip_ranges": [],
            "to_port": 65535,
            "source_groups": [
                {"GroupId": EmrManagedSecurityGroup.Kind.MASTER},
                {"GroupId": EmrManagedSecurityGroup.Kind.SLAVE},
            ],
        },
        {
            "group_name_or_id": EmrManagedSecurityGroup.Kind.MASTER,
            "from_port": -1,
            "ip_protocol": "icmp",
            "ip_ranges": [],
            "to_port": -1,
            "source_groups": [
                {"GroupId": EmrManagedSecurityGroup.Kind.MASTER},
                {"GroupId": EmrManagedSecurityGroup.Kind.SLAVE},
            ],
        },
        {
            "group_name_or_id": EmrManagedSecurityGroup.Kind.SLAVE,
            "from_port": 0,
            "ip_protocol": "tcp",
            "ip_ranges": [],
            "to_port": 65535,
            "source_groups": [
                {"GroupId": EmrManagedSecurityGroup.Kind.MASTER},
                {"GroupId": EmrManagedSecurityGroup.Kind.SLAVE},
            ],
        },
        {
            "group_name_or_id": EmrManagedSecurityGroup.Kind.SLAVE,
            "from_port": 8443,
            "ip_protocol": "tcp",
            "ip_ranges": [],
            "to_port": 8443,
            "source_groups": [{"GroupId": EmrManagedSecurityGroup.Kind.SERVICE}],
        },
        {
            "group_name_or_id": EmrManagedSecurityGroup.Kind.SLAVE,
            "from_port": 0,
            "ip_protocol": "udp",
            "ip_ranges": [],
            "to_port": 65535,
            "source_groups": [
                {"GroupId": EmrManagedSecurityGroup.Kind.MASTER},
                {"GroupId": EmrManagedSecurityGroup.Kind.SLAVE},
            ],
        },
        {
            "group_name_or_id": EmrManagedSecurityGroup.Kind.SLAVE,
            "from_port": -1,
            "ip_protocol": "icmp",
            "ip_ranges": [],
            "to_port": -1,
            "source_groups": [
                {"GroupId": EmrManagedSecurityGroup.Kind.MASTER},
                {"GroupId": EmrManagedSecurityGroup.Kind.SLAVE},
            ],
        },
        {
            "group_name_or_id": EmrManagedSecurityGroup.Kind.SERVICE,
            "from_port": 9443,
            "ip_protocol": "tcp",
            "ip_ranges": [],
            "to_port": 9443,
            "source_groups": [{"GroupId": EmrManagedSecurityGroup.Kind.MASTER}],
        },
    ]

    def __init__(self, ec2_backend: Any, vpc_id: str):
        self.ec2 = ec2_backend
        self.vpc_id = vpc_id

    def manage_security_groups(
        self,
        master_security_group: str,
        slave_security_group: str,
        service_access_security_group: str,
    ) -> Tuple[Any, Any, Any]:
        group_metadata = [
            (
                master_security_group,
                EmrManagedSecurityGroup.Kind.MASTER,
                EmrManagedMasterSecurityGroup,
            ),
            (
                slave_security_group,
                EmrManagedSecurityGroup.Kind.SLAVE,
                EmrManagedSlaveSecurityGroup,
            ),
            (
                service_access_security_group,
                EmrManagedSecurityGroup.Kind.SERVICE,
                EmrManagedServiceAccessSecurityGroup,
            ),
        ]
        managed_groups: Dict[str, Any] = {}
        for name, kind, defaults in group_metadata:
            managed_groups[kind] = self._get_or_create_sg(name, defaults)
        self._add_rules_to(managed_groups)
        return (
            managed_groups[EmrManagedSecurityGroup.Kind.MASTER],
            managed_groups[EmrManagedSecurityGroup.Kind.SLAVE],
            managed_groups[EmrManagedSecurityGroup.Kind.SERVICE],
        )

    def _get_or_create_sg(self, sg_id: str, defaults: Any) -> Any:
        find_sg = self.ec2.get_security_group_by_name_or_id
        create_sg = self.ec2.create_security_group
        group_id_or_name = sg_id or defaults.group_name
        group = find_sg(group_id_or_name, self.vpc_id)
        if group is None:
            if group_id_or_name != defaults.group_name:
                raise ValueError(
                    f"The security group '{group_id_or_name}' does not exist"
                )
            group = create_sg(defaults.group_name, defaults.description(), self.vpc_id)
        return group

    def _add_rules_to(self, managed_groups: Dict[str, Any]) -> None:
        rules_metadata = [
            (self.MANAGED_RULES_EGRESS, self.ec2.authorize_security_group_egress),
            (self.MANAGED_RULES_INGRESS, self.ec2.authorize_security_group_ingress),
        ]
        for rules, add_rule in rules_metadata:
            rendered_rules = self._render_rules(rules, managed_groups)
            for rule in rendered_rules:
                from moto.ec2.exceptions import InvalidPermissionDuplicateError

                try:
                    add_rule(vpc_id=self.vpc_id, **rule)
                except InvalidPermissionDuplicateError:
                    # If the rule already exists, we can just move on.
                    pass

    @staticmethod
    def _render_rules(  # type: ignore[misc]
        rules: Any, managed_groups: Dict[str, Any]
    ) -> List[Dict[str, Any]]:
        rendered_rules = copy.deepcopy(rules)
        for rule in rendered_rules:
            rule["group_name_or_id"] = managed_groups[rule["group_name_or_id"]].id
            rule["source_groups"] = [
                {"GroupId": managed_groups[group.get("GroupId")].id}
                for group in rule["source_groups"]
            ]
        return rendered_rules
