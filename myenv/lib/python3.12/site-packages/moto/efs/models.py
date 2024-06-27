"""Implement models for EFS resources.

See AWS docs for details:
https://docs.aws.amazon.com/efs/latest/ug/whatisefs.html
"""

import json
import time
from copy import deepcopy
from typing import Any, Dict, Iterator, List, Optional, Set, Tuple, Union

from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.common_models import CloudFormationModel
from moto.core.utils import camelcase_to_underscores, underscores_to_camelcase
from moto.ec2 import ec2_backends
from moto.ec2.exceptions import InvalidSubnetIdError
from moto.ec2.models.elastic_network_interfaces import NetworkInterface
from moto.ec2.models.subnets import Subnet
from moto.efs.exceptions import (
    AccessPointNotFound,
    BadRequest,
    FileSystemAlreadyExists,
    FileSystemInUse,
    FileSystemNotFound,
    MountTargetConflict,
    MountTargetNotFound,
    PolicyNotFound,
    SecurityGroupLimitExceeded,
    SecurityGroupNotFound,
    SubnetNotFound,
)
from moto.moto_api._internal import mock_random
from moto.utilities.tagging_service import TaggingService
from moto.utilities.utils import get_partition, md5_hash


def _lookup_az_id(account_id: str, az_name: str) -> Optional[str]:
    """Find the Availability zone ID given the AZ name."""
    ec2 = ec2_backends[account_id][az_name[:-1]]
    for zone in ec2.describe_availability_zones():
        if zone.name == az_name:
            return zone.zone_id
    return None


class AccessPoint(CloudFormationModel):
    def __init__(
        self,
        account_id: str,
        region_name: str,
        client_token: str,
        file_system_id: str,
        name: Optional[str],
        posix_user: Dict[str, Any],
        root_directory: Dict[str, str],
        context: "EFSBackend",
    ):
        self.access_point_id = f"fsap-{mock_random.get_random_hex(8)}"
        self.access_point_arn = f"arn:{get_partition(region_name)}:elasticfilesystem:{region_name}:{account_id}:access-point/{self.access_point_id}"
        self.client_token = client_token
        self.file_system_id = file_system_id
        self.name = name
        self.posix_user = posix_user
        self.account_id = account_id

        if not root_directory:
            root_directory = {"Path": "/"}

        self.root_directory = root_directory
        self.context = context

    def info_json(self) -> Dict[str, Any]:
        tags = self.context.list_tags_for_resource(self.access_point_id)
        return {
            "ClientToken": self.client_token,
            "Name": self.name,
            "Tags": tags,
            "AccessPointId": self.access_point_id,
            "AccessPointArn": self.access_point_arn,
            "FileSystemId": self.file_system_id,
            "PosixUser": self.posix_user,
            "RootDirectory": self.root_directory,
            "OwnerId": self.account_id,
            "LifeCycleState": "available",
        }

    @staticmethod
    def cloudformation_type() -> str:
        return "AWS::EFS::AccessPoint"

    @classmethod
    def create_from_cloudformation_json(  # type: ignore[misc]
        cls,
        resource_name: str,
        cloudformation_json: Any,
        account_id: str,
        region_name: str,
        **kwargs: Any,
    ) -> "AccessPoint":
        props = cloudformation_json["Properties"]
        file_system_id = props["FileSystemId"]
        posix_user = props.get("PosixUser", {})
        root_directory = props.get("RootDirectory", {})
        tags = props.get("AccessPointTags", [])

        backend: EFSBackend = efs_backends[account_id][region_name]
        return backend.create_access_point(
            resource_name,
            file_system_id=file_system_id,
            posix_user=posix_user,
            root_directory=root_directory,
            tags=tags,
        )

    def delete(self, account_id: str, region_name: str) -> None:
        backend: EFSBackend = efs_backends[account_id][region_name]
        backend.delete_access_point(self.access_point_id)


class FileSystem(CloudFormationModel):
    """A model for an EFS File System Volume."""

    def __init__(
        self,
        account_id: str,
        region_name: str,
        creation_token: str,
        file_system_id: str,
        context: "EFSBackend",
        performance_mode: str,
        encrypted: bool,
        kms_key_id: str,
        throughput_mode: str,
        provisioned_throughput_in_mibps: int,
        availability_zone_name: str,
        backup: bool,
    ):
        if availability_zone_name:
            backup = True
        if kms_key_id and not encrypted:
            raise BadRequest('If kms_key_id given, "encrypted" must be True.')

        # Save given parameters
        self.creation_token = creation_token
        self.performance_mode = performance_mode or "generalPurpose"
        self.encrypted = encrypted or False
        self.kms_key_id = kms_key_id
        self.throughput_mode = throughput_mode or "bursting"
        self.provisioned_throughput_in_mibps = provisioned_throughput_in_mibps
        self.availability_zone_name = availability_zone_name
        self.availability_zone_id = None
        if self.availability_zone_name:
            self.availability_zone_id = _lookup_az_id(
                account_id, self.availability_zone_name
            )
        self._backup = backup
        self.lifecycle_policies: List[Dict[str, str]] = []
        self.file_system_policy: Optional[str] = None

        self._context = context

        # Generate AWS-assigned parameters
        self.file_system_id = file_system_id
        self.file_system_arn = f"arn:{get_partition(region_name)}:elasticfilesystem:{region_name}:{account_id}:file-system/{self.file_system_id}"
        self.creation_time = time.time()
        self.owner_id = account_id

        # Initialize some state parameters
        self.life_cycle_state = "available"
        self._mount_targets: Dict[str, MountTarget] = {}
        self._size_value = 0

    @property
    def size_in_bytes(self) -> Dict[str, Any]:  # type: ignore[misc]
        return {
            "Value": self._size_value,
            "ValueInIA": 0,
            "ValueInStandard": self._size_value,
            "Timestamp": time.time(),
        }

    @property
    def physical_resource_id(self) -> str:
        return self.file_system_id

    @property
    def number_of_mount_targets(self) -> int:
        return len(self._mount_targets)

    @property
    def backup_policy(self) -> Optional[Dict[str, str]]:
        if self._backup:
            return {"Status": "ENABLED"}
        else:
            return None

    def info_json(self) -> Dict[str, Any]:
        ret = {
            underscores_to_camelcase(k.capitalize()): v
            for k, v in self.__dict__.items()
            if not k.startswith("_")
        }
        tags = self._context.list_tags_for_resource(self.file_system_id)
        name = ""
        for tag in tags:
            if tag["Key"] == "Name":
                name = tag["Value"]
                break

        ret.update(
            Tags=tags,
            SizeInBytes=self.size_in_bytes,
            NumberOfMountTargets=self.number_of_mount_targets,
            Name=name,
        )
        return ret

    def add_mount_target(self, subnet: Subnet, mount_target: "MountTarget") -> None:
        # Check that the mount target doesn't violate constraints.
        for other_mount_target in self._mount_targets.values():
            if other_mount_target.subnet_vpc_id != subnet.vpc_id:
                raise MountTargetConflict(
                    "requested subnet for new mount target is not in the same VPC as existing mount targets"
                )

        if subnet.availability_zone in self._mount_targets:
            raise MountTargetConflict("mount target already exists in this AZ")

        self._mount_targets[subnet.availability_zone] = mount_target

    def has_mount_target(self, subnet: Subnet) -> bool:
        return subnet.availability_zone in self._mount_targets

    def iter_mount_targets(self) -> Iterator["MountTarget"]:
        for mt in self._mount_targets.values():
            yield mt

    def remove_mount_target(self, subnet: Subnet) -> None:
        del self._mount_targets[subnet.availability_zone]

    @staticmethod
    def cloudformation_name_type() -> str:
        return ""

    @staticmethod
    def cloudformation_type() -> str:
        return "AWS::EFS::FileSystem"

    @classmethod
    def create_from_cloudformation_json(  # type: ignore[misc]
        cls,
        resource_name: str,
        cloudformation_json: Any,
        account_id: str,
        region_name: str,
        **kwargs: Any,
    ) -> "FileSystem":
        props = cloudformation_json["Properties"]
        performance_mode = props.get("PerformanceMode", "generalPurpose")
        encrypted = props.get("Encrypted", False)
        kms_key_id = props.get("KmsKeyId")
        throughput_mode = props.get("ThroughputMode", "bursting")
        provisioned_throughput_in_mibps = props.get("ProvisionedThroughputInMibps", 0)
        availability_zone_name = props.get("AvailabilityZoneName")
        backup = props.get("BackupPolicy", {}).get("Status") == "ENABLED"
        tags = props.get("FileSystemTags", [])

        lifecycle_policies = props.get("LifecyclePolicies", [])

        backend: EFSBackend = efs_backends[account_id][region_name]
        fs = backend.create_file_system(
            resource_name,
            performance_mode=performance_mode,
            encrypted=encrypted,
            kms_key_id=kms_key_id,
            throughput_mode=throughput_mode,
            provisioned_throughput_in_mibps=provisioned_throughput_in_mibps,
            availability_zone_name=availability_zone_name,
            backup=backup,
            tags=tags,
        )

        if lifecycle_policies:
            backend.put_lifecycle_configuration(
                file_system_id=fs.file_system_id, policies=lifecycle_policies
            )
        return fs

    @classmethod
    def delete_from_cloudformation_json(  # type: ignore[misc]
        cls,
        resource_name: str,
        cloudformation_json: Any,
        account_id: str,
        region_name: str,
    ) -> None:
        efs_backends[account_id][region_name].delete_file_system(resource_name)


class MountTarget(CloudFormationModel):
    """A model for an EFS Mount Target."""

    def __init__(
        self,
        account_id: str,
        file_system: FileSystem,
        subnet: Subnet,
        ip_address: Optional[str],
        security_groups: Optional[List[str]],
    ):
        # Set the simple given parameters.
        self.file_system_id = file_system.file_system_id
        self._file_system = file_system
        self._file_system.add_mount_target(subnet, self)
        self.subnet_id = subnet.id
        self._subnet = subnet
        self.vpc_id = subnet.vpc_id
        self.security_groups = security_groups

        # Check the number of security groups.
        if self.security_groups is not None and len(self.security_groups) > 5:
            raise SecurityGroupLimitExceeded(
                "The maximum number of security groups per interface has been reached."
            )

        # Get an IP address if needed, otherwise validate the one we're given.
        if ip_address is None:
            ip_address = subnet.get_available_subnet_ip(self)  # type: ignore[arg-type]
        else:
            try:
                subnet.request_ip(ip_address, self)  # type: ignore[arg-type]
            except Exception as e:
                if "IP" in str(e) and "CIDR" in str(e):
                    raise BadRequest(
                        "Address does not fall within the subnet's address range"
                    )
                else:
                    raise e
        self.ip_address = ip_address

        # Init non-user-assigned values.
        self.owner_id = account_id
        self.mount_target_id = f"fsmt-{mock_random.get_random_hex()}"
        self.life_cycle_state = "available"
        self.network_interface_id: Optional[str] = None
        self.availability_zone_id = subnet.availability_zone_id
        self.availability_zone_name = subnet.availability_zone

    def clean_up(self) -> None:
        self._file_system.remove_mount_target(self._subnet)
        self._subnet.del_subnet_ip(self.ip_address)

    def set_network_interface(self, network_interface: NetworkInterface) -> None:
        self.network_interface_id = network_interface.id

    def info_json(self) -> Dict[str, Any]:
        return {
            underscores_to_camelcase(k.capitalize()): v
            for k, v in self.__dict__.items()
            if not k.startswith("_")
        }

    @property
    def physical_resource_id(self) -> str:
        return self.mount_target_id

    @property
    def subnet_vpc_id(self) -> str:
        return self._subnet.vpc_id

    @staticmethod
    def cloudformation_name_type() -> str:
        return ""

    @staticmethod
    def cloudformation_type() -> str:
        return "AWS::EFS::MountTarget"

    @classmethod
    def create_from_cloudformation_json(  # type: ignore[misc]
        cls,
        resource_name: str,
        cloudformation_json: Any,
        account_id: str,
        region_name: str,
        **kwargs: Any,
    ) -> "MountTarget":
        props = deepcopy(cloudformation_json["Properties"])
        props = {camelcase_to_underscores(k): v for k, v in props.items()}
        return efs_backends[account_id][region_name].create_mount_target(**props)

    @classmethod
    def update_from_cloudformation_json(  # type: ignore[misc]
        cls,
        original_resource: Any,
        new_resource_name: str,
        cloudformation_json: Any,
        account_id: str,
        region_name: str,
    ) -> None:
        raise NotImplementedError(
            "Updates of EFS Mount Target via cloudformation are not yet implemented."
        )

    @classmethod
    def delete_from_cloudformation_json(  # type: ignore[misc]
        cls,
        resource_name: str,
        cloudformation_json: Any,
        account_id: str,
        region_name: str,
    ) -> None:
        efs_backends[account_id][region_name].delete_mount_target(resource_name)


class EFSBackend(BaseBackend):
    """The backend manager of EFS resources.

    This is the state-machine for each region, tracking the file systems, mount targets,
    and eventually access points that are deployed. Creating, updating, and destroying
    such resources should always go through this class.
    """

    def __init__(self, region_name: str, account_id: str):
        super().__init__(region_name, account_id)
        self.creation_tokens: Set[str] = set()
        self.access_points: Dict[str, AccessPoint] = dict()
        self.file_systems_by_id: Dict[str, FileSystem] = {}
        self.mount_targets_by_id: Dict[str, MountTarget] = {}
        self.next_markers: Dict[str, Union[List[MountTarget], List[FileSystem]]] = {}
        self.tagging_service = TaggingService()

    def _mark_description(
        self, corpus: Union[List[MountTarget], List[FileSystem]], max_items: int
    ) -> Optional[str]:
        if max_items < len(corpus):
            new_corpus = corpus[max_items:]
            new_corpus_dict = [c.info_json() for c in new_corpus]
            new_hash = md5_hash(json.dumps(new_corpus_dict).encode("utf-8"))
            next_marker = new_hash.hexdigest()
            self.next_markers[next_marker] = new_corpus
        else:
            next_marker = None
        return next_marker

    @property
    def ec2_backend(self) -> Any:  # type: ignore[misc]
        return ec2_backends[self.account_id][self.region_name]

    def create_file_system(
        self,
        creation_token: str,
        performance_mode: str,
        encrypted: bool,
        kms_key_id: str,
        throughput_mode: str,
        provisioned_throughput_in_mibps: int,
        availability_zone_name: str,
        backup: bool,
        tags: List[Dict[str, str]],
    ) -> FileSystem:
        """Create a new EFS File System Volume.

        https://docs.aws.amazon.com/efs/latest/ug/API_CreateFileSystem.html
        """
        if not creation_token:
            raise ValueError("No creation token given.")
        if creation_token in self.creation_tokens:
            raise FileSystemAlreadyExists(creation_token)

        # Create a new file system ID:
        def make_id() -> str:
            return f"fs-{mock_random.get_random_hex()}"

        fsid = make_id()
        while fsid in self.file_systems_by_id:
            fsid = make_id()
        self.file_systems_by_id[fsid] = FileSystem(
            self.account_id,
            self.region_name,
            creation_token,
            fsid,
            context=self,
            performance_mode=performance_mode,
            encrypted=encrypted,
            kms_key_id=kms_key_id,
            throughput_mode=throughput_mode,
            provisioned_throughput_in_mibps=provisioned_throughput_in_mibps,
            availability_zone_name=availability_zone_name,
            backup=backup,
        )
        self.tag_resource(fsid, tags)
        self.creation_tokens.add(creation_token)
        return self.file_systems_by_id[fsid]

    def describe_file_systems(
        self,
        marker: Optional[str] = None,
        max_items: int = 10,
        creation_token: Optional[str] = None,
        file_system_id: Optional[str] = None,
    ) -> Tuple[Optional[str], List[FileSystem]]:
        """Describe all the EFS File Systems, or specific File Systems.

        https://docs.aws.amazon.com/efs/latest/ug/API_DescribeFileSystems.html
        """
        # Restrict the possible corpus of results based on inputs.
        if creation_token and file_system_id:
            raise BadRequest(
                "Request cannot contain both a file system ID and a creation token."
            )
        elif creation_token:
            # Handle the creation token case.
            corpus = []
            for fs in self.file_systems_by_id.values():
                if fs.creation_token == creation_token:
                    corpus.append(fs)
        elif file_system_id:
            # Handle the case that a file_system_id is given.
            if file_system_id not in self.file_systems_by_id:
                raise FileSystemNotFound(file_system_id)
            corpus = [self.file_systems_by_id[file_system_id]]
        elif marker is not None:
            # Handle the case that a marker is given.
            if marker not in self.next_markers:
                raise BadRequest("Invalid Marker")
            corpus = self.next_markers[marker]  # type: ignore[assignment]
        else:
            # Handle the vanilla case.
            corpus = [fs for fs in self.file_systems_by_id.values()]

        # Handle the max_items parameter.
        file_systems = corpus[:max_items]
        next_marker = self._mark_description(corpus, max_items)
        return next_marker, file_systems

    def create_mount_target(
        self,
        file_system_id: str,
        subnet_id: str,
        ip_address: Optional[str] = None,
        security_groups: Optional[List[str]] = None,
    ) -> MountTarget:
        """Create a new EFS Mount Target for a given File System to a given subnet.

        Note that you can only create one mount target for each availability zone
        (which is implied by the subnet ID).

        https://docs.aws.amazon.com/efs/latest/ug/API_CreateMountTarget.html
        """
        # Get the relevant existing resources
        try:
            subnet = self.ec2_backend.get_subnet(subnet_id)
        except InvalidSubnetIdError:
            raise SubnetNotFound(subnet_id)
        if file_system_id not in self.file_systems_by_id:
            raise FileSystemNotFound(file_system_id)
        file_system = self.file_systems_by_id[file_system_id]

        # Validate the security groups.
        if security_groups:
            sg_lookup = {sg.id for sg in self.ec2_backend.describe_security_groups()}
            for sg_id in security_groups:
                if sg_id not in sg_lookup:
                    raise SecurityGroupNotFound(sg_id)

        # Create the new mount target
        mount_target = MountTarget(
            self.account_id, file_system, subnet, ip_address, security_groups
        )

        # Establish the network interface.
        network_interface = self.ec2_backend.create_network_interface(
            subnet, [mount_target.ip_address], group_ids=security_groups
        )
        mount_target.set_network_interface(network_interface)

        # Record the new mount target
        self.mount_targets_by_id[mount_target.mount_target_id] = mount_target
        return mount_target

    def describe_mount_targets(
        self,
        max_items: int,
        file_system_id: Optional[str],
        mount_target_id: Optional[str],
        access_point_id: Optional[str],
        marker: Optional[str],
    ) -> Tuple[Optional[str], List[MountTarget]]:
        """Describe the mount targets given an access point ID, mount target ID or a file system ID.

        https://docs.aws.amazon.com/efs/latest/ug/API_DescribeMountTargets.html
        """
        # Restrict the possible corpus of results based on inputs.
        if not (bool(file_system_id) ^ bool(mount_target_id) ^ bool(access_point_id)):
            raise BadRequest("Must specify exactly one mutually exclusive parameter.")

        if access_point_id:
            file_system_id = self.access_points[access_point_id].file_system_id

        if file_system_id:
            # Handle the case that a file_system_id is given.
            if file_system_id not in self.file_systems_by_id:
                raise FileSystemNotFound(file_system_id)
            corpus = [
                mt
                for mt in self.file_systems_by_id[file_system_id].iter_mount_targets()
            ]
        elif mount_target_id:
            if mount_target_id not in self.mount_targets_by_id:
                raise MountTargetNotFound(mount_target_id)
            # Handle mount target specification case.
            corpus = [self.mount_targets_by_id[mount_target_id]]

        # Handle the case that a marker is given. Note that the handling is quite
        # different from that in describe_file_systems.
        if marker is not None:
            if marker not in self.next_markers:
                raise BadRequest("Invalid Marker")
            corpus_mtids = {m.mount_target_id for m in corpus}
            marked_mtids = {m.mount_target_id for m in self.next_markers[marker]}  # type: ignore[union-attr]
            mt_ids = corpus_mtids & marked_mtids
            corpus = [self.mount_targets_by_id[mt_id] for mt_id in mt_ids]

        # Handle the max_items parameter.
        mount_targets = corpus[:max_items]
        next_marker = self._mark_description(corpus, max_items)
        return next_marker, mount_targets

    def delete_file_system(self, file_system_id: str) -> None:
        """Delete the file system specified by the given file_system_id.

        Note that mount targets must be deleted first.

        https://docs.aws.amazon.com/efs/latest/ug/API_DeleteFileSystem.html
        """
        if file_system_id not in self.file_systems_by_id:
            raise FileSystemNotFound(file_system_id)

        file_system = self.file_systems_by_id[file_system_id]
        if file_system.number_of_mount_targets > 0:
            raise FileSystemInUse(
                "Must delete all mount targets before deleting file system."
            )

        del self.file_systems_by_id[file_system_id]
        self.creation_tokens.remove(file_system.creation_token)

    def delete_mount_target(self, mount_target_id: str) -> None:
        """Delete a mount target specified by the given mount_target_id.

        Note that this will also delete a network interface.

        https://docs.aws.amazon.com/efs/latest/ug/API_DeleteMountTarget.html
        """
        if mount_target_id not in self.mount_targets_by_id:
            raise MountTargetNotFound(mount_target_id)

        mount_target = self.mount_targets_by_id[mount_target_id]
        self.ec2_backend.delete_network_interface(mount_target.network_interface_id)
        del self.mount_targets_by_id[mount_target_id]
        mount_target.clean_up()

    def describe_backup_policy(self, file_system_id: str) -> Dict[str, str]:
        backup_policy = self.file_systems_by_id[file_system_id].backup_policy
        if not backup_policy:
            raise PolicyNotFound("None")
        return backup_policy

    def put_lifecycle_configuration(
        self, file_system_id: str, policies: List[Dict[str, str]]
    ) -> None:
        _, fss = self.describe_file_systems(file_system_id=file_system_id)
        file_system = fss[0]
        file_system.lifecycle_policies = policies

    def describe_lifecycle_configuration(
        self, file_system_id: str
    ) -> List[Dict[str, str]]:
        _, fss = self.describe_file_systems(file_system_id=file_system_id)
        file_system = fss[0]
        return file_system.lifecycle_policies

    def describe_mount_target_security_groups(
        self, mount_target_id: str
    ) -> Optional[List[str]]:
        if mount_target_id not in self.mount_targets_by_id:
            raise MountTargetNotFound(mount_target_id)

        mount_target = self.mount_targets_by_id[mount_target_id]
        return mount_target.security_groups

    def modify_mount_target_security_groups(
        self, mount_target_id: str, security_groups: List[str]
    ) -> None:
        if mount_target_id not in self.mount_targets_by_id:
            raise MountTargetNotFound(mount_target_id)

        mount_target = self.mount_targets_by_id[mount_target_id]
        mount_target.security_groups = security_groups

        self.ec2_backend.modify_network_interface_attribute(
            eni_id=mount_target.network_interface_id, group_ids=security_groups
        )

    def create_access_point(
        self,
        client_token: str,
        tags: List[Dict[str, str]],
        file_system_id: str,
        posix_user: Dict[str, Any],
        root_directory: Dict[str, Any],
    ) -> AccessPoint:
        name = next((tag["Value"] for tag in tags if tag["Key"] == "Name"), None)
        access_point = AccessPoint(
            self.account_id,
            self.region_name,
            client_token,
            file_system_id,
            name,
            posix_user,
            root_directory,
            context=self,
        )
        self.tagging_service.tag_resource(access_point.access_point_id, tags)
        self.access_points[access_point.access_point_id] = access_point
        return access_point

    def describe_access_points(
        self,
        access_point_id: Optional[str],
        file_system_id: Optional[str],
    ) -> List[AccessPoint]:
        """
        Pagination is not yet implemented
        """
        if access_point_id:
            if access_point_id not in self.access_points:
                raise AccessPointNotFound(access_point_id)
            return [self.access_points[access_point_id]]
        elif file_system_id:
            return [
                access_point
                for access_point in self.access_points.values()
                if access_point.file_system_id == file_system_id
            ]
        return list(self.access_points.values())

    def delete_access_point(self, access_point_id: str) -> None:
        self.access_points.pop(access_point_id, None)

    def list_tags_for_resource(self, resource_id: str) -> List[Dict[str, str]]:
        return self.tagging_service.list_tags_for_resource(resource_id)["Tags"]

    def tag_resource(self, resource_id: str, tags: List[Dict[str, str]]) -> None:
        self.tagging_service.tag_resource(resource_id, tags)

    def untag_resource(self, resource_id: str, tag_keys: List[str]) -> None:
        self.tagging_service.untag_resource_using_names(resource_id, tag_keys)


efs_backends = BackendDict(EFSBackend, "efs")
