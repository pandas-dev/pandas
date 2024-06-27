from typing import Any, Dict, Iterable, List, Optional, Set

from moto.core.common_models import CloudFormationModel
from moto.packages.boto.ec2.blockdevicemapping import BlockDeviceType

from ..exceptions import (
    InvalidAMIAttributeItemValueError,
    InvalidParameterDependency,
    InvalidParameterValueError,
    InvalidSnapshotIdError,
    InvalidSnapshotInUse,
    InvalidVolumeAttachmentError,
    InvalidVolumeDetachmentError,
    InvalidVolumeIdError,
    VolumeInUseError,
)
from ..utils import (
    generic_filter,
    random_snapshot_id,
    random_volume_id,
    utc_date_and_time,
)
from .core import TaggedEC2Resource

IOPS_REQUIRED_VOLUME_TYPES = ["io1", "io2"]
IOPS_SUPPORTED_VOLUME_TYPES = ["gp3", "io1", "io2"]
THROUGHPUT_SUPPORTED_VOLUME_TYPES = ["gp3"]
GP3_DEFAULT_IOPS = 3000


class VolumeModification:
    def __init__(
        self,
        volume: "Volume",
        target_size: Optional[int] = None,
        target_volume_type: Optional[str] = None,
    ):
        if not any([target_size, target_volume_type]):
            raise InvalidParameterValueError(
                "Invalid input: Must specify at least one of size or type"
            )

        self.volume = volume
        self.original_size = volume.size
        self.original_volume_type = volume.volume_type
        self.target_size = target_size or volume.size
        self.target_volume_type = target_volume_type or volume.volume_type

        self.start_time = utc_date_and_time()
        self.end_time = utc_date_and_time()

    def get_filter_value(self, filter_name: str) -> Any:
        if filter_name == "original-size":
            return self.original_size
        elif filter_name == "original-volume-type":
            return self.original_volume_type
        elif filter_name == "target-size":
            return self.target_size
        elif filter_name == "target-volume-type":
            return self.target_volume_type
        elif filter_name == "volume-id":
            return self.volume.id


class VolumeAttachment(CloudFormationModel):
    def __init__(self, volume: "Volume", instance: Any, device: str, status: str):
        self.volume = volume
        self.attach_time = utc_date_and_time()
        self.instance = instance
        self.device = device
        self.status = status

    @staticmethod
    def cloudformation_name_type() -> str:
        return ""

    @staticmethod
    def cloudformation_type() -> str:
        # https://docs.aws.amazon.com/AWSCloudFormation/latest/UserGuide/aws-resource-ec2-volumeattachment.html
        return "AWS::EC2::VolumeAttachment"

    @classmethod
    def create_from_cloudformation_json(  # type: ignore[misc]
        cls,
        resource_name: str,
        cloudformation_json: Any,
        account_id: str,
        region_name: str,
        **kwargs: Any,
    ) -> "VolumeAttachment":
        from ..models import ec2_backends

        properties = cloudformation_json["Properties"]

        instance_id = properties["InstanceId"]
        volume_id = properties["VolumeId"]

        ec2_backend = ec2_backends[account_id][region_name]
        return ec2_backend.attach_volume(  # type: ignore[return-value]
            volume_id=volume_id,
            instance_id=instance_id,
            device_path=properties["Device"],
        )


class Volume(TaggedEC2Resource, CloudFormationModel):
    def __init__(
        self,
        ec2_backend: Any,
        volume_id: str,
        size: int,
        zone: Any,
        snapshot_id: Optional[str] = None,
        encrypted: bool = False,
        kms_key_id: Optional[str] = None,
        volume_type: Optional[str] = None,
        iops: Optional[int] = None,
        throughput: Optional[int] = None,
    ):
        self.id = volume_id
        self.volume_type = volume_type or "gp2"
        self.size = size
        self.zone = zone
        self.create_time = utc_date_and_time()
        self.attachment: Optional[VolumeAttachment] = None
        self.snapshot_id = snapshot_id
        self.ec2_backend = ec2_backend
        self.encrypted = encrypted
        self.kms_key_id = kms_key_id
        self.modifications: List[VolumeModification] = []
        self.iops = iops
        self.throughput = throughput

    def modify(
        self,
        target_size: Optional[int] = None,
        target_volume_type: Optional[str] = None,
    ) -> None:
        modification = VolumeModification(
            volume=self, target_size=target_size, target_volume_type=target_volume_type
        )
        self.modifications.append(modification)

        self.size = modification.target_size
        self.volume_type = modification.target_volume_type

    @staticmethod
    def cloudformation_name_type() -> str:
        return ""

    @staticmethod
    def cloudformation_type() -> str:
        # https://docs.aws.amazon.com/AWSCloudFormation/latest/UserGuide/aws-resource-ec2-volume.html
        return "AWS::EC2::Volume"

    @classmethod
    def create_from_cloudformation_json(  # type: ignore[misc]
        cls,
        resource_name: str,
        cloudformation_json: Any,
        account_id: str,
        region_name: str,
        **kwargs: Any,
    ) -> "Volume":
        from ..models import ec2_backends

        properties = cloudformation_json["Properties"]

        ec2_backend = ec2_backends[account_id][region_name]
        volume = ec2_backend.create_volume(
            size=properties.get("Size"), zone_name=properties.get("AvailabilityZone")
        )
        return volume

    @property
    def physical_resource_id(self) -> str:
        return self.id

    @property
    def status(self) -> str:
        if self.attachment:
            return "in-use"
        else:
            return "available"

    def get_filter_value(
        self, filter_name: str, method_name: Optional[str] = None
    ) -> Any:
        if filter_name.startswith("attachment") and not self.attachment:
            return None
        elif filter_name == "attachment.attach-time":
            return self.attachment.attach_time  # type: ignore[union-attr]
        elif filter_name == "attachment.device":
            return self.attachment.device  # type: ignore[union-attr]
        elif filter_name == "attachment.instance-id":
            return self.attachment.instance.id  # type: ignore[union-attr]
        elif filter_name == "attachment.status":
            return self.attachment.status  # type: ignore[union-attr]
        elif filter_name == "create-time":
            return self.create_time
        elif filter_name == "size":
            return self.size
        elif filter_name == "snapshot-id":
            return self.snapshot_id
        elif filter_name == "status":
            return self.status
        elif filter_name == "volume-id":
            return self.id
        elif filter_name == "encrypted":
            return str(self.encrypted).lower()
        elif filter_name == "availability-zone":
            return self.zone.name if self.zone else None
        else:
            return super().get_filter_value(filter_name, "DescribeVolumes")


class Snapshot(TaggedEC2Resource):
    def __init__(
        self,
        ec2_backend: Any,
        snapshot_id: str,
        volume: Any,
        description: str,
        encrypted: bool = False,
        owner_id: Optional[str] = None,
        from_ami: Optional[str] = None,
    ):
        self.id = snapshot_id
        self.volume = volume
        self.description = description
        self.start_time = utc_date_and_time()
        self.create_volume_permission_groups: Set[str] = set()
        self.create_volume_permission_userids: Set[str] = set()
        self.ec2_backend = ec2_backend
        self.status = "completed"
        self.encrypted = encrypted
        self.owner_id = owner_id or ec2_backend.account_id
        self.from_ami = from_ami

    def get_filter_value(
        self, filter_name: str, method_name: Optional[str] = None
    ) -> Any:
        if filter_name == "description":
            return self.description
        elif filter_name == "snapshot-id":
            return self.id
        elif filter_name == "start-time":
            return self.start_time
        elif filter_name == "volume-id":
            return self.volume.id
        elif filter_name == "volume-size":
            return self.volume.size
        elif filter_name == "encrypted":
            return str(self.encrypted).lower()
        elif filter_name == "status":
            return self.status
        elif filter_name == "owner-id":
            return self.owner_id
        else:
            return super().get_filter_value(filter_name, "DescribeSnapshots")


class EBSBackend:
    def __init__(self) -> None:
        self.volumes: Dict[str, Volume] = {}
        self.attachments: Dict[str, VolumeAttachment] = {}
        self.snapshots: Dict[str, Snapshot] = {}

    def create_volume(
        self,
        size: int,
        zone_name: str,
        snapshot_id: Optional[str] = None,
        encrypted: bool = False,
        kms_key_id: Optional[str] = None,
        volume_type: Optional[str] = None,
        iops: Optional[int] = None,
        throughput: Optional[int] = None,
    ) -> Volume:
        if kms_key_id and not encrypted:
            raise InvalidParameterDependency("KmsKeyId", "Encrypted")
        if encrypted and not kms_key_id:
            kms_key_id = self._get_default_encryption_key()
        if volume_type in IOPS_REQUIRED_VOLUME_TYPES and not iops:
            raise InvalidParameterDependency("VolumeType", "Iops")
        elif volume_type == "gp3" and not iops:
            iops = GP3_DEFAULT_IOPS
        elif volume_type not in IOPS_SUPPORTED_VOLUME_TYPES and iops:
            raise InvalidParameterDependency("VolumeType", "Iops")
        if volume_type not in THROUGHPUT_SUPPORTED_VOLUME_TYPES and throughput:
            raise InvalidParameterDependency("VolumeType", "Throughput")

        volume_id = random_volume_id()
        zone = self.get_zone_by_name(zone_name)  # type: ignore[attr-defined]
        if snapshot_id:
            snapshot = self.get_snapshot(snapshot_id)
            if size is None:
                size = snapshot.volume.size
            if snapshot.encrypted:
                encrypted = snapshot.encrypted
        volume = Volume(
            self,
            volume_id=volume_id,
            size=size,
            zone=zone,
            snapshot_id=snapshot_id,
            encrypted=encrypted,
            kms_key_id=kms_key_id,
            volume_type=volume_type,
            iops=iops,
            throughput=throughput,
        )
        self.volumes[volume_id] = volume
        return volume

    def describe_volumes(
        self, volume_ids: Optional[List[str]] = None, filters: Any = None
    ) -> List[Volume]:
        matches = list(self.volumes.values())
        if volume_ids:
            matches = [vol for vol in matches if vol.id in volume_ids]
            if len(volume_ids) > len(matches):
                unknown_ids = set(volume_ids) - set(matches)  # type: ignore[arg-type]
                raise InvalidVolumeIdError(unknown_ids)
        if filters:
            matches = generic_filter(filters, matches)
        return matches

    def modify_volume(
        self,
        volume_id: str,
        target_size: Optional[int] = None,
        target_volume_type: Optional[str] = None,
    ) -> Volume:
        volume = self.get_volume(volume_id)
        volume.modify(target_size=target_size, target_volume_type=target_volume_type)
        return volume

    def describe_volumes_modifications(
        self, volume_ids: Optional[List[str]] = None, filters: Any = None
    ) -> List[VolumeModification]:
        volumes = self.describe_volumes(volume_ids)
        modifications = []
        for volume in volumes:
            modifications.extend(volume.modifications)
        if filters:
            modifications = generic_filter(filters, modifications)
        return modifications

    def get_volume(self, volume_id: str) -> Volume:
        volume = self.volumes.get(volume_id, None)
        if not volume:
            raise InvalidVolumeIdError(volume_id)
        return volume

    def delete_volume(self, volume_id: str) -> Volume:
        if volume_id in self.volumes:
            volume = self.volumes[volume_id]
            if volume.attachment:
                raise VolumeInUseError(volume_id, volume.attachment.instance.id)
            return self.volumes.pop(volume_id)
        raise InvalidVolumeIdError(volume_id)

    def attach_volume(
        self,
        volume_id: str,
        instance_id: str,
        device_path: str,
        delete_on_termination: bool = False,
    ) -> Optional[VolumeAttachment]:
        volume = self.get_volume(volume_id)
        instance = self.get_instance(instance_id)  # type: ignore[attr-defined]

        if not volume or not instance:
            return None

        volume.attachment = VolumeAttachment(volume, instance, device_path, "attached")
        # Modify instance to capture mount of block device.
        bdt = BlockDeviceType(
            volume_id=volume_id,
            status=volume.status,
            size=volume.size,
            attach_time=utc_date_and_time(),
            delete_on_termination=delete_on_termination,
        )
        instance.block_device_mapping[device_path] = bdt
        return volume.attachment

    def detach_volume(
        self, volume_id: str, instance_id: str, device_path: str
    ) -> VolumeAttachment:
        volume = self.get_volume(volume_id)
        instance = self.get_instance(instance_id)  # type: ignore[attr-defined]

        old_attachment = volume.attachment
        if not old_attachment:
            raise InvalidVolumeAttachmentError(volume_id, instance_id)
        device_path = device_path or old_attachment.device

        try:
            del instance.block_device_mapping[device_path]
        except KeyError:
            raise InvalidVolumeDetachmentError(volume_id, instance_id, device_path)

        old_attachment.status = "detached"

        volume.attachment = None
        return old_attachment

    def create_snapshot(
        self,
        volume_id: str,
        description: str,
        owner_id: Optional[str] = None,
        from_ami: Optional[str] = None,
    ) -> Snapshot:
        snapshot_id = random_snapshot_id()
        volume = self.get_volume(volume_id)
        params = [self, snapshot_id, volume, description, volume.encrypted]
        if owner_id:
            params.append(owner_id)
        if from_ami:
            params.append(from_ami)
        snapshot = Snapshot(*params)  # type: ignore[arg-type]
        self.snapshots[snapshot_id] = snapshot
        return snapshot

    def create_snapshots(
        self, instance_spec: Dict[str, Any], description: str, tags: Dict[str, str]
    ) -> List[Snapshot]:
        """
        The CopyTagsFromSource-parameter is not yet implemented.
        """
        instance = self.get_instance(instance_spec["InstanceId"])  # type: ignore[attr-defined]
        block_device_mappings = instance.block_device_mapping

        if str(instance_spec.get("ExcludeBootVolume", False)).lower() == "true":
            volumes = [
                m.volume_id
                for k, m in block_device_mappings.items()
                if k != instance.root_device_name
            ]
        else:
            volumes = [m.volume_id for m in block_device_mappings.values()]

        snapshots = [
            self.create_snapshot(v_id, description=description) for v_id in volumes
        ]
        for snapshot in snapshots:
            snapshot.add_tags(tags)
        return snapshots

    def describe_snapshots(
        self, snapshot_ids: Optional[List[str]] = None, filters: Any = None
    ) -> List[Snapshot]:
        matches = list(self.snapshots.values())
        if snapshot_ids:
            matches = [snap for snap in matches if snap.id in snapshot_ids]
            if len(snapshot_ids) > len(matches):
                raise InvalidSnapshotIdError()
        if filters:
            matches = generic_filter(filters, matches)
        return matches

    def copy_snapshot(
        self, source_snapshot_id: str, source_region: str, description: str
    ) -> Snapshot:
        from ..models import ec2_backends

        backend = ec2_backends[self.account_id][source_region]  # type: ignore[attr-defined]
        source_snapshot = backend.describe_snapshots(snapshot_ids=[source_snapshot_id])[
            0
        ]
        snapshot_id = random_snapshot_id()
        snapshot = Snapshot(
            self,
            snapshot_id,
            volume=source_snapshot.volume,
            description=description,
            encrypted=source_snapshot.encrypted,
        )
        self.snapshots[snapshot_id] = snapshot
        return snapshot

    def get_snapshot(self, snapshot_id: str) -> Snapshot:
        snapshot = self.snapshots.get(snapshot_id, None)
        if not snapshot:
            raise InvalidSnapshotIdError()
        return snapshot

    def delete_snapshot(self, snapshot_id: str) -> Snapshot:
        if snapshot_id in self.snapshots:
            snapshot = self.snapshots[snapshot_id]
            if snapshot.from_ami and snapshot.from_ami in self.amis:  # type: ignore[attr-defined]
                raise InvalidSnapshotInUse(snapshot_id, snapshot.from_ami)
            return self.snapshots.pop(snapshot_id)
        raise InvalidSnapshotIdError()

    def get_create_volume_permission_groups(self, snapshot_id: str) -> Set[str]:
        snapshot = self.get_snapshot(snapshot_id)
        return snapshot.create_volume_permission_groups

    def get_create_volume_permission_userids(self, snapshot_id: str) -> Set[str]:
        snapshot = self.get_snapshot(snapshot_id)
        return snapshot.create_volume_permission_userids

    def add_create_volume_permission(
        self, snapshot_id: str, user_ids: List[str], groups: List[str]
    ) -> None:
        snapshot = self.get_snapshot(snapshot_id)
        if user_ids:
            snapshot.create_volume_permission_userids.update(user_ids)

        if groups and groups != ["all"]:
            raise InvalidAMIAttributeItemValueError("UserGroup", groups)
        else:
            snapshot.create_volume_permission_groups.update(groups)

    def remove_create_volume_permission(
        self,
        snapshot_id: str,
        user_ids: Optional[List[str]] = None,
        groups: Optional[Iterable[str]] = None,
    ) -> None:
        snapshot = self.get_snapshot(snapshot_id)
        if user_ids:
            snapshot.create_volume_permission_userids.difference_update(user_ids)

        if groups and groups != ["all"]:
            raise InvalidAMIAttributeItemValueError("UserGroup", groups)
        else:
            snapshot.create_volume_permission_groups.difference_update(groups)  # type: ignore[arg-type]

    def _get_default_encryption_key(self) -> str:
        # https://aws.amazon.com/kms/features/#AWS_Service_Integration
        # An AWS managed CMK is created automatically when you first create
        # an encrypted resource using an AWS service integrated with KMS.
        from moto.kms import kms_backends

        kms = kms_backends[self.account_id][self.region_name]  # type: ignore[attr-defined]
        ebs_alias = "alias/aws/ebs"
        if not kms.alias_exists(ebs_alias):
            key = kms.create_key(
                policy="",
                key_usage="ENCRYPT_DECRYPT",
                key_spec="SYMMETRIC_DEFAULT",
                description="Default master key that protects my EBS volumes when no other key is defined",
                tags=None,
            )
            kms.add_alias(key.id, ebs_alias)
        ebs_key = kms.describe_key(ebs_alias)
        return ebs_key.arn
