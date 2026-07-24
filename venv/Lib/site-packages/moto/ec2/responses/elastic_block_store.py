from moto.core.responses import ActionResult, EmptyResult

from ._base_response import EC2BaseResponse


class ElasticBlockStore(EC2BaseResponse):
    def attach_volume(self) -> ActionResult:
        volume_id = self._get_param("VolumeId")
        instance_id = self._get_param("InstanceId")
        device_path = self._get_param("Device")

        self.error_on_dryrun()

        attachment = self.ec2_backend.attach_volume(volume_id, instance_id, device_path)
        return ActionResult(
            {
                "VolumeId": attachment.volume.id,
                "InstanceId": attachment.instance.id,
                "Device": attachment.device,
                "State": "attaching",
                "AttachTime": attachment.attach_time,
            }
        )

    def copy_snapshot(self) -> ActionResult:
        source_snapshot_id = self._get_param("SourceSnapshotId")
        source_region = self._get_param("SourceRegion")
        description = self._get_param("Description")
        tags = self._parse_tag_specification()
        snapshot_tags = tags.get("snapshot", {})
        kms_key_id = self._get_param("KmsKeyId")

        self.error_on_dryrun()

        snapshot = self.ec2_backend.copy_snapshot(
            source_snapshot_id=source_snapshot_id,
            source_region=source_region,
            description=description,
            kms_key_id=kms_key_id,
        )
        snapshot.add_tags(snapshot_tags)
        return ActionResult({"SnapshotId": snapshot.id, "Tags": snapshot.tag_set})

    def create_snapshot(self) -> ActionResult:
        volume_id = self._get_param("VolumeId")
        description = self._get_param("Description")
        tags = self._parse_tag_specification()
        snapshot_tags = tags.get("snapshot", {})

        self.error_on_dryrun()

        snapshot = self.ec2_backend.create_snapshot(volume_id, description)
        snapshot.add_tags(snapshot_tags)
        return ActionResult(snapshot)

    def create_snapshots(self) -> ActionResult:
        instance_spec = self._get_param("InstanceSpecification", {})
        description = self._get_param("Description", "")
        tags = self._parse_tag_specification()
        snapshot_tags = tags.get("snapshot", {})

        self.error_on_dryrun()

        snapshots = self.ec2_backend.create_snapshots(
            instance_spec, description, snapshot_tags
        )
        return ActionResult({"Snapshots": snapshots})

    def create_volume(self) -> ActionResult:
        size = self._get_param("Size")
        zone = self._get_param("AvailabilityZone")
        snapshot_id = self._get_param("SnapshotId")
        volume_type = self._get_param("VolumeType")
        tags = self._parse_tag_specification()
        volume_tags = tags.get("volume", {})
        encrypted = self._get_bool_param("Encrypted", False)
        kms_key_id = self._get_param("KmsKeyId")
        iops = self._get_param("Iops")
        throughput = self._get_param("Throughput")
        multi_attach_enabled = self._get_param("MultiAttachEnabled")

        self.error_on_dryrun()

        volume = self.ec2_backend.create_volume(
            size=size,
            zone_name=zone,
            snapshot_id=snapshot_id,
            encrypted=encrypted,
            kms_key_id=kms_key_id,
            volume_type=volume_type,
            iops=iops,
            throughput=throughput,
            multi_attach_enabled=multi_attach_enabled,
        )
        volume.add_tags(volume_tags)
        return ActionResult(volume)

    def modify_volume(self) -> ActionResult:
        volume_id = self._get_param("VolumeId")
        target_size = self._get_param("Size")
        target_volume_type = self._get_param("VolumeType")
        target_iops = self._get_param("Iops")
        target_throughput = self._get_param("Throughput")
        target_multi_attach_enabled = self._get_param("MultiAttachEnabled")

        self.error_on_dryrun()

        volume = self.ec2_backend.modify_volume(
            volume_id,
            target_size,
            target_volume_type,
            target_iops,
            target_throughput,
            target_multi_attach_enabled,
        )
        modification = volume.modifications[-1]
        return ActionResult({"VolumeModification": modification})

    def describe_volumes_modifications(self) -> ActionResult:
        filters = self._filters_from_querystring()
        volume_ids = self._get_param("VolumeIds", [])
        modifications = self.ec2_backend.describe_volumes_modifications(
            volume_ids=volume_ids, filters=filters
        )
        return ActionResult({"VolumesModifications": modifications})

    def delete_snapshot(self) -> ActionResult:
        snapshot_id = self._get_param("SnapshotId")

        self.error_on_dryrun()

        self.ec2_backend.delete_snapshot(snapshot_id)
        return EmptyResult()

    def delete_volume(self) -> ActionResult:
        volume_id = self._get_param("VolumeId")

        self.error_on_dryrun()

        self.ec2_backend.delete_volume(volume_id)
        return EmptyResult()

    def describe_snapshots(self) -> ActionResult:
        filters = self._filters_from_querystring()
        snapshot_ids = self._get_param("SnapshotIds", [])
        snapshots = self.ec2_backend.describe_snapshots(
            snapshot_ids=snapshot_ids, filters=filters
        )
        return ActionResult({"Snapshots": snapshots})

    def describe_volumes(self) -> ActionResult:
        filters = self._filters_from_querystring()
        volume_ids = self._get_param("VolumeIds", [])
        volumes = self.ec2_backend.describe_volumes(
            volume_ids=volume_ids, filters=filters
        )
        return ActionResult({"Volumes": volumes})

    def describe_volume_attribute(self) -> str:
        raise NotImplementedError(
            "ElasticBlockStore.describe_volume_attribute is not yet implemented"
        )

    def describe_volume_status(self) -> str:
        raise NotImplementedError(
            "ElasticBlockStore.describe_volume_status is not yet implemented"
        )

    def detach_volume(self) -> ActionResult:
        volume_id = self._get_param("VolumeId")
        instance_id = self._get_param("InstanceId")
        device_path = self._get_param("Device")

        self.error_on_dryrun()

        attachment = self.ec2_backend.detach_volume(volume_id, instance_id, device_path)
        return ActionResult(
            {
                "VolumeId": attachment.volume.id,
                "InstanceId": attachment.instance.id,
                "Device": attachment.device,
                "State": "detaching",
                "AttachTime": attachment.attach_time,
            }
        )

    def enable_volume_io(self) -> str:
        self.error_on_dryrun()

        raise NotImplementedError(
            "ElasticBlockStore.enable_volume_io is not yet implemented"
        )

    def import_volume(self) -> str:
        self.error_on_dryrun()

        raise NotImplementedError(
            "ElasticBlockStore.import_volume is not yet implemented"
        )

    def describe_snapshot_attribute(self) -> ActionResult:
        snapshot_id = self._get_param("SnapshotId")
        groups = self.ec2_backend.get_create_volume_permission_groups(snapshot_id)
        user_ids = self.ec2_backend.get_create_volume_permission_userids(snapshot_id)
        permissions = [{"Group": g} for g in groups] + [{"UserId": u} for u in user_ids]
        return ActionResult(
            {
                "SnapshotId": snapshot_id,
                "CreateVolumePermissions": permissions,
            }
        )

    def modify_snapshot_attribute(self) -> ActionResult:
        snapshot_id = self._get_param("SnapshotId")
        operation_type = self._get_param("OperationType")
        groups = self._get_param("GroupNames", [])
        user_ids = self._get_param("UserIds", [])

        self.error_on_dryrun()

        if operation_type == "add":
            self.ec2_backend.add_create_volume_permission(
                snapshot_id, user_ids=user_ids, groups=groups
            )
        elif operation_type == "remove":
            self.ec2_backend.remove_create_volume_permission(
                snapshot_id, user_ids=user_ids, groups=groups
            )
        return EmptyResult()

    def modify_volume_attribute(self) -> str:
        self.error_on_dryrun()

        raise NotImplementedError(
            "ElasticBlockStore.modify_volume_attribute is not yet implemented"
        )

    def reset_snapshot_attribute(self) -> str:
        self.error_on_dryrun()

        raise NotImplementedError(
            "ElasticBlockStore.reset_snapshot_attribute is not yet implemented"
        )

    def modify_ebs_default_kms_key_id(self) -> ActionResult:
        self.error_on_dryrun()
        kms_key_id = self._get_param("KmsKeyId")
        new_default_kms_arn = self.ec2_backend.modify_ebs_default_kms_key_id(kms_key_id)
        return ActionResult({"KmsKeyId": new_default_kms_arn})
