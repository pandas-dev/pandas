from ._base_response import EC2BaseResponse


class ElasticBlockStore(EC2BaseResponse):
    def attach_volume(self) -> str:
        volume_id = self._get_param("VolumeId")
        instance_id = self._get_param("InstanceId")
        device_path = self._get_param("Device")

        self.error_on_dryrun()

        attachment = self.ec2_backend.attach_volume(volume_id, instance_id, device_path)
        template = self.response_template(ATTACHED_VOLUME_RESPONSE)
        return template.render(attachment=attachment)

    def copy_snapshot(self) -> str:
        source_snapshot_id = self._get_param("SourceSnapshotId")
        source_region = self._get_param("SourceRegion")
        description = self._get_param("Description")
        tags = self._parse_tag_specification()
        snapshot_tags = tags.get("snapshot", {})

        self.error_on_dryrun()

        snapshot = self.ec2_backend.copy_snapshot(
            source_snapshot_id, source_region, description
        )
        snapshot.add_tags(snapshot_tags)
        template = self.response_template(COPY_SNAPSHOT_RESPONSE)
        return template.render(snapshot=snapshot)

    def create_snapshot(self) -> str:
        volume_id = self._get_param("VolumeId")
        description = self._get_param("Description")
        tags = self._parse_tag_specification()
        snapshot_tags = tags.get("snapshot", {})

        self.error_on_dryrun()

        snapshot = self.ec2_backend.create_snapshot(volume_id, description)
        snapshot.add_tags(snapshot_tags)
        template = self.response_template(CREATE_SNAPSHOT_RESPONSE)
        return template.render(snapshot=snapshot)

    def create_snapshots(self) -> str:
        params = self._get_params()
        instance_spec = params.get("InstanceSpecification")
        description = params.get("Description", "")
        tags = self._parse_tag_specification()
        snapshot_tags = tags.get("snapshot", {})

        self.error_on_dryrun()

        snapshots = self.ec2_backend.create_snapshots(
            instance_spec, description, snapshot_tags
        )
        template = self.response_template(CREATE_SNAPSHOTS_RESPONSE)
        return template.render(snapshots=snapshots)

    def create_volume(self) -> str:
        size = self._get_param("Size")
        zone = self._get_param("AvailabilityZone")
        snapshot_id = self._get_param("SnapshotId")
        volume_type = self._get_param("VolumeType")
        tags = self._parse_tag_specification()
        volume_tags = tags.get("volume", {})
        encrypted = self._get_bool_param("Encrypted", if_none=False)
        kms_key_id = self._get_param("KmsKeyId")
        iops = self._get_param("Iops")
        throughput = self._get_param("Throughput")

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
        )
        volume.add_tags(volume_tags)
        template = self.response_template(CREATE_VOLUME_RESPONSE)
        return template.render(volume=volume)

    def modify_volume(self) -> str:
        volume_id = self._get_param("VolumeId")
        target_size = self._get_param("Size")
        target_volume_type = self._get_param("VolumeType")

        self.error_on_dryrun()

        volume = self.ec2_backend.modify_volume(
            volume_id, target_size, target_volume_type
        )
        template = self.response_template(MODIFY_VOLUME_RESPONSE)
        return template.render(volume=volume)

    def describe_volumes_modifications(self) -> str:
        filters = self._filters_from_querystring()
        volume_ids = self._get_multi_param("VolumeId")
        modifications = self.ec2_backend.describe_volumes_modifications(
            volume_ids=volume_ids, filters=filters
        )
        template = self.response_template(DESCRIBE_VOLUMES_MODIFICATIONS_RESPONSE)
        return template.render(modifications=modifications)

    def delete_snapshot(self) -> str:
        snapshot_id = self._get_param("SnapshotId")

        self.error_on_dryrun()

        self.ec2_backend.delete_snapshot(snapshot_id)
        return DELETE_SNAPSHOT_RESPONSE

    def delete_volume(self) -> str:
        volume_id = self._get_param("VolumeId")

        self.error_on_dryrun()

        self.ec2_backend.delete_volume(volume_id)
        return DELETE_VOLUME_RESPONSE

    def describe_snapshots(self) -> str:
        filters = self._filters_from_querystring()
        snapshot_ids = self._get_multi_param("SnapshotId")
        snapshots = self.ec2_backend.describe_snapshots(
            snapshot_ids=snapshot_ids, filters=filters
        )
        template = self.response_template(DESCRIBE_SNAPSHOTS_RESPONSE)
        return template.render(snapshots=snapshots)

    def describe_volumes(self) -> str:
        filters = self._filters_from_querystring()
        volume_ids = self._get_multi_param("VolumeId")
        volumes = self.ec2_backend.describe_volumes(
            volume_ids=volume_ids, filters=filters
        )
        template = self.response_template(DESCRIBE_VOLUMES_RESPONSE)
        return template.render(volumes=volumes)

    def describe_volume_attribute(self) -> str:
        raise NotImplementedError(
            "ElasticBlockStore.describe_volume_attribute is not yet implemented"
        )

    def describe_volume_status(self) -> str:
        raise NotImplementedError(
            "ElasticBlockStore.describe_volume_status is not yet implemented"
        )

    def detach_volume(self) -> str:
        volume_id = self._get_param("VolumeId")
        instance_id = self._get_param("InstanceId")
        device_path = self._get_param("Device")

        self.error_on_dryrun()

        attachment = self.ec2_backend.detach_volume(volume_id, instance_id, device_path)
        template = self.response_template(DETATCH_VOLUME_RESPONSE)
        return template.render(attachment=attachment)

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

    def describe_snapshot_attribute(self) -> str:
        snapshot_id = self._get_param("SnapshotId")
        groups = self.ec2_backend.get_create_volume_permission_groups(snapshot_id)
        user_ids = self.ec2_backend.get_create_volume_permission_userids(snapshot_id)
        template = self.response_template(DESCRIBE_SNAPSHOT_ATTRIBUTES_RESPONSE)
        return template.render(snapshot_id=snapshot_id, groups=groups, userIds=user_ids)

    def modify_snapshot_attribute(self) -> str:
        snapshot_id = self._get_param("SnapshotId")
        operation_type = self._get_param("OperationType")
        groups = self._get_multi_param("UserGroup")
        user_ids = self._get_multi_param("UserId")

        self.error_on_dryrun()

        if operation_type == "add":
            self.ec2_backend.add_create_volume_permission(
                snapshot_id, user_ids=user_ids, groups=groups
            )
        elif operation_type == "remove":
            self.ec2_backend.remove_create_volume_permission(
                snapshot_id, user_ids=user_ids, groups=groups
            )
        return MODIFY_SNAPSHOT_ATTRIBUTE_RESPONSE

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


CREATE_VOLUME_RESPONSE = """<CreateVolumeResponse xmlns="http://ec2.amazonaws.com/doc/2013-10-15/">
  <requestId>59dbff89-35bd-4eac-99ed-be587EXAMPLE</requestId>
  <volumeId>{{ volume.id }}</volumeId>
  <size>{{ volume.size }}</size>
  {% if volume.snapshot_id %}
    <snapshotId>{{ volume.snapshot_id }}</snapshotId>
  {% else %}
    <snapshotId/>
  {% endif %}
  <encrypted>{{ 'true' if volume.encrypted else 'false' }}</encrypted>
  {% if volume.encrypted %}
  <kmsKeyId>{{ volume.kms_key_id }}</kmsKeyId>
  {% endif %}
  <availabilityZone>{{ volume.zone.name }}</availabilityZone>
  <status>creating</status>
  <createTime>{{ volume.create_time }}</createTime>
  {% if volume.get_tags() %}
    <tagSet>
      {% for tag in volume.get_tags() %}
          <item>
          <resourceId>{{ tag.resource_id }}</resourceId>
          <resourceType>{{ tag.resource_type }}</resourceType>
          <key>{{ tag.key }}</key>
          <value>{{ tag.value }}</value>
          </item>
      {% endfor %}
    </tagSet>
  {% endif %}
  <volumeType>{{ volume.volume_type }}</volumeType>
  {% if volume.iops %}
    <iops>{{ volume.iops }}</iops>
  {% endif %}
  {% if volume.throughput %}
    <throughput>{{ volume.throughput }}</throughput>
  {% endif %}
</CreateVolumeResponse>"""

DESCRIBE_VOLUMES_RESPONSE = """<DescribeVolumesResponse xmlns="http://ec2.amazonaws.com/doc/2013-10-15/">
   <requestId>59dbff89-35bd-4eac-99ed-be587EXAMPLE</requestId>
   <volumeSet>
      {% for volume in volumes %}
          <item>
             <volumeId>{{ volume.id }}</volumeId>
             <size>{{ volume.size }}</size>
             {% if volume.snapshot_id %}
               <snapshotId>{{ volume.snapshot_id }}</snapshotId>
             {% else %}
               <snapshotId/>
             {% endif %}
             <encrypted>{{ 'true' if volume.encrypted else 'false' }}</encrypted>
             {% if volume.encrypted %}
             <kmsKeyId>{{ volume.kms_key_id }}</kmsKeyId>
             {% endif %}
             <availabilityZone>{{ volume.zone.name }}</availabilityZone>
             <status>{{ volume.status }}</status>
             <createTime>{{ volume.create_time }}</createTime>
             <attachmentSet>
                {% if volume.attachment %}
                    <item>
                       <volumeId>{{ volume.id }}</volumeId>
                       <instanceId>{{ volume.attachment.instance.id }}</instanceId>
                       <device>{{ volume.attachment.device }}</device>
                       <status>attached</status>
                       <attachTime>{{volume.attachment.attach_time}}</attachTime>
                       <deleteOnTermination>false</deleteOnTermination>
                    </item>
                {% endif %}
             </attachmentSet>
             {% if volume.get_tags() %}
               <tagSet>
                 {% for tag in volume.get_tags() %}
                   <item>
                     <resourceId>{{ tag.resource_id }}</resourceId>
                     <resourceType>{{ tag.resource_type }}</resourceType>
                     <key>{{ tag.key }}</key>
                     <value>{{ tag.value }}</value>
                   </item>
                 {% endfor %}
               </tagSet>
             {% endif %}
             <volumeType>{{ volume.volume_type }}</volumeType>
             {% if volume.iops %}
               <iops>{{ volume.iops }}</iops>
             {% endif %}
             {% if volume.throughput %}
               <throughput>{{ volume.throughput }}</throughput>
             {% endif %}
          </item>
      {% endfor %}
   </volumeSet>
</DescribeVolumesResponse>"""

DELETE_VOLUME_RESPONSE = """<DeleteVolumeResponse xmlns="http://ec2.amazonaws.com/doc/2013-10-15/">
  <requestId>59dbff89-35bd-4eac-99ed-be587EXAMPLE</requestId>
  <return>true</return>
</DeleteVolumeResponse>"""

ATTACHED_VOLUME_RESPONSE = """<AttachVolumeResponse xmlns="http://ec2.amazonaws.com/doc/2013-10-15/">
  <requestId>59dbff89-35bd-4eac-99ed-be587EXAMPLE</requestId>
  <volumeId>{{ attachment.volume.id }}</volumeId>
  <instanceId>{{ attachment.instance.id }}</instanceId>
  <device>{{ attachment.device }}</device>
  <status>attaching</status>
  <attachTime>{{attachment.attach_time}}</attachTime>
</AttachVolumeResponse>"""

DETATCH_VOLUME_RESPONSE = """<DetachVolumeResponse xmlns="http://ec2.amazonaws.com/doc/2013-10-15/">
   <requestId>59dbff89-35bd-4eac-99ed-be587EXAMPLE</requestId>
   <volumeId>{{ attachment.volume.id }}</volumeId>
   <instanceId>{{ attachment.instance.id }}</instanceId>
   <device>{{ attachment.device }}</device>
   <status>detaching</status>
   <attachTime>2013-10-04T17:38:53.000Z</attachTime>
</DetachVolumeResponse>"""

CREATE_SNAPSHOT_RESPONSE = """<CreateSnapshotResponse xmlns="http://ec2.amazonaws.com/doc/2013-10-15/">
  <requestId>59dbff89-35bd-4eac-99ed-be587EXAMPLE</requestId>
  <snapshotId>{{ snapshot.id }}</snapshotId>
  <volumeId>{{ snapshot.volume.id }}</volumeId>
  <status>pending</status>
  <startTime>{{ snapshot.start_time}}</startTime>
  <progress>60%</progress>
  <ownerId>{{ snapshot.owner_id }}</ownerId>
  <volumeSize>{{ snapshot.volume.size }}</volumeSize>
  <description>{{ snapshot.description }}</description>
  <encrypted>{{ 'true' if snapshot.encrypted else 'false' }}</encrypted>
  <tagSet>
    {% for tag in snapshot.get_tags() %}
      <item>
      <resourceId>{{ tag.resource_id }}</resourceId>
      <resourceType>{{ tag.resource_type }}</resourceType>
      <key>{{ tag.key }}</key>
      <value>{{ tag.value }}</value>
      </item>
    {% endfor %}
  </tagSet>
</CreateSnapshotResponse>"""

CREATE_SNAPSHOTS_RESPONSE = """<CreateSnapshotsResponse xmlns="http://ec2.amazonaws.com/doc/2013-10-15/">
  <requestId>59dbff89-35bd-4eac-99ed-be587EXAMPLE</requestId>
  <snapshotSet>
      {% for snapshot in snapshots %}
      <item>
          <snapshotId>{{ snapshot.id }}</snapshotId>
          <volumeId>{{ snapshot.volume.id }}</volumeId>
          <status>pending</status>
          <startTime>{{ snapshot.start_time}}</startTime>
          <progress>60%</progress>
          <ownerId>{{ snapshot.owner_id }}</ownerId>
          <volumeSize>{{ snapshot.volume.size }}</volumeSize>
          <description>{{ snapshot.description }}</description>
          <encrypted>{{ 'true' if snapshot.encrypted else 'false' }}</encrypted>
          <tagSet>
            {% for tag in snapshot.get_tags() %}
              <item>
              <resourceId>{{ tag.resource_id }}</resourceId>
              <resourceType>{{ tag.resource_type }}</resourceType>
              <key>{{ tag.key }}</key>
              <value>{{ tag.value }}</value>
              </item>
            {% endfor %}
          </tagSet>
      </item>
      {% endfor %}
  </snapshotSet>
</CreateSnapshotsResponse>"""

COPY_SNAPSHOT_RESPONSE = """<CopySnapshotResponse xmlns="http://ec2.amazonaws.com/doc/2016-11-15/">
  <requestId>59dbff89-35bd-4eac-99ed-be587EXAMPLE</requestId>
  <snapshotId>{{ snapshot.id }}</snapshotId>
  <tagSet>
    {% for tag in snapshot.get_tags() %}
      <item>
      <resourceId>{{ tag.resource_id }}</resourceId>
      <resourceType>{{ tag.resource_type }}</resourceType>
      <key>{{ tag.key }}</key>
      <value>{{ tag.value }}</value>
      </item>
    {% endfor %}
  </tagSet>
</CopySnapshotResponse>"""

DESCRIBE_SNAPSHOTS_RESPONSE = """<DescribeSnapshotsResponse xmlns="http://ec2.amazonaws.com/doc/2013-10-15/">
   <requestId>59dbff89-35bd-4eac-99ed-be587EXAMPLE</requestId>
   <snapshotSet>
      {% for snapshot in snapshots %}
          <item>
             <snapshotId>{{ snapshot.id }}</snapshotId>
            <volumeId>{{ snapshot.volume.id }}</volumeId>
             <status>{{ snapshot.status }}</status>
             <startTime>{{ snapshot.start_time}}</startTime>
             <progress>100%</progress>
             <ownerId>{{ snapshot.owner_id }}</ownerId>
            <volumeSize>{{ snapshot.volume.size }}</volumeSize>
             <description>{{ snapshot.description }}</description>
             <encrypted>{{ 'true' if snapshot.encrypted else 'false' }}</encrypted>
             <tagSet>
               {% for tag in snapshot.get_tags() %}
                 <item>
                   <resourceId>{{ tag.resource_id }}</resourceId>
                   <resourceType>{{ tag.resource_type }}</resourceType>
                   <key>{{ tag.key }}</key>
                   <value>{{ tag.value }}</value>
                 </item>
               {% endfor %}
             </tagSet>
          </item>
      {% endfor %}
   </snapshotSet>
</DescribeSnapshotsResponse>"""

DELETE_SNAPSHOT_RESPONSE = """<DeleteSnapshotResponse xmlns="http://ec2.amazonaws.com/doc/2013-10-15/">
  <requestId>59dbff89-35bd-4eac-99ed-be587EXAMPLE</requestId>
  <return>true</return>
</DeleteSnapshotResponse>"""

DESCRIBE_SNAPSHOT_ATTRIBUTES_RESPONSE = """
<DescribeSnapshotAttributeResponse xmlns="http://ec2.amazonaws.com/doc/2013-10-15/">
    <requestId>a9540c9f-161a-45d8-9cc1-1182b89ad69f</requestId>
    <snapshotId>snap-a0332ee0</snapshotId>
    <createVolumePermission>
       {% for group in groups %}
          <item>
             <group>{{ group }}</group>
          </item>
       {% endfor %}
       {% for userId in userIds %}
          <item>
             <userId>{{ userId }}</userId>
          </item>
       {% endfor %}
    </createVolumePermission>
</DescribeSnapshotAttributeResponse>
"""

MODIFY_SNAPSHOT_ATTRIBUTE_RESPONSE = """
<ModifySnapshotAttributeResponse xmlns="http://ec2.amazonaws.com/doc/2013-10-15/">
    <requestId>666d2944-9276-4d6a-be12-1f4ada972fd8</requestId>
    <return>true</return>
</ModifySnapshotAttributeResponse>
"""

MODIFY_VOLUME_RESPONSE = """<?xml version="1.0" encoding="UTF-8"?>
<ModifyVolumeResponse xmlns="http://ec2.amazonaws.com/doc/2016-11-15/">
    <requestId>59dbff89-35bd-4eac-99ed-be587EXAMPLE</requestId>
    <volumeModification>
        {% set volume_modification = volume.modifications[-1] %}
        <modificationState>modifying</modificationState>
        <originalSize>{{ volume_modification.original_size }}</originalSize>
        <originalVolumeType>{{ volume_modification.original_volume_type }}</originalVolumeType>
        <progress>0</progress>
        <startTime>{{ volume_modification.start_time }}</startTime>
        <targetSize>{{ volume_modification.target_size }}</targetSize>
        <targetVolumeType>{{ volume_modification.target_volume_type }}</targetVolumeType>
        <volumeId>{{ volume.id }}</volumeId>
    </volumeModification>
</ModifyVolumeResponse>"""

DESCRIBE_VOLUMES_MODIFICATIONS_RESPONSE = """
<?xml version="1.0" encoding="UTF-8"?>
<DescribeVolumesModificationsResponse xmlns="http://ec2.amazonaws.com/doc/2016-11-15/">
    <requestId>59dbff89-35bd-4eac-99ed-be587EXAMPLE</requestId>
    <volumeModificationSet>
      {% for modification in modifications %}
        <item>
            <endTime>{{ modification.end_time }}</endTime>
            <modificationState>completed</modificationState>
            <originalSize>{{ modification.original_size }}</originalSize>
            <originalVolumeType>{{ modification.original_volume_type }}</originalVolumeType>
            <progress>100</progress>
            <startTime>{{ modification.start_time }}</startTime>
            <targetSize>{{ modification.target_size }}</targetSize>
            <targetVolumeType>{{ modification.target_volume_type }}</targetVolumeType>
            <volumeId>{{ modification.volume.id }}</volumeId>
        </item>
      {% endfor %}
    </volumeModificationSet>
</DescribeVolumesModificationsResponse>"""
