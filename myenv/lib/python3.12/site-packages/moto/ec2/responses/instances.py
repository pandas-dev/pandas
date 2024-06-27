from copy import deepcopy
from typing import Any, Dict, List, Optional

from moto.core.utils import camelcase_to_underscores
from moto.ec2.exceptions import (
    InvalidParameterCombination,
    InvalidRequest,
    MissingParameterError,
)
from moto.ec2.utils import filter_iam_instance_profiles

from ._base_response import EC2BaseResponse


class InstanceResponse(EC2BaseResponse):
    def describe_instances(self) -> str:
        self.error_on_dryrun()
        # https://boto3.amazonaws.com/v1/documentation/api/latest/reference/services/ec2/client/describe_instances.html
        # You cannot specify this(MaxResults) parameter and the instance IDs parameter in the same request.
        if "InstanceId.1" in self.data and "MaxResults" in self.data:
            raise InvalidParameterCombination(
                "The parameter instancesSet cannot be used with the parameter maxResults"
            )
        filter_dict = self._filters_from_querystring()
        instance_ids = self._get_multi_param("InstanceId")
        token = self._get_param("NextToken")
        if instance_ids:
            reservations = self.ec2_backend.get_reservations_by_instance_ids(
                instance_ids, filters=filter_dict
            )
        else:
            reservations = self.ec2_backend.describe_instances(filters=filter_dict)

        reservation_ids = [reservation.id for reservation in reservations]
        if token:
            start = reservation_ids.index(token) + 1
        else:
            start = 0
        max_results = int(self._get_param("MaxResults", 100))
        reservations_resp = reservations[start : start + max_results]
        next_token = None
        if max_results and len(reservations) > (start + max_results):
            next_token = reservations_resp[-1].id
        template = self.response_template(EC2_DESCRIBE_INSTANCES)
        return (
            template.render(
                account_id=self.current_account,
                reservations=reservations_resp,
                next_token=next_token,
                run_instances=False,
            )
            .replace("True", "true")
            .replace("False", "false")
        )

    def run_instances(self) -> str:
        min_count = int(self._get_param("MinCount", if_none="1"))
        image_id = self._get_param("ImageId")
        owner_id = self._get_param("OwnerId")
        user_data = self._get_param("UserData")
        security_group_names = self._get_multi_param("SecurityGroup")
        kwargs = {
            "instance_type": self._get_param("InstanceType", if_none="m1.small"),
            "is_instance_type_default": not self._get_param("InstanceType"),
            "placement": self._get_param("Placement.AvailabilityZone"),
            "placement_hostid": self._get_param("Placement.HostId"),
            "region_name": self.region,
            "subnet_id": self._get_param("SubnetId"),
            "owner_id": owner_id,
            "key_name": self._get_param("KeyName"),
            "security_group_ids": self._get_multi_param("SecurityGroupId"),
            "nics": self._get_multi_param("NetworkInterface."),
            "private_ip": self._get_param("PrivateIpAddress"),
            "associate_public_ip": self._get_param("AssociatePublicIpAddress"),
            "tags": self._parse_tag_specification(),
            "ebs_optimized": self._get_param("EbsOptimized") or False,
            "instance_market_options": self._get_param(
                "InstanceMarketOptions.MarketType"
            )
            or {},
            "instance_initiated_shutdown_behavior": self._get_param(
                "InstanceInitiatedShutdownBehavior"
            ),
            "launch_template": self._get_multi_param_dict("LaunchTemplate"),
            "hibernation_options": self._get_multi_param_dict("HibernationOptions"),
            "iam_instance_profile_name": self._get_param("IamInstanceProfile.Name")
            or None,
            "iam_instance_profile_arn": self._get_param("IamInstanceProfile.Arn")
            or None,
            "monitoring_state": "enabled"
            if self._get_param("Monitoring.Enabled") == "true"
            else "disabled",
        }
        if len(kwargs["nics"]) and kwargs["subnet_id"]:
            raise InvalidParameterCombination(
                msg="Network interfaces and an instance-level subnet ID may not be specified on the same request"
            )

        mappings = self._parse_block_device_mapping()
        if mappings:
            kwargs["block_device_mappings"] = mappings

        iam_instance_profile_name = kwargs.get("iam_instance_profile_name")
        iam_instance_profile_arn = kwargs.get("iam_instance_profile_arn")
        if iam_instance_profile_arn or iam_instance_profile_name:
            # Validate the profile exists, before we error_on_dryrun and run_instances
            filter_iam_instance_profiles(
                self.current_account,
                partition=self.partition,
                iam_instance_profile_arn=iam_instance_profile_arn,
                iam_instance_profile_name=iam_instance_profile_name,
            )

        self.error_on_dryrun()

        new_reservation = self.ec2_backend.run_instances(
            image_id, min_count, user_data, security_group_names, **kwargs
        )
        if iam_instance_profile_name:
            self.ec2_backend.associate_iam_instance_profile(
                instance_id=new_reservation.instances[0].id,
                iam_instance_profile_name=iam_instance_profile_name,
            )

        if iam_instance_profile_arn:
            self.ec2_backend.associate_iam_instance_profile(
                instance_id=new_reservation.instances[0].id,
                iam_instance_profile_arn=iam_instance_profile_arn,
            )

        template = self.response_template(EC2_RUN_INSTANCES)
        return template.render(
            account_id=self.current_account,
            reservation=new_reservation,
            run_instances=True,
        )

    def terminate_instances(self) -> str:
        instance_ids = self._get_multi_param("InstanceId")

        self.error_on_dryrun()

        instances = self.ec2_backend.terminate_instances(instance_ids)
        from moto.autoscaling import autoscaling_backends
        from moto.elbv2 import elbv2_backends

        autoscaling_backends[self.current_account][
            self.region
        ].notify_terminate_instances(instance_ids)
        elbv2_backends[self.current_account][self.region].notify_terminate_instances(
            instance_ids
        )
        template = self.response_template(EC2_TERMINATE_INSTANCES)
        return template.render(instances=instances)

    def reboot_instances(self) -> str:
        instance_ids = self._get_multi_param("InstanceId")

        self.error_on_dryrun()

        instances = self.ec2_backend.reboot_instances(instance_ids)
        template = self.response_template(EC2_REBOOT_INSTANCES)
        return template.render(instances=instances)

    def stop_instances(self) -> str:
        instance_ids = self._get_multi_param("InstanceId")

        self.error_on_dryrun()

        instances = self.ec2_backend.stop_instances(instance_ids)
        template = self.response_template(EC2_STOP_INSTANCES)
        return template.render(instances=instances)

    def start_instances(self) -> str:
        instance_ids = self._get_multi_param("InstanceId")
        self.error_on_dryrun()

        instances = self.ec2_backend.start_instances(instance_ids)
        template = self.response_template(EC2_START_INSTANCES)
        return template.render(instances=instances)

    def _get_list_of_dict_params(
        self, param_prefix: str, _dct: Dict[str, Any]
    ) -> List[Any]:
        """
        Simplified version of _get_dict_param
        Allows you to pass in a custom dict instead of using self.querystring by default
        """
        params = []
        for key, value in _dct.items():
            if key.startswith(param_prefix):
                params.append(value)
        return params

    def describe_instance_status(self) -> str:
        instance_ids = self._get_multi_param("InstanceId")
        include_all_instances = self._get_param("IncludeAllInstances") == "true"
        filters = self._get_list_prefix("Filter")
        filters = [
            {"name": f["name"], "values": self._get_list_of_dict_params("value.", f)}
            for f in filters
        ]

        instances = self.ec2_backend.describe_instance_status(
            instance_ids, include_all_instances, filters
        )

        template = self.response_template(EC2_INSTANCE_STATUS)
        return template.render(instances=instances)

    def describe_instance_types(self) -> str:
        instance_type_filters = self._get_multi_param("InstanceType")
        filter_dict = self._filters_from_querystring()
        instance_types = self.ec2_backend.describe_instance_types(
            instance_type_filters, filter_dict
        )
        template = self.response_template(EC2_DESCRIBE_INSTANCE_TYPES)
        return template.render(instance_types=instance_types)

    def describe_instance_type_offerings(self) -> str:
        location_type_filters = self._get_param("LocationType")
        filter_dict = self._filters_from_querystring()
        offerings = self.ec2_backend.describe_instance_type_offerings(
            location_type_filters, filter_dict
        )
        template = self.response_template(EC2_DESCRIBE_INSTANCE_TYPE_OFFERINGS)
        return template.render(instance_type_offerings=offerings)

    def describe_instance_attribute(self) -> str:
        # TODO this and modify below should raise IncorrectInstanceState if
        # instance not in stopped state
        attribute = self._get_param("Attribute")
        instance_id = self._get_param("InstanceId")
        instance, value = self.ec2_backend.describe_instance_attribute(
            instance_id, attribute
        )

        if attribute == "groupSet":
            template = self.response_template(EC2_DESCRIBE_INSTANCE_GROUPSET_ATTRIBUTE)
        else:
            template = self.response_template(EC2_DESCRIBE_INSTANCE_ATTRIBUTE)

        return template.render(instance=instance, attribute=attribute, value=value)

    def describe_instance_credit_specifications(self) -> str:
        instance_ids = self._get_multi_param("InstanceId")
        instance = self.ec2_backend.describe_instance_credit_specifications(
            instance_ids
        )
        template = self.response_template(EC2_DESCRIBE_INSTANCE_CREDIT_SPECIFICATIONS)
        return template.render(instances=instance)

    def modify_instance_attribute(self) -> str:
        handlers = [
            self._attribute_value_handler,
            self._dot_value_instance_attribute_handler,
            self._block_device_mapping_handler,
            self._security_grp_instance_attribute_handler,
        ]

        for handler in handlers:
            success = handler()
            if success:
                return success

        msg = (
            "This specific call to ModifyInstanceAttribute has not been"
            " implemented in Moto yet. Feel free to open an issue at"
            " https://github.com/getmoto/moto/issues"
        )
        raise NotImplementedError(msg)

    def _block_device_mapping_handler(self) -> Optional[str]:
        """
        Handles requests which are generated by code similar to:

            instance.modify_attribute(
                BlockDeviceMappings=[{
                    'DeviceName': '/dev/sda1',
                    'Ebs': {'DeleteOnTermination': True}
                }]
            )

        The querystring contains information similar to:

            BlockDeviceMapping.1.Ebs.DeleteOnTermination : ['true']
            BlockDeviceMapping.1.DeviceName : ['/dev/sda1']

        For now we only support the "BlockDeviceMapping.1.Ebs.DeleteOnTermination"
        configuration, but it should be trivial to add anything else.
        """
        mapping_counter = 1
        mapping_device_name_fmt = "BlockDeviceMapping.%s.DeviceName"
        mapping_del_on_term_fmt = "BlockDeviceMapping.%s.Ebs.DeleteOnTermination"
        while True:
            mapping_device_name = mapping_device_name_fmt % mapping_counter
            if mapping_device_name not in self.querystring.keys():
                break

            mapping_del_on_term = mapping_del_on_term_fmt % mapping_counter
            del_on_term_value_str = self.querystring[mapping_del_on_term][0]
            del_on_term_value = True if "true" == del_on_term_value_str else False
            device_name_value = self.querystring[mapping_device_name][0]

            instance_id = self._get_param("InstanceId")
            instance = self.ec2_backend.get_instance(instance_id)

            self.error_on_dryrun()

            block_device_type = instance.block_device_mapping[device_name_value]
            block_device_type.delete_on_termination = del_on_term_value

            # +1 for the next device
            mapping_counter += 1

        if mapping_counter > 1:
            return EC2_MODIFY_INSTANCE_ATTRIBUTE
        return None

    def _dot_value_instance_attribute_handler(self) -> Optional[str]:
        attribute_key = None
        for key, value in self.querystring.items():
            if ".Value" in key:
                attribute_key = key
                break

        if not attribute_key:
            return None

        self.error_on_dryrun()

        value = self.querystring.get(attribute_key)[0]  # type: ignore
        normalized_attribute = camelcase_to_underscores(attribute_key.split(".")[0])
        instance_id = self._get_param("InstanceId")
        self.ec2_backend.modify_instance_attribute(
            instance_id, normalized_attribute, value
        )
        return EC2_MODIFY_INSTANCE_ATTRIBUTE

    def _attribute_value_handler(self) -> Optional[str]:
        attribute_key = self._get_param("Attribute")

        if attribute_key is None:
            return None

        self.error_on_dryrun()

        value = self._get_param("Value")
        normalized_attribute = camelcase_to_underscores(attribute_key)
        instance_id = self._get_param("InstanceId")
        self.ec2_backend.modify_instance_attribute(
            instance_id, normalized_attribute, value
        )
        return EC2_MODIFY_INSTANCE_ATTRIBUTE

    def _security_grp_instance_attribute_handler(self) -> str:
        new_security_grp_list = []
        for key in self.querystring:
            if "GroupId." in key:
                new_security_grp_list.append(self.querystring.get(key)[0])  # type: ignore

        instance_id = self._get_param("InstanceId")
        self.error_on_dryrun()

        self.ec2_backend.modify_instance_security_groups(
            instance_id, new_security_grp_list
        )
        return EC2_MODIFY_INSTANCE_ATTRIBUTE

    def _parse_block_device_mapping(self) -> List[Dict[str, Any]]:
        device_mappings = self._get_list_prefix("BlockDeviceMapping")
        mappings = []
        for device_mapping in device_mappings:
            self._validate_block_device_mapping(device_mapping)
            device_template: Dict[str, Any] = deepcopy(BLOCK_DEVICE_MAPPING_TEMPLATE)
            device_template["VirtualName"] = device_mapping.get("virtual_name")
            device_template["DeviceName"] = device_mapping.get("device_name")
            device_template["Ebs"]["SnapshotId"] = device_mapping.get(
                "ebs._snapshot_id"
            )
            device_template["Ebs"]["VolumeSize"] = device_mapping.get(
                "ebs._volume_size"
            )
            device_template["Ebs"]["DeleteOnTermination"] = self._convert_to_bool(
                device_mapping.get("ebs._delete_on_termination", False)
            )
            device_template["Ebs"]["VolumeType"] = device_mapping.get(
                "ebs._volume_type"
            )
            device_template["Ebs"]["Iops"] = device_mapping.get("ebs._iops")
            device_template["Ebs"]["Encrypted"] = self._convert_to_bool(
                device_mapping.get("ebs._encrypted", False)
            )
            device_template["Ebs"]["KmsKeyId"] = device_mapping.get("ebs._kms_key_id")
            device_template["NoDevice"] = device_mapping.get("no_device")
            mappings.append(device_template)

        return mappings

    @staticmethod
    def _validate_block_device_mapping(device_mapping: Dict[str, Any]) -> None:  # type: ignore[misc]
        from botocore import __version__ as botocore_version

        if "no_device" in device_mapping:
            assert isinstance(
                device_mapping["no_device"], str
            ), f"botocore {botocore_version} isn't limiting NoDevice to str type anymore, it is type:{type(device_mapping['no_device'])}"
            if device_mapping["no_device"] == "":
                # the only legit value it can have is empty string
                # and none of the other checks here matter if NoDevice
                # is being used
                return
            else:
                raise InvalidRequest()

        if not any(mapping for mapping in device_mapping if mapping.startswith("ebs.")):
            raise MissingParameterError("ebs")
        if (
            "ebs._volume_size" not in device_mapping
            and "ebs._snapshot_id" not in device_mapping
        ):
            raise MissingParameterError("size or snapshotId")

    @staticmethod
    def _convert_to_bool(bool_str: Any) -> bool:  # type: ignore[misc]
        if isinstance(bool_str, bool):
            return bool_str

        if isinstance(bool_str, str):
            return str(bool_str).lower() == "true"

        return False


BLOCK_DEVICE_MAPPING_TEMPLATE = {
    "VirtualName": None,
    "DeviceName": None,
    "NoDevice": None,
    "Ebs": {
        "SnapshotId": None,
        "VolumeSize": None,
        "DeleteOnTermination": None,
        "VolumeType": None,
        "Iops": None,
        "Encrypted": None,
    },
}

INSTANCE_TEMPLATE = """<item>
          <instanceId>{{ instance.id }}</instanceId>
          <imageId>{{ instance.image_id }}</imageId>
          {% if run_instances %}
          <instanceState>
            <code>0</code>
            <name>pending</name>
         </instanceState>
          {% else %}
          <instanceState>
            <code>{{ instance._state.code }}</code>
            <name>{{ instance._state.name }}</name>
         </instanceState>
         {% endif %}
          <privateDnsName>{{ instance.private_dns }}</privateDnsName>
          <publicDnsName>{{ instance.public_dns }}</publicDnsName>
          <dnsName>{{ instance.public_dns }}</dnsName>
          <reason>{{ instance._reason }}</reason>
          {% if instance.key_name is not none %}
             <keyName>{{ instance.key_name }}</keyName>
          {% endif %}
          <ebsOptimized>{{ instance.ebs_optimized }}</ebsOptimized>
          <amiLaunchIndex>{{ instance.ami_launch_index }}</amiLaunchIndex>
          <instanceType>{{ instance.instance_type }}</instanceType>
          {% if instance.iam_instance_profile %}
          <iamInstanceProfile>
            <arn>{{ instance.iam_instance_profile['Arn'] }}</arn>
            <id>{{ instance.iam_instance_profile['Id'] }}</id>
          </iamInstanceProfile>
          {% endif %}
          <launchTime>{{ instance.launch_time }}</launchTime>
          {% if instance.lifecycle %}
          <instanceLifecycle>{{ instance.lifecycle }}</instanceLifecycle>
          {% endif %}
          <placement>
            {% if instance.placement_hostid %}<hostId>{{ instance.placement_hostid }}</hostId>{% endif %}
            <availabilityZone>{{ instance.placement}}</availabilityZone>
            <groupName/>
            <tenancy>default</tenancy>
          </placement>
          <monitoring>
            <state> {{ instance.monitoring_state }} </state>
          </monitoring>
          {% if instance.subnet_id %}
            <subnetId>{{ instance.subnet_id }}</subnetId>
          {% elif instance.nics[0].subnet.id %}
            <subnetId>{{ instance.nics[0].subnet.id }}</subnetId>
          {% endif %}
          {% if instance.vpc_id %}
            <vpcId>{{ instance.vpc_id }}</vpcId>
          {% elif instance.nics[0].subnet.vpc_id %}
            <vpcId>{{ instance.nics[0].subnet.vpc_id }}</vpcId>
          {% endif %}
          <privateIpAddress>{{ instance.private_ip }}</privateIpAddress>
          {% if instance.nics[0].public_ip %}
              <ipAddress>{{ instance.nics[0].public_ip }}</ipAddress>
          {% endif %}
          <sourceDestCheck>{{ instance.source_dest_check }}</sourceDestCheck>
          <groupSet>
             {% for group in instance.dynamic_group_list %}
             <item>
               {% if group.id %}
                <groupId>{{ group.id }}</groupId>
                <groupName>{{ group.name }}</groupName>
                {% else %}
                  <groupId>{{ group }}</groupId>
                {% endif %}
             </item>
             {% endfor %}
          </groupSet>
          {% if instance.platform %}
          <platform>{{ instance.platform }}</platform>
          {% endif %}
          <virtualizationType>{{ instance.virtualization_type }}</virtualizationType>
          <stateReason>
              <code>{{ instance._state_reason.code }}</code>
              <message>{{ instance._state_reason.message }}</message>
          </stateReason>
          <architecture>{{ instance.architecture }}</architecture>
          <kernelId>{{ instance.kernel }}</kernelId>
          <rootDeviceType>ebs</rootDeviceType>
          <rootDeviceName>/dev/sda1</rootDeviceName>
          <blockDeviceMapping>
              {% for device_name,deviceobject in instance.get_block_device_mapping %}
              <item>
                 <deviceName>{{ device_name }}</deviceName>
                  <ebs>
                     <volumeId>{{ deviceobject.volume_id }}</volumeId>
                     <status>{{ deviceobject.status }}</status>
                     <attachTime>{{ deviceobject.attach_time }}</attachTime>
                     <deleteOnTermination>{{ deviceobject.delete_on_termination }}</deleteOnTermination>
                     <size>{{deviceobject.size}}</size>
                </ebs>
              </item>
             {% endfor %}
          </blockDeviceMapping>
          <clientToken>ABCDE{{ account_id }}3</clientToken>
          <hypervisor>xen</hypervisor>
          {% if instance.hibernation_options %}
          <hibernationOptions>
            <configured>{{ instance.hibernation_options.get("Configured") }}</configured>
          </hibernationOptions>
          {% endif %}
          {% if instance.get_tags() %}
          <tagSet>
            {% for tag in instance.get_tags() %}
              <item>
                <resourceId>{{ tag.resource_id }}</resourceId>
                <resourceType>{{ tag.resource_type }}</resourceType>
                <key>{{ tag.key }}</key>
                <value>{{ tag.value }}</value>
              </item>
            {% endfor %}
          </tagSet>
          {% endif %}
          <networkInterfaceSet>
            {% for nic in instance.nics.values() %}
              <item>
                <networkInterfaceId>{{ nic.id }}</networkInterfaceId>
                {% if nic.subnet %}
                  <subnetId>{{ nic.subnet.id }}</subnetId>
                  <vpcId>{{ nic.subnet.vpc_id }}</vpcId>
                {% endif %}
                <description>Primary network interface</description>
                <ownerId>{{ account_id }}</ownerId>
                <status>in-use</status>
                <macAddress>1b:2b:3c:4d:5e:6f</macAddress>
                <privateIpAddress>{{ nic.private_ip_address }}</privateIpAddress>
                <sourceDestCheck>{{ instance.source_dest_check }}</sourceDestCheck>
                <groupSet>
                  {% for group in nic.group_set %}
                  <item>
                    {% if group.id %}
                      <groupId>{{ group.id }}</groupId>
                      <groupName>{{ group.name }}</groupName>
                    {% else %}
                      <groupId>{{ group }}</groupId>
                    {% endif %}
                  </item>
                  {% endfor %}
                </groupSet>
                <attachment>
                  <attachmentId>{{ nic.attachment_id }}</attachmentId>
                  <deviceIndex>{{ nic.device_index }}</deviceIndex>
                  <status>attached</status>
                  <attachTime>2015-01-01T00:00:00Z</attachTime>
                  <deleteOnTermination>true</deleteOnTermination>
                </attachment>
                {% if nic.public_ip %}
                  <association>
                    <publicIp>{{ nic.public_ip }}</publicIp>
                    <ipOwnerId>{{ account_id }}</ipOwnerId>
                  </association>
                {% endif %}
                <privateIpAddressesSet>
                  <item>
                    <privateIpAddress>{{ nic.private_ip_address }}</privateIpAddress>
                    <primary>true</primary>
                    {% if nic.public_ip %}
                      <association>
                        <publicIp>{{ nic.public_ip }}</publicIp>
                        <ipOwnerId>{{ account_id }}</ipOwnerId>
                      </association>
                    {% endif %}
                  </item>
                </privateIpAddressesSet>
              </item>
            {% endfor %}
          </networkInterfaceSet>
        </item>"""

EC2_RUN_INSTANCES = (
    """<RunInstancesResponse xmlns="http://ec2.amazonaws.com/doc/2013-10-15/">
  <requestId>59dbff89-35bd-4eac-99ed-be587EXAMPLE</requestId>
  <reservationId>{{ reservation.id }}</reservationId>
  <ownerId>{{ account_id }}</ownerId>
  <groupSet>
    <item>
      <groupId>sg-245f6a01</groupId>
      <groupName>default</groupName>
    </item>
  </groupSet>
  <instancesSet>
    {% for instance in reservation.instances %}
        """
    + INSTANCE_TEMPLATE
    + """
    {% endfor %}
  </instancesSet>
  </RunInstancesResponse>"""
)

EC2_DESCRIBE_INSTANCES = (
    """<DescribeInstancesResponse xmlns="http://ec2.amazonaws.com/doc/2013-10-15/">
  <requestId>fdcdcab1-ae5c-489e-9c33-4637c5dda355</requestId>
      <reservationSet>
        {% for reservation in reservations %}
          <item>
            <reservationId>{{ reservation.id }}</reservationId>
            <ownerId>{{ account_id }}</ownerId>
            <groupSet>
              {% for group in reservation.dynamic_group_list %}
              <item>
      {% if group.id %}
                <groupId>{{ group.id }}</groupId>
                <groupName>{{ group.name }}</groupName>
                {% else %}
                <groupId>{{ group }}</groupId>
                {% endif %}
              </item>
              {% endfor %}
            </groupSet>
            <instancesSet>
                {% for instance in reservation.instances %}
                  """
    + INSTANCE_TEMPLATE
    + """
                {% endfor %}
            </instancesSet>
          </item>
        {% endfor %}
      </reservationSet>
      {% if next_token %}
      <nextToken>{{ next_token }}</nextToken>
      {% endif %}
</DescribeInstancesResponse>"""
)

EC2_TERMINATE_INSTANCES = """
<TerminateInstancesResponse xmlns="http://ec2.amazonaws.com/doc/2013-10-15/">
  <requestId>59dbff89-35bd-4eac-99ed-be587EXAMPLE</requestId>
  <instancesSet>
    {% for instance, previous_state in instances %}
      <item>
        <instanceId>{{ instance.id }}</instanceId>
        <previousState>
          <code>{{ previous_state.code }}</code>
          <name>{{ previous_state.name }}</name>
        </previousState>
        <currentState>
          <code>32</code>
          <name>shutting-down</name>
        </currentState>
      </item>
    {% endfor %}
  </instancesSet>
</TerminateInstancesResponse>"""

EC2_STOP_INSTANCES = """
<StopInstancesResponse xmlns="http://ec2.amazonaws.com/doc/2013-10-15/">
  <requestId>59dbff89-35bd-4eac-99ed-be587EXAMPLE</requestId>
  <instancesSet>
    {% for instance, previous_state in instances %}
      <item>
        <instanceId>{{ instance.id }}</instanceId>
        <previousState>
          <code>{{ previous_state.code }}</code>
          <name>{{ previous_state.name }}</name>
        </previousState>
        <currentState>
          <code>64</code>
          <name>stopping</name>
        </currentState>
      </item>
    {% endfor %}
  </instancesSet>
</StopInstancesResponse>"""

EC2_START_INSTANCES = """
<StartInstancesResponse xmlns="http://ec2.amazonaws.com/doc/2013-10-15/">
  <requestId>59dbff89-35bd-4eac-99ed-be587EXAMPLE</requestId>
  <instancesSet>
    {% for instance, previous_state in instances %}
      <item>
        <instanceId>{{ instance.id }}</instanceId>
        <previousState>
          <code>{{ previous_state.code }}</code>
          <name>{{ previous_state.name }}</name>
        </previousState>
        <currentState>
          <code>0</code>
          <name>pending</name>
        </currentState>
      </item>
    {% endfor %}
  </instancesSet>
</StartInstancesResponse>"""

EC2_REBOOT_INSTANCES = """<RebootInstancesResponse xmlns="http://ec2.amazonaws.com/doc/2013-10-15/">
  <requestId>59dbff89-35bd-4eac-99ed-be587EXAMPLE</requestId>
  <return>true</return>
</RebootInstancesResponse>"""

EC2_DESCRIBE_INSTANCE_ATTRIBUTE = """<DescribeInstanceAttributeResponse xmlns="http://ec2.amazonaws.com/doc/2013-10-15/">
  <requestId>59dbff89-35bd-4eac-99ed-be587EXAMPLE</requestId>
  <instanceId>{{ instance.id }}</instanceId>
  <{{ attribute }}>
    {% if value is not none %}
    <value>{{ value }}</value>
    {% endif %}
  </{{ attribute }}>
</DescribeInstanceAttributeResponse>"""

EC2_DESCRIBE_INSTANCE_CREDIT_SPECIFICATIONS = """<DescribeInstanceCreditSpecificationsResponse xmlns="http://ec2.amazonaws.com/doc/2016-11-15/">
    <requestId>1b234b5c-d6ef-7gh8-90i1-j2345678901</requestId>
    <instanceCreditSpecificationSet>
       {% for instance in instances %}
      <item>
        <instanceId>{{ instance.id }}</instanceId>
        <cpuCredits>standard</cpuCredits>
      </item>
    {% endfor %}
    </instanceCreditSpecificationSet>
</DescribeInstanceCreditSpecificationsResponse>"""

EC2_DESCRIBE_INSTANCE_GROUPSET_ATTRIBUTE = """<DescribeInstanceAttributeResponse xmlns="http://ec2.amazonaws.com/doc/2013-10-15/">
  <requestId>59dbff89-35bd-4eac-99ed-be587EXAMPLE</requestId>
  <instanceId>{{ instance.id }}</instanceId>
  <{{ attribute }}>
    {% for sg in value %}
      <item>
        <groupId>{{ sg.id }}</groupId>
      </item>
    {% endfor %}
  </{{ attribute }}>
</DescribeInstanceAttributeResponse>"""

EC2_MODIFY_INSTANCE_ATTRIBUTE = """<ModifyInstanceAttributeResponse xmlns="http://ec2.amazonaws.com/doc/2013-10-15/">
  <requestId>59dbff89-35bd-4eac-99ed-be587EXAMPLE</requestId>
  <return>true</return>
</ModifyInstanceAttributeResponse>"""

EC2_INSTANCE_STATUS = """<?xml version="1.0" encoding="UTF-8"?>
<DescribeInstanceStatusResponse xmlns="http://ec2.amazonaws.com/doc/2013-10-15/">
    <requestId>59dbff89-35bd-4eac-99ed-be587EXAMPLE</requestId>
    <instanceStatusSet>
      {% for instance in instances %}
        <item>
            <instanceId>{{ instance.id }}</instanceId>
            <availabilityZone>{{ instance.placement }}</availabilityZone>
            <instanceState>
                <code>{{ instance.state_code }}</code>
                <name>{{ instance.state }}</name>
            </instanceState>
            {% if instance.state_code == 16 %}
              <systemStatus>
                  <status>ok</status>
                  <details>
                      <item>
                          <name>reachability</name>
                          <status>passed</status>
                      </item>
                  </details>
              </systemStatus>
              <instanceStatus>
                  <status>ok</status>
                  <details>
                      <item>
                          <name>reachability</name>
                          <status>passed</status>
                      </item>
                  </details>
              </instanceStatus>
            {% else %}
              <systemStatus>
                  <status>not-applicable</status>
              </systemStatus>
              <instanceStatus>
                  <status>not-applicable</status>
              </instanceStatus>
            {% endif %}
        </item>
      {% endfor %}
    </instanceStatusSet>
</DescribeInstanceStatusResponse>"""

EC2_DESCRIBE_INSTANCE_TYPES = """<?xml version="1.0" encoding="UTF-8"?>
<DescribeInstanceTypesResponse xmlns="http://ec2.amazonaws.com/doc/2016-11-15/">
    <requestId>f8b86168-d034-4e65-b48d-3b84c78e64af</requestId>
    <instanceTypeSet>
    {% for instance_type in instance_types %}
        <item>
            <autoRecoverySupported>{{ instance_type.AutoRecoverySupported }}</autoRecoverySupported>
            <bareMetal>{{ instance_type.BareMetal }}</bareMetal>
            <burstablePerformanceSupported>{{ instance_type.BurstablePerformanceSupported }}</burstablePerformanceSupported>
            <currentGeneration>{{ instance_type.CurrentGeneration }}</currentGeneration>
            <dedicatedHostsSupported>{{ instance_type.DedicatedHostsSupported }}</dedicatedHostsSupported>
            <ebsInfo>
                <ebsOptimizedInfo>
                    <baselineBandwidthInMbps>{{ instance_type.get('EbsInfo', {}).get('EbsOptimizedInfo', {}).get('BaselineBandwidthInMbps', 0) | int }}</baselineBandwidthInMbps>
                    <baselineIops>{{  instance_type.get('EbsInfo', {}).get('EbsOptimizedInfo', {}).get('BaselineIops', 0) | int }}</baselineIops>
                    <baselineThroughputInMBps>{{  instance_type.get('EbsInfo', {}).get('EbsOptimizedInfo', {}).get('BaselineThroughputInMBps', 0.0) | float }}</baselineThroughputInMBps>
                    <maximumBandwidthInMbps>{{  instance_type.get('EbsInfo', {}).get('EbsOptimizedInfo', {}).get('MaximumBandwidthInMbps', 0) | int }}</maximumBandwidthInMbps>
                    <maximumIops>{{  instance_type.get('EbsInfo', {}).get('EbsOptimizedInfo', {}).get('MaximumIops', 0) | int }}</maximumIops>
                    <maximumThroughputInMBps>{{ instance_type.get('EbsInfo', {}).get('EbsOptimizedInfo', {}).get('MaximumThroughputInMBps', 0.0) | float }}</maximumThroughputInMBps>
                </ebsOptimizedInfo>
                <ebsOptimizedSupport>{{ instance_type.get('EbsInfo', {}).get('EbsOptimizedSupport', 'default') }}</ebsOptimizedSupport>
                <encryptionSupport>{{ instance_type.get('EbsInfo', {}).get('EncryptionSupport', 'supported') }}</encryptionSupport>
                <nvmeSupport>{{ instance_type.get('EbsInfo', {}).get('NvmeSupport', 'required') }}</nvmeSupport>
            </ebsInfo>
            <networkInfo>
                <defaultNetworkCardIndex>{{ instance_type.get('NetworkInfo', {}).get('DefaultNetworkCardIndex', 0) | int }}</defaultNetworkCardIndex>
                <efaSupported>{{ instance_type.get('NetworkInfo', {}).get('EfaSupported', False) }}</efaSupported>
                <enaSrdSupported>{{ instance_type.get('NetworkInfo', {}).get('EnaSrdSupported', False) }}</enaSrdSupported>
                <enaSupport>{{ instance_type.get('NetworkInfo', {}).get('EnaSupport', 'unsupported') }}</enaSupport>
                <encryptionInTransitSupported>{{ instance_type.get('NetworkInfo', {}).get('EncryptionInTransitSupported', False) }}</encryptionInTransitSupported>
                <ipv4AddressesPerInterface>{{ instance_type.get('NetworkInfo', {}).get('Ipv4AddressesPerInterface', 0) | int }}</ipv4AddressesPerInterface>
                <ipv6AddressesPerInterface>{{ instance_type.get('NetworkInfo', {}).get('Ipv6AddressesPerInterface', 0) | int }}</ipv6AddressesPerInterface>
                <ipv6Supported>{{ instance_type.get('NetworkInfo', {}).get('Ipv6Supported', False) }}</ipv6Supported>
                <maximumNetworkCards>{{ instance_type.get('NetworkInfo', {}).get('MaximumNetworkCards', 0) | int }}</maximumNetworkCards>
                <maximumNetworkInterfaces>{{ instance_type.get('NetworkInfo', {}).get('MaximumNetworkInterfaces', 0) | int }}</maximumNetworkInterfaces>
                <networkCards>
                  {% for network_card in instance_type.get('NetworkInfo', {}).get('NetworkCards', []) %}
                    <item>
                        <baselineBandwidthInGbps>{{ network_card.get('BaselineBandwidthInGbps', 0.0) | float }}</baselineBandwidthInGbps>
                        <maximumNetworkInterfaces>{{ network_card.get('MaximumNetworkInterfaces', 0) | int }}</maximumNetworkInterfaces>
                        <networkCardIndex>{{ network_card.get('NetworkCardIndex', 0) | int }}</networkCardIndex>
                        <networkPerformance>{{ network_card.get('NetworkPerformance', 'Up to 25 Schmeckles') }}</networkPerformance>
                        <peakBandwidthInGbps>{{ network_card.get('PeakBandwidthInGbps', 0.0) | float }}</peakBandwidthInGbps>
                    </item>
                  {% endfor %}
                </networkCards>
                <networkPerformance>{{ instance_type.get('NetworkInfo', {}).get('NetworkPerformance', 'Up to 25 Schmeckles') }}</networkPerformance>
            </networkInfo>
            <freeTierEligible>{{ instance_type.FreeTierEligible }}</freeTierEligible>
            <hibernationSupported>{{ instance_type.HibernationSupported }}</hibernationSupported>
            <hypervisor>{{ instance_type.get('Hypervisor', 'motovisor') }}</hypervisor>
            <instanceStorageSupported>{{ instance_type.InstanceStorageSupported }}</instanceStorageSupported>
            <placementGroupInfo>
                <supportedStrategies>
                  {% for strategy in instance_type.get('PlacementGroupInfo', {}).get('SupportedStrategies', []) %}
                    <item>{{ strategy }}</item>
                  {% endfor %}
                </supportedStrategies>
            </placementGroupInfo>
            <supportedRootDeviceTypes>
              {% for dev_type in instance_type.get('SupportedRootDeviceTypes', []) %}
                <item>{{ dev_type }}</item>
              {% endfor %}
            </supportedRootDeviceTypes>
            <supportedUsageClasses>
              {% for usage_class in instance_type.get('SupportedUsageClasses', []) %}
                <item>{{ usage_class }}</item>
              {% endfor %}
            </supportedUsageClasses>
            <supportedVirtualizationTypes>
              {% for supported_vtype in instance_type.get('SupportedVirtualizationTypes', []) %}
                <item>{{ supported_vtype }}</item>
              {% endfor %}
            </supportedVirtualizationTypes>
            <instanceType>{{ instance_type.InstanceType }}</instanceType>
            <vCpuInfo>
                <defaultVCpus>{{ instance_type.get('VCpuInfo', {}).get('DefaultVCpus', 0)|int }}</defaultVCpus>
                <defaultCores>{{ instance_type.get('VCpuInfo', {}).get('DefaultCores', 0)|int }}</defaultCores>
                <defaultThreadsPerCore>{{ instance_type.get('VCpuInfo').get('DefaultThreadsPerCore', 0)|int }}</defaultThreadsPerCore>
                <validCores>
                  {% for valid_core in instance_type.get("VCpuInfo", {}).get('ValidCores', []) %}
                    <item>{{ valid_core }}</item>
                  {% endfor %}
                </validCores>
                <validThreadsPerCore>
                  {% for threads_per_core in instance_type.get("VCpuInfo", {}).get('ValidThreadsPerCore', []) %}
                    <item>{{ threads_per_core }}</item>
                  {% endfor %}
                </validThreadsPerCore>
            </vCpuInfo>
            <memoryInfo>
                <sizeInMiB>{{ instance_type.get('MemoryInfo', {}).get('SizeInMiB', 0)|int }}</sizeInMiB>
            </memoryInfo>
            <instanceStorageInfo>
                <totalSizeInGB>{{ instance_type.get('InstanceStorageInfo', {}).get('TotalSizeInGB', 0)|int }}</totalSizeInGB>
            </instanceStorageInfo>
            <processorInfo>
                <supportedArchitectures>
                    {% for arch in instance_type.get('ProcessorInfo', {}).get('SupportedArchitectures', []) %}
                    <item>
                        {{ arch }}
                    </item>
                    {% endfor %}
                </supportedArchitectures>
                <sustainedClockSpeedInGhz>{{ instance_type.get('ProcessorInfo', {}).get('SustainedClockSpeedInGhz', 0.0) | float }}</sustainedClockSpeedInGhz>
            </processorInfo>
            {% if instance_type.get('GpuInfo', {})|length > 0 %}
            <gpuInfo>
                <gpus>
                    {% for gpu in instance_type.get('GpuInfo').get('Gpus') %}
                    <item>
                        <count>{{ gpu['Count']|int }}</count>
                        <manufacturer>{{ gpu['Manufacturer'] }}</manufacturer>
                        <memoryInfo>
                            <sizeInMiB>{{ gpu['MemoryInfo']['SizeInMiB']|int }}</sizeInMiB>
                        </memoryInfo>
                        <name>{{ gpu['Name'] }}</name>
                    </item>
                    {% endfor %}
                </gpus>
                <totalGpuMemoryInMiB>{{ instance_type['GpuInfo']['TotalGpuMemoryInMiB']|int }}</totalGpuMemoryInMiB>
            </gpuInfo>
            {% endif %}
        </item>
    {% endfor %}
    </instanceTypeSet>
</DescribeInstanceTypesResponse>"""


EC2_DESCRIBE_INSTANCE_TYPE_OFFERINGS = """<?xml version="1.0" encoding="UTF-8"?>
<DescribeInstanceTypeOfferingsResponse xmlns="http://api.outscale.com/wsdl/fcuext/2014-04-15/">
    <requestId>f8b86168-d034-4e65-b48d-3b84c78e64af</requestId>
    <instanceTypeOfferingSet>
    {% for offering in instance_type_offerings %}
        <item>
            <instanceType>{{ offering.InstanceType }}</instanceType>
            <location>{{ offering.Location }}</location>
            <locationType>{{ offering.LocationType }}</locationType>
        </item>
    {% endfor %}
    </instanceTypeOfferingSet>
</DescribeInstanceTypeOfferingsResponse>"""
