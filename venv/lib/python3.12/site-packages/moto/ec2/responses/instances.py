from copy import deepcopy
from typing import Any

from moto.core.responses import ActionResult, EmptyResult
from moto.core.types import Base64EncodedString
from moto.core.utils import camelcase_to_underscores
from moto.ec2.exceptions import (
    InvalidParameterCombination,
    InvalidRequest,
    MissingParameterError,
)
from moto.ec2.utils import filter_iam_instance_profiles, parse_user_data

from ._base_response import EC2BaseResponse


class InstanceResponse(EC2BaseResponse):
    def describe_instances(self) -> ActionResult:
        self.error_on_dryrun()
        # https://boto3.amazonaws.com/v1/documentation/api/latest/reference/services/ec2/client/describe_instances.html
        # You cannot specify this(MaxResults) parameter and the instance IDs parameter in the same request.
        if "InstanceId.1" in self.data and "MaxResults" in self.data:
            raise InvalidParameterCombination(
                "The parameter instancesSet cannot be used with the parameter maxResults"
            )
        filter_dict = self._filters_from_querystring()
        instance_ids = self._get_param("InstanceIds", [])
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
        result = {"Reservations": reservations_resp, "NextToken": next_token}
        return ActionResult(result)

    def run_instances(self) -> ActionResult:
        min_count = int(self._get_param("MinCount", if_none="1"))
        image_id = self._get_param("ImageId")
        user_data = parse_user_data(self._get_param("UserData"))
        security_group_names = self._get_param("SecurityGroups", [])
        kwargs = {
            "instance_type": self._get_param("InstanceType", "m1.small"),
            "is_instance_type_default": not self._get_param("InstanceType"),
            "placement": self._get_param("Placement.AvailabilityZone"),
            "placement_hostid": self._get_param("Placement.HostId"),
            "region_name": self.region,
            "subnet_id": self._get_param("SubnetId"),
            "key_name": self._get_param("KeyName"),
            "security_group_ids": self._get_param("SecurityGroupIds", []),
            "nics": self._get_param("NetworkInterfaces", []),
            "private_ip": self._get_param("PrivateIpAddress"),
            "associate_public_ip": self._get_param("AssociatePublicIpAddress"),
            "tags": self._parse_tag_specification(),
            "ebs_optimized": self._get_param("EbsOptimized") or False,
            "disable_api_stop": self._get_param("DisableApiStop") or False,
            "instance_market_options": self._get_param(
                "InstanceMarketOptions.MarketType", {}
            ),
            "instance_initiated_shutdown_behavior": self._get_param(
                "InstanceInitiatedShutdownBehavior"
            ),
            "launch_template": self._get_param("LaunchTemplate", {}),
            "hibernation_options": self._get_param("HibernationOptions", {}),
            "iam_instance_profile_name": self._get_param("IamInstanceProfile.Name"),
            "iam_instance_profile_arn": self._get_param("IamInstanceProfile.Arn"),
            "monitoring_state": "enabled"
            if self._get_param("Monitoring.Enabled")
            else "disabled",
            "ipv6_address_count": self._get_int_param("Ipv6AddressCount"),
            "metadata_options": self._get_param("MetadataOptions", {}),
            "client_token": self._get_param("ClientToken", "ABCDE"),
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
        return ActionResult(new_reservation)

    def terminate_instances(self) -> ActionResult:
        instance_ids = self._get_param("InstanceIds", [])

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
        result = {
            "TerminatingInstances": [
                {
                    "InstanceId": instance.id,
                    "CurrentState": {"Code": 32, "Name": "shutting-down"},
                    "PreviousState": previous_state,
                }
                for instance, previous_state in instances
            ]
        }
        return ActionResult(result)

    def reboot_instances(self) -> ActionResult:
        instance_ids = self._get_param("InstanceIds", [])
        self.error_on_dryrun()
        self.ec2_backend.reboot_instances(instance_ids)
        return EmptyResult()

    def stop_instances(self) -> ActionResult:
        instance_ids = self._get_param("InstanceIds", [])
        self.error_on_dryrun()
        instances = self.ec2_backend.stop_instances(instance_ids)
        result = {
            "StoppingInstances": [
                {
                    "InstanceId": instance.id,
                    "CurrentState": {"Code": 64, "Name": "stopping"},
                    "PreviousState": previous_state,
                }
                for instance, previous_state in instances
            ]
        }
        return ActionResult(result)

    def start_instances(self) -> ActionResult:
        instance_ids = self._get_param("InstanceIds", [])
        self.error_on_dryrun()
        instances = self.ec2_backend.start_instances(instance_ids)
        result = {
            "StartingInstances": [
                {
                    "InstanceId": instance.id,
                    "CurrentState": {"Code": 0, "Name": "pending"},
                    "PreviousState": previous_state,
                }
                for instance, previous_state in instances
            ]
        }
        return ActionResult(result)

    def describe_instance_status(self) -> ActionResult:
        instance_ids = self._get_param("InstanceIds", [])
        include_all_instances = self._get_bool_param("IncludeAllInstances", False)
        filters = [
            {"name": k, "values": v}
            for k, v in self._filters_from_querystring().items()
        ]
        instances = self.ec2_backend.describe_instance_status(
            instance_ids, include_all_instances, filters
        )
        result = {"InstanceStatuses": instances}
        return ActionResult(result)

    def describe_instance_types(self) -> ActionResult:
        instance_type_filters = self._get_param("InstanceTypes", [])
        filter_dict = self._filters_from_querystring()
        instance_types = self.ec2_backend.describe_instance_types(
            instance_type_filters, filter_dict
        )
        result = {"InstanceTypes": instance_types}
        return ActionResult(result)

    def describe_instance_type_offerings(self) -> ActionResult:
        location_type_filters = self._get_param("LocationType")
        filter_dict = self._filters_from_querystring()
        offerings = self.ec2_backend.describe_instance_type_offerings(
            location_type_filters, filter_dict
        )
        result = {"InstanceTypeOfferings": offerings}
        return ActionResult(result)

    def describe_instance_attribute(self) -> ActionResult:
        # TODO this and modify below should raise IncorrectInstanceState if
        # instance not in stopped state
        attribute = self._get_param("Attribute")
        instance_id = self._get_param("InstanceId")
        instance, value = self.ec2_backend.describe_instance_attribute(
            instance_id, attribute
        )
        attribute_name = attribute[0].upper() + attribute[1:]
        attribute_value: Any = {"Value": value}
        if attribute_name == "GroupSet":
            attribute_name = "Groups"
            attribute_value = [{"GroupId": group.id} for group in value]
        elif attribute_name == "UserData" and value is None:
            attribute_value = ""
        result = {"InstanceId": instance.id, attribute_name: attribute_value}
        return ActionResult(result)

    def describe_instance_credit_specifications(self) -> ActionResult:
        instance_ids = self._get_param("InstanceIds", [])
        instances = self.ec2_backend.describe_instance_credit_specifications(
            instance_ids
        )
        result = {
            "InstanceCreditSpecifications": [
                {"InstanceId": instance.id, "CpuCredits": "standard"}
                for instance in instances
            ],
        }
        return ActionResult(result)

    def modify_instance_attribute(self) -> ActionResult:
        handlers = [
            self._attribute_value_handler,
            self._dot_value_instance_attribute_handler,
            self._block_device_mapping_handler,
            self._security_grp_instance_attribute_handler,
        ]

        for handler in handlers:
            success = handler()
            if success:
                return EmptyResult()

        msg = (
            "This specific call to ModifyInstanceAttribute has not been"
            " implemented in Moto yet. Feel free to open an issue at"
            " https://github.com/getmoto/moto/issues"
        )
        raise NotImplementedError(msg)

    def modify_instance_metadata_options(self) -> ActionResult:
        instance_id = self._get_param("InstanceId")
        tokens = self._get_param("HttpTokens")
        hop_limit = self._get_int_param("HttpPutResponseHopLimit")
        endpoint = self._get_param("HttpEndpoint")
        dry_run = False
        http_protocol = self._get_param("HttpProtocolIpv6")
        metadata_tags = self._get_param("InstanceMetadataTags")
        options = self.ec2_backend.modify_instance_metadata_options(
            instance_id=instance_id,
            http_tokens=tokens,
            hop_limit=hop_limit,
            http_endpoint=endpoint,
            dry_run=dry_run,
            http_protocol=http_protocol,
            metadata_tags=metadata_tags,
        )
        result = {"InstanceId": instance_id, "InstanceMetadataOptions": options}
        return ActionResult(result)

    def get_instance_uefi_data(self) -> ActionResult:
        instance_id = self._get_param("InstanceId")
        uefi_data = self.ec2_backend.get_instance_uefi_data(instance_id)
        result = {
            "InstanceId": instance_id,
            "UefiData": Base64EncodedString.from_encoded_bytes(uefi_data),
        }
        return ActionResult(result)

    def _block_device_mapping_handler(self) -> bool:
        """
        Handles requests which are generated by code similar to:

            instance.modify_attribute(
                BlockDeviceMappings=[{
                    'DeviceName': '/dev/sda1',
                    'Ebs': {'DeleteOnTermination': True}
                }]
            )

        For now we only support the "BlockDeviceMapping.1.Ebs.DeleteOnTermination"
        configuration, but it should be trivial to add anything else.
        """
        attribute_modified = False
        for mapping in self._get_param("BlockDeviceMappings", []):
            device_name = mapping.get("DeviceName")
            del_on_term_value = mapping.get("Ebs", {}).get("DeleteOnTermination")

            instance_id = self._get_param("InstanceId")
            instance = self.ec2_backend.get_instance(instance_id)

            self.error_on_dryrun()

            if del_on_term_value is None:
                continue
            block_device_type = instance.block_device_mapping[device_name]
            block_device_type.delete_on_termination = del_on_term_value

            attribute_modified = True

        if attribute_modified:
            return True
        return False

    def _dot_value_instance_attribute_handler(self) -> bool:
        attribute_modified = False
        for attribute in (
            "InstanceType",
            "Kernel",
            "UserData",
            "DisableApiTermination",
            "SourceDestCheck",
            "EbsOptimized",
            "DisableApiStop",
        ):
            if self._get_param(f"{attribute}.Value") is not None:
                self.error_on_dryrun()
                instance_id = self._get_param("InstanceId")
                attr_name = camelcase_to_underscores(attribute)
                attr_value = self._get_param(f"{attribute}.Value")
                if attribute == "UserData" and attr_value:
                    attr_value = parse_user_data(attr_value)
                self.ec2_backend.modify_instance_attribute(
                    instance_id, attr_name, attr_value
                )
                attribute_modified = True
                break
        if attribute_modified:
            return True
        return False

    def _attribute_value_handler(self) -> bool:
        attribute_key = self._get_param("Attribute")

        if attribute_key is None:
            return False

        self.error_on_dryrun()

        value = self._get_param("Value")
        normalized_attribute = camelcase_to_underscores(attribute_key)
        instance_id = self._get_param("InstanceId")
        self.ec2_backend.modify_instance_attribute(
            instance_id, normalized_attribute, value
        )
        return True

    def _security_grp_instance_attribute_handler(self) -> bool:
        new_security_grp_list = self._get_param("Groups", [])
        self.error_on_dryrun()
        instance_id = self._get_param("InstanceId")
        self.ec2_backend.modify_instance_security_groups(
            instance_id, new_security_grp_list
        )
        return True

    def _parse_block_device_mapping(self) -> list[dict[str, Any]]:
        device_mappings = self._get_param("BlockDeviceMappings", [])
        mappings = []
        for device_mapping in device_mappings:
            self._validate_block_device_mapping(device_mapping)
            device_template: dict[str, Any] = deepcopy(BLOCK_DEVICE_MAPPING_TEMPLATE)
            device_template["VirtualName"] = device_mapping.get("VirtualName")
            device_template["DeviceName"] = device_mapping.get("DeviceName")
            device_template["Ebs"]["SnapshotId"] = device_mapping.get("Ebs", {}).get(
                "SnapshotId"
            )
            device_template["Ebs"]["VolumeSize"] = device_mapping.get("Ebs", {}).get(
                "VolumeSize"
            )
            device_template["Ebs"]["DeleteOnTermination"] = device_mapping.get(
                "Ebs", {}
            ).get("DeleteOnTermination", False)
            device_template["Ebs"]["VolumeType"] = device_mapping.get("Ebs", {}).get(
                "VolumeType"
            )
            device_template["Ebs"]["Iops"] = device_mapping.get("Ebs", {}).get("Iops")
            device_template["Ebs"]["Encrypted"] = device_mapping.get("Ebs", {}).get(
                "VolumeSize", False
            )
            device_template["Ebs"]["KmsKeyId"] = device_mapping.get("Ebs", {}).get(
                "KmsKeyId"
            )
            device_template["NoDevice"] = device_mapping.get("NoDevice")
            mappings.append(device_template)

        return mappings

    @staticmethod
    def _validate_block_device_mapping(device_mapping: dict[str, Any]) -> None:  # type: ignore[misc]
        from botocore import __version__ as botocore_version

        if "NoDevice" in device_mapping:
            assert isinstance(device_mapping["NoDevice"], str), (
                f"botocore {botocore_version} isn't limiting NoDevice to str type anymore, it is type:{type(device_mapping['no_device'])}"
            )
            if device_mapping["NoDevice"] == "":
                # the only legit value it can have is empty string
                # and none of the other checks here matter if NoDevice
                # is being used
                return
            else:
                raise InvalidRequest()

        if "Ebs" not in device_mapping:
            raise MissingParameterError("ebs")
        if not device_mapping.get("Ebs", {}).get(
            "VolumeSize"
        ) and not device_mapping.get("Ebs", {}).get("SnapshotId"):
            raise MissingParameterError("size or snapshotId")


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
