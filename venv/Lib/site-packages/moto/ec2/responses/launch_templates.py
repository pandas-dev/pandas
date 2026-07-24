from moto.core.responses import ActionResult
from moto.ec2.exceptions import InvalidParameterCombination
from moto.ec2.utils import parse_user_data

from ._base_response import EC2BaseResponse


class LaunchTemplates(EC2BaseResponse):
    def create_launch_template(self) -> ActionResult:
        name = self._get_param("LaunchTemplateName")
        version_description = self._get_param("VersionDescription")
        tag_spec = self._parse_tag_specification()

        parsed_template_data = self._get_param("LaunchTemplateData", {})
        parsed_template_data["UserData"] = parse_user_data(
            self._get_param("LaunchTemplateData.UserData")
        )
        self.error_on_dryrun()

        if tag_spec:
            if "TagSpecifications" not in parsed_template_data:
                parsed_template_data["TagSpecifications"] = []
            converted_tag_spec = []
            for resource_type, tags in tag_spec.items():
                converted_tag_spec.append(
                    {
                        "ResourceType": resource_type,
                        "Tags": [
                            {"Key": key, "Value": value} for key, value in tags.items()
                        ],
                    }
                )

            parsed_template_data["TagSpecifications"].extend(converted_tag_spec)

        template = self.ec2_backend.create_launch_template(
            name, version_description, parsed_template_data, tag_spec
        )
        version = template.default_version()

        result = {
            "LaunchTemplate": {
                "CreateTime": version.create_time,
                "CreatedBy": f"arn:{self.partition}:iam::{self.current_account}:root",
                "DefaultVersionNumber": template.default_version_number,
                "LatestVersionNumber": version.number,
                "LaunchTemplateId": template.id,
                "LaunchTemplateName": template.name,
                "Tags": template.tags,
            },
        }

        return ActionResult(result)

    def create_launch_template_version(self) -> ActionResult:
        name = self._get_param("LaunchTemplateName")
        tmpl_id = self._get_param("LaunchTemplateId")
        if name:
            template = self.ec2_backend.get_launch_template_by_name(name)
        if tmpl_id:
            template = self.ec2_backend.get_launch_template(tmpl_id)

        version_description = self._get_param("VersionDescription")

        template_data = self._get_param("LaunchTemplateData", {})

        self.error_on_dryrun()

        version = template.create_version(template_data, version_description)

        result = {
            "LaunchTemplateVersion": {
                "CreateTime": version.create_time,
                "CreatedBy": f"arn:{self.partition}:iam::{self.current_account}:root",
                "DefaultVersion": template.is_default(version),
                "LaunchTemplateData": version.data,
                "LaunchTemplateId": template.id,
                "LaunchTemplateName": template.name,
                "VersionDescription": version.description,
                "VersionNumber": version.number,
            },
        }
        return ActionResult(result)

    def delete_launch_template(self) -> ActionResult:
        name = self._get_param("LaunchTemplateName")
        tid = self._get_param("LaunchTemplateId")

        self.error_on_dryrun()

        template = self.ec2_backend.delete_launch_template(name, tid)

        result = {
            "LaunchTemplate": {
                "DefaultVersionNumber": template.default_version_number,
                "LaunchTemplateId": template.id,
                "LaunchTemplateName": template.name,
            },
        }

        return ActionResult(result)

    def describe_launch_template_versions(self) -> ActionResult:
        name = self._get_param("LaunchTemplateName")
        template_id = self._get_param("LaunchTemplateId")

        max_results = self._get_int_param("MaxResults", 15)
        versions = self._get_param("Versions", [])
        min_version = self._get_int_param("MinVersion")
        max_version = self._get_int_param("MaxVersion")

        self.error_on_dryrun()

        ret_versions = self.ec2_backend.describe_launch_template_versions(
            template_name=name,
            template_id=template_id,
            versions=versions,
            min_version=min_version,
            max_version=max_version,
            max_results=max_results,
        )

        result = {
            "LaunchTemplateVersions": [
                {
                    "CreateTime": version.create_time,
                    "CreatedBy": f"arn:{self.partition}:iam::{self.current_account}:root",
                    "DefaultVersion": template.is_default(version),
                    "LaunchTemplateData": version.data,
                    "LaunchTemplateId": template.id,
                    "LaunchTemplateName": template.name,
                    "VersionDescription": version.description,
                    "VersionNumber": version.number,
                }
                for template, version in ret_versions
            ]
        }

        return ActionResult(result)

    def describe_launch_templates(self) -> ActionResult:
        max_results = self._get_int_param("MaxResults", 15)
        template_names = self._get_param("LaunchTemplateNames", [])
        template_ids = self._get_param("LaunchTemplateIds", [])
        filters = self._filters_from_querystring()

        self.error_on_dryrun()

        templates = self.ec2_backend.describe_launch_templates(
            template_names=template_names,
            template_ids=template_ids,
            filters=filters,
        )

        templates = templates[:max_results]

        result = {
            "LaunchTemplates": [
                {
                    "CreateTime": template.create_time,
                    "CreatedBy": f"arn:{self.partition}:iam::{self.current_account}:root",
                    "DefaultVersionNumber": template.default_version_number,
                    "LatestVersionNumber": template.latest_version_number,
                    "LaunchTemplateId": template.id,
                    "LaunchTemplateName": template.name,
                    "Tags": template.tags,
                }
                for template in templates
            ]
        }

        return ActionResult(result)

    def get_launch_template_data(self) -> ActionResult:
        instance_id = self._get_param("InstanceId")
        instance = self.ec2_backend.get_launch_template_data(instance_id)
        # Result is based on Moto's original XML jinja template, including hardcoded values.
        result = {
            "LaunchTemplateData": {
                "BlockDeviceMappings": [
                    {
                        "DeviceName": "string",
                        "Ebs": {
                            "Encrypted": device.encrypted,
                            "DeleteOnTermination": device.delete_on_termination,
                            "SnapshotId": device.snapshot_id,
                            "VolumeSize": device.size,
                            "VolumeType": device.volume_type,
                        },
                    }
                    for device_name, device in instance.block_device_mapping.items()
                ],
                "CapacityReservationSpecification": {
                    "CapacityReservationPreference": "open",
                },
                "CreditSpecification": {"CpuCredits": "standard"},
                "DisableApiStop": instance.disable_api_stop,
                "DisableApiTermination": instance.disable_api_termination,
                "EbsOptimized": instance.ebs_optimized,
                "EnclaveOptions": {"Enabled": False},
                "HibernationOptions": {"Configured": False},
                "ImageId": instance.image_id,
                "InstanceInitiatedShutdownBehavior": instance.instance_initiated_shutdown_behavior,
                "InstanceType": instance.instance_type,
                "KeyName": instance.key_name,
                "MaintenanceOptions": {"AutoRecovery": "default"},
                "MetadataOptions": {
                    "HttpTokens": "optional",
                    "HttpPutResponseHopLimit": 1,
                    "HttpEndpoint": "enabled",
                    "HttpProtocolIpv6": "disabled",
                    "InstanceMetadataTags": "disabled",
                },
                "Monitoring": {"Enabled": instance.monitored},
                "NetworkInterfaces": [
                    {
                        "AssociatePublicIpAddress": True,
                        "DeleteOnTermination": nic.delete_on_termination,
                        "Description": "",
                        "DeviceIndex": nic.device_index,
                        "Groups": [
                            {"GroupId": group.group_id}
                            for group in nic.group_set
                            if nic.group_set
                        ],
                        "InterfaceType": nic.interface_type,
                        "PrivateIpAddresses": [
                            {
                                "Primary": addr["Primary"],
                                "PrivateIpAddress": addr["PrivateIpAddress"],
                            }
                            for addr in nic.private_ip_addresses
                        ],
                        "SubnetId": nic.subnet.id,
                        "NetworkCardIndex": nic_index,
                    }
                    for nic_index, nic in instance.nics.items()
                ],
                "Placement": {
                    "AvailabilityZone": instance.availability_zone,
                    "GroupName": "",
                    "Tenancy": "default",
                },
                "PrivateDnsNameOptions": {
                    "HostnameType": "ip-name",
                    "EnableResourceNameDnsARecord": True,
                    "EnableResourceNameDnsAAAARecord": False,
                },
                "TagSpecifications": [
                    {
                        "ResourceType": "instance",
                        "Tags": [
                            {
                                "Key": tag["key"],
                                "Value": tag["value"],
                            }
                            for tag in instance.get_tags()
                        ],
                    }
                ]
                if instance.get_tags()
                else [],
            }
        }
        return ActionResult(result)

    def modify_launch_template(self) -> ActionResult:
        template_name = self._get_param("LaunchTemplateName")
        template_id = self._get_param("LaunchTemplateId")
        default_version = self._get_param("DefaultVersion")

        if template_name and template_id:
            raise InvalidParameterCombination(
                "Either provide launch template ID or launch template name to modify the template."
            )

        self.error_on_dryrun()

        template = self.ec2_backend.modify_launch_template(
            template_name=template_name,
            template_id=template_id,
            default_version=default_version,
        )

        result = {
            "LaunchTemplate": {
                "CreateTime": template.create_time,
                "CreatedBy": f"arn:{self.partition}:iam::{self.current_account}:root",
                "DefaultVersionNumber": template.default_version_number,
                "LatestVersionNumber": template.latest_version_number,
                "LaunchTemplateId": template.id,
                "LaunchTemplateName": template.name,
                "Tags": template.tags,
            },
        }
        return ActionResult(result)
