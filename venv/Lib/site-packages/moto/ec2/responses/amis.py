from moto.core.responses import ActionResult, EmptyResult

from ..exceptions import AuthFailureRestricted, InvalidRequest
from ._base_response import EC2BaseResponse


class AmisResponse(EC2BaseResponse):
    def create_image(self) -> ActionResult:
        name = self._get_param("Name")
        description = self._get_param("Description", "")
        instance_id = self._get_param("InstanceId")
        tag_specifications = self._get_param("TagSpecifications", [])

        self.error_on_dryrun()

        image = self.ec2_backend.create_image(
            instance_id,
            name,
            description,
            tag_specifications=tag_specifications,
        )
        result = {"ImageId": image.id}
        return ActionResult(result)

    def copy_image(self) -> ActionResult:
        source_image_id = self._get_param("SourceImageId")
        source_region = self._get_param("SourceRegion")
        name = self._get_param("Name")
        description = self._get_param("Description")

        self.error_on_dryrun()

        image = self.ec2_backend.copy_image(
            source_image_id, source_region, name, description
        )
        result = {"ImageId": image.id}
        return ActionResult(result)

    def deregister_image(self) -> ActionResult:
        ami_id = self._get_param("ImageId")

        self.error_on_dryrun()

        self.ec2_backend.deregister_image(ami_id)
        result = {"Return": True}
        return ActionResult(result)

    def describe_images(self) -> ActionResult:
        self.error_on_dryrun()
        ami_ids = self._get_param("ImageIds", [])
        filters = self._filters_from_querystring()
        owners = self._get_param("Owners", [])
        exec_users = self._get_param("ExecutableUsers", [])
        images = self.ec2_backend.describe_images(
            ami_ids=ami_ids, filters=filters, exec_users=exec_users, owners=owners
        )
        result = {"Images": images}
        return ActionResult(result)

    def describe_image_attribute(self) -> ActionResult:
        ami_id = self._get_param("ImageId")
        attribute_name = self._get_param("Attribute")

        # only valid attributes as per
        # https://boto3.amazonaws.com/v1/documentation/api/latest/reference/services/ec2/client/describe_image_attribute.html
        valid_attributes_list = {
            "description": "description",
            "kernel": "kernel_id",
            "ramdisk": "ramdisk",
            "launchPermission": "launch_permissions",
            "productCodes": "product_codes",
            "blockDeviceMapping": "bdm",
            "sriovNetSupport": "sriov",
            "bootMode": "boot_mode",
            "tpmSupport": "tmp",
            "uefiData": "uefi",
            "lastLaunchedTime": "lld",
            "imdsSupport": "imds",
        }
        if attribute_name not in valid_attributes_list:
            raise InvalidRequest
        elif attribute_name == "blockDeviceMapping":
            # replicate real aws behaviour and throw and error
            # https://github.com/aws/aws-cli/issues/1083
            raise AuthFailureRestricted

        attribute_value = None
        launch_permissions = None
        if attribute_name == "launchPermission":
            launch_permissions = self.ec2_backend.describe_image_attribute(
                ami_id, valid_attributes_list[attribute_name]
            )
        else:
            attribute_value = self.ec2_backend.describe_image_attribute(
                ami_id, valid_attributes_list[attribute_name]
            )

        result = {"ImageId": ami_id}
        if attribute_name == "productCodes":
            result["ProductCodes"] = [
                {"ProductCodeId": code, "ProductCodeType": "marketplace"}
                for code in attribute_value  # type: ignore[union-attr]
            ]
        elif attribute_name == "launchPermission":
            result["LaunchPermissions"] = launch_permissions
        else:
            result[attribute_name] = {"Value": attribute_value}

        return ActionResult(result)

    def modify_image_attribute(self) -> ActionResult:
        ami_id = self._get_param("ImageId")
        launch_permissions_to_add = self._get_param("LaunchPermission", {}).get(
            "Add", []
        )
        launch_permissions_to_remove = self._get_param("LaunchPermission", {}).get(
            "Remove", []
        )
        # If only one OperationType is added, the other attributes are submitted as different variables
        operation_type = self._get_param("OperationType")
        if operation_type in ["add", "remove"]:
            group = self._get_param("UserGroups", [])
            lp = (
                launch_permissions_to_add
                if operation_type == "add"
                else launch_permissions_to_remove
            )
            if group:
                lp.append({"Group": group[0]})

            for user_id in self._get_param("UserIds", []):
                lp.append({"UserId": user_id})

            org_arn = self._get_param("OrganizationArns", [])
            if org_arn:
                lp.append({"OrganizationArn": org_arn[0]})

            ou_arn = self._get_param("OrganizationalUnitArns", [])
            if ou_arn:
                lp.append({"OrganizationalUnitArn": ou_arn[0]})

        self.error_on_dryrun()

        self.ec2_backend.modify_image_attribute(
            ami_id=ami_id,
            launch_permissions_to_add=launch_permissions_to_add,
            launch_permissions_to_remove=launch_permissions_to_remove,
        )
        return EmptyResult()

    def register_image(self) -> ActionResult:
        name = self._get_param("Name")
        description = self._get_param("Description", "")

        self.error_on_dryrun()

        image = self.ec2_backend.register_image(name, description)
        result = {"ImageId": image.id}
        return ActionResult(result)

    def reset_image_attribute(self) -> str:
        self.error_on_dryrun()

        raise NotImplementedError("AMIs.reset_image_attribute is not yet implemented")
