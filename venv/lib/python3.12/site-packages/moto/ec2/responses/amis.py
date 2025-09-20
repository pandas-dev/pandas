from moto.core.responses import ActionResult, EmptyResult

from ..exceptions import AuthFailureRestricted, InvalidRequest
from ._base_response import EC2BaseResponse


class AmisResponse(EC2BaseResponse):
    def create_image(self) -> ActionResult:
        name = self.querystring.get("Name")[0]  # type: ignore[index]
        description = self._get_param("Description", if_none="")
        instance_id = self._get_param("InstanceId")
        tag_specifications = self._get_multi_param("TagSpecification")

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
        ami_ids = self._get_multi_param("ImageId")
        filters = self._filters_from_querystring()
        owners = self._get_multi_param("Owner")
        exec_users = self._get_multi_param("ExecutableBy")
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
        launch_permissions_to_add = list(
            self._get_params().get("LaunchPermission", {}).get("Add", {}).values()
        )
        launch_permissions_to_remove = list(
            self._get_params().get("LaunchPermission", {}).get("Remove", {}).values()
        )
        # If only one OperationType is added, the other attributes are submitted as different variables
        operation_type = self._get_param("OperationType")
        if operation_type in ["add", "remove"]:
            group = self._get_param("UserGroup.1")
            lp = (
                launch_permissions_to_add
                if operation_type == "add"
                else launch_permissions_to_remove
            )
            if group:
                lp.append({"Group": group})

            for user_id in self._get_multi_param("UserId"):
                lp.append({"UserId": user_id})

            org_arn = self._get_param("OrganizationArn.1")
            if org_arn:
                lp.append({"OrganizationArn": org_arn})

            ou_arn = self._get_param("OrganizationalUnitArn.1")
            if ou_arn:
                lp.append({"OrganizationalUnitArn": ou_arn})

        self.error_on_dryrun()

        self.ec2_backend.modify_image_attribute(
            ami_id=ami_id,
            launch_permissions_to_add=launch_permissions_to_add,
            launch_permissions_to_remove=launch_permissions_to_remove,
        )
        return EmptyResult()

    def register_image(self) -> ActionResult:
        name = self.querystring.get("Name")[0]  # type: ignore[index]
        description = self._get_param("Description", if_none="")

        self.error_on_dryrun()

        image = self.ec2_backend.register_image(name, description)
        result = {"ImageId": image.id}
        return ActionResult(result)

    def reset_image_attribute(self) -> str:
        self.error_on_dryrun()

        raise NotImplementedError("AMIs.reset_image_attribute is not yet implemented")
