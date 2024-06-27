import json
import os
import re
from os import environ
from typing import Any, Dict, List, Optional, Set, cast

from moto import settings
from moto.utilities.utils import load_resource

from ..exceptions import (
    InvalidAMIAttributeItemValueError,
    InvalidAMIIdError,
    InvalidTaggableResourceType,
    MalformedAMIIdError,
    UnvailableAMIIdError,
)
from ..utils import (
    generic_filter,
    random_ami_id,
    utc_date_and_time,
)
from .core import TaggedEC2Resource
from .instances import Instance

if "MOTO_AMIS_PATH" in environ:
    with open(environ["MOTO_AMIS_PATH"], "r", encoding="utf-8") as f:
        AMIS: List[Dict[str, Any]] = json.load(f)
else:
    AMIS = load_resource(__name__, "../resources/amis.json")


class Ami(TaggedEC2Resource):
    def __init__(  # pylint: disable=dangerous-default-value
        self,
        ec2_backend: Any,
        ami_id: str,
        instance: Optional[Instance] = None,
        source_ami: Optional["Ami"] = None,
        name: Optional[str] = None,
        description: Optional[str] = None,
        owner_id: Optional[str] = None,
        owner_alias: Optional[str] = None,
        public: bool = False,
        virtualization_type: Optional[str] = None,
        architecture: Optional[str] = None,
        state: str = "available",
        creation_date: Optional[str] = None,
        platform: Optional[str] = None,
        image_type: str = "machine",
        image_location: Optional[str] = None,
        hypervisor: Optional[str] = None,
        root_device_type: str = "standard",
        root_device_name: str = "/dev/sda1",
        sriov: str = "simple",
        region_name: str = "us-east-1a",
        snapshot_description: Optional[str] = None,
        product_codes: Set[str] = set(),
        boot_mode: str = "uefi",
    ):
        self.ec2_backend = ec2_backend
        self.id = ami_id
        self.state = state
        self.name = name
        self.image_type = image_type
        self.image_location = image_location
        self.owner_id = owner_id or ec2_backend.account_id
        self.owner_alias = owner_alias
        self.description = description
        self.virtualization_type = virtualization_type
        self.architecture = architecture
        self.kernel_id = None
        self.platform = platform
        self.hypervisor = hypervisor
        self.root_device_name = root_device_name
        self.root_device_type = root_device_type
        self.sriov = sriov
        self.creation_date = creation_date or utc_date_and_time()
        self.product_codes = product_codes
        self.boot_mode = boot_mode

        if instance:
            self.instance = instance
            self.instance_id = instance.id
            self.virtualization_type = instance.virtualization_type
            self.architecture = instance.architecture
            self.kernel_id = instance.kernel
            self.platform = instance.platform

        elif source_ami:
            """
            http://docs.aws.amazon.com/AWSEC2/latest/UserGuide/CopyingAMIs.html
            "We don't copy launch permissions, user-defined tags, or Amazon S3 bucket permissions from the source AMI to the new AMI."
            ~ 2014.09.29
            """
            self.virtualization_type = source_ami.virtualization_type
            self.architecture = source_ami.architecture
            self.kernel_id = source_ami.kernel_id
            self.platform = source_ami.platform
            if not name:
                self.name = source_ami.name
            if not description:
                self.description = source_ami.description

        self.launch_permissions: List[Dict[str, str]] = []

        if public:
            self.launch_permissions.append({"Group": "all"})

        # AWS auto-creates these, we should reflect the same.
        volume = self.ec2_backend.create_volume(size=15, zone_name=region_name)
        snapshot_description = (
            snapshot_description or f"Auto-created snapshot for AMI {self.id}"
        )
        self.ebs_snapshot = self.ec2_backend.create_snapshot(
            volume.id, snapshot_description, self.owner_id, from_ami=ami_id
        )
        self.ec2_backend.delete_volume(volume.id)

    @property
    def is_public(self) -> bool:
        return {"Group": "all"} in self.launch_permissions

    @property
    def is_public_string(self) -> str:
        return str(self.is_public).lower()

    def get_filter_value(
        self, filter_name: str, method_name: Optional[str] = None
    ) -> Any:
        if filter_name == "virtualization-type":
            return self.virtualization_type
        elif filter_name == "kernel-id":
            return self.kernel_id
        elif filter_name in ["architecture", "platform"]:
            return getattr(self, filter_name)
        elif filter_name == "image-id":
            return self.id
        elif filter_name == "is-public":
            return self.is_public_string
        elif filter_name == "state":
            return self.state
        elif filter_name == "name":
            return self.name
        elif filter_name == "owner-id":
            return self.owner_id
        elif filter_name == "owner-alias":
            return self.owner_alias
        elif filter_name == "product-code":
            return self.product_codes
        elif filter_name == "product-code.type":
            return "marketplace"  # devpay is not (yet?) supported
        else:
            return super().get_filter_value(filter_name, "DescribeImages")


class AmiBackend:
    AMI_REGEX = re.compile("ami-[a-z0-9]+")

    def __init__(self) -> None:
        self.amis: Dict[str, Ami] = {}
        self.deleted_amis: List[str] = list()
        self._load_amis()

    def _load_amis(self) -> None:
        if "MOTO_AMIS_PATH" not in os.environ and not settings.ec2_load_default_amis():
            return
        for ami in AMIS:
            ami_id = ami["ami_id"]
            # we are assuming the default loaded amis are owned by amazon
            # owner_alias is required for terraform owner filters
            ami["owner_alias"] = "amazon"
            self.amis[ami_id] = Ami(self, **ami)
        if "MOTO_AMIS_PATH" not in environ:
            for path in ["latest_amis", "ecs/optimized_amis"]:
                try:
                    latest_amis = cast(
                        List[Dict[str, Any]],
                        load_resource(
                            __name__,
                            f"../resources/{path}/{self.region_name}.json",  # type: ignore[attr-defined]
                        ),
                    )
                    for ami in latest_amis:
                        ami_id = ami["ami_id"]
                        ami["owner_alias"] = "amazon"
                        self.amis[ami_id] = Ami(self, **ami)
                except FileNotFoundError:
                    # Will error on unknown (new) regions - just return an empty list here
                    pass

    def create_image(
        self,
        instance_id: str,
        name: str,
        description: str,
        tag_specifications: List[Dict[str, Any]],
    ) -> Ami:
        # TODO: check that instance exists and pull info from it.
        ami_id = random_ami_id()
        instance = self.get_instance(instance_id)  # type: ignore[attr-defined]
        tags = []
        for tag_specification in tag_specifications:
            resource_type = tag_specification["ResourceType"]
            if resource_type == "image":
                tags += tag_specification["Tag"]
            elif resource_type == "snapshot":
                raise NotImplementedError()
            else:
                raise InvalidTaggableResourceType(resource_type)

        ami = Ami(
            self,
            ami_id,
            instance=instance,
            source_ami=None,
            name=name,
            description=description,
            owner_id=None,
            snapshot_description=f"Created by CreateImage({instance_id}) for {ami_id}",
        )
        for tag in tags:
            ami.add_tag(tag["Key"], tag["Value"])
        self.amis[ami_id] = ami
        return ami

    def copy_image(
        self,
        source_image_id: str,
        source_region: str,
        name: Optional[str] = None,
        description: Optional[str] = None,
    ) -> Ami:
        from ..models import ec2_backends

        source_backend = ec2_backends[self.account_id][source_region]  # type: ignore[attr-defined]
        source_ami = source_backend.describe_images(ami_ids=[source_image_id])[0]
        ami_id = random_ami_id()
        ami = Ami(
            self,
            ami_id,
            instance=None,
            source_ami=source_ami,
            name=name,
            description=description,
        )
        self.amis[ami_id] = ami
        return ami

    def describe_images(
        self,
        ami_ids: Optional[List[str]] = None,
        filters: Optional[Dict[str, Any]] = None,
        exec_users: Optional[List[str]] = None,
        owners: Optional[List[str]] = None,
    ) -> List[Ami]:
        images = list(self.amis.copy().values())

        if ami_ids and len(ami_ids):
            # boto3 seems to default to just searching based on ami ids if that parameter is passed
            # and if no images are found, it raises an errors
            # Note that we can search for images that have been previously deleted, without raising any errors
            malformed_ami_ids = [
                ami_id for ami_id in ami_ids if not ami_id.startswith("ami-")
            ]
            if malformed_ami_ids:
                raise MalformedAMIIdError(malformed_ami_ids)

            images = [ami for ami in images if ami.id in ami_ids]
            deleted_images = [
                ami_id for ami_id in ami_ids if ami_id in self.deleted_amis
            ]
            if len(images) + len(deleted_images) == 0:
                raise InvalidAMIIdError(ami_ids)
        else:
            # Limit images by launch permissions
            if exec_users:
                tmp_images = []
                for ami in images:
                    for user_id in exec_users:
                        for lp in ami.launch_permissions:
                            if lp.get("UserId") == user_id:
                                tmp_images.append(ami)
                images = tmp_images

            # Limit by owner ids
            if owners:
                # support filtering by Owners=['self']
                if "self" in owners:
                    owners = list(
                        map(lambda o: self.account_id if o == "self" else o, owners)  # type: ignore[attr-defined]
                    )
                images = [
                    ami
                    for ami in images
                    if ami.owner_id in owners or ami.owner_alias in owners
                ]

            # Generic filters
            if filters:
                return generic_filter(filters, images)

        return images

    def deregister_image(self, ami_id: str) -> None:
        if ami_id in self.amis:
            self.amis.pop(ami_id)
            self.deleted_amis.append(ami_id)
        elif ami_id in self.deleted_amis:
            raise UnvailableAMIIdError(ami_id)
        else:
            raise InvalidAMIIdError(ami_id)

    def validate_permission_targets(self, permissions: List[Dict[str, str]]) -> None:
        for perm in permissions:
            # If anything is invalid, nothing is added. (No partial success.)
            # AWS docs:
            # The AWS account ID is a 12-digit number, such as 123456789012, that you use to construct Amazon Resource Names (ARNs)."
            # http://docs.aws.amazon.com/general/latest/gr/acct-identifiers.html

            if "UserId" in perm and (
                len(perm["UserId"]) != 12 or not perm["UserId"].isdigit()
            ):
                raise InvalidAMIAttributeItemValueError("userId", perm["UserId"])

            if "Group" in perm and perm["Group"] != "all":
                raise InvalidAMIAttributeItemValueError("UserGroup", perm["Group"])

    def modify_image_attribute(
        self,
        ami_id: str,
        launch_permissions_to_add: List[Dict[str, str]],
        launch_permissions_to_remove: List[Dict[str, str]],
    ) -> None:
        ami = self.describe_images(ami_ids=[ami_id])[0]
        self.validate_permission_targets(launch_permissions_to_add)
        self.validate_permission_targets(launch_permissions_to_remove)
        for lp in launch_permissions_to_add:
            if lp not in ami.launch_permissions:
                ami.launch_permissions.append(lp)
        for lp in launch_permissions_to_remove:
            try:
                ami.launch_permissions.remove(lp)
            except ValueError:
                # The LaunchPermission may not exist
                pass

    def register_image(
        self, name: Optional[str] = None, description: Optional[str] = None
    ) -> Ami:
        ami_id = random_ami_id()
        ami = Ami(
            self,
            ami_id,
            instance=None,
            source_ami=None,
            name=name,
            description=description,
        )
        self.amis[ami_id] = ami
        return ami

    def describe_image_attribute(self, ami_id: str, attribute_name: str) -> Any:
        return self.amis[ami_id].__getattribute__(attribute_name)
