"""DirectoryServiceBackend class with methods for supported APIs."""

from datetime import datetime, timezone
from typing import Any, Dict, List, Optional, Tuple

from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.common_models import BaseModel
from moto.ds.exceptions import (
    ClientException,
    DirectoryLimitExceededException,
    EntityAlreadyExistsException,
    EntityDoesNotExistException,
    InvalidParameterException,
    TagLimitExceededException,
    ValidationException,
)
from moto.ds.utils import PAGINATION_MODEL
from moto.ds.validations import validate_args
from moto.ec2 import ec2_backends
from moto.ec2.exceptions import InvalidSubnetIdError
from moto.moto_api._internal import mock_random
from moto.utilities.paginator import paginate
from moto.utilities.tagging_service import TaggingService


class Directory(BaseModel):  # pylint: disable=too-many-instance-attributes
    """Representation of a Simple AD Directory.

    When the "create" API for a Simple AD or a Microsoft AD directory is
    invoked, two domain controllers and a DNS server are supposed to be
    created.  That is NOT done for the fake directories.

    However, the DnsIpAddrs attribute is supposed to contain the IP addresses
    of the DNS servers.  For a AD Connecter, the DnsIpAddrs are provided when
    the directory is created, but the ConnectSettings.ConnectIps values should
    contain the IP addresses of the DNS servers or domain controllers in the
    directory to which the AD connector is connected.

    Instead, the dns_ip_addrs attribute or ConnectIPs attribute for the fake
    directories will contain IPs picked from the subnets' CIDR blocks.
    """

    # The assumption here is that the limits are the same for all regions.
    CLOUDONLY_DIRECTORIES_LIMIT = 10
    CLOUDONLY_MICROSOFT_AD_LIMIT = 20
    CONNECTED_DIRECTORIES_LIMIT = 10

    MAX_TAGS_PER_DIRECTORY = 50

    def __init__(
        self,
        account_id: str,
        region: str,
        name: str,
        password: str,
        directory_type: str,
        size: Optional[str] = None,
        vpc_settings: Optional[Dict[str, Any]] = None,
        connect_settings: Optional[Dict[str, Any]] = None,
        short_name: Optional[str] = None,
        description: Optional[str] = None,
        edition: Optional[str] = None,
    ):  # pylint: disable=too-many-arguments
        self.account_id = account_id
        self.region = region
        self.name = name
        self.password = password
        self.directory_type = directory_type
        self.size = size
        self.vpc_settings = vpc_settings
        self.connect_settings = connect_settings
        self.short_name = short_name
        self.description = description
        self.edition = edition

        # Calculated or default values for the directory attributes.
        self.directory_id = f"d-{mock_random.get_random_hex(10)}"
        self.access_url = f"{self.directory_id}.awsapps.com"
        self.alias = self.directory_id
        self.desired_number_of_domain_controllers = 0
        self.sso_enabled = False
        self.stage = "Active"
        self.launch_time = datetime.now(timezone.utc).isoformat()
        self.stage_last_updated_date_time = datetime.now(timezone.utc).isoformat()

        if self.directory_type == "ADConnector":
            self.security_group_id = self.create_security_group(
                self.connect_settings["VpcId"]  # type: ignore[index]
            )
            self.eni_ids, self.subnet_ips = self.create_eni(
                self.security_group_id,
                self.connect_settings["SubnetIds"],  # type: ignore[index]
            )
            self.connect_settings["SecurityGroupId"] = self.security_group_id  # type: ignore[index]
            self.connect_settings["ConnectIps"] = self.subnet_ips  # type: ignore[index]
            self.dns_ip_addrs = self.connect_settings["CustomerDnsIps"]  # type: ignore[index]

        else:
            self.security_group_id = self.create_security_group(
                self.vpc_settings["VpcId"]  # type: ignore[index]
            )
            self.eni_ids, self.subnet_ips = self.create_eni(
                self.security_group_id,
                self.vpc_settings["SubnetIds"],  # type: ignore[index]
            )
            self.vpc_settings["SecurityGroupId"] = self.security_group_id  # type: ignore[index]
            self.dns_ip_addrs = self.subnet_ips

    def create_security_group(self, vpc_id: str) -> str:
        """Create security group for the network interface."""
        security_group_info = ec2_backends[self.account_id][
            self.region
        ].create_security_group(
            name=f"{self.directory_id}_controllers",
            description=(
                f"AWS created security group for {self.directory_id} "
                f"directory controllers"
            ),
            vpc_id=vpc_id,
        )
        return security_group_info.id

    def delete_security_group(self) -> None:
        """Delete the given security group."""
        ec2_backends[self.account_id][self.region].delete_security_group(
            group_id=self.security_group_id
        )

    def create_eni(
        self, security_group_id: str, subnet_ids: List[str]
    ) -> Tuple[List[str], List[str]]:
        """Return ENI ids and primary addresses created for each subnet."""
        eni_ids = []
        subnet_ips = []
        for subnet_id in subnet_ids:
            eni_info = ec2_backends[self.account_id][
                self.region
            ].create_network_interface(
                subnet=subnet_id,
                private_ip_address=None,  # type: ignore[arg-type]
                group_ids=[security_group_id],
                description=f"AWS created network interface for {self.directory_id}",
            )
            eni_ids.append(eni_info.id)
            subnet_ips.append(eni_info.private_ip_address)
        return eni_ids, subnet_ips  # type: ignore[return-value]

    def delete_eni(self) -> None:
        """Delete ENI for each subnet and the security group."""
        for eni_id in self.eni_ids:
            ec2_backends[self.account_id][self.region].delete_network_interface(eni_id)

    def update_alias(self, alias: str) -> None:
        """Change default alias to given alias."""
        self.alias = alias
        self.access_url = f"{alias}.awsapps.com"

    def enable_sso(self, new_state: bool) -> None:
        """Enable/disable sso based on whether new_state is True or False."""
        self.sso_enabled = new_state

    def to_dict(self) -> Dict[str, Any]:
        """Create a dictionary of attributes for Directory."""
        attributes = {
            "AccessUrl": self.access_url,
            "Alias": self.alias,
            "DirectoryId": self.directory_id,
            "DesiredNumberOfDomainControllers": self.desired_number_of_domain_controllers,
            "DnsIpAddrs": self.dns_ip_addrs,
            "LaunchTime": self.launch_time,
            "Name": self.name,
            "SsoEnabled": self.sso_enabled,
            "Stage": self.stage,
            "StageLastUpdatedDateTime": self.stage_last_updated_date_time,
            "Type": self.directory_type,
        }

        if self.edition:
            attributes["Edition"] = self.edition
        if self.size:
            attributes["Size"] = self.size
        if self.short_name:
            attributes["ShortName"] = self.short_name
        if self.description:
            attributes["Description"] = self.description

        if self.vpc_settings:
            attributes["VpcSettings"] = self.vpc_settings
        else:
            attributes["ConnectSettings"] = self.connect_settings
            attributes["ConnectSettings"]["CustomerDnsIps"] = None
        return attributes


class DirectoryServiceBackend(BaseBackend):
    """Implementation of DirectoryService APIs."""

    def __init__(self, region_name: str, account_id: str):
        super().__init__(region_name, account_id)
        self.directories: Dict[str, Directory] = {}
        self.tagger = TaggingService()

    def _verify_subnets(self, region: str, vpc_settings: Dict[str, Any]) -> None:
        """Verify subnets are valid, else raise an exception.

        If settings are valid, add AvailabilityZones to vpc_settings.
        """
        if len(vpc_settings["SubnetIds"]) != 2:
            raise InvalidParameterException(
                "Invalid subnet ID(s). They must correspond to two subnets "
                "in different Availability Zones."
            )

        # Subnet IDs are checked before the VPC ID.  The Subnet IDs must
        # be valid and in different availability zones.
        try:
            subnets = ec2_backends[self.account_id][region].describe_subnets(
                subnet_ids=vpc_settings["SubnetIds"]
            )
        except InvalidSubnetIdError as exc:
            raise InvalidParameterException(
                "Invalid subnet ID(s). They must correspond to two subnets "
                "in different Availability Zones."
            ) from exc

        regions = [subnet.availability_zone for subnet in subnets]
        if regions[0] == regions[1]:
            raise ClientException(
                "Invalid subnet ID(s). The two subnets must be in "
                "different Availability Zones."
            )

        vpcs = ec2_backends[self.account_id][region].describe_vpcs()
        if vpc_settings["VpcId"] not in [x.id for x in vpcs]:
            raise ClientException("Invalid VPC ID.")
        vpc_settings["AvailabilityZones"] = regions

    def connect_directory(
        self,
        region: str,
        name: str,
        short_name: str,
        password: str,
        description: str,
        size: str,
        connect_settings: Dict[str, Any],
        tags: List[Dict[str, str]],
    ) -> str:  # pylint: disable=too-many-arguments
        """Create a fake AD Connector."""
        if len(self.directories) > Directory.CONNECTED_DIRECTORIES_LIMIT:
            raise DirectoryLimitExceededException(
                f"Directory limit exceeded. A maximum of "
                f"{Directory.CONNECTED_DIRECTORIES_LIMIT} directories may be created"
            )

        validate_args(
            [
                ("password", password),
                ("size", size),
                ("name", name),
                ("description", description),
                ("shortName", short_name),
                (
                    "connectSettings.vpcSettings.subnetIds",
                    connect_settings["SubnetIds"],
                ),
                (
                    "connectSettings.customerUserName",
                    connect_settings["CustomerUserName"],
                ),
                ("connectSettings.customerDnsIps", connect_settings["CustomerDnsIps"]),
            ]
        )
        # ConnectSettings and VpcSettings both have a VpcId and Subnets.
        self._verify_subnets(region, connect_settings)

        errmsg = self.tagger.validate_tags(tags or [])
        if errmsg:
            raise ValidationException(errmsg)
        if len(tags) > Directory.MAX_TAGS_PER_DIRECTORY:
            raise DirectoryLimitExceededException("Tag Limit is exceeding")

        directory = Directory(
            self.account_id,
            region,
            name,
            password,
            "ADConnector",
            size=size,
            connect_settings=connect_settings,
            short_name=short_name,
            description=description,
        )
        self.directories[directory.directory_id] = directory
        self.tagger.tag_resource(directory.directory_id, tags or [])
        return directory.directory_id

    def create_directory(
        self,
        region: str,
        name: str,
        short_name: str,
        password: str,
        description: str,
        size: str,
        vpc_settings: Dict[str, Any],
        tags: List[Dict[str, str]],
    ) -> str:  # pylint: disable=too-many-arguments
        """Create a fake Simple Ad Directory."""
        if len(self.directories) > Directory.CLOUDONLY_DIRECTORIES_LIMIT:
            raise DirectoryLimitExceededException(
                f"Directory limit exceeded. A maximum of "
                f"{Directory.CLOUDONLY_DIRECTORIES_LIMIT} directories may be created"
            )

        # botocore doesn't look for missing vpc_settings, but boto3 does.
        if not vpc_settings:
            raise InvalidParameterException("VpcSettings must be specified.")
        validate_args(
            [
                ("password", password),
                ("size", size),
                ("name", name),
                ("description", description),
                ("shortName", short_name),
                ("vpcSettings.subnetIds", vpc_settings["SubnetIds"]),
            ]
        )
        self._verify_subnets(region, vpc_settings)

        errmsg = self.tagger.validate_tags(tags or [])
        if errmsg:
            raise ValidationException(errmsg)
        if len(tags) > Directory.MAX_TAGS_PER_DIRECTORY:
            raise DirectoryLimitExceededException("Tag Limit is exceeding")

        directory = Directory(
            self.account_id,
            region,
            name,
            password,
            "SimpleAD",
            size=size,
            vpc_settings=vpc_settings,
            short_name=short_name,
            description=description,
        )
        self.directories[directory.directory_id] = directory
        self.tagger.tag_resource(directory.directory_id, tags or [])
        return directory.directory_id

    def _validate_directory_id(self, directory_id: str) -> None:
        """Raise an exception if the directory id is invalid or unknown."""
        # Validation of ID takes precedence over a check for its existence.
        validate_args([("directoryId", directory_id)])
        if directory_id not in self.directories:
            raise EntityDoesNotExistException(
                f"Directory {directory_id} does not exist"
            )

    def create_alias(self, directory_id: str, alias: str) -> Dict[str, str]:
        """Create and assign an alias to a directory."""
        self._validate_directory_id(directory_id)

        # The default alias name is the same as the directory name.  Check
        # whether this directory was already given an alias.
        directory = self.directories[directory_id]
        if directory.alias != directory_id:
            raise InvalidParameterException(
                "The directory in the request already has an alias. That "
                "alias must be deleted before a new alias can be created."
            )

        # Is the alias already in use?
        if alias in [x.alias for x in self.directories.values()]:
            raise EntityAlreadyExistsException(f"Alias '{alias}' already exists.")
        validate_args([("alias", alias)])

        directory.update_alias(alias)
        return {"DirectoryId": directory_id, "Alias": alias}

    def create_microsoft_ad(
        self,
        region: str,
        name: str,
        short_name: str,
        password: str,
        description: str,
        vpc_settings: Dict[str, Any],
        edition: str,
        tags: List[Dict[str, str]],
    ) -> str:  # pylint: disable=too-many-arguments
        """Create a fake Microsoft Ad Directory."""
        if len(self.directories) > Directory.CLOUDONLY_MICROSOFT_AD_LIMIT:
            raise DirectoryLimitExceededException(
                f"Directory limit exceeded. A maximum of "
                f"{Directory.CLOUDONLY_MICROSOFT_AD_LIMIT} directories may be created"
            )

        # boto3 looks for missing vpc_settings for create_microsoft_ad().
        validate_args(
            [
                ("password", password),
                ("edition", edition),
                ("name", name),
                ("description", description),
                ("shortName", short_name),
                ("vpcSettings.subnetIds", vpc_settings["SubnetIds"]),
            ]
        )
        self._verify_subnets(region, vpc_settings)

        errmsg = self.tagger.validate_tags(tags or [])
        if errmsg:
            raise ValidationException(errmsg)
        if len(tags) > Directory.MAX_TAGS_PER_DIRECTORY:
            raise DirectoryLimitExceededException("Tag Limit is exceeding")

        directory = Directory(
            self.account_id,
            region,
            name,
            password,
            "MicrosoftAD",
            vpc_settings=vpc_settings,
            short_name=short_name,
            description=description,
            edition=edition,
        )
        self.directories[directory.directory_id] = directory
        self.tagger.tag_resource(directory.directory_id, tags or [])
        return directory.directory_id

    def delete_directory(self, directory_id: str) -> str:
        """Delete directory with the matching ID."""
        self._validate_directory_id(directory_id)
        self.directories[directory_id].delete_eni()
        self.directories[directory_id].delete_security_group()
        self.tagger.delete_all_tags_for_resource(directory_id)
        self.directories.pop(directory_id)
        return directory_id

    def disable_sso(
        self,
        directory_id: str,
        username: Optional[str] = None,
        password: Optional[str] = None,
    ) -> None:
        """Disable single-sign on for a directory."""
        self._validate_directory_id(directory_id)
        validate_args([("ssoPassword", password), ("userName", username)])
        directory = self.directories[directory_id]
        directory.enable_sso(False)

    def enable_sso(
        self,
        directory_id: str,
        username: Optional[str] = None,
        password: Optional[str] = None,
    ) -> None:
        """Enable single-sign on for a directory."""
        self._validate_directory_id(directory_id)
        validate_args([("ssoPassword", password), ("userName", username)])

        directory = self.directories[directory_id]
        if directory.alias == directory_id:
            raise ClientException(
                f"An alias is required before enabling SSO. DomainId={directory_id}"
            )

        directory = self.directories[directory_id]
        directory.enable_sso(True)

    @paginate(pagination_model=PAGINATION_MODEL)
    def describe_directories(
        self, directory_ids: Optional[List[str]] = None
    ) -> List[Directory]:
        """Return info on all directories or directories with matching IDs."""
        for directory_id in directory_ids or self.directories:
            self._validate_directory_id(directory_id)

        directories = list(self.directories.values())
        if directory_ids:
            directories = [x for x in directories if x.directory_id in directory_ids]
        return sorted(directories, key=lambda x: x.launch_time)

    def get_directory_limits(self) -> Dict[str, Any]:
        """Return hard-coded limits for the directories."""
        counts = {"SimpleAD": 0, "MicrosoftAD": 0, "ConnectedAD": 0}
        for directory in self.directories.values():
            if directory.directory_type == "SimpleAD":
                counts["SimpleAD"] += 1
            elif directory.directory_type in ["MicrosoftAD", "SharedMicrosoftAD"]:
                counts["MicrosoftAD"] += 1
            elif directory.directory_type == "ADConnector":
                counts["ConnectedAD"] += 1

        return {
            "CloudOnlyDirectoriesLimit": Directory.CLOUDONLY_DIRECTORIES_LIMIT,
            "CloudOnlyDirectoriesCurrentCount": counts["SimpleAD"],
            "CloudOnlyDirectoriesLimitReached": counts["SimpleAD"]
            == Directory.CLOUDONLY_DIRECTORIES_LIMIT,
            "CloudOnlyMicrosoftADLimit": Directory.CLOUDONLY_MICROSOFT_AD_LIMIT,
            "CloudOnlyMicrosoftADCurrentCount": counts["MicrosoftAD"],
            "CloudOnlyMicrosoftADLimitReached": counts["MicrosoftAD"]
            == Directory.CLOUDONLY_MICROSOFT_AD_LIMIT,
            "ConnectedDirectoriesLimit": Directory.CONNECTED_DIRECTORIES_LIMIT,
            "ConnectedDirectoriesCurrentCount": counts["ConnectedAD"],
            "ConnectedDirectoriesLimitReached": counts["ConnectedAD"]
            == Directory.CONNECTED_DIRECTORIES_LIMIT,
        }

    def add_tags_to_resource(
        self, resource_id: str, tags: List[Dict[str, str]]
    ) -> None:
        """Add or overwrite one or more tags for specified directory."""
        self._validate_directory_id(resource_id)
        errmsg = self.tagger.validate_tags(tags)
        if errmsg:
            raise ValidationException(errmsg)
        if len(tags) > Directory.MAX_TAGS_PER_DIRECTORY:
            raise TagLimitExceededException("Tag limit exceeded")
        self.tagger.tag_resource(resource_id, tags)

    def remove_tags_from_resource(self, resource_id: str, tag_keys: List[str]) -> None:
        """Removes tags from a directory."""
        self._validate_directory_id(resource_id)
        self.tagger.untag_resource_using_names(resource_id, tag_keys)

    @paginate(pagination_model=PAGINATION_MODEL)
    def list_tags_for_resource(self, resource_id: str) -> List[Dict[str, str]]:
        """List all tags on a directory."""
        self._validate_directory_id(resource_id)
        return self.tagger.list_tags_for_resource(resource_id).get("Tags")  # type: ignore[return-value]


ds_backends = BackendDict(DirectoryServiceBackend, service_name="ds")
