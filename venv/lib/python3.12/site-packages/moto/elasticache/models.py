import copy
import random
import string
from re import compile as re_compile
from typing import Any, Dict, List, Optional

from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.common_models import BaseModel
from moto.core.utils import utcnow
from moto.utilities.paginator import paginate
from moto.utilities.utils import get_partition

from ..moto_api._internal import mock_random
from .exceptions import (
    CacheClusterAlreadyExists,
    CacheClusterNotFound,
    CacheSubnetGroupAlreadyExists,
    CacheSubnetGroupNotFound,
    InvalidARNFault,
    InvalidParameterCombinationException,
    InvalidParameterValueException,
    InvalidSubnet,
    ReplicationGroupAlreadyExists,
    ReplicationGroupNotFound,
    UserAlreadyExists,
    UserNotFound,
)
from .utils import PAGINATION_MODEL, AuthenticationTypes


class User(BaseModel):
    def __init__(
        self,
        account_id: str,
        region: str,
        user_id: str,
        user_name: str,
        access_string: str,
        engine: str,
        no_password_required: bool,
        passwords: Optional[List[str]] = None,
        authentication_type: Optional[str] = None,
        tags: Optional[List[Dict[str, str]]] = None,
    ):
        self.id = user_id
        self.name = user_name
        self.engine = engine

        self.passwords = passwords or []
        self.access_string = access_string
        self.no_password_required = no_password_required
        self.status = "active"
        self.minimum_engine_version = "6.0"
        self.user_group_ids: List[str] = []
        self.region = region
        self.arn = f"arn:{get_partition(self.region)}:elasticache:{self.region}:{account_id}:user:{self.id}"
        self.authentication_type = authentication_type
        self.tags = tags or []

    @property
    def authentication(self) -> dict[str, Any]:
        return {
            "Type": "no-password"
            if self.no_password_required
            else self.authentication_type,
            "PasswordCount": len(self.passwords) if self.passwords else None,
        }

    def get_tags(self) -> List[Dict[str, str]]:
        return self.tags


class CacheCluster(BaseModel):
    def __init__(
        self,
        account_id: str,
        region_name: str,
        cache_cluster_id: str,
        replication_group_id: Optional[str],
        az_mode: Optional[str],
        preferred_availability_zone: Optional[str],
        num_cache_nodes: Optional[int],
        cache_node_type: Optional[str],
        engine: Optional[str],
        engine_version: Optional[str],
        cache_parameter_group_name: Optional[str],
        cache_subnet_group_name: Optional[str],
        transit_encryption_enabled: Optional[bool],
        network_type: Optional[str],
        ip_discovery: Optional[str],
        snapshot_name: Optional[str],
        preferred_maintenance_window: Optional[str],
        port: Optional[int],
        notification_topic_arn: Optional[str],
        auto_minor_version_upgrade: Optional[bool],
        snapshot_retention_limit: Optional[int],
        snapshot_window: Optional[str],
        auth_token: Optional[str],
        outpost_mode: Optional[str],
        preferred_outpost_arn: Optional[str],
        preferred_availability_zones: Optional[List[str]],
        cache_security_group_names: Optional[List[str]],
        security_group_ids: Optional[List[str]],
        tags: Optional[List[Dict[str, str]]],
        snapshot_arns: Optional[List[str]],
        preferred_outpost_arns: Optional[List[str]],
        log_delivery_configurations: List[Dict[str, Any]],
        cache_node_ids_to_remove: Optional[List[str]],
        cache_node_ids_to_reboot: Optional[List[str]],
    ):
        if tags is None:
            tags = []
        self.cache_cluster_id = cache_cluster_id
        self.az_mode = az_mode
        self.preferred_availability_zone = preferred_availability_zone
        self.preferred_availability_zones = preferred_availability_zones or []
        self.engine = engine or "redis"
        self.engine_version = engine_version
        if engine == "redis":
            self.num_cache_nodes = 1
            self.replication_group_id = replication_group_id
            self.snapshot_arns = snapshot_arns or []
            self.snapshot_name = snapshot_name
            self.snapshot_window = snapshot_window
        if engine == "memcached":
            if num_cache_nodes is None:
                self.num_cache_nodes = 1
            elif 1 <= num_cache_nodes <= 40:
                self.num_cache_nodes = num_cache_nodes
        self.cache_node_type = cache_node_type
        self.cache_parameter_group_name = cache_parameter_group_name
        self.cache_subnet_group_name = cache_subnet_group_name
        self.cache_security_group_names = cache_security_group_names or []
        self.security_group_ids = security_group_ids or []
        self.tags = tags
        self.preferred_maintenance_window = preferred_maintenance_window
        self.port = port or 6379
        self.notification_topic_arn = notification_topic_arn
        self.auto_minor_version_upgrade = auto_minor_version_upgrade
        self.snapshot_retention_limit = snapshot_retention_limit or 0
        self.auth_token = auth_token
        self.outpost_mode = outpost_mode
        self.preferred_outpost_arn = preferred_outpost_arn
        self.preferred_outpost_arns = preferred_outpost_arns or []
        self.log_delivery_configurations = log_delivery_configurations or []
        self.transit_encryption_enabled = transit_encryption_enabled
        self.network_type = network_type
        self.ip_discovery = ip_discovery
        self.cache_node_ids_to_remove = cache_node_ids_to_remove
        self.cache_node_ids_to_reboot = cache_node_ids_to_reboot

        self.cache_cluster_create_time = utcnow()
        self.auth_token_last_modified_date = utcnow()
        self.cache_cluster_status = "available"
        self.arn = f"arn:{get_partition(region_name)}:elasticache:{region_name}:{account_id}:cluster:{cache_cluster_id}"
        self.cache_node_id = str(mock_random.uuid4())
        self.status = "available"

    def get_tags(self) -> List[Dict[str, str]]:
        return self.tags


class CacheSubnetGroup(BaseModel):
    def __init__(
        self,
        account_id: str,
        region_name: str,
        cache_subnet_group_name: str,
        cache_subnet_group_description: str,
        subnet_ids: List[str],
        tags: Optional[List[Dict[str, str]]],
    ):
        if tags is None:
            tags = []
        self.cache_subnet_group_name = cache_subnet_group_name
        self.cache_subnet_group_description = cache_subnet_group_description
        self.subnet_ids = subnet_ids
        self.tags = tags

        # Only import ec2_backends if necessary
        from moto.ec2.models import ec2_backends

        ec2_backend = ec2_backends[account_id][region_name]
        self.supported_network_types = []
        self.subnets = []
        vpc_exists = False
        try:
            # Get VPC details from provided subnet IDs
            subnets = ec2_backend.describe_subnets(subnet_ids=subnet_ids)
            vpc_exists = True
        except Exception as e:
            # Should raise InvalidSubnet if subnet_ids are invalid
            if "InvalidSubnet" in str(e):
                for subnet_id in subnet_ids:
                    subnet_response: Dict[str, Any] = {}
                    subnet_response["subnet_identifier"] = subnet_id
                    subnet_response["Subnet_availability_zone"] = {"Name": "us-east-1a"}
                    subnet_response["supported_network_types"] = ["ipv4"]
                    self.subnets.append(subnet_response)
                vpcs = ["vpc-0123456789abcdef0"]
                self.supported_network_types = ["ipv4"]

        if vpc_exists:
            vpcs = []
            for subnet in subnets:
                subnet_response = {}
                vpcs.append(subnet.vpc_id)
                subnet_response["subnet_identifier"] = subnet.id
                subnet_response["subnet_availability_zone"] = subnet.availability_zone
                subnet_response["supported_network_types"] = []
                if subnet.vpc_id != vpcs[0]:
                    raise InvalidSubnet(subnet_id=subnet.id)

                # ipv6 native subnets only appends ipv6
                # You can't mix ipv6 native subnets with other types of subnets
                if subnet.ipv6_native:
                    self.supported_network_types.append("ipv6")
                    subnet_response["supported_network_types"].append("ipv6")

                # ipv4 only and dual_stack subnets both append ipv4
                elif subnet.cidr_block:
                    self.supported_network_types.append("ipv4")
                    subnet_response["supported_network_types"].append("ipv4")

                if subnet.ipv6_cidr_block_associations and subnet.cidr_block:
                    self.supported_network_types.append("dual_stack")
                    subnet_response["supported_network_types"].append("dual_stack")

                self.subnets.append(subnet_response)

            if self.supported_network_types:
                self.supported_network_types = list(set(self.supported_network_types))

        self.arn = f"arn:aws:elasticache:{region_name}:{account_id}:subnetgroup:{cache_subnet_group_name}"
        self.vpc_id = vpcs[0] if vpcs else None

    def get_tags(self) -> List[Dict[str, str]]:
        return self.tags


class ReplicationGroup(BaseModel):
    def __init__(
        self,
        account_id: str,
        region_name: str,
        replication_group_id: str,
        preferred_cache_cluster_azs: Optional[List[str]],
        num_cache_clusters: Optional[int],
        num_node_groups: Optional[int],
        replicas_per_node_group: Optional[int],
        node_group_configuration: List[Dict[str, Any]],
        preferred_maintenance_window: Optional[str],
        replication_group_description: str,
        primary_cluster_id: Optional[str],
        automatic_failover_enabled: bool,
        global_replication_group_id: Optional[str],
        multi_az_enabled: Optional[bool],
        port: Optional[int],
        snapshot_retention_limit: int,
        snapshot_window: Optional[str],
        cluster_mode: Optional[str],
        cache_node_type: Optional[str],
        auth_token: Optional[str],
        transit_encryption_enabled: Optional[bool],
        at_rest_encryption_enabled: Optional[bool],
        kms_key_id: Optional[str],
        user_group_ids: List[str],
        log_delivery_configurations: List[Dict[str, Any]],
        data_tiering_enabled: Optional[bool],
        auto_minor_version_upgrade: Optional[bool],
        engine: Optional[str],
        network_type: Optional[str],
        ip_discovery: Optional[str],
        transit_encryption_mode: Optional[str],
        cache_security_group_names: Optional[List[str]],
        cache_subnet_group_name: Optional[str],
        security_group_ids: Optional[List[str]],
        tags: Optional[List[Dict[str, str]]],
        notification_topic_arn: Optional[str],
        serverless_cache_snapshot_name: Optional[str],
        cache_parameter_group_name: Optional[str],
        engine_version: Optional[str],
        snapshot_arns: Optional[List[str]],
        snapshot_name: Optional[str],
    ):
        tags = tags or []
        self.cluster_mode = cluster_mode or "disabled"
        self.replication_group_id = replication_group_id
        self.description = replication_group_description
        self.port = port or 6379
        self.status = "available"

        preferred_cache_cluster_azs = preferred_cache_cluster_azs or []
        replicas_per_node_group = replicas_per_node_group or 0

        # This field is only present if we are creating a
        # secondary replication group and associating it with an existing
        # global replication group and primary cluster.
        if global_replication_group_id:
            self.global_replication_group_id = global_replication_group_id
            self.global_replication_group_member_role = "SECONDARY"
        else:
            self.global_replication_group_id = ""

        # Leaving off PendingModifiedValues for now
        # until global_replication_group is implemented.
        self.primary_cluster_id = primary_cluster_id

        self.node_groups = []
        random_str = "".join(
            random.choices(string.ascii_lowercase + string.digits, k=6)
        )
        region_short = "use1"
        replication_group_domain = (
            f"{replication_group_id}.{random_str}.{region_short}.cache.amazonaws.com"
        )
        self.member_clusters: List[str] = []
        self.member_clusters_outpost_arns: List[str] = []
        if not num_node_groups:
            num_node_groups = 1

        for i in range(num_node_groups):
            node_group: Dict[str, Any] = {}
            node_group["node_group_id"] = f"{i + 1:0>4}"
            node_group["status"] = "available"
            node_group["node_group_members"] = []

            if self.cluster_mode == "enabled":
                if i == 0:
                    node_group["slots"] = "0-5461"
                else:
                    # slots cannot overlap among node groups
                    node_group["slots"] = (
                        f"{(5460 * (i)) + (i + 1)}-{(5460 * (i + 1)) + (i + 1)}"
                    )

                replicas_per_node_group = replicas_per_node_group or 0
                self._set_node_members_clusters_enabled(
                    member_clusters=self.member_clusters,
                    replication_group_id=replication_group_id,
                    replicas_per_node_group=replicas_per_node_group,
                    node_group_configuration=node_group_configuration,
                    node_group=node_group,
                    member_clusters_outpost_arns=self.member_clusters_outpost_arns,
                )

                self.configuration_endpoint: Dict[str, Any] = {}
                self.configuration_endpoint["address"] = (
                    f"clustercfg.{replication_group_domain}"
                )
                self.configuration_endpoint["port"] = self.port

            # self.cluster_mode is disabled or not set
            else:
                node_group["primary_endpoint_address"] = (
                    f"master.{replication_group_domain}"
                )
                node_group["reader_endpoint_address"] = (
                    f"replica.{replication_group_domain}"
                )
                node_group["port"] = self.port
                node_group["primary_endpoint"] = {
                    "address": node_group["primary_endpoint_address"],
                    "port": self.port,
                }

                num_cache_clusters = num_cache_clusters or 1

                self._set_node_members_clusters_disabled(
                    member_clusters=self.member_clusters,
                    replication_group_id=replication_group_id,
                    replication_group_domain=replication_group_domain,
                    node_group=node_group,
                    preferred_cache_cluster_azs=preferred_cache_cluster_azs,
                    num_cache_clusters=num_cache_clusters,
                    replicas_per_node_group=replicas_per_node_group,
                )
            self.node_groups.append(node_group)

        self.automatic_failover = (
            "enabled" if automatic_failover_enabled else "disabled"
        )
        self.multi_az = "enabled" if multi_az_enabled else "disabled"
        self.snapshot_retention_limit = snapshot_retention_limit or 0
        self.snapshot_window = snapshot_window or "05:00-6:00"
        self.snapshotting_cluster_id = None
        if self.snapshot_retention_limit > 0:
            if len(self.member_clusters) > 1:
                self.snapshotting_cluster_id = self.member_clusters[1]
            else:
                self.snapshotting_cluster_id = self.member_clusters[0]
        self.cluster_enabled = True if self.cluster_mode == "enabled" else False
        self.cache_node_type = cache_node_type or "cache.t4g.micro"
        self.auth_token_enabled = bool(auth_token)
        self.auth_token_last_modified_date = (
            utcnow() if self.auth_token_enabled else None
        )
        self.transit_encryption_enabled = transit_encryption_enabled or False
        self.at_rest_encryption_enabled = at_rest_encryption_enabled or False
        self.kms_key_id = kms_key_id
        self.arn = f"arn:aws:elasticache:{region_name}:{account_id}:replicationgroup:{replication_group_id}"
        self.user_group_ids = user_group_ids or []
        self.log_delivery_configurations = self._get_log_delivery_configurations(
            log_delivery_configurations
        )
        self.replication_group_create_time = utcnow()
        self.data_tiering = "enabled" if data_tiering_enabled else "disabled"
        self.auto_minor_version_upgrade = auto_minor_version_upgrade
        self.network_type = network_type or "ipv4"
        self.ip_discovery = ip_discovery or "ipv4"
        self.transit_encryption_mode = transit_encryption_mode
        self.engine = engine or "redis"

        # Parameters not returned by create or describe replication group
        self.cache_security_group_names = cache_security_group_names
        self.cache_subnet_group_name = cache_subnet_group_name
        self.security_group_ids = security_group_ids
        self.tags = tags
        self.preferred_maintenance_window = preferred_maintenance_window
        self.notification_topic_arn = notification_topic_arn
        self.serverless_cache_snapshot_name = serverless_cache_snapshot_name
        self.cache_parameter_group_name = cache_parameter_group_name
        self.engine_version = engine_version
        self.snapshot_arns = snapshot_arns
        self.snapshot_name = snapshot_name

    @property
    def global_replication_group_info(self) -> dict[str, str]:
        return {
            "GlobalReplicationGroupId": self.global_replication_group_id,
            "GlobalReplicationGroupMemberRole": self.global_replication_group_member_role,
        }

    def get_tags(self) -> List[Dict[str, str]]:
        return self.tags

    def _get_log_delivery_configurations(
        self, log_delivery_configurations: List[Dict[str, Any]]
    ) -> List[Dict[str, Any]]:
        log_delivery_configurations_resp = []
        if log_delivery_configurations:
            for log_delivery_configuration in log_delivery_configurations:
                log_config_resp = copy.copy(log_delivery_configuration)
                log_config_resp["status"] = "active"
                log_config_resp["message"] = "Log delivery configuration is active."

                log_delivery_configurations_resp.append(log_config_resp)
        return log_delivery_configurations_resp

    # This method will set populate the node_group with node_group_members
    # and member_clusters with the cache_cluster_ids of the primary and replicas.
    def _set_node_members_clusters_enabled(
        self,
        member_clusters: List[str],
        replication_group_id: str,
        replicas_per_node_group: int,
        node_group_configuration: List[Dict[str, Any]],
        node_group: Dict[str, Any],
        member_clusters_outpost_arns: List[str],
    ) -> None:
        replica_count = replicas_per_node_group
        primary_node = {}
        primary_node["cache_cluster_id"] = (
            f"{replication_group_id}-{node_group['node_group_id']}-001"
        )
        primary_node["cache_node_id"] = node_group["node_group_id"]

        # Set the PreferredAvailabilityZone,
        # default to us-east-1a for moto only
        primary_node["preferred_availability_zone"] = "us-east-1a"

        # With node_group_configurations
        node_config = {}
        if node_group_configuration:
            for config in node_group_configuration:
                if node_group["node_group_id"] == config.get("NodeGroupId"):
                    node_config = config
                    # overwrite slots if part of node config
                    if config.get("Slots"):
                        node_group["slots"] = config.get("Slots")

                    replica_count_val = config.get("ReplicaCount")
                    if replica_count_val is not None:
                        replica_count = int(replica_count_val)
                    else:
                        replica_count = replicas_per_node_group

                    primary_node["preferred_availability_zone"] = config.get(
                        "PrimaryAvailabilityZone", "us-east-1a"
                    )
                    primary_outpost_arn = config.get("PrimaryOutpostArn", "")
                    if primary_outpost_arn:
                        primary_node["preferred_outpost_arn"] = primary_outpost_arn
                        member_clusters_outpost_arns.append(primary_outpost_arn)

        node_group["node_group_members"].append(primary_node)
        member_clusters.append(primary_node["cache_cluster_id"])

        #  Loop through the number of replicas per node group
        for r in range(1, replica_count + 1):
            replica_node = {}
            # r + 1 because primary is 001
            replica_node["cache_cluster_id"] = (
                f"{replication_group_id}-{node_group['node_group_id']}-{r + 1:0>3}"
            )
            replica_node["cache_node_id"] = node_group["node_group_id"]

            # Set the ReplicaAvailabiilityZone,
            # default to us-east-1b for moto only
            replica_node["preferred_availability_zone"] = "us-east-1b"

            if node_config:
                replica_az = node_config.get("ReplicaAvailabilityZones")
                if replica_az:
                    if len(replica_az) >= r:
                        replica_node["preferred_availability_zone"] = replica_az[r - 1]
                    elif len(replica_az) >= 1:
                        replica_node["preferred_availability_zone"] = replica_az[1]
                replica_outpost_arns = node_config.get("ReplicaOutpostArns", [])
                if replica_outpost_arns:
                    if len(replica_outpost_arns) >= r:
                        replica_node["preferred_outpost_arn"] = replica_outpost_arns[
                            r - 1
                        ]
                    elif len(replica_outpost_arns) >= 1:
                        replica_node["preferred_outpost_arn"] = replica_outpost_arns[1]

                    member_clusters_outpost_arns.append(
                        replica_node["preferred_outpost_arn"]
                    )

            node_group["node_group_members"].append(replica_node)
            member_clusters.append(replica_node["cache_cluster_id"])

    def _set_node_members_clusters_disabled(
        self,
        member_clusters: List[str],
        replication_group_id: str,
        replication_group_domain: str,
        node_group: Dict[str, Any],
        preferred_cache_cluster_azs: List[str],
        num_cache_clusters: int,
        replicas_per_node_group: int,
    ) -> None:
        primary_node: Dict[str, Any] = {}
        primary_node["cache_cluster_id"] = f"{replication_group_id}-001"
        primary_node["cache_node_id"] = node_group["node_group_id"]
        primary_node["read_endpoint"] = {}
        primary_node["read_endpoint"]["address"] = (
            f"{primary_node['cache_cluster_id']}.{replication_group_domain}"
        )
        primary_node["read_endpoint"]["port"] = self.port

        if preferred_cache_cluster_azs:
            primary_node["preferred_availability_zone"] = preferred_cache_cluster_azs[0]
        else:
            primary_node["preferred_availability_zone"] = "us-east-1a"

        primary_node["current_role"] = "primary"
        node_group["node_group_members"].append(primary_node)
        member_clusters.append(primary_node["cache_cluster_id"])

        if replicas_per_node_group:
            num_cache_clusters = replicas_per_node_group + 1
        else:
            num_cache_clusters = 1

        # Use num_cache_clusters when cluster_mode is disabled
        for r in range(1, num_cache_clusters):
            replica_node: Dict[str, Any] = {}
            # r + 1 because primary is 001
            replica_node["cache_cluster_id"] = f"{replication_group_id}-{r + 1:0>3}"
            replica_node["cache_node_id"] = node_group["node_group_id"]
            replica_node["read_endpoint"] = {}
            replica_node["read_endpoint"]["address"] = (
                f"{replica_node['cache_cluster_id']}.{replication_group_domain}"
            )
            replica_node["read_endpoint"]["port"] = self.port

            if preferred_cache_cluster_azs:
                if len(preferred_cache_cluster_azs) > r:
                    replica_az = preferred_cache_cluster_azs[r]
                elif len(preferred_cache_cluster_azs) >= 2:
                    replica_az = preferred_cache_cluster_azs[1]
                else:
                    replica_az = preferred_cache_cluster_azs[0]
            else:
                replica_az = "us-east-1b"

            replica_node["preferred_availability_zone"] = replica_az
            replica_node["current_role"] = "replica"

            node_group["node_group_members"].append(replica_node)
            member_clusters.append(replica_node["cache_cluster_id"])


class ElastiCacheBackend(BaseBackend):
    """Implementation of ElastiCache APIs."""

    def __init__(self, region_name: str, account_id: str):
        super().__init__(region_name, account_id)
        self.arn_regex = re_compile(
            r"^arn:aws:elasticache:.*:[0-9]*:(cluster|snapshot|subnetgroup|replicationgroup|user):.*$"
        )
        self.users = dict()
        self.users["default"] = User(
            account_id=self.account_id,
            region=self.region_name,
            user_id="default",
            user_name="default",
            engine="redis",
            access_string="on ~* +@all",
            no_password_required=True,
        )

        self.cache_clusters: Dict[str, Any] = dict()
        self.cache_subnet_groups: Dict[str, CacheSubnetGroup] = dict()
        self.replication_groups: Dict[str, ReplicationGroup] = dict()

    def create_user(
        self,
        user_id: str,
        user_name: str,
        engine: str,
        passwords: List[str],
        access_string: str,
        no_password_required: bool,
        authentication_type: str,  # contain it to the str in the enums TODO
        tags: Optional[List[Dict[str, str]]] = None,
    ) -> User:
        if user_id in self.users:
            raise UserAlreadyExists

        if authentication_type not in AuthenticationTypes._value2member_map_:
            raise InvalidParameterValueException(
                f"Input Authentication type: {authentication_type} is not in the allowed list: [password,no-password-required,iam]"
            )

        if (
            no_password_required
            and authentication_type != AuthenticationTypes.NOPASSWORD
        ):
            raise InvalidParameterCombinationException(
                f"No password required flag is true but provided authentication type is {authentication_type}"
            )

        if passwords and authentication_type != AuthenticationTypes.PASSWORD:
            raise InvalidParameterCombinationException(
                f"Password field is not allowed with authentication type: {authentication_type}"
            )

        if not passwords and authentication_type == AuthenticationTypes.PASSWORD:
            raise InvalidParameterCombinationException(
                "A user with Authentication Mode: password, must have at least one password"
            )

        user = User(
            account_id=self.account_id,
            region=self.region_name,
            user_id=user_id,
            user_name=user_name,
            engine=engine,
            passwords=passwords,
            access_string=access_string,
            no_password_required=no_password_required,
            authentication_type=authentication_type,
            tags=tags,
        )
        self.users[user_id] = user
        return user

    def delete_user(self, user_id: str) -> User:
        if user_id in self.users:
            user = self.users[user_id]
            if user.status == "active":
                user.status = "deleting"
            return user
        raise UserNotFound(user_id)

    def describe_users(self, user_id: Optional[str]) -> List[User]:
        """
        Only the `user_id` parameter is currently supported.
        Pagination is not yet implemented.
        """
        if user_id:
            if user_id in self.users:
                user = self.users[user_id]
                if user.status == "deleting":
                    self.users.pop(user_id)
                return [user]
            else:
                raise UserNotFound(user_id)
        return list(self.users.values())

    def create_cache_cluster(
        self,
        cache_cluster_id: str,
        replication_group_id: str,
        az_mode: str,
        preferred_availability_zone: str,
        num_cache_nodes: int,
        cache_node_type: str,
        engine: str,
        engine_version: str,
        cache_parameter_group_name: str,
        cache_subnet_group_name: str,
        transit_encryption_enabled: bool,
        network_type: str,
        ip_discovery: str,
        snapshot_name: str,
        preferred_maintenance_window: str,
        port: int,
        notification_topic_arn: str,
        auto_minor_version_upgrade: bool,
        snapshot_retention_limit: int,
        snapshot_window: str,
        auth_token: str,
        outpost_mode: str,
        preferred_outpost_arn: str,
        preferred_availability_zones: List[str],
        cache_security_group_names: List[str],
        security_group_ids: List[str],
        tags: List[Dict[str, str]],
        snapshot_arns: List[str],
        preferred_outpost_arns: List[str],
        log_delivery_configurations: List[Dict[str, Any]],
        cache_node_ids_to_remove: List[str],
        cache_node_ids_to_reboot: List[str],
    ) -> CacheCluster:
        if cache_cluster_id in self.cache_clusters:
            raise CacheClusterAlreadyExists(cache_cluster_id)
        cache_cluster = CacheCluster(
            account_id=self.account_id,
            region_name=self.region_name,
            cache_cluster_id=cache_cluster_id,
            replication_group_id=replication_group_id,
            az_mode=az_mode,
            preferred_availability_zone=preferred_availability_zone,
            preferred_availability_zones=preferred_availability_zones,
            num_cache_nodes=num_cache_nodes,
            cache_node_type=cache_node_type,
            engine=engine,
            engine_version=engine_version,
            cache_parameter_group_name=cache_parameter_group_name,
            cache_subnet_group_name=cache_subnet_group_name,
            cache_security_group_names=cache_security_group_names,
            security_group_ids=security_group_ids,
            tags=tags,
            snapshot_arns=snapshot_arns,
            snapshot_name=snapshot_name,
            preferred_maintenance_window=preferred_maintenance_window,
            port=port,
            notification_topic_arn=notification_topic_arn,
            auto_minor_version_upgrade=auto_minor_version_upgrade,
            snapshot_retention_limit=snapshot_retention_limit,
            snapshot_window=snapshot_window,
            auth_token=auth_token,
            outpost_mode=outpost_mode,
            preferred_outpost_arn=preferred_outpost_arn,
            preferred_outpost_arns=preferred_outpost_arns,
            log_delivery_configurations=log_delivery_configurations,
            transit_encryption_enabled=transit_encryption_enabled,
            network_type=network_type,
            ip_discovery=ip_discovery,
            cache_node_ids_to_remove=cache_node_ids_to_remove,
            cache_node_ids_to_reboot=cache_node_ids_to_reboot,
        )
        self.cache_clusters[cache_cluster_id] = cache_cluster
        return cache_cluster

    @paginate(PAGINATION_MODEL)
    def describe_cache_clusters(
        self,
        cache_cluster_id: str,
        max_records: int,
        marker: str,
    ) -> List[CacheCluster]:
        if max_records is None:
            max_records = 100
        if cache_cluster_id:
            if cache_cluster_id in self.cache_clusters:
                cache_cluster = self.cache_clusters[cache_cluster_id]
                return list([cache_cluster])
            else:
                raise CacheClusterNotFound(cache_cluster_id)
        cache_clusters = list(self.cache_clusters.values())[:max_records]

        return cache_clusters

    def delete_cache_cluster(self, cache_cluster_id: str) -> CacheCluster:
        if cache_cluster_id:
            if cache_cluster_id in self.cache_clusters:
                cache_cluster = self.cache_clusters[cache_cluster_id]
                cache_cluster.cache_cluster_status = "deleting"
                return cache_cluster
        raise CacheClusterNotFound(cache_cluster_id)

    def create_cache_subnet_group(
        self,
        cache_subnet_group_name: str,
        cache_subnet_group_description: str,
        subnet_ids: List[str],
        tags: Optional[List[Dict[str, str]]],
    ) -> CacheSubnetGroup:
        if cache_subnet_group_name in self.cache_subnet_groups:
            raise CacheSubnetGroupAlreadyExists(cache_subnet_group_name)

        cache_subnet_group = CacheSubnetGroup(
            account_id=self.account_id,
            region_name=self.region_name,
            cache_subnet_group_name=cache_subnet_group_name,
            cache_subnet_group_description=cache_subnet_group_description,
            subnet_ids=subnet_ids,
            tags=tags,
        )
        self.cache_subnet_groups[cache_subnet_group_name] = cache_subnet_group
        return cache_subnet_group

    @paginate(PAGINATION_MODEL)
    def describe_cache_subnet_groups(
        self,
        cache_subnet_group_name: str,
    ) -> List[CacheSubnetGroup]:
        if cache_subnet_group_name:
            if cache_subnet_group_name in self.cache_subnet_groups:
                cache_subnet_group = self.cache_subnet_groups[cache_subnet_group_name]
                return list([cache_subnet_group])
            else:
                raise CacheSubnetGroupNotFound(cache_subnet_group_name)
        cache_subnet_groups = list(self.cache_subnet_groups.values())
        return cache_subnet_groups

    def list_tags_for_resource(self, arn: str) -> List[Dict[str, str]]:
        if self.arn_regex.match(arn):
            arn_breakdown = arn.split(":")
            resource_type = arn_breakdown[len(arn_breakdown) - 2]
            resource_name = arn_breakdown[len(arn_breakdown) - 1]
            if resource_type == "cluster":
                if resource_name in self.cache_clusters:
                    return self.cache_clusters[resource_name].get_tags()
            elif resource_type == "subnetgroup":
                if resource_name in self.cache_subnet_groups:
                    return self.cache_subnet_groups[resource_name].get_tags()
            elif resource_type == "replicationgroup":
                if resource_name in self.replication_groups:
                    return self.replication_groups[resource_name].get_tags()
            elif resource_type == "user":
                if resource_name in self.users:
                    return self.users[resource_name].get_tags()
            else:
                return []
        else:
            raise InvalidARNFault(arn)
        return []

    def create_replication_group(
        self,
        replication_group_id: str,
        replication_group_description: str,
        global_replication_group_id: str,
        primary_cluster_id: str,
        automatic_failover_enabled: bool,
        multi_az_enabled: bool,
        num_cache_clusters: int,
        preferred_cache_cluster_azs: List[str],
        num_node_groups: int,
        replicas_per_node_group: int,
        node_group_configuration: List[Dict[str, Any]],
        cache_node_type: str,
        engine: str,
        engine_version: str,
        cache_parameter_group_name: str,
        cache_subnet_group_name: str,
        cache_security_group_names: List[str],
        security_group_ids: List[str],
        tags: List[Dict[str, str]],
        snapshot_arns: List[str],
        snapshot_name: str,
        preferred_maintenance_window: str,
        port: int,
        notification_topic_arn: str,
        auto_minor_version_upgrade: bool,
        snapshot_retention_limit: int,
        snapshot_window: str,
        auth_token: str,
        transit_encryption_enabled: bool,
        at_rest_encryption_enabled: bool,
        kms_key_id: str,
        user_group_ids: List[str],
        log_delivery_configurations: List[Dict[str, Any]],
        data_tiering_enabled: bool,
        network_type: str,
        ip_discovery: str,
        transit_encryption_mode: str,
        cluster_mode: str,
        serverless_cache_snapshot_name: str,
    ) -> ReplicationGroup:
        if replication_group_id in self.replication_groups:
            raise ReplicationGroupAlreadyExists(replication_group_id)
        replication_group = ReplicationGroup(
            account_id=self.account_id,
            region_name=self.region_name,
            replication_group_id=replication_group_id,
            replication_group_description=replication_group_description,
            global_replication_group_id=global_replication_group_id,
            primary_cluster_id=primary_cluster_id,
            automatic_failover_enabled=automatic_failover_enabled,
            multi_az_enabled=multi_az_enabled,
            num_cache_clusters=num_cache_clusters,
            preferred_cache_cluster_azs=preferred_cache_cluster_azs,
            num_node_groups=num_node_groups,
            replicas_per_node_group=replicas_per_node_group,
            node_group_configuration=node_group_configuration,
            cache_node_type=cache_node_type,
            engine=engine,
            engine_version=engine_version,
            cache_parameter_group_name=cache_parameter_group_name,
            cache_subnet_group_name=cache_subnet_group_name,
            cache_security_group_names=cache_security_group_names,
            security_group_ids=security_group_ids,
            tags=tags,
            snapshot_arns=snapshot_arns,
            snapshot_name=snapshot_name,
            preferred_maintenance_window=preferred_maintenance_window,
            port=port,
            notification_topic_arn=notification_topic_arn,
            auto_minor_version_upgrade=auto_minor_version_upgrade,
            snapshot_retention_limit=snapshot_retention_limit,
            snapshot_window=snapshot_window,
            auth_token=auth_token,
            transit_encryption_enabled=transit_encryption_enabled,
            at_rest_encryption_enabled=at_rest_encryption_enabled,
            kms_key_id=kms_key_id,
            user_group_ids=user_group_ids,
            log_delivery_configurations=log_delivery_configurations,
            data_tiering_enabled=data_tiering_enabled,
            network_type=network_type,
            ip_discovery=ip_discovery,
            transit_encryption_mode=transit_encryption_mode,
            cluster_mode=cluster_mode,
            serverless_cache_snapshot_name=serverless_cache_snapshot_name,
        )

        self.replication_groups[replication_group_id] = replication_group
        return replication_group

    @paginate(PAGINATION_MODEL)
    def describe_replication_groups(
        self, replication_group_id: str
    ) -> List[ReplicationGroup]:
        if replication_group_id:
            if replication_group_id in self.replication_groups:
                replication_group = self.replication_groups[replication_group_id]
                return list([replication_group])
            else:
                raise ReplicationGroupNotFound(replication_group_id)
        replication_groups = list(self.replication_groups.values())
        return replication_groups

    def delete_replication_group(
        self, replication_group_id: str, retain_primary_cluster: Optional[bool]
    ) -> ReplicationGroup:
        if replication_group_id:
            if replication_group_id in self.replication_groups:
                replication_group = self.replication_groups[replication_group_id]
                replication_group.status = "deleting"
            if not retain_primary_cluster:
                replication_group = self.replication_groups[replication_group_id]
                primary_id = replication_group.primary_cluster_id
                if primary_id and primary_id in replication_group.member_clusters:
                    replication_group.member_clusters.remove(primary_id)
                replication_group.primary_cluster_id = ""
            return replication_group
        raise ReplicationGroupNotFound(replication_group_id)


elasticache_backends = BackendDict(ElastiCacheBackend, "elasticache")
