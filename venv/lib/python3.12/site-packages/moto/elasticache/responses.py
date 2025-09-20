from moto.core.responses import ActionResult, BaseResponse

from .exceptions import (
    InvalidParameterCombinationException,
    InvalidParameterValueException,
    PasswordTooShort,
)
from .models import ElastiCacheBackend, elasticache_backends
from .utils import VALID_AUTH_MODE_KEYS, VALID_ENGINE_TYPES, AuthenticationTypes


class ElastiCacheResponse(BaseResponse):
    """Handler for ElastiCache requests and responses."""

    def __init__(self) -> None:
        super().__init__(service_name="elasticache")
        self.automated_parameter_parsing = True

    @property
    def elasticache_backend(self) -> ElastiCacheBackend:
        """Return backend instance specific for this region."""
        return elasticache_backends[self.current_account][self.region]

    def create_user(self) -> ActionResult:
        params = self._get_params()
        user_id = params.get("UserId")
        user_name = params.get("UserName")
        engine = params.get("Engine", "").lower()
        passwords = params.get("Passwords", [])
        no_password_required = self._get_bool_param("NoPasswordRequired")
        authentication_mode = params.get("AuthenticationMode")
        authentication_type = "null"
        tags = params.get("Tags")

        if no_password_required is not None:
            authentication_type = (
                AuthenticationTypes.NOPASSWORD.value
                if no_password_required
                else AuthenticationTypes.PASSWORD.value
            )

        if passwords:
            authentication_type = AuthenticationTypes.PASSWORD.value

        if engine not in VALID_ENGINE_TYPES:
            raise InvalidParameterValueException(
                f'Unknown parameter for Engine: "{engine}", must be one of: {", ".join(VALID_ENGINE_TYPES)}'
            )

        if authentication_mode:
            for key, _ in authentication_mode.original_items():
                if key not in VALID_AUTH_MODE_KEYS:
                    raise InvalidParameterValueException(
                        f'Unknown parameter in AuthenticationMode: "{key}", must be one of: {", ".join(VALID_AUTH_MODE_KEYS)}'
                    )

            authentication_type = authentication_mode.get("Type")
            authentication_passwords = authentication_mode.get("Passwords", [])

            if passwords and authentication_passwords:
                raise InvalidParameterCombinationException(
                    "Passwords provided via multiple arguments. Use only one argument"
                )

            # if passwords is empty, then we can use the authentication_passwords
            passwords = passwords if passwords else authentication_passwords

        if any([len(p) < 16 for p in passwords]):
            raise PasswordTooShort

        access_string = params.get("AccessString")
        user = self.elasticache_backend.create_user(
            user_id=user_id,  # type: ignore[arg-type]
            user_name=user_name,  # type: ignore[arg-type]
            engine=engine,  # type: ignore[arg-type]
            passwords=passwords,
            access_string=access_string,  # type: ignore[arg-type]
            no_password_required=no_password_required,
            authentication_type=authentication_type,
            tags=tags,
        )
        return ActionResult(user)

    def delete_user(self) -> ActionResult:
        params = self._get_params()
        user_id = params.get("UserId")
        user = self.elasticache_backend.delete_user(user_id=user_id)  # type: ignore[arg-type]
        return ActionResult(user)

    def describe_users(self) -> ActionResult:
        params = self._get_params()
        user_id = params.get("UserId")
        users = self.elasticache_backend.describe_users(user_id=user_id)
        return ActionResult({"Users": users})

    def create_cache_cluster(self) -> ActionResult:
        cache_cluster_id = self._get_param("CacheClusterId")
        replication_group_id = self._get_param("ReplicationGroupId")
        az_mode = self._get_param("AZMode")
        preferred_availability_zone = self._get_param("PreferredAvailabilityZone")
        preferred_availability_zones = self._get_param("PreferredAvailabilityZones")
        num_cache_nodes = self._get_int_param("NumCacheNodes")
        cache_node_type = self._get_param("CacheNodeType")
        engine = self._get_param("Engine")
        engine_version = self._get_param("EngineVersion")
        cache_parameter_group_name = self._get_param("CacheParameterGroupName")
        cache_subnet_group_name = self._get_param("CacheSubnetGroupName")
        cache_security_group_names = self._get_param("CacheSecurityGroupNames")
        security_group_ids = self._get_param("SecurityGroupIds")
        tags = self._get_param("Tags", [])
        snapshot_arns = self._get_param("SnapshotArns")
        snapshot_name = self._get_param("SnapshotName")
        preferred_maintenance_window = self._get_param("PreferredMaintenanceWindow")
        port = self._get_param("Port")
        notification_topic_arn = self._get_param("NotificationTopicArn")
        auto_minor_version_upgrade = self._get_bool_param("AutoMinorVersionUpgrade")
        snapshot_retention_limit = self._get_int_param("SnapshotRetentionLimit")
        snapshot_window = self._get_param("SnapshotWindow")
        auth_token = self._get_param("AuthToken")
        outpost_mode = self._get_param("OutpostMode")
        preferred_outpost_arn = self._get_param("PreferredOutpostArn")
        preferred_outpost_arns = self._get_param("PreferredOutpostArns")
        log_delivery_configurations = self._get_param("LogDeliveryConfigurations")
        transit_encryption_enabled = self._get_bool_param("TransitEncryptionEnabled")
        network_type = self._get_param("NetworkType")
        ip_discovery = self._get_param("IpDiscovery")
        # Define the following attributes as they're included in the response even during creation of a cache cluster
        cache_node_ids_to_remove = self._get_param("CacheNodeIdsToRemove", [])
        cache_node_ids_to_reboot = self._get_param("CacheNodeIdsToReboot", [])
        cache_cluster = self.elasticache_backend.create_cache_cluster(
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
        return ActionResult({"CacheCluster": cache_cluster})

    def describe_cache_clusters(self) -> ActionResult:
        cache_cluster_id = self._get_param("CacheClusterId")
        max_records = self._get_int_param("MaxRecords")
        marker = self._get_param("Marker")
        cache_clusters, marker = self.elasticache_backend.describe_cache_clusters(
            cache_cluster_id=cache_cluster_id,
            marker=marker,
            max_records=max_records,
        )
        result = {"CacheClusters": cache_clusters, "Marker": marker}
        return ActionResult(result)

    def delete_cache_cluster(self) -> ActionResult:
        cache_cluster_id = self._get_param("CacheClusterId")
        cache_cluster = self.elasticache_backend.delete_cache_cluster(
            cache_cluster_id=cache_cluster_id,
        )
        return ActionResult({"CacheCluster": cache_cluster})

    def list_tags_for_resource(self) -> ActionResult:
        arn = self._get_param("ResourceName")
        tags = self.elasticache_backend.list_tags_for_resource(arn)
        return ActionResult({"TagList": tags})

    def create_cache_subnet_group(self) -> ActionResult:
        cache_subnet_group_name = self._get_param("CacheSubnetGroupName")
        cache_subnet_group_description = self._get_param("CacheSubnetGroupDescription")
        subnet_ids = self._get_param("SubnetIds", [])
        tags = self._get_param("Tags", [])
        cache_subnet_group = self.elasticache_backend.create_cache_subnet_group(
            cache_subnet_group_name=cache_subnet_group_name,
            cache_subnet_group_description=cache_subnet_group_description,
            subnet_ids=subnet_ids,
            tags=tags,
        )
        return ActionResult({"CacheSubnetGroup": cache_subnet_group})

    def describe_cache_subnet_groups(self) -> ActionResult:
        cache_subnet_group_name = self._get_param("CacheSubnetGroupName")
        max_records = self._get_param("MaxRecords")
        marker = self._get_param("Marker")
        cache_subnet_groups, marker = (
            self.elasticache_backend.describe_cache_subnet_groups(
                cache_subnet_group_name=cache_subnet_group_name,
                marker=marker,
                max_records=max_records,
            )
        )
        result = {"CacheSubnetGroups": cache_subnet_groups, "Marker": marker}
        return ActionResult(result)

    def create_replication_group(self) -> ActionResult:
        params = self._get_params()
        replication_group_id = self._get_param("ReplicationGroupId")
        replication_group_description = self._get_param("ReplicationGroupDescription")
        global_replication_group_id = self._get_param("GlobalReplicationGroupId")
        primary_cluster_id = self._get_param("PrimaryClusterId")
        automatic_failover_enabled = self._get_bool_param("AutomaticFailoverEnabled")
        multi_az_enabled = self._get_bool_param("MultiAZEnabled")
        num_cache_clusters = self._get_int_param("NumCacheClusters")
        preferred_cache_cluster_azs = self._get_param("PreferredCacheClusterAZs", [])
        num_node_groups = self._get_int_param("NumNodeGroups")
        replicas_per_node_group = self._get_int_param("ReplicasPerNodeGroup")
        node_group_configuration = self._get_param("NodeGroupConfiguration", [])
        cache_node_type = self._get_param("CacheNodeType")
        engine = self._get_param("Engine")
        engine_version = self._get_param("EngineVersion")
        cache_parameter_group_name = self._get_param("CacheParameterGroupName")
        cache_subnet_group_name = self._get_param("CacheSubnetGroupName")
        cache_security_group_names = self._get_param("CacheSecurityGroupNames")
        security_group_ids = self._get_param("SecurityGroupIds")
        tags = self._get_param("Tags", [])
        snapshot_arns = self._get_param("SnapshotArns")
        snapshot_name = self._get_param("SnapshotName")
        preferred_maintenance_window = self._get_param("PreferredMaintenanceWindow")
        port = self._get_param("Port")
        notification_topic_arn = self._get_param("NotificationTopicArn")
        auto_minor_version_upgrade = self._get_param("AutoMinorVersionUpgrade")
        snapshot_retention_limit = self._get_int_param("SnapshotRetentionLimit")
        snapshot_window = self._get_param("SnapshotWindow")
        auth_token = self._get_param("AuthToken")
        transit_encryption_enabled = self._get_bool_param("TransitEncryptionEnabled")
        at_rest_encryption_enabled = self._get_bool_param("AtRestEncryptionEnabled")
        kms_key_id = self._get_param("KmsKeyId")
        user_group_ids = self._get_param("UserGroupIds")
        log_delivery_configurations = params.get("LogDeliveryConfigurations", [])
        data_tiering_enabled = self._get_param("DataTieringEnabled")
        network_type = self._get_param("NetworkType")
        ip_discovery = self._get_param("IpDiscovery")
        transit_encryption_mode = self._get_param("TransitEncryptionMode")
        cluster_mode = self._get_param("ClusterMode")
        serverless_cache_snapshot_name = self._get_param("ServerlessCacheSnapshotName")
        replication_group = self.elasticache_backend.create_replication_group(
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
        result = {"ReplicationGroup": replication_group}
        return ActionResult(result)

    def describe_replication_groups(self) -> ActionResult:
        replication_group_id = self._get_param("ReplicationGroupId")
        max_records = self._get_param("MaxRecords")
        marker = self._get_param("Marker")
        replication_groups, marker = (
            self.elasticache_backend.describe_replication_groups(
                replication_group_id=replication_group_id,
                marker=marker,
                max_records=max_records,
            )
        )
        result = {"ReplicationGroups": replication_groups, "Marker": marker}
        return ActionResult(result)

    def delete_replication_group(self) -> ActionResult:
        replication_group_id = self._get_param("ReplicationGroupId")
        retain_primary_cluster = self._get_bool_param("RetainPrimaryCluster")
        replication_group = self.elasticache_backend.delete_replication_group(
            replication_group_id=replication_group_id,
            retain_primary_cluster=retain_primary_cluster,
        )
        return ActionResult({"ReplicationGroup": replication_group})
