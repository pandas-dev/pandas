from moto.core.responses import BaseResponse

from .exceptions import PasswordRequired, PasswordTooShort
from .models import ElastiCacheBackend, elasticache_backends


class ElastiCacheResponse(BaseResponse):
    """Handler for ElastiCache requests and responses."""

    def __init__(self) -> None:
        super().__init__(service_name="elasticache")

    @property
    def elasticache_backend(self) -> ElastiCacheBackend:
        """Return backend instance specific for this region."""
        return elasticache_backends[self.current_account][self.region]

    def create_user(self) -> str:
        params = self._get_params()
        user_id = params.get("UserId")
        user_name = params.get("UserName")
        engine = params.get("Engine")
        passwords = params.get("Passwords", [])
        no_password_required = self._get_bool_param("NoPasswordRequired", False)
        password_required = not no_password_required
        if password_required and not passwords:
            raise PasswordRequired
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
        )
        template = self.response_template(CREATE_USER_TEMPLATE)
        return template.render(user=user)

    def delete_user(self) -> str:
        params = self._get_params()
        user_id = params.get("UserId")
        user = self.elasticache_backend.delete_user(user_id=user_id)  # type: ignore[arg-type]
        template = self.response_template(DELETE_USER_TEMPLATE)
        return template.render(user=user)

    def describe_users(self) -> str:
        params = self._get_params()
        user_id = params.get("UserId")
        users = self.elasticache_backend.describe_users(user_id=user_id)
        template = self.response_template(DESCRIBE_USERS_TEMPLATE)
        return template.render(users=users)

    def create_cache_cluster(self) -> str:
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
        tags = (self._get_multi_param_dict("Tags") or {}).get("Tag", [])
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
        template = self.response_template(CREATE_CACHE_CLUSTER_TEMPLATE)
        return template.render(cache_cluster=cache_cluster)

    def describe_cache_clusters(self) -> str:
        cache_cluster_id = self._get_param("CacheClusterId")
        max_records = self._get_int_param("MaxRecords")
        marker = self._get_param("Marker")

        cache_clusters, marker = self.elasticache_backend.describe_cache_clusters(
            cache_cluster_id=cache_cluster_id,
            marker=marker,
            max_records=max_records,
        )
        template = self.response_template(DESCRIBE_CACHE_CLUSTERS_TEMPLATE)
        return template.render(marker=marker, cache_clusters=cache_clusters)

    def delete_cache_cluster(self) -> str:
        cache_cluster_id = self._get_param("CacheClusterId")
        cache_cluster = self.elasticache_backend.delete_cache_cluster(
            cache_cluster_id=cache_cluster_id,
        )
        template = self.response_template(DELETE_CACHE_CLUSTER_TEMPLATE)
        return template.render(cache_cluster=cache_cluster)

    def list_tags_for_resource(self) -> str:
        arn = self._get_param("ResourceName")
        template = self.response_template(LIST_TAGS_FOR_RESOURCE_TEMPLATE)
        tags = self.elasticache_backend.list_tags_for_resource(arn)
        return template.render(tags=tags)


USER_TEMPLATE = """<UserId>{{ user.id }}</UserId>
    <UserName>{{ user.name }}</UserName>
    <Status>{{ user.status }}</Status>
    <Engine>{{ user.engine }}</Engine>
    <MinimumEngineVersion>{{ user.minimum_engine_version }}</MinimumEngineVersion>
    <AccessString>{{ user.access_string }}</AccessString>
    <UserGroupIds>
{% for usergroupid in user.usergroupids %}
      <member>{{ usergroupid }}</member>
{% endfor %}
    </UserGroupIds>
    <Authentication>
      {% if user.no_password_required %}
      <Type>no-password</Type>
      {% else %}
      <Type>password</Type>
      <PasswordCount>{{ user.passwords|length }}</PasswordCount>
      {% endif %}
    </Authentication>
    <ARN>{{ user.arn }}</ARN>"""

CREATE_USER_TEMPLATE = (
    """<CreateUserResponse xmlns="http://elasticache.amazonaws.com/doc/2015-02-02/">
          <ResponseMetadata>
            <RequestId>1549581b-12b7-11e3-895e-1334aEXAMPLE</RequestId>
          </ResponseMetadata>
          <CreateUserResult>
            """
    + USER_TEMPLATE
    + """
  </CreateUserResult>
</CreateUserResponse>"""
)

DELETE_USER_TEMPLATE = (
    """<DeleteUserResponse xmlns="http://elasticache.amazonaws.com/doc/2015-02-02/">
          <ResponseMetadata>
            <RequestId>1549581b-12b7-11e3-895e-1334aEXAMPLE</RequestId>
          </ResponseMetadata>
          <DeleteUserResult>
            """
    + USER_TEMPLATE
    + """
  </DeleteUserResult>
</DeleteUserResponse>"""
)

DESCRIBE_USERS_TEMPLATE = (
    """<DescribeUsersResponse xmlns="http://elasticache.amazonaws.com/doc/2015-02-02/">
          <ResponseMetadata>
            <RequestId>1549581b-12b7-11e3-895e-1334aEXAMPLE</RequestId>
          </ResponseMetadata>
          <DescribeUsersResult>
            <Users>
        {% for user in users %}
              <member>
                """
    + USER_TEMPLATE
    + """
      </member>
{% endfor %}
    </Users>
    <Marker></Marker>
  </DescribeUsersResult>
</DescribeUsersResponse>"""
)

CREATE_CACHE_CLUSTER_TEMPLATE = """<CreateCacheClusterResponse xmlns="http://elasticache.amazonaws.com/doc/2015-02-02/">
  <ResponseMetadata>
    <RequestId>1549581b-12b7-11e3-895e-1334aEXAMPLE</RequestId>
  </ResponseMetadata>
  <CreateCacheClusterResult>
    <CacheCluster>
  <CacheClusterId>{{ cache_cluster.cache_cluster_id }}</CacheClusterId>
  <ConfigurationEndpoint>
    <Address>example.cache.amazonaws.com</Address>
    <Port>{{ cache_cluster.port }}</Port>
  </ConfigurationEndpoint>
  <ClientDownloadLandingPage></ClientDownloadLandingPage>
  <CacheNodeType>{{ cache_cluster.cache_node_type }}</CacheNodeType>
  <Engine>{{ cache_cluster.engine }}</Engine>
  <EngineVersion>{{ cache_cluster.engine_version }}</EngineVersion>
  <CacheClusterStatus>available</CacheClusterStatus>
  <NumCacheNodes>{{ cache_cluster.num_cache_nodes }}</NumCacheNodes>
  <PreferredAvailabilityZone>{{ cache_cluster.preferred_availability_zone }}</PreferredAvailabilityZone>
  <PreferredOutpostArn>{{ cache_cluster.preferred_outpost_arn }}</PreferredOutpostArn>
  <CacheClusterCreateTime>{{ cache_cluster.cache_cluster_create_time }}</CacheClusterCreateTime>
  <PreferredMaintenanceWindow>{{ cache_cluster.preferred_maintenance_window }}</PreferredMaintenanceWindow>
  {% if cache_cluster.cache_node_ids_to_remove != [] %}
  <PendingModifiedValues>
    <NumCacheNodes>{{ cache_cluster.num_cache_nodes }}</NumCacheNodes>
    {% for cache_node_id_to_remove in cache_cluster.cache_node_ids_to_remove %}
    <CacheNodeIdsToRemove>{{ cache_node_id_to_remove }}</CacheNodeIdsToRemove>
    {% endfor %}
    <EngineVersion>{{ cache_cluster.engine_version }}</EngineVersion>
    <CacheNodeType>{{ cache_cluster.cache_node_type }}</CacheNodeType>
    <AuthTokenStatus>SETTING</AuthTokenStatus>
    <LogDeliveryConfigurations>
    {% for log_delivery_configuration in cache_cluster.log_delivery_configurations %}
      <LogType>{{ log_delivery_configuration.LogType }}</LogType>
      <DestinationType>{{ log_delivery_configuration.DestinationType }}</DestinationType>
      <DestinationDetails>
        <CloudWatchLogsDetails>
          <LogGroup>{{ log_delivery_configuration.LogGroup }}</LogGroup>
        </CloudWatchLogsDetails>
        <KinesisFirehoseDetails>
          <DeliveryStream>{{ log_delivery_configuration.DeliveryStream }}</DeliveryStream>
        </KinesisFirehoseDetails>
      </DestinationDetails>
      <LogFormat>{{ log_delivery_configuration.LogFormat }}</LogFormat>
      {% endfor %}
    </LogDeliveryConfigurations>
    <TransitEncryptionEnabled>{{ cache_cluster.transit_encryption_enabled }}</TransitEncryptionEnabled>
    <TransitEncryptionMode>preferred</TransitEncryptionMode>
  </PendingModifiedValues>
  {% endif %}
  <NotificationConfiguration>
    <TopicArn>{{ cache_cluster.notification_topic_arn }}</TopicArn>
    <TopicStatus>active</TopicStatus>
  </NotificationConfiguration>
  <CacheSecurityGroups>
  {% for cache_security_group_name in cache_cluster.cache_security_group_names %}
    <CacheSecurityGroupName>{{ cache_security_group_name }}</CacheSecurityGroupName>
    {% endfor %}
    <Status>active</Status>
  </CacheSecurityGroups>
  <CacheParameterGroup>
    <CacheParameterGroupName>{{ cache_cluster.cache_parameter_group_name }}</CacheParameterGroupName>
    <ParameterApplyStatus>active</ParameterApplyStatus>
    {% for cache_node_id_to_reboot in cache_cluster.cache_node_ids_to_reboot %}
    <CacheNodeIdsToReboot>
    {{ cache_node_id_to_reboot }}
    </CacheNodeIdsToReboot>
    {% endfor %}
  </CacheParameterGroup>
  <CacheSubnetGroupName>{{ cache_cluster.cache_subnet_group_name }}</CacheSubnetGroupName>
  <CacheNodes>
    <CacheNodeId>{{ cache_cluster.cache_node_id }}</CacheNodeId>
    <CacheNodeStatus>{{ cache_cluster.cache_node_status }}</CacheNodeStatus>
    <CacheNodeCreateTime>{{ cache_cluster.cache_cluster_create_time }}</CacheNodeCreateTime>
    <Endpoint>
      <Address>{{ cache_cluster.address }}</Address>
      <Port>{{ cache_cluster.port }}</Port>
    </Endpoint>
    <ParameterGroupStatus>active</ParameterGroupStatus>
    <SourceCacheNodeId>{{ cache_cluster.cache_node_id }}</SourceCacheNodeId>
    <CustomerAvailabilityZone>{{ cache_cluster.preferred_availability_zone }}</CustomerAvailabilityZone>
    <CustomerOutpostArn>{{ cache_cluster.preferred_output_arn }}</CustomerOutpostArn>
  </CacheNodes>
  <AutoMinorVersionUpgrade>{{ cache_cluster.auto_minor_version_upgrade }}</AutoMinorVersionUpgrade>
  <SecurityGroups>
  {% for security_group_id in cache_cluster.security_group_ids %}
    <SecurityGroupId>{{ security_group_id }}</SecurityGroupId>
    <Status>active</Status>
    {% endfor %}
  </SecurityGroups>
  {% if cache_cluster.engine == "redis" %}
  <ReplicationGroupId>{{ cache_cluster.replication_group_id }}</ReplicationGroupId>
  <SnapshotRetentionLimit>{{ cache_cluster.snapshot_retention_limit }}</SnapshotRetentionLimit>
  <SnapshotWindow>{{ cache_cluster.snapshot_window }}</SnapshotWindow>
  {% endif %}
  <AuthTokenEnabled>true</AuthTokenEnabled>
  <AuthTokenLastModifiedDate>{{ cache_cluster.cache_cluster_create_time }}</AuthTokenLastModifiedDate>
  <TransitEncryptionEnabled>{{ cache_cluster.transit_encryption_enabled }}</TransitEncryptionEnabled>
  <AtRestEncryptionEnabled>true</AtRestEncryptionEnabled>
  <ARN>{{ cache_cluster.arn }}</ARN>
  <ReplicationGroupLogDeliveryEnabled>true</ReplicationGroupLogDeliveryEnabled>
  <LogDeliveryConfigurations>
  {% for log_delivery_configuration in cache_cluster.log_delivery_configurations %}
      <LogType>{{ log_delivery_configuration.LogType }}</LogType>
      <DestinationType>{{ log_delivery_configuration.DestinationType }}</DestinationType>
      <DestinationDetails>
        <CloudWatchLogsDetails>
          <LogGroup>{{ log_delivery_configuration.LogGroup }}</LogGroup>
        </CloudWatchLogsDetails>
        <KinesisFirehoseDetails>
          <DeliveryStream>{{ log_delivery_configuration.DeliveryStream }}</DeliveryStream>
        </KinesisFirehoseDetails>
      </DestinationDetails>
      <LogFormat>{{ log_delivery_configuration.LogFormat }}</LogFormat>
      <Status>active</Status>
      <Message></Message>
  {% endfor %}
  </LogDeliveryConfigurations>
  <NetworkType>{{ cache_cluster.network_type }}</NetworkType>
  <IpDiscovery>{{ cache_cluster.ip_discovery }}</IpDiscovery>
  <TransitEncryptionMode>preferred</TransitEncryptionMode>
</CacheCluster>
  </CreateCacheClusterResult>
</CreateCacheClusterResponse>"""

DESCRIBE_CACHE_CLUSTERS_TEMPLATE = """<DescribeCacheClustersResponse xmlns="http://elasticache.amazonaws.com/doc/2015-02-02/">
  <ResponseMetadata>
    <RequestId>1549581b-12b7-11e3-895e-1334aEXAMPLE</RequestId>
  </ResponseMetadata>
  <DescribeCacheClustersResult>
    {% if marker %}<Marker>{{ marker }}</Marker>{% endif %}
    <CacheClusters>
{% for cache_cluster in cache_clusters %}
      <member>
        <CacheClusterId>{{ cache_cluster.cache_cluster_id }}</CacheClusterId>
        <ConfigurationEndpoint>{{ cache_cluster.configuration_endpoint }}</ConfigurationEndpoint>
        <ClientDownloadLandingPage>{{ cache_cluster.client_download_landing_page }}</ClientDownloadLandingPage>
        <CacheNodeType>{{ cache_cluster.cache_node_type }}</CacheNodeType>
        <Engine>{{ cache_cluster.engine }}</Engine>
        <EngineVersion>{{ cache_cluster.engine_version }}</EngineVersion>
        <CacheClusterStatus>{{ cache_cluster.cache_cluster_status }}</CacheClusterStatus>
        <NumCacheNodes>{{ cache_cluster.num_cache_nodes }}</NumCacheNodes>
        <PreferredAvailabilityZone>{{ cache_cluster.preferred_availability_zone }}</PreferredAvailabilityZone>
        <PreferredOutpostArn>{{ cache_cluster.preferred_outpost_arn }}</PreferredOutpostArn>
        <CacheClusterCreateTime>{{ cache_cluster.cache_cluster_create_time }}</CacheClusterCreateTime>
        <PreferredMaintenanceWindow>{{ cache_cluster.preferred_maintenance_window }}</PreferredMaintenanceWindow>
        <PendingModifiedValues>{{ cache_cluster.pending_modified_values }}</PendingModifiedValues>
        <NotificationConfiguration>{{ cache_cluster.notification_configuration }}</NotificationConfiguration>
        <CacheSecurityGroups>
{% for cache_security_group in cache_cluster.cache_security_groups %}
          <member>
            <CacheSecurityGroupName>{{ cache_security_group.cache_security_group_name }}</CacheSecurityGroupName>
            <Status>{{ cache_security_group.status }}</Status>
          </member>
{% endfor %}
        </CacheSecurityGroups>
        <CacheParameterGroup>{{ cache_cluster.cache_parameter_group }}</CacheParameterGroup>
        <CacheSubnetGroupName>{{ cache_cluster.cache_subnet_group_name }}</CacheSubnetGroupName>
        <CacheNodes>
{% for cache_node in cache_cluster.cache_nodes %}
          <member>
            <CacheNodeId>{{ cache_node.cache_node_id }}</CacheNodeId>
            <CacheNodeStatus>{{ cache_node.cache_node_status }}</CacheNodeStatus>
            <CacheNodeCreateTime>{{ cache_node.cache_node_create_time }}</CacheNodeCreateTime>
            <Endpoint>{{ cache_node.endpoint }}</Endpoint>
            <ParameterGroupStatus>{{ cache_node.parameter_group_status }}</ParameterGroupStatus>
            <SourceCacheNodeId>{{ cache_node.source_cache_node_id }}</SourceCacheNodeId>
            <CustomerAvailabilityZone>{{ cache_node.customer_availability_zone }}</CustomerAvailabilityZone>
            <CustomerOutpostArn>{{ cache_node.customer_outpost_arn }}</CustomerOutpostArn>
          </member>
{% endfor %}
        </CacheNodes>
        <AutoMinorVersionUpgrade>{{ cache_cluster.auto_minor_version_upgrade }}</AutoMinorVersionUpgrade>
        <SecurityGroups>
{% for security_group in cache_cluster.security_groups %}
          <member>
            <SecurityGroupId>{{ security_group.security_group_id }}</SecurityGroupId>
            <Status>{{ security_group.status }}</Status>
          </member>
{% endfor %}
        </SecurityGroups>
        <ReplicationGroupId>{{ cache_cluster.replication_group_id }}</ReplicationGroupId>
        <SnapshotRetentionLimit>{{ cache_cluster.snapshot_retention_limit }}</SnapshotRetentionLimit>
        <SnapshotWindow>{{ cache_cluster.snapshot_window }}</SnapshotWindow>
        <AuthTokenEnabled>{{ cache_cluster.auth_token_enabled }}</AuthTokenEnabled>
        <AuthTokenLastModifiedDate>{{ cache_cluster.auth_token_last_modified_date }}</AuthTokenLastModifiedDate>
        <TransitEncryptionEnabled>{{ cache_cluster.transit_encryption_enabled }}</TransitEncryptionEnabled>
        <AtRestEncryptionEnabled>{{ cache_cluster.at_rest_encryption_enabled }}</AtRestEncryptionEnabled>
        <ARN>{{ cache_cluster.arn }}</ARN>
        <ReplicationGroupLogDeliveryEnabled>{{ cache_cluster.replication_group_log_delivery_enabled }}</ReplicationGroupLogDeliveryEnabled>
        <LogDeliveryConfigurations>
{% for log_delivery_configuration in cache_cluster.log_delivery_configurations %}
          <member>
            <LogType>{{ log_delivery_configuration.log_type }}</LogType>
            <DestinationType>{{ log_delivery_configuration.destination_type }}</DestinationType>
            <DestinationDetails>{{ log_delivery_configuration.destination_details }}</DestinationDetails>
            <LogFormat>{{ log_delivery_configuration.log_format }}</LogFormat>
            <Status>{{ log_delivery_configuration.status }}</Status>
            <Message>{{ log_delivery_configuration.message }}</Message>
          </member>
{% endfor %}
        </LogDeliveryConfigurations>
        <NetworkType>{{ cache_cluster.network_type }}</NetworkType>
        <IpDiscovery>{{ cache_cluster.ip_discovery }}</IpDiscovery>
        <TransitEncryptionMode>{{ cache_cluster.transit_encryption_mode }}</TransitEncryptionMode>
      </member>
{% endfor %}
    </CacheClusters>
  </DescribeCacheClustersResult>
</DescribeCacheClustersResponse>"""

DELETE_CACHE_CLUSTER_TEMPLATE = """<DeleteCacheClusterResponse xmlns="http://elasticache.amazonaws.com/doc/2015-02-02/">
  <ResponseMetadata>
    <RequestId>1549581b-12b7-11e3-895e-1334aEXAMPLE</RequestId>
  </ResponseMetadata>
  <DeleteCacheClusterResult>
    <CacheCluster>
  <CacheClusterId>{{ cache_cluster.cache_cluster_id }}</CacheClusterId>
  <ConfigurationEndpoint>
    <Address>example.cache.amazonaws.com</Address>
    <Port>{{ cache_cluster.port }}</Port>
  </ConfigurationEndpoint>
  <ClientDownloadLandingPage></ClientDownloadLandingPage>
  <CacheNodeType>{{ cache_cluster.cache_node_type }}</CacheNodeType>
  <Engine>{{ cache_cluster.engine }}</Engine>
  <EngineVersion>{{ cache_cluster.engine_version }}</EngineVersion>
  <CacheClusterStatus>available</CacheClusterStatus>
  <NumCacheNodes>{{ cache_cluster.num_cache_nodes }}</NumCacheNodes>
  <PreferredAvailabilityZone>{{ cache_cluster.preferred_availability_zone }}</PreferredAvailabilityZone>
  <PreferredOutpostArn>{{ cache_cluster.preferred_outpost_arn }}</PreferredOutpostArn>
  <CacheClusterCreateTime>{{ cache_cluster.cache_cluster_create_time }}</CacheClusterCreateTime>
  <PreferredMaintenanceWindow>{{ cache_cluster.preferred_maintenance_window }}</PreferredMaintenanceWindow>
  {% if cache_cluster.cache_node_ids_to_remove != [] %}
  <PendingModifiedValues>
    <NumCacheNodes>{{ cache_cluster.num_cache_nodes }}</NumCacheNodes>
    {% for cache_node_id_to_remove in cache_cluster.cache_node_ids_to_remove %}
    <CacheNodeIdsToRemove>{{ cache_node_id_to_remove }}</CacheNodeIdsToRemove>
    {% endfor %}
    <EngineVersion>{{ cache_cluster.engine_version }}</EngineVersion>
    <CacheNodeType>{{ cache_cluster.cache_node_type }}</CacheNodeType>
    <AuthTokenStatus>SETTING</AuthTokenStatus>
    <LogDeliveryConfigurations>
    {% for log_delivery_configuration in cache_cluster.log_delivery_configurations %}
      <LogType>{{ log_delivery_configuration.LogType }}</LogType>
      <DestinationType>{{ log_delivery_configuration.DestinationType }}</DestinationType>
      <DestinationDetails>
        <CloudWatchLogsDetails>
          <LogGroup>{{ log_delivery_configuration.LogGroup }}</LogGroup>
        </CloudWatchLogsDetails>
        <KinesisFirehoseDetails>
          <DeliveryStream>{{ log_delivery_configuration.DeliveryStream }}</DeliveryStream>
        </KinesisFirehoseDetails>
      </DestinationDetails>
      <LogFormat>{{ log_delivery_configuration.LogFormat }}</LogFormat>
      {% endfor %}
    </LogDeliveryConfigurations>
    <TransitEncryptionEnabled>{{ cache_cluster.transit_encryption_enabled }}</TransitEncryptionEnabled>
    <TransitEncryptionMode>preferred</TransitEncryptionMode>
  </PendingModifiedValues>
  {% endif %}
  <NotificationConfiguration>
    <TopicArn>{{ cache_cluster.notification_topic_arn }}</TopicArn>
    <TopicStatus>active</TopicStatus>
  </NotificationConfiguration>
  <CacheSecurityGroups>
  {% for cache_security_group_name in cache_cluster.cache_security_group_names %}
    <CacheSecurityGroupName>{{ cache_security_group_name }}</CacheSecurityGroupName>
    {% endfor %}
    <Status>active</Status>
  </CacheSecurityGroups>
  <CacheParameterGroup>
    <CacheParameterGroupName>{{ cache_cluster.cache_parameter_group_name }}</CacheParameterGroupName>
    <ParameterApplyStatus>active</ParameterApplyStatus>
    {% for cache_node_id_to_reboot in cache_cluster.cache_node_ids_to_reboot %}
    <CacheNodeIdsToReboot>
    {{ cache_node_id_to_reboot }}
    </CacheNodeIdsToReboot>
    {% endfor %}
  </CacheParameterGroup>
  <CacheSubnetGroupName>{{ cache_cluster.cache_subnet_group_name }}</CacheSubnetGroupName>
  <CacheNodes>
    <CacheNodeId>{{ cache_cluster.cache_node_id }}</CacheNodeId>
    <CacheNodeStatus>{{ cache_cluster.cache_node_status }}</CacheNodeStatus>
    <CacheNodeCreateTime>{{ cache_cluster.cache_cluster_create_time }}</CacheNodeCreateTime>
    <Endpoint>
      <Address>{{ cache_cluster.address }}</Address>
      <Port>{{ cache_cluster.port }}</Port>
    </Endpoint>
    <ParameterGroupStatus>active</ParameterGroupStatus>
    <SourceCacheNodeId>{{ cache_cluster.cache_node_id }}</SourceCacheNodeId>
    <CustomerAvailabilityZone>{{ cache_cluster.preferred_availability_zone }}</CustomerAvailabilityZone>
    <CustomerOutpostArn>{{ cache_cluster.preferred_output_arn }}</CustomerOutpostArn>
  </CacheNodes>
  <AutoMinorVersionUpgrade>{{ cache_cluster.auto_minor_version_upgrade }}</AutoMinorVersionUpgrade>
  <SecurityGroups>
  {% for security_group_id in cache_cluster.security_group_ids %}
    <SecurityGroupId>{{ security_group_id }}</SecurityGroupId>
    <Status>active</Status>
    {% endfor %}
  </SecurityGroups>
  {% if cache_cluster.engine == "redis" %}
  <ReplicationGroupId>{{ cache_cluster.replication_group_id }}</ReplicationGroupId>
  <SnapshotRetentionLimit>{{ cache_cluster.snapshot_retention_limit }}</SnapshotRetentionLimit>
  <SnapshotWindow>{{ cache_cluster.snapshot_window }}</SnapshotWindow>
  {% endif %}
  <AuthTokenEnabled>true</AuthTokenEnabled>
  <AuthTokenLastModifiedDate>{{ cache_cluster.cache_cluster_create_time }}</AuthTokenLastModifiedDate>
  <TransitEncryptionEnabled>{{ cache_cluster.transit_encryption_enabled }}</TransitEncryptionEnabled>
  <AtRestEncryptionEnabled>true</AtRestEncryptionEnabled>
  <ARN>{{ cache_cluster.arn }}</ARN>
  <ReplicationGroupLogDeliveryEnabled>true</ReplicationGroupLogDeliveryEnabled>
  <LogDeliveryConfigurations>
  {% for log_delivery_configuration in cache_cluster.log_delivery_configurations %}
      <LogType>{{ log_delivery_configuration.LogType }}</LogType>
      <DestinationType>{{ log_delivery_configuration.DestinationType }}</DestinationType>
      <DestinationDetails>
        <CloudWatchLogsDetails>
          <LogGroup>{{ log_delivery_configuration.LogGroup }}</LogGroup>
        </CloudWatchLogsDetails>
        <KinesisFirehoseDetails>
          <DeliveryStream>{{ log_delivery_configuration.DeliveryStream }}</DeliveryStream>
        </KinesisFirehoseDetails>
      </DestinationDetails>
      <LogFormat>{{ log_delivery_configuration.LogFormat }}</LogFormat>
      <Status>active</Status>
      <Message></Message>
  {% endfor %}
  </LogDeliveryConfigurations>
  <NetworkType>{{ cache_cluster.network_type }}</NetworkType>
  <IpDiscovery>{{ cache_cluster.ip_discovery }}</IpDiscovery>
  <TransitEncryptionMode>preferred</TransitEncryptionMode>
</CacheCluster>
  </DeleteCacheClusterResult>
</DeleteCacheClusterResponse>"""

LIST_TAGS_FOR_RESOURCE_TEMPLATE = """<ListTagsForResourceResponse xmlns="http://elasticache.amazonaws.com/doc/2015-02-02/">
  <ListTagsForResourceResult>
    <TagList>
    {%- for tag in tags -%}
      <Tag>
        <Key>{{ tag['Key'] }}</Key>
        <Value>{{ tag['Value'] }}</Value>
      </Tag>
    {%- endfor -%}
    </TagList>
  </ListTagsForResourceResult>
  <ResponseMetadata>
    <RequestId>8c21ba39-a598-11e4-b688-194eaf8658fa</RequestId>
  </ResponseMetadata>
</ListTagsForResourceResponse>"""
