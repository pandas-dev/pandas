from ._base_response import EC2BaseResponse


class SpotInstances(EC2BaseResponse):
    def cancel_spot_instance_requests(self) -> str:
        request_ids = self._get_multi_param("SpotInstanceRequestId")

        self.error_on_dryrun()

        requests = self.ec2_backend.cancel_spot_instance_requests(request_ids)
        template = self.response_template(CANCEL_SPOT_INSTANCES_TEMPLATE)
        return template.render(requests=requests)

    def create_spot_datafeed_subscription(self) -> None:
        self.error_on_dryrun()

        raise NotImplementedError(
            "SpotInstances.create_spot_datafeed_subscription is not yet implemented"
        )

    def delete_spot_datafeed_subscription(self) -> None:
        self.error_on_dryrun()

        raise NotImplementedError(
            "SpotInstances.delete_spot_datafeed_subscription is not yet implemented"
        )

    def describe_spot_datafeed_subscription(self) -> None:
        raise NotImplementedError(
            "SpotInstances.describe_spot_datafeed_subscription is not yet implemented"
        )

    def describe_spot_instance_requests(self) -> str:
        spot_instance_ids = self._get_multi_param("SpotInstanceRequestId")
        filters = self._filters_from_querystring()
        requests = self.ec2_backend.describe_spot_instance_requests(
            filters=filters, spot_instance_ids=spot_instance_ids
        )
        template = self.response_template(DESCRIBE_SPOT_INSTANCES_TEMPLATE)
        return template.render(requests=requests)

    def describe_spot_price_history(self) -> str:
        instance_types_filters = self._get_multi_param("InstanceType")
        filter_dict = self._filters_from_querystring()
        prices = self.ec2_backend.describe_spot_price_history(
            instance_types_filters, filter_dict
        )
        template = self.response_template(DESCRIBE_SPOT_PRICE_HISTORY_TEMPLATE)
        return template.render(prices=prices)

    def request_spot_instances(self) -> str:
        price = self._get_param("SpotPrice")
        image_id = self._get_param("LaunchSpecification.ImageId")
        count = self._get_int_param("InstanceCount", 1)
        spot_instance_type = self._get_param("Type", "one-time")
        valid_from = self._get_param("ValidFrom")
        valid_until = self._get_param("ValidUntil")
        launch_group = self._get_param("LaunchGroup")
        availability_zone_group = self._get_param("AvailabilityZoneGroup")
        key_name = self._get_param("LaunchSpecification.KeyName")
        security_groups = self._get_multi_param("LaunchSpecification.SecurityGroup")
        user_data = self._get_param("LaunchSpecification.UserData")
        instance_type = self._get_param("LaunchSpecification.InstanceType", "m1.small")
        placement = self._get_param("LaunchSpecification.Placement.AvailabilityZone")
        kernel_id = self._get_param("LaunchSpecification.KernelId")
        ramdisk_id = self._get_param("LaunchSpecification.RamdiskId")
        monitoring_enabled = self._get_param("LaunchSpecification.Monitoring.Enabled")
        subnet_id = self._get_param("LaunchSpecification.SubnetId")
        instance_interruption_behaviour = self._get_param(
            "InstanceInterruptionBehavior"
        )
        tags = self._parse_tag_specification()

        self.error_on_dryrun()

        requests = self.ec2_backend.request_spot_instances(
            price=price,
            image_id=image_id,
            count=count,
            spot_instance_type=spot_instance_type,
            valid_from=valid_from,
            valid_until=valid_until,
            launch_group=launch_group,
            availability_zone_group=availability_zone_group,
            key_name=key_name,
            security_groups=security_groups,
            user_data=user_data,
            instance_type=instance_type,
            placement=placement,
            kernel_id=kernel_id,
            ramdisk_id=ramdisk_id,
            monitoring_enabled=monitoring_enabled,
            subnet_id=subnet_id,
            instance_interruption_behaviour=instance_interruption_behaviour,
            tags=tags,
        )

        template = self.response_template(REQUEST_SPOT_INSTANCES_TEMPLATE)
        return template.render(requests=requests)


REQUEST_SPOT_INSTANCES_TEMPLATE = """<RequestSpotInstancesResponse xmlns="http://ec2.amazonaws.com/doc/2013-10-15/">
  <requestId>59dbff89-35bd-4eac-99ed-be587EXAMPLE</requestId>
  <spotInstanceRequestSet>
    {% for request in requests %}
    <item>
      <spotInstanceRequestId>{{ request.id }}</spotInstanceRequestId>
      <spotPrice>{{ request.price }}</spotPrice>
      <type>{{ request.type }}</type>
      <state>{{ request.state }}</state>
      <status>
        <code>{{ request.status }}</code>
        <updateTime>2015-01-01T00:00:00.000Z</updateTime>
        <message>{{ request.status_message }}</message>
      </status>
      <instanceId>{{ request.instance_id }}</instanceId>
      <availabilityZoneGroup>{{ request.availability_zone_group }}</availabilityZoneGroup>
      <launchSpecification>
        <imageId>{{ request.launch_specification.image_id }}</imageId>
        <keyName>{{ request.launch_specification.key_name }}</keyName>
        <groupSet>
          {% for group in request.launch_specification.groups %}
          <item>
            <groupId>{{ group.id }}</groupId>
            <groupName>{{ group.name }}</groupName>
          </item>
          {% endfor %}
        </groupSet>
        <kernelId>{{ request.launch_specification.kernel }}</kernelId>
        <ramdiskId>{{ request.launch_specification.ramdisk }}</ramdiskId>
        <subnetId>{{ request.launch_specification.subnet_id }}</subnetId>
        <instanceType>{{ request.launch_specification.instance_type }}</instanceType>
        <blockDeviceMapping/>
        <monitoring>
          <enabled>{{ request.launch_specification.monitored }}</enabled>
        </monitoring>
        <ebsOptimized>{{ request.launch_specification.ebs_optimized }}</ebsOptimized>
        <PlacementRequestType>
          <availabilityZone>{{ request.launch_specification.placement }}</availabilityZone>
          <groupName></groupName>
        </PlacementRequestType>
      </launchSpecification>
      <launchGroup>{{ request.launch_group }}</launchGroup>
      <createTime>2015-01-01T00:00:00.000Z</createTime>
      {% if request.valid_from %}
      <validFrom>{{ request.valid_from }}</validFrom>
      {% endif %}
      {% if request.valid_until %}
      <validUntil>{{ request.valid_until }}</validUntil>
      {% endif %}
      <productDescription>Linux/UNIX</productDescription>
    </item>
    {% endfor %}
 </spotInstanceRequestSet>
</RequestSpotInstancesResponse>"""

DESCRIBE_SPOT_INSTANCES_TEMPLATE = """<DescribeSpotInstanceRequestsResponse xmlns="http://ec2.amazonaws.com/doc/2013-10-15/">
  <requestId>59dbff89-35bd-4eac-99ed-be587EXAMPLE</requestId>
  <spotInstanceRequestSet>
    {% for request in requests %}
    <item>
      <spotInstanceRequestId>{{ request.id }}</spotInstanceRequestId>
      <spotPrice>{{ request.price }}</spotPrice>
      <type>{{ request.type }}</type>
      <state>{{ request.state }}</state>
      <status>
        <code>{{ request.status }}</code>
        <updateTime>2015-01-01T00:00:00.000Z</updateTime>
        <message>{{ request.status_message }}</message>
      </status>
      <instanceId>{{ request.instance.id }}</instanceId>
      {% if request.availability_zone_group %}
        <availabilityZoneGroup>{{ request.availability_zone_group }}</availabilityZoneGroup>
      {% endif %}
      <launchSpecification>
        <imageId>{{ request.launch_specification.image_id }}</imageId>
        {% if request.launch_specification.key_name %}
          <keyName>{{ request.launch_specification.key_name }}</keyName>
        {% endif %}
        <groupSet>
          {% for group in request.launch_specification.groups %}
          <item>
            <groupId>{{ group.id }}</groupId>
            <groupName>{{ group.name }}</groupName>
          </item>
          {% endfor %}
        </groupSet>
        {% if request.launch_specification.kernel %}
        <kernelId>{{ request.launch_specification.kernel }}</kernelId>
        {% endif %}
        {% if request.launch_specification.ramdisk %}
        <ramdiskId>{{ request.launch_specification.ramdisk }}</ramdiskId>
        {% endif %}
        {% if request.launch_specification.subnet_id %}
        <subnetId>{{ request.launch_specification.subnet_id }}</subnetId>
        {% endif %}
        <instanceType>{{ request.launch_specification.instance_type }}</instanceType>
        <blockDeviceMapping/>
        <monitoring>
          <enabled>{{ request.launch_specification.monitored }}</enabled>
        </monitoring>
        <ebsOptimized>{{ request.launch_specification.ebs_optimized }}</ebsOptimized>
        {% if request.launch_specification.placement %}
          <PlacementRequestType>
            <availabilityZone>{{ request.launch_specification.placement }}</availabilityZone>
            <groupName></groupName>
          </PlacementRequestType>
        {% endif %}
      </launchSpecification>
      <tagSet>
        {% for tag in request.get_tags() %}
          <item>
            <resourceId>{{ tag.resource_id }}</resourceId>
            <resourceType>{{ tag.resource_type }}</resourceType>
            <key>{{ tag.key }}</key>
            <value>{{ tag.value }}</value>
          </item>
        {% endfor %}
      </tagSet>
      {% if request.launch_group %}
        <launchGroup>{{ request.launch_group }}</launchGroup>
      {% endif %}
        <createTime>2015-01-01T00:00:00.000Z</createTime>
      {% if request.valid_from %}
        <validFrom>{{ request.valid_from }}</validFrom>
      {% endif %}
      {% if request.valid_until %}
        <validUntil>{{ request.valid_until }}</validUntil>
      {% endif %}
      <productDescription>Linux/UNIX</productDescription>
      <instanceInterruptionBehavior>{{ request.instance_interruption_behaviour }}</instanceInterruptionBehavior>
    </item>
    {% endfor %}
  </spotInstanceRequestSet>
</DescribeSpotInstanceRequestsResponse>"""

CANCEL_SPOT_INSTANCES_TEMPLATE = """<CancelSpotInstanceRequestsResponse xmlns="http://ec2.amazonaws.com/doc/2013-10-15/">
  <requestId>59dbff89-35bd-4eac-99ed-be587EXAMPLE</requestId>
  <spotInstanceRequestSet>
    {% for request in requests %}
    <item>
      <spotInstanceRequestId>{{ request.id }}</spotInstanceRequestId>
      <state>cancelled</state>
    </item>
    {% endfor %}
  </spotInstanceRequestSet>
</CancelSpotInstanceRequestsResponse>"""

DESCRIBE_SPOT_PRICE_HISTORY_TEMPLATE = """<DescribeSpotPriceHistoryResponse xmlns="http://ec2.amazonaws.com/doc/2013-10-15/">
  <requestId>59dbff89-35bd-4eac-99ed-be587EXAMPLE</requestId>
  <spotPriceHistorySet>
    {% for price in prices %}
    <item>
      <instanceType>{{ price.InstanceType }}</instanceType>
      <productDescription>Linux/UNIX (Amazon VPC)</productDescription>
      <spotPrice>0.00001</spotPrice>
      <availabilityZone>{{ price.Location }}</availabilityZone>
      <timestamp>2006-01-02T15:04:05.999999999Z</timestamp>
    </item>
    {% endfor %}
  </spotPriceHistorySet>
  </DescribeSpotPriceHistoryResponse>"""
