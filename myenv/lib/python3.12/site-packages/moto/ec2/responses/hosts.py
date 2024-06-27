from ._base_response import EC2BaseResponse


class HostsResponse(EC2BaseResponse):
    def allocate_hosts(self) -> str:
        params = self._get_params()
        quantity = int(params.get("Quantity"))  # type: ignore
        host_recovery = params.get("HostRecovery")
        zone = params.get("AvailabilityZone")
        instance_type = params.get("InstanceType")
        instance_family = params.get("InstanceFamily")
        auto_placement = params.get("AutoPlacement")
        tags = self._parse_tag_specification()
        host_tags = tags.get("dedicated-host", {})
        host_ids = self.ec2_backend.allocate_hosts(
            quantity,
            host_recovery,
            zone,
            instance_type,
            instance_family,
            auto_placement,
            host_tags,
        )
        template = self.response_template(EC2_ALLOCATE_HOSTS)
        return template.render(host_ids=host_ids)

    def describe_hosts(self) -> str:
        host_ids = list(self._get_params().get("HostId", {}).values())
        filters = self._filters_from_querystring()
        hosts = self.ec2_backend.describe_hosts(host_ids, filters)
        template = self.response_template(EC2_DESCRIBE_HOSTS)
        return template.render(account_id=self.current_account, hosts=hosts)

    def modify_hosts(self) -> str:
        params = self._get_params()
        host_ids = list(self._get_params().get("HostId", {}).values())
        auto_placement = params.get("AutoPlacement")
        host_recovery = params.get("HostRecovery")
        instance_type = params.get("InstanceType")
        instance_family = params.get("InstanceFamily")
        self.ec2_backend.modify_hosts(
            host_ids, auto_placement, host_recovery, instance_type, instance_family
        )
        template = self.response_template(EC2_MODIFY_HOSTS)
        return template.render(host_ids=host_ids)

    def release_hosts(self) -> str:
        host_ids = list(self._get_params().get("HostId", {}).values())
        self.ec2_backend.release_hosts(host_ids)
        template = self.response_template(EC2_RELEASE_HOSTS)
        return template.render(host_ids=host_ids)


EC2_ALLOCATE_HOSTS = """<AllocateHostsResult xmlns="http://ec2.amazonaws.com/doc/2013-10-15/">
<requestId>fdcdcab1-ae5c-489e-9c33-4637c5dda355</requestId>
    <hostIdSet>
        {% for host_id in host_ids %}
            <item>{{ host_id }}</item>
        {% endfor %}
    </hostIdSet>
</AllocateHostsResult>"""


EC2_DESCRIBE_HOSTS = """<DescribeHostsResult xmlns="http://ec2.amazonaws.com/doc/2013-10-15/">
<requestId>fdcdcab1-ae5c-489e-9c33-4637c5dda355</requestId>
    <hostSet>
        {% for host in hosts %}
            <item>
                <allocationTime>{{ host.allocation_time }}</allocationTime>
                <autoPlacement>{{ host.auto_placement }}</autoPlacement>
                <availabilityZone>{{ host.zone }}</availabilityZone>
                <availableCapacity></availableCapacity>
                <hostId>{{ host.id }}</hostId>
                <state>{{ host.state }}</state>
                <hostProperties>
                  {% if host.instance_type %}
                    <instanceType>{{ host.instance_type }}</instanceType>
                  {% endif %}
                  {% if host.instance_family %}
                    <instanceFamily>{{ host.instance_family }}</instanceFamily>
                  {% endif %}
                </hostProperties>
                <hostReservationId>reserv_id</hostReservationId>
                <instances>
                </instances>
                <ownerId>{{ account_id }}</ownerId>
                <hostRecovery>{{ host.host_recovery }}</hostRecovery>
                <tagSet>
                  {% for tag in host.get_tags() %}
                    <item>
                      <key>{{ tag.key }}</key>
                      <value>{{ tag.value }}</value>
                    </item>
                  {% endfor %}
                </tagSet>
            </item>
        {% endfor %}
    </hostSet>
</DescribeHostsResult>"""


EC2_MODIFY_HOSTS = """<ModifyHostsResult xmlns="http://ec2.amazonaws.com/doc/2013-10-15/">
<requestId>fdcdcab1-ae5c-489e-9c33-4637c5dda355</requestId>
    <successful>
        {% for host_id in host_ids %}
            <item>{{ host_id }}</item>
        {% endfor %}
    </successful>
</ModifyHostsResult>"""


EC2_RELEASE_HOSTS = """<ReleaseHostsResult xmlns="http://ec2.amazonaws.com/doc/2013-10-15/">
<requestId>fdcdcab1-ae5c-489e-9c33-4637c5dda355</requestId>
    <successful>
        {% for host_id in host_ids %}
            <item>{{ host_id }}</item>
        {% endfor %}
    </successful>
</ReleaseHostsResult>"""
