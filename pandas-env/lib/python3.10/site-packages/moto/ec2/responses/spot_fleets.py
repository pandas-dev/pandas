from ._base_response import EC2BaseResponse


class SpotFleets(EC2BaseResponse):
    def cancel_spot_fleet_requests(self) -> str:
        spot_fleet_request_ids = self._get_multi_param("SpotFleetRequestId.")
        terminate_instances = self._get_bool_param("TerminateInstances")
        spot_fleets = self.ec2_backend.cancel_spot_fleet_requests(
            spot_fleet_request_ids, terminate_instances
        )
        template = self.response_template(CANCEL_SPOT_FLEETS_TEMPLATE)
        return template.render(spot_fleets=spot_fleets)

    def describe_spot_fleet_instances(self) -> str:
        spot_fleet_request_id = self._get_param("SpotFleetRequestId")

        spot_requests = self.ec2_backend.describe_spot_fleet_instances(
            spot_fleet_request_id
        )
        template = self.response_template(DESCRIBE_SPOT_FLEET_INSTANCES_TEMPLATE)
        return template.render(
            spot_request_id=spot_fleet_request_id, spot_requests=spot_requests
        )

    def describe_spot_fleet_requests(self) -> str:
        spot_fleet_request_ids = self._get_multi_param("SpotFleetRequestId.")

        requests = self.ec2_backend.describe_spot_fleet_requests(spot_fleet_request_ids)
        template = self.response_template(DESCRIBE_SPOT_FLEET_TEMPLATE)
        return template.render(requests=requests)

    def modify_spot_fleet_request(self) -> str:
        spot_fleet_request_id = self._get_param("SpotFleetRequestId")
        target_capacity = self._get_int_param("TargetCapacity")
        terminate_instances = self._get_param(
            "ExcessCapacityTerminationPolicy", if_none="Default"
        )
        self.ec2_backend.modify_spot_fleet_request(
            spot_fleet_request_id, target_capacity, terminate_instances
        )
        return self.response_template(MODIFY_SPOT_FLEET_REQUEST_TEMPLATE).render()

    def request_spot_fleet(self) -> str:
        spot_config = self._get_multi_param_dict("SpotFleetRequestConfig")
        spot_price = spot_config.get("SpotPrice")
        target_capacity = spot_config["TargetCapacity"]
        iam_fleet_role = spot_config["IamFleetRole"]
        allocation_strategy = spot_config["AllocationStrategy"]
        instance_interruption_behaviour = spot_config.get(
            "InstanceInterruptionBehavior"
        )

        launch_specs = spot_config.get("LaunchSpecifications")
        launch_template_config = list(
            self._get_params()
            .get("SpotFleetRequestConfig", {})
            .get("LaunchTemplateConfigs", {})
            .values()
        )
        tag_specifications = spot_config.get("TagSpecification")

        request = self.ec2_backend.request_spot_fleet(
            spot_price=spot_price,
            target_capacity=target_capacity,
            iam_fleet_role=iam_fleet_role,
            allocation_strategy=allocation_strategy,
            launch_specs=launch_specs,
            launch_template_config=launch_template_config,
            instance_interruption_behaviour=instance_interruption_behaviour,
            tag_specifications=tag_specifications,
        )

        template = self.response_template(REQUEST_SPOT_FLEET_TEMPLATE)
        return template.render(request=request)


REQUEST_SPOT_FLEET_TEMPLATE = """<RequestSpotFleetResponse xmlns="http://ec2.amazonaws.com/doc/2016-09-15/">
    <requestId>60262cc5-2bd4-4c8d-98ed-example</requestId>
    <spotFleetRequestId>{{ request.id }}</spotFleetRequestId>
</RequestSpotFleetResponse>"""

MODIFY_SPOT_FLEET_REQUEST_TEMPLATE = """<ModifySpotFleetRequestResponse xmlns="http://ec2.amazonaws.com/doc/2016-09-15/">
    <requestId>21681fea-9987-aef3-2121-example</requestId>
    <return>true</return>
</ModifySpotFleetRequestResponse>"""

DESCRIBE_SPOT_FLEET_TEMPLATE = """<DescribeSpotFleetRequestsResponse xmlns="http://ec2.amazonaws.com/doc/2016-09-15/">
    <requestId>4d68a6cc-8f2e-4be1-b425-example</requestId>
    <spotFleetRequestConfigSet>
        {% for request in requests %}
        <item>
            <spotFleetRequestId>{{ request.id }}</spotFleetRequestId>
            <spotFleetRequestState>{{ request.state }}</spotFleetRequestState>
            <tagSet>
                {% for key, value in request.tags.get('spot-fleet-request', {}).items() %}
                <item>
                    <key>{{ key }}</key>
                    <value>{{ value }}</value>
                </item>
                {% endfor %}
            </tagSet>
            <spotFleetRequestConfig>
                {% if request.spot_price %}
                <spotPrice>{{ request.spot_price }}</spotPrice>
                {% endif %}
                <targetCapacity>{{ request.target_capacity }}</targetCapacity>
                <iamFleetRole>{{ request.iam_fleet_role }}</iamFleetRole>
                <allocationStrategy>{{ request.allocation_strategy }}</allocationStrategy>
                <fulfilledCapacity>{{ request.fulfilled_capacity }}</fulfilledCapacity>
                <launchSpecifications>
                    {% for launch_spec in request.launch_specs %}
                    <item>
                        <subnetId>{{ launch_spec.subnet_id }}</subnetId>
                        <ebsOptimized>{{ launch_spec.ebs_optimized }}</ebsOptimized>
                        <imageId>{{ launch_spec.image_id }}</imageId>
                        <instanceType>{{ launch_spec.instance_type }}</instanceType>
                        <iamInstanceProfile><arn>{{ launch_spec.iam_instance_profile }}</arn></iamInstanceProfile>
                        <keyName>{{ launch_spec.key_name }}</keyName>
                        <monitoring><enabled>{{ launch_spec.monitoring }}</enabled></monitoring>
                        {% if launch_spec.spot_price %}
                        <spotPrice>{{ launch_spec.spot_price }}</spotPrice>
                        {% endif %}
                        <userData>{{ launch_spec.user_data }}</userData>
                        <weightedCapacity>{{ launch_spec.weighted_capacity }}</weightedCapacity>
                        <groupSet>
                            {% for group in launch_spec.group_set %}
                            <item>
                                <groupId>{{ group }}</groupId>
                            </item>
                            {% endfor %}
                        </groupSet>
                        <tagSpecificationSet>
                            {% for resource_type in launch_spec.tag_specifications %}
                            <item>
                                <resourceType>{{ resource_type }}</resourceType>
                                <tag>
                                {% for key, value in launch_spec.tag_specifications[resource_type].items() %}
                                    <item>
                                        <key>{{ key }}</key>
                                        <value>{{ value }}</value>
                                    </item>
                                {% endfor %}
                                </tag>
                            </item>
                            {% endfor %}
                        </tagSpecificationSet>
                    </item>
                    {% endfor %}
                </launchSpecifications>
            </spotFleetRequestConfig>
        </item>
        {% endfor %}
    </spotFleetRequestConfigSet>
</DescribeSpotFleetRequestsResponse>"""

DESCRIBE_SPOT_FLEET_INSTANCES_TEMPLATE = """<DescribeSpotFleetInstancesResponse xmlns="http://ec2.amazonaws.com/doc/2016-09-15/">
    <requestId>cfb09950-45e2-472d-a6a9-example</requestId>
    <spotFleetRequestId>{{ spot_request_id }}</spotFleetRequestId>
    <activeInstanceSet>
        {% for spot_request in spot_requests %}
        <item>
            <instanceId>{{ spot_request.instance.id }}</instanceId>
            <spotInstanceRequestId>{{ spot_request.id }}</spotInstanceRequestId>
            <instanceType>{{ spot_request.instance.instance_type }}</instanceType>
        </item>
        {% endfor %}
    </activeInstanceSet>
</DescribeSpotFleetInstancesResponse>
"""

CANCEL_SPOT_FLEETS_TEMPLATE = """<CancelSpotFleetRequestsResponse xmlns="http://ec2.amazonaws.com/doc/2016-09-15/">
    <requestId>e12d2fe5-6503-4b4b-911c-example</requestId>
    <unsuccessfulFleetRequestSet/>
    <successfulFleetRequestSet>
        {% for spot_fleet in spot_fleets %}
        <item>
            <spotFleetRequestId>{{ spot_fleet.id }}</spotFleetRequestId>
            <currentSpotFleetRequestState>cancelled_terminating</currentSpotFleetRequestState>
            <previousSpotFleetRequestState>active</previousSpotFleetRequestState>
        </item>
        {% endfor %}
    </successfulFleetRequestSet>
</CancelSpotFleetRequestsResponse>"""
