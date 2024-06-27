from ._base_response import EC2BaseResponse


class Fleets(EC2BaseResponse):
    def delete_fleets(self) -> str:
        fleet_ids = self._get_multi_param("FleetId.")
        terminate_instances = self._get_bool_param("TerminateInstances")
        fleets = self.ec2_backend.delete_fleets(fleet_ids, terminate_instances)
        template = self.response_template(DELETE_FLEETS_TEMPLATE)
        return template.render(fleets=fleets)

    def describe_fleet_instances(self) -> str:
        fleet_id = self._get_param("FleetId")

        instances = self.ec2_backend.describe_fleet_instances(fleet_id)
        template = self.response_template(DESCRIBE_FLEET_INSTANCES_TEMPLATE)
        return template.render(fleet_id=fleet_id, instances=instances)

    def describe_fleets(self) -> str:
        fleet_ids = self._get_multi_param("FleetId.")

        requests = self.ec2_backend.describe_fleets(fleet_ids)
        template = self.response_template(DESCRIBE_FLEETS_TEMPLATE)
        rend = template.render(requests=requests)
        return rend

    def create_fleet(self) -> str:
        on_demand_options = self._get_multi_param_dict("OnDemandOptions")
        spot_options = self._get_multi_param_dict("SpotOptions")
        target_capacity_specification = self._get_multi_param_dict(
            "TargetCapacitySpecification"
        )
        launch_template_configs = self._get_multi_param(
            param_prefix="LaunchTemplateConfigs"
        )

        excess_capacity_termination_policy = self._get_param(
            "ExcessCapacityTerminationPolicy"
        )
        replace_unhealthy_instances = self._get_param("ReplaceUnhealthyInstances")
        terminate_instances_with_expiration = self._get_param(
            "TerminateInstancesWithExpiration", if_none=True
        )
        fleet_type = self._get_param("Type", if_none="maintain")
        valid_from = self._get_param("ValidFrom")
        valid_until = self._get_param("ValidUntil")

        tag_specifications = self._get_multi_param("TagSpecification")

        request = self.ec2_backend.create_fleet(
            on_demand_options=on_demand_options,
            spot_options=spot_options,
            target_capacity_specification=target_capacity_specification,
            launch_template_configs=launch_template_configs,
            excess_capacity_termination_policy=excess_capacity_termination_policy,
            replace_unhealthy_instances=replace_unhealthy_instances,
            terminate_instances_with_expiration=terminate_instances_with_expiration,
            fleet_type=fleet_type,
            valid_from=valid_from,
            valid_until=valid_until,
            tag_specifications=tag_specifications,
        )

        template = self.response_template(CREATE_FLEET_TEMPLATE)
        return template.render(request=request)


CREATE_FLEET_TEMPLATE = """<CreateFleetResponse xmlns="http://ec2.amazonaws.com/doc/2016-11-15/">
    <requestId>60262cc5-2bd4-4c8d-98ed-example</requestId>
    <fleetId>{{ request.id }}</fleetId>
    {% if request.fleet_type == "instant" %}
    <fleetInstanceSet>
        {% for instance in request.on_demand_instances %}
        <item>
            <instanceType>{{ instance["instance"].instance_type }}</instanceType>
            <lifecycle>on-demand</lifecycle>
            <instanceIds>
                <item>{{ instance["instance"].id }}</item>
            </instanceIds>
        </item>
        {% endfor %}
        {% for instance in request.spot_requests %}
        <item>
            <instanceType>{{ instance.instance.instance_type }}</instanceType>
            <lifecycle>spot</lifecycle>
            <instanceIds>
                <item>{{ instance.instance.id }}</item>
            </instanceIds>
        </item>
        {% endfor %}
    </fleetInstanceSet>
    {% endif %}
</CreateFleetResponse>"""

DESCRIBE_FLEETS_TEMPLATE = """<DescribeFleetsResponse xmlns="http://ec2.amazonaws.com/doc/2016-09-15/">
    <requestId>4d68a6cc-8f2e-4be1-b425-example</requestId>
    <fleetSet>
        {% for request in requests %}
        <item>
            <fleetId>{{ request.id }}</fleetId>
            <fleetState>{{ request.state }}</fleetState>
            <excessCapacityTerminationPolicy>{{ request.excess_capacity_termination_policy }}</excessCapacityTerminationPolicy>
            <fulfilledCapacity>{{ request.fulfilled_capacity }}</fulfilledCapacity>
            <fulfilledOnDemandCapacity>{{ request.fulfilled_on_demand_capacity }}</fulfilledOnDemandCapacity>
            <launchTemplateConfigs>
                {% for config in request.launch_template_configs %}
                <item>
                    <launchTemplateSpecification>
                        <launchTemplateId>{{ config.LaunchTemplateSpecification.LaunchTemplateId }}</launchTemplateId>
                        <version>{{ config.LaunchTemplateSpecification.Version }}</version>
                    </launchTemplateSpecification>
                    {% if config.Overrides %}
                    <overrides>
                        {% for override in config.Overrides %}
                        <item>
                            {% if override.AvailabilityZone %}
                            <availabilityZone>{{ override.AvailabilityZone }}</availabilityZone>
                            {% endif %}
                            {% if override.InstanceType %}
                            <instanceType>{{ override.InstanceType }}</instanceType>
                            {% endif %}
                            {% if override.InstanceRequirements %}
                            <instanceRequirements>
                                {% if override.InstanceRequirements.AcceleratorCount %}
                                <acceleratorCount>
                                    {% if override.InstanceRequirements.AcceleratorCount.Max %}
                                    <max>{{ override.InstanceRequirements.AcceleratorCount.Max }}</max>
                                    {% endif %}
                                    {% if override.InstanceRequirements.AcceleratorCount.Min %}
                                    <min>{{ override.InstanceRequirements.AcceleratorCount.Min }}</min>
                                    {% endif %}
                                </acceleratorCount>
                                {% endif %}
                                {% if override.InstanceRequirements.AcceleratorManufacturer %}
                                <acceleratorManufacturerSet>
                                    {% for manufacturer in override.InstanceRequirements.AcceleratorManufacturer %}
                                    <item>{{ manufacturer }}</item>
                                    {% endfor %}
                                </acceleratorManufacturerSet>
                                {% endif %}
                                {% if override.InstanceRequirements.AcceleratorName %}
                                <acceleratorNameSet>
                                    {% for name in override.InstanceRequirements.AcceleratorName %}
                                    <item>{{ name }}</item>
                                    {% endfor %}
                                </acceleratorNameSet>
                                {% endif %}
                                {% if override.InstanceRequirements.AcceleratorTotalMemoryMiB %}
                                <acceleratorTotalMemoryMiB>
                                    {% if override.InstanceRequirements.AcceleratorTotalMemoryMiB.Max %}
                                    <max>{{ override.InstanceRequirements.AcceleratorTotalMemoryMiB.Max }}</max>
                                    {% endif %}
                                    {% if override.InstanceRequirements.AcceleratorTotalMemoryMiB.Min %}
                                    <min>{{ override.InstanceRequirements.AcceleratorTotalMemoryMiB.Min }}</min>
                                    {% endif %}
                                </acceleratorTotalMemoryMiB>
                                {% endif %}
                                {% if override.InstanceRequirements.AcceleratorType %}
                                <acceleratorTypeSet>
                                    {% for type in override.InstanceRequirements.AcceleratorType %}
                                    <item>{{ type }}</item>
                                    {% endfor %}
                                </acceleratorTypeSet>
                                {% endif %}
                                {% if override.InstanceRequirements.BareMetal %}
                                <bareMetal>{{ override.InstanceRequirements.BareMetal }}</bareMetal>
                                {% endif %}
                                {% if override.InstanceRequirements.BaselineEbsBandwidthMbps %}
                                <baselineEbsBandwidthMbps>
                                    {% if override.InstanceRequirements.BaselineEbsBandwidthMbps.Min %}
                                    <min>{{ override.InstanceRequirements.BaselineEbsBandwidthMbps.Min }}</min>
                                    {% endif %}
                                    {% if override.InstanceRequirements.BaselineEbsBandwidthMbps.Max %}
                                    <max>{{ override.InstanceRequirements.BaselineEbsBandwidthMbps.Max }}</max>
                                    {% endif %}
                                </baselineEbsBandwidthMbps>
                                {% endif %}
                                {% if override.InstanceRequirements.BurstablePerformance %}
                                <burstablePerformance>{{ override.InstanceRequirements.BurstablePerformance }}</burstablePerformance>
                                {% endif %}
                                {% if override.InstanceRequirements.CpuManufacturer %}
                                <cpuManufacturerSet>
                                    {% for manufacturer in override.InstanceRequirements.CpuManufacturer %}
                                    <item>{{ manufacturer }}</item>
                                    {% endfor %}
                                </cpuManufacturerSet>
                                {% endif %}
                                {% if override.InstanceRequirements.ExcludedInstanceType %}
                                <excludedInstanceTypeSet>
                                    {% for type in override.InstanceRequirements.ExcludedInstanceType %}
                                    <item>{{ type }}</item>
                                    {% endfor %}
                                </excludedInstanceTypeSet>
                                {% endif %}
                                {% if override.InstanceRequirements.InstanceGeneration %}
                                <instanceGenerationSet>
                                    {% for generation in override.InstanceRequirements.InstanceGeneration %}
                                    <item>{{ generation }}</item>
                                    {% endfor %}
                                </instanceGenerationSet>
                                {% endif %}
                                {% if override.InstanceRequirements.LocalStorage %}
                                <localStorage>{{ override.InstanceRequirements.LocalStorage }}</localStorage>
                                {% endif %}
                                {% if override.InstanceRequirements.LocalStorageType %}
                                <localStorageTypeSet>
                                    {% for type in override.InstanceRequirements.LocalStorageType %}
                                    <item>{{ type }}</item>
                                    {% endfor %}
                                </localStorageTypeSet>
                                {% endif %}
                                {% if override.InstanceRequirements.MemoryGiBPerVCpu %}
                                <memoryGiBPerVCpu>
                                    {% if override.InstanceRequirements.MemoryGiBPerVCpu.Min %}
                                    <min>{{ override.InstanceRequirements.MemoryGiBPerVCpu.Min }}</min>
                                    {% endif %}
                                    {% if override.InstanceRequirements.MemoryGiBPerVCpu.Max %}
                                    <max>{{ override.InstanceRequirements.MemoryGiBPerVCpu.Max }}</max>
                                    {% endif %}
                                </memoryGiBPerVCpu>
                                {% endif %}
                                {% if override.InstanceRequirements.MemoryMiB %}
                                <memoryMiB>
                                    {% if override.InstanceRequirements.MemoryMiB.Min %}
                                    <min>{{ override.InstanceRequirements.MemoryMiB.Min }}</min>
                                    {% endif %}
                                    {% if override.InstanceRequirements.MemoryMiB.Max %}
                                    <max>{{ override.InstanceRequirements.MemoryMiB.Max }}</max>
                                    {% endif %}
                                </memoryMiB>
                                {% endif %}
                                {% if override.InstanceRequirements.NetworkInterfaceCount %}
                                <networkInterfaceCount>
                                    {% if override.InstanceRequirements.NetworkInterfaceCount.Max %}
                                    <max>{{ override.InstanceRequirements.NetworkInterfaceCount.Max }}</max>
                                    {% endif %}
                                    {% if override.InstanceRequirements.NetworkInterfaceCount.Min %}
                                    <min>{{ override.InstanceRequirements.NetworkInterfaceCount.Min }}</min>
                                    {% endif %}
                                </networkInterfaceCount>
                                {% endif %}
                                {% if override.InstanceRequirements.OnDemandMaxPricePercentageOverLowestPrice %}
                                <onDemandMaxPricePercentageOverLowestPrice>{{ override.InstanceRequirements.OnDemandMaxPricePercentageOverLowestPrice }}</onDemandMaxPricePercentageOverLowestPrice>
                                {% endif %}
                                {% if override.InstanceRequirements.RequireHibernateSupport %}
                                <requireHibernateSupport>{{ override.InstanceRequirements.RequireHibernateSupport }}</requireHibernateSupport>
                                {% endif %}
                                {% if override.InstanceRequirements.SpotMaxPricePercentageOverLowestPrice %}
                                <spotMaxPricePercentageOverLowestPrice>{{ override.InstanceRequirements.SpotMaxPricePercentageOverLowestPrice }}</spotMaxPricePercentageOverLowestPrice>
                                {% endif %}
                                {% if override.InstanceRequirements.TotalLocalStorageGB %}
                                <totalLocalStorageGB>
                                    {% if override.InstanceRequirements.TotalLocalStorageGB.Min %}
                                    <min>{{ override.InstanceRequirements.TotalLocalStorageGB.Min }}</min>
                                    {% endif %}
                                    {% if override.InstanceRequirements.TotalLocalStorageGB.Max %}
                                    <max>{{ override.InstanceRequirements.TotalLocalStorageGB.Max }}</max>
                                    {% endif %}
                                </totalLocalStorageGB>
                                {% endif %}
                                {% if override.InstanceRequirements.VCpuCount %}
                                <vCpuCount>
                                    {% if override.InstanceRequirements.VCpuCount.Min %}
                                    <min>{{ override.InstanceRequirements.VCpuCount.Min }}</min>
                                    {% endif %}
                                    {% if override.InstanceRequirements.VCpuCount.Max %}
                                    <max>{{ override.InstanceRequirements.VCpuCount.Max }}</max>
                                    {% endif %}
                                </vCpuCount>
                                {% endif %}
                            </instanceRequirements>
                            {% endif %}
                            {% if override.MaxPrice %}
                            <maxPrice>{{ override.MaxPrice }}</maxPrice>
                            {% endif %}
                            {% if override.Placement %}
                            <placement>
                                {% if override.Placement.GroupName %}
                                <groupName>{{ override.Placement.GroupName }}</groupName>
                                {% endif %}
                            </placement>
                            {% endif %}
                            {% if override.Priority %}
                            <priority>{{ override.Priority }}</priority>
                            {% endif %}
                            {% if override.SubnetId %}
                            <subnetId>{{ override.SubnetId }}</subnetId>
                            {% endif %}
                            {% if override.WeightedCapacity %}
                            <weightedCapacity>{{ override.WeightedCapacity }}</weightedCapacity>
                            {% endif %}
                        </item>
                        {% endfor %}
                    </overrides>
                    {% endif %}
                </item>
                {% endfor %}
            </launchTemplateConfigs>
            <targetCapacitySpecification>
                <totalTargetCapacity>{{ request.target_capacity }}</totalTargetCapacity>
                {% if request.on_demand_target_capacity %}
                <onDemandTargetCapacity>{{ request.on_demand_target_capacity }}</onDemandTargetCapacity>
                {% endif %}
                {% if request.spot_target_capacity %}
                <spotTargetCapacity>{{ request.spot_target_capacity }}</spotTargetCapacity>
                {% endif %}
                <defaultTargetCapacityType>{{ request.target_capacity_specification.DefaultTargetCapacityType }}</defaultTargetCapacityType>
            </targetCapacitySpecification>
            {% if request.spot_options %}
            <spotOptions>
                {% if request.spot_options.AllocationStrategy %}
                <allocationStrategy>{{ request.spot_options.AllocationStrategy }}</allocationStrategy>
                {% endif %}
                {% if request.spot_options.InstanceInterruptionBehavior %}
                <instanceInterruptionBehavior>{{ request.spot_options.InstanceInterruptionBehavior }}</instanceInterruptionBehavior>
                {% endif %}
                {% if request.spot_options.InstancePoolsToUseCount %}
                <instancePoolsToUseCount>{{ request.spot_options.InstancePoolsToUseCount }}</instancePoolsToUseCount>
                {% endif %}
                {% if request.spot_options.MaintenanceStrategies %}
                <maintenanceStrategies>
                    {% if request.spot_options.MaintenanceStrategies.CapacityRebalance %}
                    <capacityRebalance>
                        {% if request.spot_options.MaintenanceStrategies.CapacityRebalance.ReplacementStrategy %}
                        <replacementStrategy>{{ request.spot_options.MaintenanceStrategies.CapacityRebalance.ReplacementStrategy }}</replacementStrategy>
                        {% endif %}
                        {% if request.spot_options.MaintenanceStrategies.CapacityRebalance.TerminationDelay %}
                        <terminationDelay>{{ request.spot_options.MaintenanceStrategies.CapacityRebalance.TerminationDelay }}</terminationDelay>
                        {% endif %}
                    </capacityRebalance>
                    {% endif %}
                </maintenanceStrategies>
                {% endif %}
                {% if request.spot_options.MaxTotalPrice %}
                <maxTotalPrice>{{ request.spot_options.MaxTotalPrice }}</maxTotalPrice>
                {% endif %}
                {% if request.spot_options.MinTargetCapacity %}
                <minTargetCapacity>{{ request.spot_options.MinTargetCapacity }}</minTargetCapacity>
                {% endif %}
                {% if request.spot_options.SingleAvailabilityZone %}
                <singleAvailabilityZone>{{ request.spot_options.SingleAvailabilityZone }}</singleAvailabilityZone>
                {% endif %}
                {% if request.spot_options.SingleInstanceType %}
                <singleInstanceType>{{ request.spot_options.SingleInstanceType }}</singleInstanceType>
                {% endif %}
            </spotOptions>
            {% endif %}
            <!-- {'AllocationStrategy': 'lowest-price', 'MaxTotalPrice': '50', 'MinTargetCapacity': 1, 'SingleAvailabilityZone': True, 'SingleInstanceType': True} -->
            {% if request.on_demand_options %}
            <onDemandOptions>
                {% if request.on_demand_options.AllocationStrategy %}
                <allocationStrategy>{{ request.on_demand_options.AllocationStrategy }}</allocationStrategy>
                {% endif %}
                {% if request.on_demand_options.MaxTotalPrice %}
                <maxTotalPrice>{{ request.on_demand_options.MaxTotalPrice }}</maxTotalPrice>
                {% endif %}
                {% if request.on_demand_options.MinTargetCapacity %}
                <minTargetCapacity>{{ request.on_demand_options.MinTargetCapacity }}</minTargetCapacity>
                {% endif %}
                {% if request.on_demand_options.SingleAvailabilityZone %}
                <singleAvailabilityZone>{{ request.on_demand_options.SingleAvailabilityZone }}</singleAvailabilityZone>
                {% endif %}
                {% if request.on_demand_options.SingleInstanceType %}
                <singleInstanceType>{{ request.on_demand_options.SingleInstanceType }}</singleInstanceType>
                {% endif %}
                {% if request.on_demand_options.CapacityReservationOptions %}
                <capacityReservationOptions>
                    {% if request.on_demand_options.CapacityReservationOptions.UsageStrategy %}
                    <usageStrategy>{{ request.on_demand_options.CapacityReservationOptions.UsageStrategy }}</usageStrategy>
                    {% endif %}
                </capacityReservationOptions>
                {% endif %}
            </onDemandOptions>
            {% endif %}
            <terminateInstancesWithExpiration>{{ request.terminate_instances_with_expiration }}</terminateInstancesWithExpiration>
            <type>{{ request.fleet_type }}</type>
            {% if request.valid_from %}
            <validFrom>{{ request.valid_from }}</validFrom>
            {% endif %}
            {% if request.valid_until %}
            <validUntil>{{ request.valid_until }}</validUntil>
            {% endif %}
            <replaceUnhealthyInstances>{{ request.replace_unhealthy_instances }}</replaceUnhealthyInstances>
            <tagSet>
                {% for tag in request.tags %}
                <item>
                    <key>{{ tag.key }}</key>
                    <value>{{ tag.value }}</value>
                </item>
                {% endfor %}
            </tagSet>
        </item>
        {% endfor %}
    </fleetSet>
</DescribeFleetsResponse>"""

DESCRIBE_FLEET_INSTANCES_TEMPLATE = """<DescribeFleetInstancesResponse xmlns="http://ec2.amazonaws.com/doc/2016-09-15/">
    <requestId>cfb09950-45e2-472d-a6a9-example</requestId>
    <fleetId>{{ fleet_id }}</fleetId>
    <activeInstanceSet>
        {% for i in instances %}
        <item>
            <instanceId>{{ i.instance.id }}</instanceId>
            {% if i.id %}
            <spotInstanceRequestId>{{ i.id }}</spotInstanceRequestId>
            {% endif %}
            <instanceType>{{ i.instance.instance_type }}</instanceType>
            <instanceHealth>healthy</instanceHealth>
        </item>
        {% endfor %}
    </activeInstanceSet>
</DescribeFleetInstancesResponse>
"""

DELETE_FLEETS_TEMPLATE = """<DeleteFleetResponse xmlns="http://ec2.amazonaws.com/doc/2016-09-15/">
    <requestId>e12d2fe5-6503-4b4b-911c-example</requestId>
    <unsuccessfulFleetDeletionSet/>
    <successfulFleetDeletionSet>
        {% for fleet in fleets %}
        <item>
            <fleetId>{{ fleet.id }}</fleetId>
            <currentFleetState>{{ fleet.state }}</currentFleetState>
            <previousFleetState>active</previousFleetState>
        </item>
        {% endfor %}
    </successfulFleetDeletionSet>
</DeleteFleetResponse>"""
