import json

from moto.core.responses import BaseResponse

from .models import OpsWorksBackend, opsworks_backends


class OpsWorksResponse(BaseResponse):
    def __init__(self) -> None:
        super().__init__(service_name="opsworks")

    @property
    def opsworks_backend(self) -> OpsWorksBackend:
        return opsworks_backends[self.current_account][self.region]

    def create_stack(self) -> str:
        kwargs = dict(
            name=self._get_param("Name"),
            region=self._get_param("Region"),
            vpcid=self._get_param("VpcId"),
            attributes=self._get_param("Attributes"),
            default_instance_profile_arn=self._get_param("DefaultInstanceProfileArn"),
            default_os=self._get_param("DefaultOs"),
            hostname_theme=self._get_param("HostnameTheme"),
            default_availability_zone=self._get_param("DefaultAvailabilityZone"),
            default_subnet_id=self._get_param("DefaultInstanceProfileArn"),
            custom_json=self._get_param("CustomJson"),
            configuration_manager=self._get_param("ConfigurationManager"),
            chef_configuration=self._get_param("ChefConfiguration"),
            use_custom_cookbooks=self._get_param("UseCustomCookbooks"),
            use_opsworks_security_groups=self._get_param("UseOpsworksSecurityGroups"),
            custom_cookbooks_source=self._get_param("CustomCookbooksSource"),
            default_ssh_keyname=self._get_param("DefaultSshKeyName"),
            default_root_device_type=self._get_param("DefaultRootDeviceType"),
            service_role_arn=self._get_param("ServiceRoleArn"),
            agent_version=self._get_param("AgentVersion"),
        )
        stack = self.opsworks_backend.create_stack(**kwargs)
        return json.dumps({"StackId": stack.id}, indent=1)

    def create_layer(self) -> str:
        kwargs = dict(
            stack_id=self._get_param("StackId"),
            layer_type=self._get_param("Type"),
            name=self._get_param("Name"),
            shortname=self._get_param("Shortname"),
            attributes=self._get_param("Attributes"),
            custom_instance_profile_arn=self._get_param("CustomInstanceProfileArn"),
            custom_json=self._get_param("CustomJson"),
            custom_security_group_ids=self._get_param("CustomSecurityGroupIds"),
            packages=self._get_param("Packages"),
            volume_configurations=self._get_param("VolumeConfigurations"),
            enable_autohealing=self._get_param("EnableAutoHealing"),
            auto_assign_elastic_ips=self._get_param("AutoAssignElasticIps"),
            auto_assign_public_ips=self._get_param("AutoAssignPublicIps"),
            custom_recipes=self._get_param("CustomRecipes"),
            install_updates_on_boot=self._get_param("InstallUpdatesOnBoot"),
            use_ebs_optimized_instances=self._get_param("UseEbsOptimizedInstances"),
            lifecycle_event_configuration=self._get_param(
                "LifecycleEventConfiguration"
            ),
        )
        layer = self.opsworks_backend.create_layer(**kwargs)
        return json.dumps({"LayerId": layer.id}, indent=1)

    def create_app(self) -> str:
        kwargs = dict(
            stack_id=self._get_param("StackId"),
            name=self._get_param("Name"),
            app_type=self._get_param("Type"),
            shortname=self._get_param("Shortname"),
            description=self._get_param("Description"),
            datasources=self._get_param("DataSources"),
            app_source=self._get_param("AppSource"),
            domains=self._get_param("Domains"),
            enable_ssl=self._get_param("EnableSsl"),
            ssl_configuration=self._get_param("SslConfiguration"),
            attributes=self._get_param("Attributes"),
            environment=self._get_param("Environment"),
        )
        app = self.opsworks_backend.create_app(**kwargs)
        return json.dumps({"AppId": app.id}, indent=1)

    def create_instance(self) -> str:
        kwargs = dict(
            stack_id=self._get_param("StackId"),
            layer_ids=self._get_param("LayerIds"),
            instance_type=self._get_param("InstanceType"),
            auto_scale_type=self._get_param("AutoScalingType"),
            hostname=self._get_param("Hostname"),
            os=self._get_param("Os"),
            ami_id=self._get_param("AmiId"),
            ssh_keyname=self._get_param("SshKeyName"),
            availability_zone=self._get_param("AvailabilityZone"),
            virtualization_type=self._get_param("VirtualizationType"),
            subnet_id=self._get_param("SubnetId"),
            architecture=self._get_param("Architecture"),
            root_device_type=self._get_param("RootDeviceType"),
            block_device_mappings=self._get_param("BlockDeviceMappings"),
            install_updates_on_boot=self._get_param("InstallUpdatesOnBoot"),
            ebs_optimized=self._get_param("EbsOptimized"),
            agent_version=self._get_param("AgentVersion"),
        )
        opsworks_instance = self.opsworks_backend.create_instance(**kwargs)
        return json.dumps({"InstanceId": opsworks_instance.id}, indent=1)

    def describe_stacks(self) -> str:
        stack_ids = self._get_param("StackIds")
        stacks = self.opsworks_backend.describe_stacks(stack_ids)
        return json.dumps({"Stacks": stacks}, indent=1)

    def describe_layers(self) -> str:
        stack_id = self._get_param("StackId")
        layer_ids = self._get_param("LayerIds")
        layers = self.opsworks_backend.describe_layers(stack_id, layer_ids)
        return json.dumps({"Layers": layers}, indent=1)

    def describe_apps(self) -> str:
        stack_id = self._get_param("StackId")
        app_ids = self._get_param("AppIds")
        apps = self.opsworks_backend.describe_apps(stack_id, app_ids)
        return json.dumps({"Apps": apps}, indent=1)

    def describe_instances(self) -> str:
        instance_ids = self._get_param("InstanceIds")
        layer_id = self._get_param("LayerId")
        stack_id = self._get_param("StackId")
        instances = self.opsworks_backend.describe_instances(
            instance_ids, layer_id, stack_id
        )
        return json.dumps({"Instances": instances}, indent=1)

    def start_instance(self) -> str:
        instance_id = self._get_param("InstanceId")
        self.opsworks_backend.start_instance(instance_id)
        return ""
