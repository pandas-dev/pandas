import re
from collections import OrderedDict
from typing import Any, Dict, Iterable, List, Optional

from botocore.exceptions import ParamValidationError

from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.common_models import BaseModel, CloudFormationModel
from moto.core.exceptions import RESTError
from moto.core.utils import iso_8601_datetime_with_milliseconds
from moto.ec2.models import EC2Backend, ec2_backends
from moto.ec2.models.subnets import Subnet
from moto.moto_api._internal import mock_random
from moto.utilities.tagging_service import TaggingService

from ..elb.models import register_certificate
from ..utilities.utils import ARN_PARTITION_REGEX
from .exceptions import (
    ActionTargetGroupNotFoundError,
    DuplicateListenerError,
    DuplicateLoadBalancerName,
    DuplicatePriorityError,
    DuplicateTargetGroupName,
    InvalidActionTypeError,
    InvalidConditionFieldError,
    InvalidConditionValueError,
    InvalidConfigurationRequest,
    InvalidDescribeRulesRequest,
    InvalidLoadBalancerActionException,
    InvalidModifyRuleArgumentsError,
    InvalidProtocol,
    InvalidProtocolValue,
    InvalidStatusCodeActionTypeError,
    InvalidTargetError,
    InvalidTargetGroupNameError,
    ListenerNotFoundError,
    LoadBalancerNotFoundError,
    PriorityInUseError,
    ResourceInUseError,
    RuleNotFoundError,
    SubnetNotFoundError,
    TargetGroupNotFoundError,
    TargetNotRunning,
    TooManyCertificatesError,
    TooManyTagsError,
    ValidationError,
)
from .utils import make_arn_for_load_balancer, make_arn_for_target_group

ALLOWED_ACTIONS = [
    "redirect",
    "authenticate-cognito",
    "authenticate-oidc",
    "fixed-response",
    "forward",
]


class FakeHealthStatus(BaseModel):
    def __init__(
        self,
        instance_id: str,
        port: str,
        health_port: Optional[str],
        status: str,
        reason: Optional[str] = None,
        description: Optional[str] = None,
    ):
        self.instance_id = instance_id
        self.port = port
        self.health_check_port = health_port
        self.status = status
        self.reason = reason
        self.description = description

    @property
    def target(self) -> dict[str, str]:
        return {
            "Id": self.instance_id,
            "Port": self.port,
        }

    @property
    def target_health(self) -> dict[str, Optional[str]]:
        return {
            "State": self.status,
            "Reason": self.reason,
            "Description": self.description,
        }


class FakeTargetGroup(CloudFormationModel):
    HTTP_CODE_REGEX = re.compile(r"(?:(?:\d+-\d+|\d+),?)+")

    def __init__(
        self,
        name: str,
        account_id: str,
        region_name: str,
        vpc_id: str,
        protocol: str,
        port: str,
        protocol_version: Optional[str] = None,
        healthcheck_protocol: Optional[str] = None,
        healthcheck_port: Optional[str] = None,
        healthcheck_path: Optional[str] = None,
        healthcheck_interval_seconds: Optional[int] = None,
        healthcheck_timeout_seconds: Optional[int] = None,
        healthcheck_enabled: Optional[str] = None,
        healthy_threshold_count: Optional[str] = None,
        unhealthy_threshold_count: Optional[str] = None,
        matcher: Optional[Dict[str, Any]] = None,
        target_type: Optional[str] = None,
        ip_address_type: Optional[str] = None,
    ):
        # TODO: default values differs when you add Network Load balancer
        self.name = name
        self.account_id = account_id
        self.region_name = region_name
        tg_id = mock_random.get_random_hex(length=16)
        self.arn = make_arn_for_target_group(
            account_id=self.account_id,
            tg_id=tg_id,
            name=name,
            region_name=self.region_name,
        )
        self.vpc_id = vpc_id
        if target_type == "lambda":
            self.protocol = None
            self.protocol_version = None
        elif target_type == "alb":
            self.protocol = "TCP"
            self.protocol_version = None
        else:
            self.protocol = protocol
            self.protocol_version = protocol_version
        self.port = port
        self.health_check_protocol = healthcheck_protocol or self.protocol
        self.health_check_port = None
        if target_type != "lambda":
            self.health_check_port = healthcheck_port or "traffic-port"
        self.health_check_path = healthcheck_path
        self.health_check_interval_seconds = healthcheck_interval_seconds or 30
        self.health_check_timeout_seconds = healthcheck_timeout_seconds or 10
        self.ip_address_type = (
            ip_address_type or "ipv4" if self.protocol != "GENEVE" else None
        )
        self.health_check_enabled = (
            healthcheck_enabled.lower() == "true"
            if healthcheck_enabled in ["true", "false"]
            else True
        )
        self.healthy_threshold_count = healthy_threshold_count or "5"
        self.unhealthy_threshold_count = unhealthy_threshold_count or "2"
        self.load_balancer_arns: List[str] = []
        if self.health_check_protocol != "TCP":
            self.matcher: Dict[str, Any] = matcher or {"HttpCode": "200"}
            self.health_check_path = self.health_check_path
            if target_type != "lambda":
                self.health_check_port = self.health_check_port or str(self.port)
        self.target_type = target_type or "instance"

        self.attributes = {
            "deregistration_delay.timeout_seconds": 300,
            "stickiness.enabled": "false",
            "load_balancing.algorithm.type": "round_robin",
            "load_balancing.cross_zone.enabled": "use_load_balancer_configuration",
            "slow_start.duration_seconds": 0,
            "waf.fail_open.enabled": "false",
        }
        if target_type == "lambda":
            self.attributes["lambda.multi_value_headers.enabled"] = "false"
        if self.protocol in ["HTTP", "HTTPS"]:
            self.attributes["stickiness.type"] = "lb_cookie"
        if self.protocol in ["TCP", "UDP", "TCP_UDP"]:
            self.attributes["stickiness.type"] = "source_ip"

        self.targets: Dict[str, Dict[str, Any]] = OrderedDict()
        self.deregistered_targets: Dict[str, Dict[str, Any]] = OrderedDict()
        self.terminated_targets: Dict[str, Dict[str, Any]] = OrderedDict()

    @property
    def physical_resource_id(self) -> str:
        return self.arn

    def register(self, targets: List[Dict[str, Any]]) -> None:
        for target in targets:
            if instance := self.ec2_backend.get_instance_by_id(target["id"]):
                if instance.state != "running":
                    raise TargetNotRunning(instance_id=target["id"])

            self.targets[target["id"]] = {
                "id": target["id"],
                "port": target.get("port", self.port),
            }
            self.deregistered_targets.pop(target["id"], None)

    def deregister(self, targets: List[Dict[str, Any]]) -> None:
        for target in targets:
            t = self.targets.pop(target["id"], None)
            if not t:
                raise InvalidTargetError()
            self.deregistered_targets[target["id"]] = t

    def deregister_terminated_instances(self, instance_ids: List[str]) -> None:
        for target_id in list(self.targets.keys()):
            if target_id in instance_ids:
                t = self.targets.pop(target_id)
                self.terminated_targets[target_id] = t

    def health_for(
        self, target: Dict[str, Any], ec2_backend: EC2Backend
    ) -> FakeHealthStatus:
        t = self.targets.get(target["id"])
        if t is None:
            port = self.port
            if "port" in target:
                port = target["port"]
            if target["id"] in self.deregistered_targets:
                return FakeHealthStatus(
                    target["id"],
                    port,
                    self.health_check_port,
                    "unused",
                    "Target.NotRegistered",
                    "Target is not registered to the target group",
                )
            if target["id"] in self.terminated_targets:
                return FakeHealthStatus(
                    target["id"],
                    port,
                    self.health_check_port,
                    "draining",
                    "Target.DeregistrationInProgress",
                    "Target deregistration is in progress",
                )

            if target["id"].startswith("i-"):  # EC2 instance ID
                return FakeHealthStatus(
                    target["id"],
                    target.get("Port", 80),
                    self.health_check_port,
                    "unused",
                    "Target.NotRegistered",
                    "Target is not registered to the target group",
                )

            return FakeHealthStatus(
                target["id"],
                port,
                self.health_check_port,
                "unavailable",
                "Target.NotRegistered",
                "Target is not registered",
            )
        if t["id"].startswith("i-"):  # EC2 instance ID
            instance = ec2_backend.get_instance_by_id(t["id"])
            if instance and instance.state == "stopped":
                return FakeHealthStatus(
                    t["id"],
                    t["port"],
                    self.health_check_port,
                    "unused",
                    "Target.InvalidState",
                    "Target is in the stopped state",
                )
        return FakeHealthStatus(t["id"], t["port"], self.health_check_port, "healthy")

    @staticmethod
    def cloudformation_name_type() -> str:
        return "Name"

    @staticmethod
    def cloudformation_type() -> str:
        # https://docs.aws.amazon.com/AWSCloudFormation/latest/UserGuide/aws-resource-elasticloadbalancingv2-targetgroup.html
        return "AWS::ElasticLoadBalancingV2::TargetGroup"

    @classmethod
    def create_from_cloudformation_json(  # type: ignore[misc]
        cls,
        resource_name: str,
        cloudformation_json: Any,
        account_id: str,
        region_name: str,
        **kwargs: Any,
    ) -> "FakeTargetGroup":
        properties = cloudformation_json["Properties"]

        elbv2_backend = elbv2_backends[account_id][region_name]

        vpc_id = properties.get("VpcId")
        protocol = properties.get("Protocol")
        port = properties.get("Port")
        healthcheck_protocol = properties.get("HealthCheckProtocol")
        healthcheck_port = properties.get("HealthCheckPort")
        healthcheck_path = properties.get("HealthCheckPath")
        healthcheck_interval_seconds = properties.get("HealthCheckIntervalSeconds")
        healthcheck_timeout_seconds = properties.get("HealthCheckTimeoutSeconds")
        healthy_threshold_count = properties.get("HealthyThresholdCount")
        unhealthy_threshold_count = properties.get("UnhealthyThresholdCount")
        matcher = properties.get("Matcher")
        target_type = properties.get("TargetType")

        target_group = elbv2_backend.create_target_group(
            name=resource_name,
            vpc_id=vpc_id,
            protocol=protocol,
            port=port,
            healthcheck_protocol=healthcheck_protocol,
            healthcheck_port=healthcheck_port,
            healthcheck_path=healthcheck_path,
            healthcheck_interval_seconds=healthcheck_interval_seconds,
            healthcheck_timeout_seconds=healthcheck_timeout_seconds,
            healthy_threshold_count=healthy_threshold_count,
            unhealthy_threshold_count=unhealthy_threshold_count,
            matcher=matcher,
            target_type=target_type,
        )
        return target_group

    @property
    def ec2_backend(self) -> EC2Backend:
        return ec2_backends[self.account_id][self.region_name]


class FakeListener(CloudFormationModel):
    def __init__(
        self,
        load_balancer_arn: str,
        arn: str,
        protocol: str,
        port: str,
        ssl_policy: str,
        certificate: Optional[str],
        default_actions: List["FakeAction"],
        alpn_policy: Optional[List[str]],
    ):
        self.load_balancer_arn = load_balancer_arn
        self.arn = arn
        self.protocol = (protocol or "").upper()
        self.port = port
        self.ssl_policy = ssl_policy
        self.certificate = certificate
        self.certificates = [certificate] if certificate is not None else []
        self.default_actions = default_actions
        self.alpn_policy = alpn_policy or []
        self._non_default_rules: Dict[str, FakeListenerRule] = OrderedDict()
        self._default_rule: Dict[int, FakeRule] = OrderedDict()
        self._default_rule[0] = FakeRule(
            listener_arn=self.arn,
            conditions=[],
            priority="default",
            actions=default_actions,
            is_default=True,
        )

        self.attributes = {"tcp.idle_timeout.seconds": "350"}

    @property
    def physical_resource_id(self) -> str:
        return self.arn

    @property
    def rules(self) -> Dict[Any, "FakeRule"]:  # type: ignore[misc]
        return OrderedDict(
            list(self._non_default_rules.items()) + list(self._default_rule.items())  # type: ignore
        )

    def remove_rule(self, arn: str) -> None:
        self._non_default_rules.pop(arn)

    def register(self, arn: str, rule: "FakeListenerRule") -> None:
        self._non_default_rules[arn] = rule
        sorted(self._non_default_rules.values(), key=lambda x: x.priority)

    @staticmethod
    def cloudformation_name_type() -> str:
        return ""

    @staticmethod
    def cloudformation_type() -> str:
        # https://docs.aws.amazon.com/AWSCloudFormation/latest/UserGuide/aws-resource-elasticloadbalancingv2-listener.html
        return "AWS::ElasticLoadBalancingV2::Listener"

    @classmethod
    def create_from_cloudformation_json(  # type: ignore[misc]
        cls,
        resource_name: str,
        cloudformation_json: Any,
        account_id: str,
        region_name: str,
        **kwargs: Any,
    ) -> "FakeListener":
        properties = cloudformation_json["Properties"]

        elbv2_backend = elbv2_backends[account_id][region_name]
        load_balancer_arn = properties.get("LoadBalancerArn")
        protocol = properties.get("Protocol")
        port = properties.get("Port")
        ssl_policy = properties.get("SslPolicy")
        certificates = properties.get("Certificates")

        default_actions = elbv2_backend.convert_and_validate_properties(properties)
        certificates = elbv2_backend.convert_and_validate_certificates(certificates)
        if certificates:
            certificate = certificates[0].get("certificate_arn")
        else:
            certificate = None
        listener = elbv2_backend.create_listener(
            load_balancer_arn, protocol, port, ssl_policy, certificate, default_actions
        )
        return listener

    @classmethod
    def update_from_cloudformation_json(  # type: ignore[misc]
        cls,
        original_resource: Any,
        new_resource_name: str,
        cloudformation_json: Any,
        account_id: str,
        region_name: str,
    ) -> "FakeListener":
        properties = cloudformation_json["Properties"]

        elbv2_backend = elbv2_backends[account_id][region_name]
        protocol = properties.get("Protocol")
        port = properties.get("Port")
        ssl_policy = properties.get("SslPolicy")
        certificates = properties.get("Certificates")

        default_actions = elbv2_backend.convert_and_validate_properties(properties)
        certificates = elbv2_backend.convert_and_validate_certificates(certificates)
        listener = elbv2_backend.modify_listener(
            original_resource.arn,
            port,
            protocol,
            ssl_policy,
            certificates,
            default_actions,
        )
        return listener


class FakeListenerRule(CloudFormationModel):
    def __init__(
        self,
        listener_arn: str,
        arn: str,
        conditions: List[Dict[str, Any]],
        priority: int,
        actions: List["FakeAction"],
    ):
        self.listener_arn = listener_arn
        self.arn = arn
        self.conditions = conditions
        self.actions = actions
        self.priority = priority
        self.is_default = False

    @property
    def physical_resource_id(self) -> str:
        return self.arn

    @staticmethod
    def cloudformation_type() -> str:
        # https://docs.aws.amazon.com/AWSCloudFormation/latest/UserGuide/aws-resource-elasticloadbalancingv2-listenerrule.html
        return "AWS::ElasticLoadBalancingV2::ListenerRule"

    @classmethod
    def create_from_cloudformation_json(  # type: ignore[misc]
        cls,
        resource_name: str,
        cloudformation_json: Any,
        account_id: str,
        region_name: str,
        **kwargs: Any,
    ) -> "FakeListenerRule":
        properties = cloudformation_json["Properties"]
        elbv2_backend = elbv2_backends[account_id][region_name]
        listener_arn = properties.get("ListenerArn")
        priority = properties.get("Priority")
        conditions = properties.get("Conditions")

        actions = elbv2_backend.convert_and_validate_action_properties(properties)
        listener_rule = elbv2_backend.create_rule(
            listener_arn, conditions, priority, actions
        )
        return listener_rule

    @classmethod
    def update_from_cloudformation_json(  # type: ignore[misc]
        cls,
        original_resource: Any,
        new_resource_name: str,
        cloudformation_json: Any,
        account_id: str,
        region_name: str,
    ) -> "FakeListenerRule":
        properties = cloudformation_json["Properties"]

        elbv2_backend = elbv2_backends[account_id][region_name]
        conditions = properties.get("Conditions")

        actions = elbv2_backend.convert_and_validate_action_properties(properties)
        listener_rule = elbv2_backend.modify_rule(
            original_resource.arn, conditions, actions
        )
        return listener_rule  # type: ignore[return-value]


class FakeRule(BaseModel):
    def __init__(
        self,
        listener_arn: str,
        conditions: List[Dict[str, Any]],
        priority: Any,
        actions: List["FakeAction"],
        is_default: bool,
    ):
        self.listener_arn = listener_arn
        self.arn = (
            listener_arn.replace(":listener/", ":listener-rule/") + f"/{id(self)}"
        )
        self.conditions = conditions
        self.priority = priority  # int or 'default'
        self.actions = actions
        self.is_default = is_default


class FakeAction(BaseModel):
    def __init__(self, data: Dict[str, Any]):
        self.data = data
        self.type = data.get("Type")

        if "ForwardConfig" in self.data:
            if "TargetGroupStickinessConfig" not in self.data["ForwardConfig"]:
                self.data["ForwardConfig"]["TargetGroupStickinessConfig"] = {
                    "Enabled": False
                }

            for target_group in self.data["ForwardConfig"].get("TargetGroups", []):
                target_group.setdefault("Weight", 1)
        # Dynamically give our Action class all properties of the source data.
        self.__dict__.update(self.data)


class FakeBackend(BaseModel):
    def __init__(self, instance_port: str):
        self.instance_port = instance_port
        self.policy_names: List[str] = []

    def __repr__(self) -> str:
        return f"FakeBackend(inp: {self.instance_port}, policies: {self.policy_names})"


class AvailabilityZone:
    def __init__(self, subnet: Subnet) -> None:
        self.subnet_id = subnet.id
        self.zone_name = subnet.availability_zone


class FakeLoadBalancer(CloudFormationModel):
    VALID_ATTRS = {
        "access_logs.s3.enabled",
        "access_logs.s3.bucket",
        "access_logs.s3.prefix",
        "connection_logs.s3.bucket",
        "connection_logs.s3.enabled",
        "connection_logs.s3.prefix",
        "client_keep_alive.seconds",
        "deletion_protection.enabled",
        "dns_record.client_routing_policy",
        "idle_timeout.timeout_seconds",
        "ipv6.deny_all_igw_traffic",
        "load_balancing.cross_zone.enabled",
        "routing.http.desync_mitigation_mode",
        "routing.http.drop_invalid_header_fields.enabled",
        "routing.http.preserve_host_header.enabled",
        "routing.http.x_amzn_tls_version_and_cipher_suite.enabled",
        "routing.http.xff_client_port.enabled",
        "routing.http.xff_header_processing.mode",
        "routing.http2.enabled",
        "waf.fail_open.enabled",
        "zonal_shift.config.enabled",
    }

    def __init__(
        self,
        name: str,
        security_groups: List[str],
        subnets: List[Subnet],
        vpc_id: str,
        arn: str,
        dns_name: str,
        state: str,
        scheme: str = "internet-facing",
        loadbalancer_type: Optional[str] = None,
    ):
        self.name = name
        self.created_time = iso_8601_datetime_with_milliseconds()
        self.scheme = scheme
        self.security_groups = security_groups
        self.subnets = subnets or []
        self.vpc_id = vpc_id
        self.listeners: Dict[str, FakeListener] = OrderedDict()
        self.tags: Dict[str, Any] = {}
        self.arn = arn
        self.dns_name = dns_name
        self.state = state
        self.type = loadbalancer_type or "application"

        self.ip_address_type = "ipv4"
        self.attrs = {
            # "access_logs.s3.enabled": "false",  # commented out for TF compatibility
            "access_logs.s3.bucket": "",
            "access_logs.s3.prefix": "",
            "deletion_protection.enabled": "false",
            # "idle_timeout.timeout_seconds": "60",  # commented out for TF compatibility
            "load_balancing.cross_zone.enabled": "false",
        }
        # TODO: This was hardcoded in the original XML template and still needs to be implemented.
        self.canonical_hosted_zone_id = "Z2P70J7EXAMPLE"

    @property
    def load_balancer_state(self) -> dict[str, str]:
        return {"Code": self.state}

    @property
    def availability_zones(self) -> list[AvailabilityZone]:
        return [AvailabilityZone(subnet) for subnet in self.subnets]

    @property
    def physical_resource_id(self) -> str:
        return self.arn

    def activate(self) -> None:
        if self.state == "provisioning":
            self.state = "active"

    def delete(self, account_id: str, region: str) -> None:
        """Not exposed as part of the ELB API - used for CloudFormation."""
        elbv2_backends[account_id][region].delete_load_balancer(self.arn)

    @staticmethod
    def cloudformation_name_type() -> str:
        return "Name"

    @staticmethod
    def cloudformation_type() -> str:
        # https://docs.aws.amazon.com/AWSCloudFormation/latest/UserGuide/aws-resource-elasticloadbalancingv2-loadbalancer.html
        return "AWS::ElasticLoadBalancingV2::LoadBalancer"

    @classmethod
    def create_from_cloudformation_json(  # type: ignore[misc]
        cls,
        resource_name: str,
        cloudformation_json: Any,
        account_id: str,
        region_name: str,
        **kwargs: Any,
    ) -> "FakeLoadBalancer":
        properties = cloudformation_json["Properties"]

        elbv2_backend = elbv2_backends[account_id][region_name]

        security_groups = properties.get("SecurityGroups")
        subnet_ids = properties.get("Subnets")
        scheme = properties.get("Scheme", "internet-facing")

        load_balancer = elbv2_backend.create_load_balancer(
            resource_name, security_groups, subnet_ids, scheme=scheme
        )
        return load_balancer

    @classmethod
    def has_cfn_attr(cls, attr: str) -> bool:
        return attr in [
            "DNSName",
            "LoadBalancerName",
            "CanonicalHostedZoneID",
            "LoadBalancerFullName",
            "SecurityGroups",
        ]

    def get_cfn_attribute(self, attribute_name: str) -> Any:
        """
        Implemented attributes:
        * DNSName
        * LoadBalancerName

        Not implemented:
        * CanonicalHostedZoneID
        * LoadBalancerFullName
        * SecurityGroups

        This method is similar to models.py:FakeLoadBalancer.get_cfn_attribute()
        """
        from moto.cloudformation.exceptions import UnformattedGetAttTemplateException

        not_implemented_yet = [
            "CanonicalHostedZoneID",
            "LoadBalancerFullName",
            "SecurityGroups",
        ]
        if attribute_name == "DNSName":
            return self.dns_name
        elif attribute_name == "LoadBalancerName":
            return self.name
        elif attribute_name in not_implemented_yet:
            raise NotImplementedError(
                f'"Fn::GetAtt" : [ "{{0}}" , "{attribute_name}" ]"'
            )
        else:
            raise UnformattedGetAttTemplateException()


class ELBv2Backend(BaseBackend):
    def __init__(self, region_name: str, account_id: str):
        super().__init__(region_name, account_id)
        self.target_groups: Dict[str, FakeTargetGroup] = OrderedDict()
        self.load_balancers: Dict[str, FakeLoadBalancer] = OrderedDict()
        self.tagging_service = TaggingService()

    @property
    def ec2_backend(self) -> EC2Backend:
        return ec2_backends[self.account_id][self.region_name]

    def create_load_balancer(
        self,
        name: str,
        security_groups: List[str],
        subnet_ids: List[str],
        subnet_mappings: Optional[List[Dict[str, str]]] = None,
        scheme: str = "internet-facing",
        loadbalancer_type: Optional[str] = None,
        tags: Optional[List[Dict[str, str]]] = None,
    ) -> FakeLoadBalancer:
        vpc_id = None
        subnets = []
        state = "provisioning"

        if not subnet_ids and not subnet_mappings:
            raise SubnetNotFoundError()
        for subnet_id in subnet_ids:
            subnet = self.ec2_backend.get_subnet(subnet_id)
            if subnet is None:
                raise SubnetNotFoundError()
            subnets.append(subnet)
        for subnet_mapping in subnet_mappings or []:
            subnet_id = subnet_mapping["SubnetId"]
            subnet = self.ec2_backend.get_subnet(subnet_id)
            if subnet is None:
                raise SubnetNotFoundError()
            subnets.append(subnet)

        vpc_id = subnets[0].vpc_id
        lb_id = mock_random.get_random_hex(length=16)
        arn = make_arn_for_load_balancer(
            account_id=self.account_id,
            lb_id=lb_id,
            name=name,
            region_name=self.region_name,
        )
        dns_name = f"{name}-1.{self.region_name}.elb.amazonaws.com"

        if arn in self.load_balancers:
            raise DuplicateLoadBalancerName()

        new_load_balancer = FakeLoadBalancer(
            name=name,
            security_groups=security_groups,
            arn=arn,
            scheme=scheme,
            subnets=subnets,
            vpc_id=vpc_id,
            dns_name=dns_name,
            state=state,
            loadbalancer_type=loadbalancer_type,
        )
        self.load_balancers[arn] = new_load_balancer
        self.tagging_service.tag_resource(arn, tags)
        return new_load_balancer

    def convert_and_validate_action_properties(
        self, properties: Dict[str, Any]
    ) -> List[Dict[str, Any]]:
        # transform Actions to confirm with the rest of the code and XML templates
        default_actions = []
        for i, action in enumerate(properties["Actions"]):
            action_type = action["Type"]
            if action_type in ALLOWED_ACTIONS:
                default_actions.append(action)
            else:
                raise InvalidActionTypeError(action_type, i + 1)
        return default_actions

    def create_rule(
        self,
        listener_arn: str,
        conditions: List[Dict[str, Any]],
        priority: int,
        actions: List[Dict[str, Any]],
        tags: Optional[List[Dict[str, str]]] = None,
    ) -> FakeListenerRule:
        fake_actions = [FakeAction(action) for action in actions]
        listeners = self.describe_listeners(None, [listener_arn])
        if not listeners:
            raise ListenerNotFoundError()
        listener = listeners[0]

        # validate conditions
        # see: https://docs.aws.amazon.com/cli/latest/reference/elbv2/create-rule.html
        self._validate_conditions(conditions)

        # TODO: check QueryStringConfig condition
        # TODO: check HttpRequestMethodConfig condition
        # TODO: check SourceIpConfig condition
        # TODO: check pattern of value for 'host-header'
        # TODO: check pattern of value for 'path-pattern'

        # validate Priority
        for rule in listener.rules.values():
            if rule.priority == priority:
                raise PriorityInUseError()

        self._validate_actions(fake_actions)
        for action in fake_actions:
            if action.type == "forward":
                found_arns = self._get_target_group_arns_from(action_data=action.data)
                for arn in found_arns:
                    target_group = self.target_groups[arn]
                    target_group.load_balancer_arns.append(listener.load_balancer_arn)

        arn = listener_arn.replace(":listener/", ":listener-rule/")
        arn += f"/{mock_random.get_random_hex(16)}"

        # TODO: check for error 'TooManyRegistrationsForTargetId'
        # TODO: check for error 'TooManyRules'

        # create rule
        listener_rule = FakeListenerRule(
            listener.arn, arn, conditions, priority, fake_actions
        )
        listener.register(arn, listener_rule)
        self.tagging_service.tag_resource(arn, tags)
        return listener_rule

    def _validate_conditions(self, conditions: List[Dict[str, Any]]) -> None:
        for condition in conditions:
            if "Field" in condition:
                field = condition["Field"]
                if field not in [
                    "host-header",
                    "http-header",
                    "http-request-method",
                    "path-pattern",
                    "query-string",
                    "source-ip",
                ]:
                    raise InvalidConditionFieldError(field)
                if "Values" in condition and field not in [
                    "host-header",
                    "path-pattern",
                ]:
                    raise InvalidConditionValueError(
                        f"The 'Values' field is not compatible with '{field}'"
                    )
                else:
                    method_name = "_validate_" + field.replace("-", "_") + "_condition"
                    func = getattr(self, method_name)
                    func(condition)

    def _validate_host_header_condition(self, condition: Dict[str, Any]) -> None:
        values = None
        if "HostHeaderConfig" in condition:
            values = condition["HostHeaderConfig"]["Values"]
        elif "Values" in condition:
            values = condition["Values"]
            if len(values) > 1:
                raise InvalidConditionValueError(
                    "The 'host-header' field contains too many values; the limit is '1'"
                )
        if values is None or len(values) == 0:
            raise InvalidConditionValueError("A condition value must be specified")
        for value in values:
            if len(value) > 128:
                raise InvalidConditionValueError(
                    "The 'host-header' value is too long; the limit is '128'"
                )

    def _validate_http_header_condition(self, condition: Dict[str, Any]) -> None:
        if "HttpHeaderConfig" in condition:
            config = condition["HttpHeaderConfig"]
            name = config.get("HttpHeaderName")
            if len(name) > 40:
                raise InvalidConditionValueError(
                    "The 'HttpHeaderName' value is too long; the limit is '40'"
                )
            values = config["Values"]
            for value in values:
                if len(value) > 128:
                    raise InvalidConditionValueError(
                        "The 'http-header' value is too long; the limit is '128'"
                    )
        else:
            raise InvalidConditionValueError(
                "A 'HttpHeaderConfig' must be specified with 'http-header'"
            )

    def _validate_http_request_method_condition(
        self, condition: Dict[str, Any]
    ) -> None:
        if "HttpRequestMethodConfig" in condition:
            for value in condition["HttpRequestMethodConfig"]["Values"]:
                if len(value) > 40:
                    raise InvalidConditionValueError(
                        "The 'http-request-method' value is too long; the limit is '40'"
                    )
                if not re.match("[A-Z_-]+", value):
                    raise InvalidConditionValueError(
                        "The 'http-request-method' value is invalid; the allowed characters are A-Z, hyphen and underscore"
                    )
        else:
            raise InvalidConditionValueError(
                "A 'HttpRequestMethodConfig' must be specified with 'http-request-method'"
            )

    def _validate_path_pattern_condition(self, condition: Dict[str, Any]) -> None:
        values = None
        if "PathPatternConfig" in condition:
            values = condition["PathPatternConfig"]["Values"]
        elif "Values" in condition:
            values = condition["Values"]
            if len(values) > 1:
                raise InvalidConditionValueError(
                    "The 'path-pattern' field contains too many values; the limit is '1'"
                )
        if values is None or len(values) == 0:
            raise InvalidConditionValueError("A condition value must be specified")
        if condition.get("Values") and condition.get("PathPatternConfig"):
            raise InvalidConditionValueError(
                "You cannot provide both Values and 'PathPatternConfig' for a condition of type 'path-pattern'"
            )
        for value in values:
            if len(value) > 128:
                raise InvalidConditionValueError(
                    "The 'path-pattern' value is too long; the limit is '128'"
                )

    def _validate_source_ip_condition(self, condition: Dict[str, Any]) -> None:
        if "SourceIpConfig" in condition:
            values = condition["SourceIpConfig"].get("Values", [])
            if len(values) == 0:
                raise InvalidConditionValueError(
                    "A 'source-ip' value must be specified"
                )
        else:
            raise InvalidConditionValueError(
                "A 'SourceIpConfig' must be specified with 'source-ip'"
            )

    def _validate_query_string_condition(self, condition: Dict[str, Any]) -> None:
        if "QueryStringConfig" in condition:
            config = condition["QueryStringConfig"]
            values = config["Values"]
            for value in values:
                if "Value" not in value:
                    raise InvalidConditionValueError(
                        "A 'Value' must be specified in 'QueryStringKeyValuePair'"
                    )
                if "Key" in value and len(value["Key"]) > 128:
                    raise InvalidConditionValueError(
                        "The 'Key' value is too long; the limit is '128'"
                    )
                if len(value["Value"]) > 128:
                    raise InvalidConditionValueError(
                        "The 'Value' value is too long; the limit is '128'"
                    )
        else:
            raise InvalidConditionValueError(
                "A 'QueryStringConfig' must be specified with 'query-string'"
            )

    def _get_target_group_arns_from(self, action_data: Dict[str, Any]) -> List[Any]:
        if "TargetGroupArn" in action_data:
            return [action_data["TargetGroupArn"]]
        elif "ForwardConfig" in action_data:
            return [
                tg["TargetGroupArn"]
                for tg in action_data["ForwardConfig"].get("TargetGroups", [])
            ]
        else:
            return []

    def _validate_actions(self, actions: List[FakeAction]) -> None:
        # validate Actions
        target_group_arns = [
            target_group.arn for target_group in self.target_groups.values()
        ]
        for i, action in enumerate(actions):
            index = i + 1
            action_type = action.type
            if action_type == "forward":
                found_arns = self._get_target_group_arns_from(action_data=action.data)
                for target_group_arn in found_arns:
                    if target_group_arn not in target_group_arns:
                        raise ActionTargetGroupNotFoundError(target_group_arn)
            elif action_type == "fixed-response":
                self._validate_fixed_response_action(action, i, index)
            elif action_type in [
                "redirect",
                "authenticate-cognito",
                "authenticate-oidc",
            ]:
                pass
            # pass if listener rule has forward_config as an Action property
            elif action_type == "forward" and "ForwardConfig" in action.data.keys():
                pass
            else:
                raise InvalidActionTypeError(action_type, index)

    def _validate_port_and_protocol(
        self, loadbalancer_type: str, port: Optional[str], protocol: Optional[str]
    ) -> None:
        if loadbalancer_type in {"application", "network"}:
            if port is None:
                raise ValidationError("A listener port must be specified")
            if protocol is None:
                raise ValidationError("A listener protocol must be specified")

            valid_application_protocols = ["HTTP", "HTTPS"]
            valid_network_protocols = ["UDP", "TCP", "TLS", "TCP_UDP"]
            all_valid_protocols = (
                valid_application_protocols + valid_network_protocols + ["GENEVE"]
            )

            if protocol not in all_valid_protocols:
                raise InvalidProtocolValue(
                    protocol, valid_application_protocols + valid_network_protocols
                )

            if loadbalancer_type == "application":
                if protocol not in valid_application_protocols:
                    raise InvalidProtocol(protocol, valid_application_protocols)
            else:
                if protocol not in valid_network_protocols:
                    raise InvalidProtocol(protocol, valid_network_protocols)
        else:
            if port is not None:
                raise ValidationError(
                    "A port cannot be specified for gateway listeners"
                )
            if protocol is not None:
                raise ValidationError(
                    "A protocol cannot be specified for gateway listeners"
                )

    def _validate_fixed_response_action(
        self, action: FakeAction, i: int, index: int
    ) -> None:
        status_code = action.data.get("FixedResponseConfig", {}).get("StatusCode")
        if status_code is None:
            raise ParamValidationError(
                report='Missing required parameter in Actions[%s].FixedResponseConfig: "StatusCode"'
                % i
            )
        expression = r"^(2|4|5)\d\d$"
        if not re.match(expression, status_code):
            raise InvalidStatusCodeActionTypeError(
                f"1 validation error detected: Value '{status_code}' at 'actions.{index}.member.fixedResponseConfig.statusCode' failed to satisfy constraint: \
Member must satisfy regular expression pattern: {expression}"
            )
        content_type = action.data["FixedResponseConfig"].get("ContentType")
        if content_type and content_type not in [
            "text/plain",
            "text/css",
            "text/html",
            "application/javascript",
            "application/json",
        ]:
            raise InvalidLoadBalancerActionException(
                "The ContentType must be one of:'text/html', 'application/json', 'application/javascript', 'text/css', 'text/plain'"
            )

    def create_target_group(self, name: str, **kwargs: Any) -> FakeTargetGroup:
        protocol = kwargs.get("protocol")
        target_type = kwargs.get("target_type")

        if len(name) > 32:
            raise InvalidTargetGroupNameError(
                f"Target group name '{name}' cannot be longer than '32' characters"
            )
        if not re.match(r"^[a-zA-Z0-9\-]+$", name):
            raise InvalidTargetGroupNameError(
                f"Target group name '{name}' can only contain characters that are alphanumeric characters or hyphens(-)"
            )

        # undocumented validation
        if not re.match(r"(?!.*--)(?!^-)(?!.*-$)^[A-Za-z0-9-]+$", name):
            raise InvalidTargetGroupNameError(
                "1 validation error detected: Value '%s' at 'targetGroup.targetGroupArn.targetGroupName' failed to satisfy constraint: Member must satisfy regular expression pattern: (?!.*--)(?!^-)(?!.*-$)^[A-Za-z0-9-]+$"
                % name
            )

        if name.startswith("-") or name.endswith("-"):
            raise InvalidTargetGroupNameError(
                f"Target group name '{name}' cannot begin or end with '-'"
            )
        for target_group in self.target_groups.values():
            if target_group.name == name:
                raise DuplicateTargetGroupName()

        valid_protocols = ["HTTPS", "HTTP", "TCP", "TLS", "UDP", "TCP_UDP", "GENEVE"]
        if (
            kwargs.get("healthcheck_protocol")
            and kwargs["healthcheck_protocol"] not in valid_protocols
        ):
            raise InvalidConditionValueError(
                f"Value {kwargs['healthcheck_protocol']} at 'healthCheckProtocol' failed to satisfy constraint: "
                f"Member must satisfy enum value set: {valid_protocols}"
            )
        if kwargs.get("protocol") and kwargs["protocol"] not in valid_protocols:
            raise InvalidConditionValueError(
                f"Value {kwargs['protocol']} at 'protocol' failed to satisfy constraint: "
                f"Member must satisfy enum value set: {valid_protocols}"
            )

        if (
            kwargs.get("matcher")
            and kwargs["matcher"].get("HttpCode")
            and FakeTargetGroup.HTTP_CODE_REGEX.match(kwargs["matcher"]["HttpCode"])
            is None
        ):
            raise RESTError(
                "InvalidParameterValue",
                "HttpCode must be like 200 | 200-399 | 200,201 ...",
            )

        if target_type in ("instance", "ip", "alb"):
            for param in ("protocol", "port", "vpc_id"):
                if not kwargs.get(param):
                    param = "VPC ID" if param == "vpc_id" else param.lower()
                    raise ValidationError(f"A {param} must be specified")

        if target_type == "lambda":
            for param in ["protocol", "port", "vpc_id"]:
                if kwargs.get(param) is not None:
                    param = "VPC ID" if param == "vpc_id" else param.capitalize()
                    raise ValidationError(
                        f"{param} cannot be specified for target groups with target type 'lambda'"
                    )

        if kwargs.get("vpc_id"):
            from moto.ec2.exceptions import InvalidVPCIdError

            try:
                self.ec2_backend.get_vpc(kwargs.get("vpc_id") or "")
            except InvalidVPCIdError:
                raise ValidationError(
                    f"The VPC ID '{kwargs.get('vpc_id')}' is not found"
                )

        kwargs_patch = {}

        conditions: Dict[str, Any] = {
            "target_lambda": {
                "healthcheck_interval_seconds": 35,
                "healthcheck_timeout_seconds": 30,
                "unhealthy_threshold_count": 2,
                "healthcheck_enabled": "false",
                "healthcheck_path": "/",
            },
            "target_alb": {
                "healthcheck_protocol": "HTTP",
                "healthcheck_path": "/",
                "healthcheck_timeout_seconds": 6,
                "matcher": {"HttpCode": "200-399"},
            },
            "protocol_GENEVE": {
                "healthcheck_interval_seconds": 10,
                "healthcheck_port": 80,
                "healthcheck_timeout_seconds": 5,
                "healthcheck_protocol": "TCP",
                "unhealthy_threshold_count": 2,
            },
            "protocol_HTTP_HTTPS": {
                "healthcheck_timeout_seconds": 5,
                "protocol_version": "HTTP1",
                "healthcheck_path": "/",
                "unhealthy_threshold_count": 2,
                "healthcheck_interval_seconds": 30,
            },
            "protocol_TCP": {
                "healthcheck_timeout_seconds": 10,
            },
            "protocol_TCP_TCP_UDP_UDP_TLS": {
                "healthcheck_protocol": "TCP",
                "unhealthy_threshold_count": 2,
                "healthcheck_interval_seconds": 30,
            },
        }

        if target_type == "lambda":
            kwargs_patch.update(
                {k: kwargs.get(k) or v for k, v in conditions["target_lambda"].items()}
            )

        if protocol == "GENEVE":
            kwargs_patch.update(
                {
                    k: kwargs.get(k) or v
                    for k, v in conditions["protocol_GENEVE"].items()
                }
            )

        if protocol in ("HTTP", "HTTPS"):
            kwargs_patch.update(
                {
                    k: kwargs.get(k) or v
                    for k, v in conditions["protocol_HTTP_HTTPS"].items()
                }
            )

        if protocol == "TCP":
            kwargs_patch.update(
                {k: kwargs.get(k) or v for k, v in conditions["protocol_TCP"].items()}
            )

        if protocol in ("TCP", "TCP_UDP", "UDP", "TLS"):
            kwargs_patch.update(
                {
                    k: kwargs.get(k) or v
                    for k, v in conditions["protocol_TCP_TCP_UDP_UDP_TLS"].items()
                }
            )

        if target_type == "alb":
            kwargs_patch.update(
                {k: kwargs.get(k) or v for k, v in conditions["target_alb"].items()}
            )

        original_kwargs = dict(kwargs)
        kwargs.update(kwargs_patch)

        healthcheck_timeout_seconds = kwargs.get("healthcheck_timeout_seconds", 10)
        healthcheck_interval_seconds = kwargs.get("healthcheck_interval_seconds", 30)

        if (
            healthcheck_timeout_seconds is not None
            and healthcheck_interval_seconds is not None
        ):
            if healthcheck_interval_seconds < healthcheck_timeout_seconds:
                message = f"Health check timeout '{healthcheck_timeout_seconds}' must be smaller than or equal to the interval '{healthcheck_interval_seconds}'"
                if protocol in ("HTTP", "HTTPS"):
                    message = f"Health check timeout '{healthcheck_timeout_seconds}' must be smaller than the interval '{healthcheck_interval_seconds}'"
                raise ValidationError(message)
            both_values_supplied = (
                original_kwargs.get("healthcheck_timeout_seconds") is not None
                and original_kwargs.get("healthcheck_interval_seconds") is not None
            )
            if (
                both_values_supplied
                and healthcheck_interval_seconds == healthcheck_timeout_seconds
                and protocol in ("HTTP", "HTTPS")
            ):
                raise ValidationError(
                    f"Health check timeout '{healthcheck_timeout_seconds}' must be smaller than the interval '{healthcheck_interval_seconds}'"
                )

        tags = kwargs.pop("tags", None)
        target_group = FakeTargetGroup(
            name=name,
            account_id=self.account_id,
            region_name=self.region_name,
            **kwargs,
        )
        self.target_groups[target_group.arn] = target_group
        if tags:
            self.add_tags(resource_arns=[target_group.arn], tags=tags)
        return target_group

    def modify_target_group_attributes(
        self, target_group_arn: str, attributes: Dict[str, Any]
    ) -> None:
        target_group = self.target_groups.get(target_group_arn)
        if not target_group:
            raise TargetGroupNotFoundError()

        deregistration_delay_timeout_seconds = attributes.get(
            "deregistration_delay.timeout_seconds"
        )
        if deregistration_delay_timeout_seconds:
            if int(deregistration_delay_timeout_seconds) not in range(0, 3600):
                raise ValidationError(
                    f"'deregistration_delay.timeout_seconds' value '{deregistration_delay_timeout_seconds}' must be between '0-3600' inclusive"
                )
            if target_group.target_type == "lambda":
                raise InvalidConfigurationRequest(
                    "A target group with target type 'lambda' does not support the attribute deregistration_delay.timeout_seconds"
                )

        stickiness_type = attributes.get("stickiness.type")
        # TODO: strict type checking for app_cookie
        stickiness_cookie_name = attributes.get("stickiness.app_cookie.cookie_name")
        if stickiness_type:
            if target_group.protocol == "GENEVE":
                if stickiness_cookie_name:
                    # TODO: generalise error message
                    raise ValidationError(
                        "Target group attribute key 'stickiness.app_cookie.cookie_name' is not recognized"
                    )
                elif stickiness_type not in [
                    "source_ip_dest_ip_proto",
                    "source_ip_dest_ip",
                ]:
                    raise ValidationError(
                        f"'{stickiness_type}' must be one of [source_ip_dest_ip_proto, source_ip_dest_ip]"
                    )
            if stickiness_type == "source_ip":
                if target_group.protocol in ["HTTP", "HTTPS"]:
                    raise InvalidConfigurationRequest(
                        f"Stickiness type 'source_ip' is not supported for target groups with the {target_group.protocol} protocol"
                    )
                elif target_group.protocol == "TLS":
                    raise InvalidConfigurationRequest(
                        "You cannot enable stickiness on target groups with the TLS protocol"
                    )
            elif stickiness_type == "lb_cookie":
                if target_group.protocol in ["TCP", "TLS", "UDP", "TCP_UDP"]:
                    raise InvalidConfigurationRequest(
                        f"Stickiness type 'lb_cookie' is not supported for target groups with the {target_group.protocol} protocol"
                    )
            elif stickiness_type == "app_cookie":
                if not stickiness_cookie_name:
                    raise InvalidConfigurationRequest(
                        "You must set an application cookie name to enable stickiness of type 'app_cookie'"
                    )
                if target_group.protocol in ["TCP", "TLS", "UDP", "TCP_UDP", "GENEVE"]:
                    raise InvalidConfigurationRequest(
                        f"Stickiness type 'app_cookie' is not supported for target groups with the {target_group.protocol} protocol"
                    )
            elif stickiness_type in ["source_ip_dest_ip", "source_ip_dest_ip_proto"]:
                if target_group.protocol in [
                    "HTTP",
                    "HTTPS",
                    "TCP",
                    "TLS",
                    "UDP",
                    "TCP_UDP",
                ]:
                    raise ValidationError(
                        "'Stickiness type' must be one of [app_cookie, lb_cookie, source_ip]"
                    )

        target_group.attributes.update(attributes)

    def convert_and_validate_certificates(
        self, certificates: List[Dict[str, Any]]
    ) -> List[Dict[str, Any]]:
        # transform default certificate to conform with the rest of the code and XML templates
        for cert in certificates or []:
            cert["certificate_arn"] = cert["CertificateArn"]

        return certificates

    def convert_and_validate_properties(
        self, properties: Dict[str, Any]
    ) -> List[Dict[str, Any]]:
        # transform default actions to confirm with the rest of the code and XML templates
        # Caller: CF create/update for type "AWS::ElasticLoadBalancingV2::Listener"
        default_actions = []
        for i, action in enumerate(properties["DefaultActions"]):
            action_type = action["Type"]
            if action_type == "forward":
                default_actions.append(
                    {"Type": action_type, "TargetGroupArn": action["TargetGroupArn"]}
                )
            elif action_type in ALLOWED_ACTIONS:
                default_actions.append(action)
            else:
                raise InvalidActionTypeError(action_type, i + 1)
        return default_actions

    def create_listener(
        self,
        load_balancer_arn: str,
        protocol: str,
        port: str,
        ssl_policy: str,
        certificate: Optional[str],
        actions: List[Dict[str, Any]],
        alpn_policy: Optional[List[str]] = None,
        tags: Optional[List[Dict[str, str]]] = None,
    ) -> FakeListener:
        default_actions = [FakeAction(action) for action in actions]
        balancer = self.load_balancers.get(load_balancer_arn)
        if balancer is None:
            raise LoadBalancerNotFoundError()
        if port in balancer.listeners:
            raise DuplicateListenerError()

        self._validate_port_and_protocol(balancer.type, port, protocol)
        self._validate_actions(default_actions)

        arn = load_balancer_arn.replace(":loadbalancer/", ":listener/")
        arn += f"/{mock_random.get_random_hex(16)}"

        listener = FakeListener(
            load_balancer_arn,
            arn,
            protocol,
            port,
            ssl_policy,
            certificate,
            default_actions,
            alpn_policy,
        )
        if certificate and not re.search(f"{ARN_PARTITION_REGEX}:iam:", certificate):
            register_certificate(
                account_id=self.account_id,
                region=self.region_name,
                arn_certificate=certificate,
                arn_user=arn,
            )

        balancer.listeners[listener.arn] = listener
        for action in default_actions:
            if action.type == "forward":
                found_arns = self._get_target_group_arns_from(action_data=action.data)
                for arn in found_arns:
                    target_group = self.target_groups[arn]
                    target_group.load_balancer_arns.append(load_balancer_arn)

        self.tagging_service.tag_resource(listener.arn, tags)

        return listener

    def describe_load_balancers(
        self, arns: Optional[List[str]], names: Optional[List[str]]
    ) -> List[FakeLoadBalancer]:
        balancers = list(self.load_balancers.values())
        arns = arns or []
        names = names or []
        if not arns and not names:
            for balancer in balancers:
                balancer.activate()
            return balancers

        matched_balancers = []
        matched_balancer = None

        for arn in arns:
            for balancer in balancers:
                balancer.activate()
                if balancer.arn == arn:
                    matched_balancer = balancer
            if matched_balancer is None:
                raise LoadBalancerNotFoundError()
            elif matched_balancer not in matched_balancers:
                matched_balancers.append(matched_balancer)

        for name in names:
            for balancer in balancers:
                balancer.activate()
                if balancer.name == name:
                    matched_balancer = balancer
            if matched_balancer is None:
                raise LoadBalancerNotFoundError()
            elif matched_balancer not in matched_balancers:
                matched_balancers.append(matched_balancer)

        return matched_balancers

    def describe_rules(
        self, listener_arn: Optional[str], rule_arns: Optional[List[str]]
    ) -> List[FakeRule]:
        if listener_arn is None and not rule_arns:
            raise InvalidDescribeRulesRequest(
                "You must specify either listener rule ARNs or a listener ARN"
            )
        if listener_arn is not None and rule_arns is not None:
            raise InvalidDescribeRulesRequest(
                "Listener rule ARNs and a listener ARN cannot be specified at the same time"
            )
        if listener_arn:
            listener = self.describe_listeners(None, [listener_arn])[0]
            return list(listener.rules.values())

        # search for rule arns
        matched_rules = []
        for load_balancer_arn in self.load_balancers:
            listeners = self.load_balancers.get(load_balancer_arn).listeners.values()  # type: ignore
            for listener in listeners:
                for rule in listener.rules.values():
                    if rule.arn in rule_arns:  # type: ignore[operator]
                        matched_rules.append(rule)
        if len(matched_rules) != len(rule_arns):  # type: ignore
            raise RuleNotFoundError("One or more rules not found")
        return matched_rules

    def describe_target_groups(
        self,
        load_balancer_arn: Optional[str],
        target_group_arns: List[str],
        names: Optional[List[str]],
    ) -> Iterable[FakeTargetGroup]:
        args = sum(bool(arg) for arg in [load_balancer_arn, target_group_arns, names])

        if args > 1:
            raise ValidationError(
                "Target group names, target group ARNs, and a load balancer ARN cannot be specified at the same time"
            )

        if load_balancer_arn:
            if load_balancer_arn not in self.load_balancers:
                raise LoadBalancerNotFoundError()
            target_groups = [
                tg
                for tg in self.target_groups.values()
                if load_balancer_arn in tg.load_balancer_arns
            ]
            if target_groups is None or len(target_groups) == 0:
                raise TargetGroupNotFoundError()
            return sorted(target_groups, key=lambda tg: tg.name)

        if target_group_arns:
            try:
                target_groups = [self.target_groups[arn] for arn in target_group_arns]
                return sorted(target_groups, key=lambda tg: tg.name)
            except KeyError:
                raise TargetGroupNotFoundError()

        if names:
            matched = []
            for name in names:
                found = None
                for target_group in self.target_groups.values():
                    if target_group.name == name:
                        found = target_group
                if not found:
                    raise TargetGroupNotFoundError()
                matched.append(found)
            return sorted(matched, key=lambda tg: tg.name)

        return sorted(self.target_groups.values(), key=lambda tg: tg.name)

    def describe_listeners(
        self, load_balancer_arn: Optional[str], listener_arns: List[str]
    ) -> List[FakeListener]:
        if load_balancer_arn:
            if load_balancer_arn not in self.load_balancers:
                raise LoadBalancerNotFoundError()
            return list(self.load_balancers.get(load_balancer_arn).listeners.values())  # type: ignore

        matched = []
        for load_balancer in self.load_balancers.values():
            for listener_arn in listener_arns:
                listener = load_balancer.listeners.get(listener_arn)
                if listener:
                    matched.append(listener)
        if listener_arns and len(matched) == 0:
            raise ListenerNotFoundError()
        return matched

    def delete_load_balancer(self, arn: str) -> None:
        self.load_balancers.pop(arn, None)

    def delete_rule(self, arn: str) -> None:
        for load_balancer_arn in self.load_balancers:
            listeners = self.load_balancers.get(load_balancer_arn).listeners.values()  # type: ignore[union-attr]
            for listener in listeners:
                for rule in listener.rules.values():
                    if rule.arn == arn:
                        listener.remove_rule(rule.arn)
                        return

        # should raise RuleNotFound Error according to the AWS API doc
        # however, boto3 does't raise error even if rule is not found

    def delete_target_group(self, target_group_arn: str) -> Optional[FakeTargetGroup]:  # type: ignore[return]
        if target_group_arn not in self.target_groups:
            raise TargetGroupNotFoundError()

        target_group = self.target_groups[target_group_arn]
        if target_group:
            if self._any_listener_using(target_group_arn):
                raise ResourceInUseError(
                    f"The target group '{target_group_arn}' is currently in use by a listener or a rule"
                )
            del self.target_groups[target_group_arn]
            return target_group

    def delete_listener(self, listener_arn: str) -> FakeListener:
        for load_balancer in self.load_balancers.values():
            listener = load_balancer.listeners.pop(listener_arn, None)
            if listener:
                return listener
        raise ListenerNotFoundError()

    def modify_rule(
        self,
        rule_arn: str,
        conditions: List[Dict[str, Any]],
        actions: List[Dict[str, Any]],
    ) -> FakeRule:
        fake_actions = [FakeAction(action) for action in actions]
        # if conditions or actions is empty list, do not update the attributes
        if not conditions and not fake_actions:
            raise InvalidModifyRuleArgumentsError()
        rules = self.describe_rules(listener_arn=None, rule_arns=[rule_arn])
        if not rules:
            raise RuleNotFoundError()
        rule = rules[0]

        self._validate_conditions(conditions)
        # TODO: check pattern of value for 'host-header'
        # TODO: check pattern of value for 'path-pattern'

        # validate Actions
        self._validate_actions(fake_actions)

        # TODO: check for error 'TooManyRegistrationsForTargetId'
        # TODO: check for error 'TooManyRules'

        # modify rule
        if conditions:
            rule.conditions = conditions
        if actions:
            rule.actions = fake_actions
        return rule

    def register_targets(self, target_group_arn: str, instances: List[Any]) -> None:
        target_group = self.target_groups.get(target_group_arn)
        if target_group is None:
            raise TargetGroupNotFoundError()
        target_group.register(instances)

    def deregister_targets(
        self, target_group_arn: str, instances: List[Dict[str, Any]]
    ) -> None:
        target_group = self.target_groups.get(target_group_arn)
        if target_group is None:
            raise TargetGroupNotFoundError()
        target_group.deregister(instances)

    def describe_target_health(
        self, target_group_arn: str, targets: List[Dict[str, Any]]
    ) -> List[FakeHealthStatus]:
        target_group = self.target_groups.get(target_group_arn)
        if target_group is None:
            raise TargetGroupNotFoundError()

        if not targets:
            targets = list(target_group.targets.values())
        return [target_group.health_for(target, self.ec2_backend) for target in targets]

    def set_rule_priorities(
        self, rule_priorities: List[Dict[str, Any]]
    ) -> List[FakeRule]:
        # validate
        priorities = [rule_priority["priority"] for rule_priority in rule_priorities]
        for priority in set(priorities):
            if priorities.count(priority) > 1:
                raise DuplicatePriorityError(priority)

        # validate
        for rule_priority in rule_priorities:
            given_rule_arn = rule_priority["rule_arn"]
            priority = rule_priority["priority"]
            _given_rules = self.describe_rules(
                listener_arn=None, rule_arns=[given_rule_arn]
            )
            if not _given_rules:
                raise RuleNotFoundError()
            given_rule = _given_rules[0]
            listeners = self.describe_listeners(None, [given_rule.listener_arn])
            listener = listeners[0]
            for rule_in_listener in listener.rules.values():
                if rule_in_listener.priority == priority:
                    raise PriorityInUseError()
        # modify
        modified_rules = []
        for rule_priority in rule_priorities:
            given_rule_arn = rule_priority["rule_arn"]
            priority = rule_priority["priority"]
            _given_rules = self.describe_rules(
                listener_arn=None, rule_arns=[given_rule_arn]
            )
            if not _given_rules:
                raise RuleNotFoundError()
            given_rule = _given_rules[0]
            given_rule.priority = priority
            modified_rules.append(given_rule)
        return modified_rules

    def set_ip_address_type(self, arn: str, ip_type: str) -> None:
        if ip_type not in ("ipv4", "dualstack"):
            raise RESTError(
                "ValidationError",
                f"1 validation error detected: Value '{ip_type}' at 'ipAddressType' failed to satisfy constraint: Member must satisfy enum value set: [ipv4, dualstack]",
            )

        balancer = self.load_balancers.get(arn)
        if balancer is None:
            raise LoadBalancerNotFoundError()

        if ip_type == "dualstack" and balancer.scheme == "internal":
            raise RESTError(
                "InvalidConfigurationRequest",
                "Internal load balancers cannot be dualstack",
            )

        balancer.ip_address_type = ip_type

    def set_security_groups(self, arn: str, sec_groups: List[str]) -> None:
        balancer = self.load_balancers.get(arn)
        if balancer is None:
            raise LoadBalancerNotFoundError()

        # Check all security groups exist
        for sec_group_id in sec_groups:
            if self.ec2_backend.get_security_group_from_id(sec_group_id) is None:
                raise RESTError(
                    "InvalidSecurityGroup",
                    f"Security group {sec_group_id} does not exist",
                )

        balancer.security_groups = sec_groups

    def set_subnets(
        self, arn: str, subnets: List[str], subnet_mappings: List[Dict[str, Any]]
    ) -> list[AvailabilityZone]:
        balancer = self.load_balancers.get(arn)
        if balancer is None:
            raise LoadBalancerNotFoundError()

        subnet_objects = []
        sub_zone_list: Dict[str, str] = {}
        for subnet_id in subnets:
            try:
                subnet = self._get_subnet(sub_zone_list, subnet_id)

                sub_zone_list[subnet.availability_zone] = subnet.id
                subnet_objects.append(subnet)
            except Exception:
                raise SubnetNotFoundError()

        for subnet_mapping in subnet_mappings:
            subnet_id = subnet_mapping["SubnetId"]
            subnet = self._get_subnet(sub_zone_list, subnet_id)

            sub_zone_list[subnet.availability_zone] = subnet.id
            subnet_objects.append(subnet)

        if len(sub_zone_list) < 2:
            raise RESTError(
                "InvalidConfigurationRequest",
                "More than 1 availability zone must be specified",
            )

        balancer.subnets = subnet_objects

        return [AvailabilityZone(subnet) for subnet in subnet_objects]

    def _get_subnet(self, sub_zone_list: Dict[str, str], subnet_id: str) -> Subnet:
        subnet = self.ec2_backend.get_subnet(subnet_id)
        if subnet.availability_zone in sub_zone_list:
            raise RESTError(
                "InvalidConfigurationRequest",
                "More than 1 subnet cannot be specified for 1 availability zone",
            )
        return subnet

    def modify_load_balancer_attributes(
        self, arn: str, attrs: Dict[str, Any]
    ) -> Dict[str, Any]:
        balancer = self.load_balancers.get(arn)
        if balancer is None:
            raise LoadBalancerNotFoundError()

        for key in attrs:
            if key not in FakeLoadBalancer.VALID_ATTRS:
                raise RESTError("InvalidConfigurationRequest", f"Key {key} not valid")

        balancer.attrs.update(attrs)
        return balancer.attrs

    def describe_load_balancer_attributes(self, arn: str) -> Dict[str, Any]:
        balancer = self.load_balancers.get(arn)
        if balancer is None:
            raise LoadBalancerNotFoundError()

        return balancer.attrs

    def modify_target_group(
        self,
        arn: str,
        health_check_proto: Optional[str] = None,
        health_check_port: Optional[str] = None,
        health_check_path: Optional[str] = None,
        health_check_interval: Optional[int] = None,
        health_check_timeout: Optional[int] = None,
        healthy_threshold_count: Optional[str] = None,
        unhealthy_threshold_count: Optional[str] = None,
        http_codes: Optional[str] = None,
        health_check_enabled: Optional[bool] = None,
    ) -> FakeTargetGroup:
        target_group = self.target_groups.get(arn)
        if target_group is None:
            raise TargetGroupNotFoundError()

        if (
            http_codes is not None
            and FakeTargetGroup.HTTP_CODE_REGEX.match(http_codes) is None
        ):
            raise RESTError(
                "InvalidParameterValue",
                "HttpCode must be like 200 | 200-399 | 200,201 ...",
            )

        if http_codes is not None and target_group.protocol in ["HTTP", "HTTPS"]:
            target_group.matcher["HttpCode"] = http_codes
        if health_check_interval is not None:
            target_group.health_check_interval_seconds = health_check_interval
        if health_check_path is not None:
            target_group.health_check_path = health_check_path
        if health_check_port is not None:
            target_group.health_check_port = health_check_port
        if health_check_proto is not None:
            target_group.health_check_protocol = health_check_proto
        if health_check_timeout is not None:
            target_group.health_check_timeout_seconds = health_check_timeout
        if health_check_enabled is not None:
            target_group.health_check_enabled = health_check_enabled
        if healthy_threshold_count is not None:
            target_group.healthy_threshold_count = healthy_threshold_count
        if unhealthy_threshold_count is not None:
            target_group.unhealthy_threshold_count = unhealthy_threshold_count

        return target_group

    def modify_listener(
        self,
        arn: str,
        port: Optional[str] = None,
        protocol: Optional[str] = None,
        ssl_policy: Optional[str] = None,
        certificates: Optional[List[Dict[str, Any]]] = None,
        actions: Optional[List[Dict[str, Any]]] = None,
    ) -> FakeListener:
        default_actions = [FakeAction(action) for action in actions]  # type: ignore
        listener = self.describe_listeners(load_balancer_arn=None, listener_arns=[arn])[
            0
        ]

        if port is not None:
            listener.port = port

        if protocol not in (None, "HTTP", "HTTPS", "TCP"):
            raise RESTError(
                "UnsupportedProtocol", f"Protocol {protocol} is not supported"
            )

        # HTTPS checks
        protocol_becomes_https = protocol == "HTTPS"
        protocol_stays_https = protocol is None and listener.protocol == "HTTPS"
        if protocol_becomes_https or protocol_stays_https:
            # Check certificates exist
            if certificates:
                default_cert = certificates[0]
                default_cert_arn = default_cert["certificate_arn"]
                if not self._certificate_exists(certificate_arn=default_cert_arn):
                    raise RESTError(
                        "CertificateNotFound",
                        f"Certificate {default_cert_arn} not found",
                    )
                listener.certificate = default_cert_arn
                # TODO: Calling describe_listener_certificates after this operation returns a wrong result
                listener.certificates = [c["certificate_arn"] for c in certificates]
            elif len(certificates) == 0 and len(listener.certificates) == 0:  # type: ignore[arg-type]
                raise RESTError(
                    "CertificateWereNotPassed",
                    "You must provide a list containing exactly one certificate if the listener protocol is HTTPS.",
                )
            # else:
            # The list was not provided, meaning we just keep the existing certificates (if any)

        if protocol is not None:
            listener.protocol = protocol

        if ssl_policy is not None:
            # Its already validated in responses.py
            listener.ssl_policy = ssl_policy

        if default_actions is not None and default_actions != []:
            # Is currently not validated
            listener.default_actions = default_actions
            listener._default_rule[0].actions = default_actions

        return listener

    def _certificate_exists(self, certificate_arn: str) -> bool:
        """
        Verify the provided certificate exists in either ACM or IAM
        """
        from moto.acm.models import CertificateNotFound, acm_backends

        try:
            acm_backend = acm_backends[self.account_id][self.region_name]
            acm_backend.get_certificate(certificate_arn)
            return True
        except CertificateNotFound:
            pass

        from moto.iam import iam_backends

        cert = iam_backends[self.account_id][self.partition].get_certificate_by_arn(
            certificate_arn
        )
        if cert is not None:
            return True

        # ACM threw an error, and IAM did not return a certificate
        # Safe to assume it doesn't exist when we get here
        return False

    def _any_listener_using(self, target_group_arn: str) -> bool:
        for load_balancer in self.load_balancers.values():
            for listener in load_balancer.listeners.values():
                for rule in listener.rules.values():
                    for action in rule.actions:
                        found_arns = self._get_target_group_arns_from(
                            action_data=action.data
                        )
                        if target_group_arn in found_arns:
                            return True
        return False

    def notify_terminate_instances(self, instance_ids: List[str]) -> None:
        for target_group in self.target_groups.values():
            target_group.deregister_terminated_instances(instance_ids)

    def add_listener_certificates(
        self, arn: str, certificates: List[Dict[str, Any]]
    ) -> List[str]:
        listener = self.describe_listeners(load_balancer_arn=None, listener_arns=[arn])[
            0
        ]
        # https://docs.aws.amazon.com/elasticloadbalancing/latest/application/load-balancer-limits.html
        if len(certificates) + len(listener.certificates) > 25:
            raise TooManyCertificatesError()
        listener.certificates.extend([c["certificate_arn"] for c in certificates])
        return listener.certificates

    def describe_listener_certificates(self, arn: str) -> List[str]:
        listener = self.describe_listeners(load_balancer_arn=None, listener_arns=[arn])[
            0
        ]
        return listener.certificates

    def remove_listener_certificates(
        self, arn: str, certificates: List[Dict[str, Any]]
    ) -> None:
        listener = self.describe_listeners(load_balancer_arn=None, listener_arns=[arn])[
            0
        ]
        cert_arns = [c["certificate_arn"] for c in certificates]
        listener.certificates = [c for c in listener.certificates if c not in cert_arns]

    def add_tags(self, resource_arns: List[str], tags: List[Dict[str, str]]) -> None:
        tag_dict = self.tagging_service.flatten_tag_list(tags)
        for arn in resource_arns:
            existing = self.tagging_service.get_tag_dict_for_resource(arn)
            for key in tag_dict:
                if len(existing) >= 10 and key not in existing:
                    raise TooManyTagsError()
            self._get_resource_by_arn(arn)
            self.tagging_service.tag_resource(arn, tags)

    def remove_tags(self, resource_arns: List[str], tag_keys: List[str]) -> None:
        for arn in resource_arns:
            self._get_resource_by_arn(arn)
            self.tagging_service.untag_resource_using_names(arn, tag_keys)

    def describe_tags(self, resource_arns: List[str]) -> Dict[str, Dict[str, str]]:
        return {
            arn: self.tagging_service.get_tag_dict_for_resource(arn)
            for arn in resource_arns
        }

    def describe_listener_attributes(self, listener_arn: str) -> Dict[str, str]:
        listener = self.describe_listeners(
            load_balancer_arn=None, listener_arns=[listener_arn]
        )[0]
        return listener.attributes

    def modify_listener_attributes(
        self, listener_arn: str, attrs: List[Dict[str, str]]
    ) -> Dict[str, str]:
        listener = self.describe_listeners(
            load_balancer_arn=None, listener_arns=[listener_arn]
        )[0]
        listener.attributes.update({a["Key"]: a["Value"] for a in attrs})
        return listener.attributes

    def _get_resource_by_arn(self, arn: str) -> Any:
        if ":targetgroup" in arn:
            resource: Any = self.target_groups.get(arn)
            if not resource:
                raise TargetGroupNotFoundError()
        elif ":loadbalancer" in arn:
            resource = self.load_balancers.get(arn)
            if not resource:
                raise LoadBalancerNotFoundError()
        elif ":listener-rule" in arn:
            lb_arn = arn.replace(":listener-rule", ":loadbalancer").rsplit("/", 2)[0]
            balancer = self.load_balancers.get(lb_arn)
            if not balancer:
                raise LoadBalancerNotFoundError()
            listener_arn = arn.replace(":listener-rule", ":listener").rsplit("/", 1)[0]
            listener = balancer.listeners.get(listener_arn)
            if not listener:
                raise ListenerNotFoundError()
            resource = listener.rules.get(arn)
            if not resource:
                raise RuleNotFoundError()
        elif ":listener" in arn:
            lb_arn, _, _ = arn.replace(":listener", ":loadbalancer").rpartition("/")
            balancer = self.load_balancers.get(lb_arn)
            if not balancer:
                raise LoadBalancerNotFoundError()
            resource = balancer.listeners.get(arn)
            if not resource:
                raise ListenerNotFoundError()
        else:
            raise LoadBalancerNotFoundError()
        return resource


elbv2_backends = BackendDict(ELBv2Backend, "elbv2")
