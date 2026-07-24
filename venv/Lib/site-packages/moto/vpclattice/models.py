import random
import uuid
from collections.abc import Iterator
from datetime import datetime, timezone
from typing import Any

from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.common_models import BaseModel
from moto.core.resource_tagging import TaggableResourcesMixin, TaggedResource
from moto.utilities.paginator import paginate
from moto.utilities.tagging_service import TaggingService
from moto.vpclattice.exceptions import (
    ResourceNotFoundException,
    ValidationException,
)


class VPCLatticeService(BaseModel):
    def __init__(
        self,
        region: str,
        account_id: str,
        auth_type: str,
        certificate_arn: str | None,
        client_token: str,
        custom_domain_name: str | None,
        name: str,
        tags: dict[str, str] | None,
    ) -> None:
        self.id: str = f"svc-{str(uuid.uuid4())[:17]}"
        self.auth_type: str = auth_type
        self.certificate_arn: str = certificate_arn or ""
        self.client_token: str = client_token
        self.custom_domain_name: str = custom_domain_name or ""
        self.dns_entry: VPCLatticeDNSEntry = VPCLatticeDNSEntry(
            region, self.id, self.custom_domain_name
        )
        self.name: str = name
        self.arn: str = f"arn:aws:vpc-lattice:{region}:{account_id}:service/{self.id}"
        self.status: str = "ACTIVE"
        self.tags: dict[str, str] = tags or {}

    def to_dict(self) -> dict[str, Any]:
        return {
            "arn": self.arn,
            "authType": self.auth_type,
            "certificateArn": self.certificate_arn,
            "customDomainName": self.custom_domain_name,
            "dnsEntry": self.dns_entry.to_dict(),
            "id": self.id,
            "name": self.name,
            "status": self.status,
        }


class VPCLatticeServiceNetwork(BaseModel):
    def __init__(
        self,
        region: str,
        account_id: str,
        auth_type: str,
        client_token: str,
        name: str,
        sharing_config: dict[str, Any] | None,
        tags: dict[str, str] | None,
    ) -> None:
        self.auth_type: str = auth_type
        self.client_token: str = client_token
        self.id: str = f"sn-{str(uuid.uuid4())[:17]}"
        self.name: str = name
        self.arn: str = (
            f"arn:aws:vpc-lattice:{region}:{account_id}:servicenetwork/{self.id}"
        )
        self.sharing_config: dict[str, Any] = sharing_config or {}
        self.tags: dict[str, str] = tags or {}

    def to_dict(self) -> dict[str, Any]:
        return {
            "arn": self.arn,
            "authType": self.auth_type,
            "id": self.id,
            "name": self.name,
            "sharingConfig": self.sharing_config,
        }


class VPCLatticeListener(BaseModel):
    def __init__(
        self,
        region: str,
        account_id: str,
        service_identifier: str,
        name: str,
        protocol: str,
        port: int | None,
        default_action: dict[str, Any],
        client_token: str,
        tags: dict[str, str] | None,
    ) -> None:
        self.id: str = f"listener-{str(uuid.uuid4())[:17]}"
        self.arn: str = f"arn:aws:vpc-lattice:{region}:{account_id}:service/{service_identifier}/listener/{self.id}"
        self.service_identifier: str = service_identifier
        self.service_arn: str = (
            f"arn:aws:vpc-lattice:{region}:{account_id}:service/{service_identifier}"
        )
        self.name: str = name
        self.protocol: str = protocol
        self.port: int = port or 443
        self.default_action: dict[str, Any] = default_action
        self.client_token: str = client_token
        self.createdAt = datetime.now(timezone.utc).isoformat()
        self.last_updated_at = datetime.now(timezone.utc).isoformat()
        self.tags: dict[str, str] = tags or {}

    def to_dict(self) -> dict[str, Any]:
        return {
            "id": self.id,
            "arn": self.arn,
            "name": self.name,
            "protocol": self.protocol,
            "port": self.port,
            "serviceArn": self.service_arn,
            "serviceId": self.service_identifier,
            "defaultAction": self.default_action,
            "createdAt": self.createdAt,
            "lastUpdatedAt": self.last_updated_at,
        }


class VPCLatticeServiceNetworkResourceAssociation(BaseModel):
    def __init__(
        self,
        region: str,
        account_id: str,
        service_network: VPCLatticeServiceNetwork,
        # resource_configuration: VPCLatticeResourceConfiguration,
        private_dns_enabled: bool,
        client_token: str,
        tags: dict[str, str] | None,
    ) -> None:
        now_iso = datetime.now(timezone.utc).isoformat()
        self.id: str = f"snra-{str(uuid.uuid4())[:17]}"
        self.arn: str = (
            f"arn:aws:vpc-lattice:{region}:{account_id}:"
            f"servicenetworkresourceassociation/{self.id}"
        )
        self.service_network_id: str = service_network.id
        self.service_network_name: str = service_network.name
        self.service_network_arn: str = service_network.arn
        # self.resource_configuration_id: str = resource_configuration.id
        # self.resource_configuration_arn: str = resource_configuration.arn
        # self.resource_configuration_name: str = resource_configuration.name
        self.created_by: str = account_id
        self.created_at: str = now_iso
        self.status: str = "ACTIVE"
        self.client_token: str = client_token
        self.last_updated_at: str = now_iso
        self.private_dns_enabled: bool = private_dns_enabled
        self.tags: dict[str, str] = tags or {}

    def to_dict(self) -> dict[str, Any]:
        return {
            "id": self.id,
            "arn": self.arn,
            "status": self.status,
            "createdBy": self.created_by,
            "createdAt": self.created_at,
            # "resourceConfigurationId": self.resource_configuration_id,
            # "resourceConfigurationArn": self.resource_configuration_arn,
            # "resourceConfigurationName": self.resource_configuration_name,
            "serviceNetworkId": self.service_network_id,
            "serviceNetworkName": self.service_network_name,
            "serviceNetworkArn": self.service_network_arn,
            "lastUpdatedAt": self.last_updated_at,
            "privateDnsEnabled": self.private_dns_enabled,
        }


class VPCLatticeServiceNetworkServiceAssociation(BaseModel):
    def __init__(
        self,
        region: str,
        account_id: str,
        service_network_identifier: str,
        service_identifier: str,
        client_token: str,
        tags: dict[str, str] | None,
    ) -> None:
        now_iso = datetime.now(timezone.utc).isoformat()
        self.id: str = f"snsa-{str(uuid.uuid4())[:17]}"
        self.arn: str = f"arn:aws:vpc-lattice:{region}:{account_id}:servicenetworkserviceassociation/{self.id}"
        self.service_network_id: str = service_network_identifier
        self.service_id: str = service_identifier
        self.created_by: str = account_id
        self.created_at: str = now_iso
        self.last_updated_at: str = now_iso
        self.status: str = "ACTIVE"
        self.client_token: str = client_token
        self.dns_entry: VPCLatticeDNSEntry = VPCLatticeDNSEntry(region, self.id)
        self.custom_domain_name: str = self.dns_entry.domain_name
        self.tags: dict[str, str] = tags or {}

    def to_dict(self) -> dict[str, Any]:
        return {
            "id": self.id,
            "arn": self.arn,
            "serviceNetworkId": self.service_network_id,
            "serviceId": self.service_id,
            "createdBy": self.created_by,
            "createdAt": self.created_at,
            "lastUpdatedAt": self.last_updated_at,
            "customDomainName": self.custom_domain_name,
            "dnsEntry": self.dns_entry.to_dict(),
            "status": self.status,
        }


class VPCLatticeServiceNetworkVpcAssociation(BaseModel):
    def __init__(
        self,
        region: str,
        account_id: str,
        client_token: str,
        security_group_ids: list[str] | None,
        service_network: VPCLatticeServiceNetwork,
        tags: dict[str, str] | None,
        vpc_identifier: str,
    ) -> None:
        now_iso = datetime.now(timezone.utc).isoformat()
        self.id: str = f"snva-{str(uuid.uuid4())[:17]}"
        self.arn: str = f"arn:aws:vpc-lattice:{region}:{account_id}:servicenetworkvpcassociation/{self.id}"
        self.created_by: str = "user"
        self.created_at: str = now_iso
        self.last_updated_at: str = now_iso
        self.security_group_ids: list[str] = security_group_ids or []
        self.service_network_id: str = service_network.id
        self.service_network_name: str = service_network.name
        self.service_network_arn: str = service_network.arn
        self.status: str = "ACTIVE"
        self.vpc_id: str = vpc_identifier
        self.tags: dict[str, str] = tags or {}

    def to_dict(self) -> dict[str, Any]:
        return {
            "id": self.id,
            "status": self.status,
            "arn": self.arn,
            "createdBy": self.created_by,
            "createdAt": self.created_at,
            "serviceNetworkId": self.service_network_id,
            "serviceNetworkName": self.service_network_name,
            "serviceNetworkArn": self.service_network_arn,
            "vpcId": self.vpc_id,
            "lastUpdatedAt": self.last_updated_at,
            "securityGroupIds": self.security_group_ids,
        }


class VPCLatticeRule(BaseModel):
    def __init__(
        self,
        region: str,
        account_id: str,
        action: dict[str, Any],
        client_token: str,
        listener_identifier: str,
        match: dict[str, Any],
        name: str,
        priority: int,
        service_identifier: str,
        tags: dict[str, str],
    ) -> None:
        self.action: dict[str, Any] = action or {}
        self.id: str = f"rule-{str(uuid.uuid4())[:17]}"
        self.arn: str = (
            f"arn:aws:vpc-lattice:{region}:{account_id}:service/{service_identifier}"
            f"/listener/listener-{listener_identifier}/rule/{self.id}"
        )
        self.client_token: str = client_token
        self.listener_identifier: str = listener_identifier
        self.match: dict[str, Any] = match or {}
        self.name: str = name
        self.priority: int = priority
        self.service_identifier: str = service_identifier
        self.tags: dict[str, str] = tags or {}

    def to_dict(self) -> dict[str, Any]:
        return {
            "arn": self.arn,
            "id": self.id,
            "name": self.name,
            "priority": self.priority,
            "action": self.action,
            "match": self.match,
            "serviceIdentifier": self.service_identifier,
            "listenerIdentifier": self.listener_identifier,
            "tags": self.tags,
        }


class VPCLatticeTargetGroup(BaseModel):
    def __init__(
        self,
        region: str,
        account_id: str,
        name: str,
        target_type: str,
        config: dict[str, Any] | None,
        client_token: str,
        tags: dict[str, str] | None,
    ) -> None:
        self.id: str = f"tg-{str(uuid.uuid4())[:17]}"
        self.arn: str = (
            f"arn:aws:vpc-lattice:{region}:{account_id}:targetgroup/{self.id}"
        )
        self.name: str = name
        self.type: str = target_type
        self.config: dict[str, Any] = config or {}
        self.client_token: str = client_token
        self.tags: dict[str, str] = tags or {}
        self.status: str = "CREATE_IN_PROGRESS"
        self.vpcIdentifer: str | None = self.config.get("vpcIdentifier")
        self.createdAt = datetime.now(timezone.utc).isoformat()
        self.last_updated_at = datetime.now(timezone.utc).isoformat()

    def to_dict(self) -> dict[str, Any]:
        return {
            "id": self.id,
            "arn": self.arn,
            "name": self.name,
            "type": self.type,
            "config": self.config,
            "createdAt": self.createdAt,
            "vpcIdentifier": self.vpcIdentifer,
            "lastUpdatedAt": self.last_updated_at,
            "status": self.status,
        }


class VPCLatticeDNSEntry:
    def __init__(
        self,
        region_name: str,
        service_id: str,
        custom_domain_name: str | None = None,
    ) -> None:
        self.domain_name: str = (
            custom_domain_name or f"{service_id}.{region_name}.vpclattice.amazonaws.com"
        )
        self.hosted_zone_id: str = f"Z{random.randint(100000, 999999)}XYZ"

    def to_dict(self) -> dict[str, str]:
        return {"domainName": self.domain_name, "hostedZoneId": self.hosted_zone_id}


class VPCLatticeAccessLogSubscription(BaseModel):
    def __init__(
        self,
        region: str,
        account_id: str,
        destinationArn: str,
        resourceArn: str,
        resourceId: str,  # resourceIdentifier
        serviceNetworkLogType: str | None,
        tags: dict[str, str] | None,
    ) -> None:
        self.id: str = f"als-{str(uuid.uuid4())[:17]}"
        self.arn: str = (
            f"arn:aws:vpc-lattice:{region}:{account_id}:accesslogsubscription/{self.id}"
        )
        self.created_at = datetime.now(timezone.utc).isoformat()
        self.destinationArn = destinationArn
        self.last_updated_at = datetime.now(timezone.utc).isoformat()
        self.resourceArn = resourceArn
        self.resourceId = resourceId
        self.serviceNetworkLogType = serviceNetworkLogType or "SERVICE"
        self.tags = tags or {}

    def to_dict(self) -> dict[str, Any]:
        return {
            "arn": self.arn,
            "createdAt": self.created_at,
            "destinationArn": self.destinationArn,
            "id": self.id,
            "lastUpdatedAt": self.last_updated_at,
            "resourceArn": self.resourceArn,
            "resourceId": self.resourceId,
            "serviceNetworkLogType": self.serviceNetworkLogType,
            "tags": self.tags,
        }


class AuthPolicy:
    def __init__(
        self,
        policy_name: str,
        policy: str,
        created_at: str,
        last_updated_at: str,
        state: str,
    ) -> None:
        self.policy_name = policy_name
        self.policy = policy
        self.created_at = created_at
        self.last_updated_at = last_updated_at
        self.state = state


class VPCLatticeBackend(BaseBackend, TaggableResourcesMixin):
    SERVICE_NAMESPACE = "vpc-lattice"

    PAGINATION_MODEL = {
        "list_services": {
            "input_token": "next_token",
            "limit_key": "max_results",
            "limit_default": 50,
            "unique_attribute": "id",
        },
        "list_service_networks": {
            "input_token": "next_token",
            "limit_key": "max_results",
            "limit_default": 50,
            "unique_attribute": "id",
        },
        "list_service_network_vpc_associations": {
            "input_token": "next_token",
            "limit_key": "max_results",
            "limit_default": 50,
            "unique_attribute": "id",
        },
        "list_listeners": {
            "input_token": "next_token",
            "limit_key": "max_results",
            "limit_default": 50,
            "unique_attribute": "id",
        },
        "list_target_groups": {
            "input_token": "next_token",
            "limit_key": "max_results",
            "limit_default": 50,
            "unique_attribute": "id",
        },
        "list_service_network_resource_associations": {
            "input_token": "next_token",
            "limit_key": "max_results",
            "limit_default": 50,
            "unique_attribute": "id",
        },
        "list_service_network_service_associations": {
            "input_token": "next_token",
            "limit_key": "max_results",
            "limit_default": 50,
            "unique_attribute": "id",
        },
    }

    def __init__(self, region_name: str, account_id: str) -> None:
        super().__init__(region_name, account_id)
        self.services: dict[str, VPCLatticeService] = {}
        self.listeners: dict[str, VPCLatticeListener] = {}
        self.service_networks: dict[str, VPCLatticeServiceNetwork] = {}
        self.service_network_resource_associations: dict[
            str, VPCLatticeServiceNetworkResourceAssociation
        ] = {}
        self.service_network_service_associations: dict[
            str, VPCLatticeServiceNetworkServiceAssociation
        ] = {}
        self.service_network_vpc_associations: dict[
            str, VPCLatticeServiceNetworkVpcAssociation
        ] = {}
        self.rules: dict[str, VPCLatticeRule] = {}
        self.tagger: TaggingService = TaggingService()
        self.target_groups: dict[str, VPCLatticeTargetGroup] = {}
        self.access_log_subscriptions: dict[str, VPCLatticeAccessLogSubscription] = {}
        self.resource_policies: dict[str, str] = {}
        self.auth_policies: dict[str, AuthPolicy] = {}

    def create_service(
        self,
        auth_type: str,
        certificate_arn: str | None,
        client_token: str,
        custom_domain_name: str | None,
        name: str,
        tags: dict[str, str] | None,
    ) -> VPCLatticeService:
        service = VPCLatticeService(
            self.region_name,
            self.account_id,
            auth_type,
            certificate_arn,
            client_token,
            custom_domain_name,
            name,
            tags,
        )
        self.services[service.id] = service
        self.tag_resource(service.arn, tags or {})
        return service

    def get_service(self, service_identifier: str) -> VPCLatticeService:
        service = self.services.get(service_identifier)
        if not service:
            raise ResourceNotFoundException(service_identifier)
        return service

    @paginate(pagination_model=PAGINATION_MODEL)
    def list_services(self) -> list[VPCLatticeService]:
        return list(self.services.values())

    def create_service_network(
        self,
        auth_type: str,
        client_token: str,
        name: str,
        sharing_config: dict[str, Any] | None,
        tags: dict[str, str] | None,
    ) -> VPCLatticeServiceNetwork:
        """
        WARNING: This method currently does NOT fail if there is a disassociation in progress.
        """
        sn = VPCLatticeServiceNetwork(
            self.region_name,
            self.account_id,
            auth_type,
            client_token,
            name,
            sharing_config,
            tags,
        )
        self.service_networks[sn.id] = sn
        self.tag_resource(sn.arn, tags or {})
        return sn

    def get_service_network(
        self, service_network_identifier: str
    ) -> VPCLatticeServiceNetwork:
        service_network = self.service_networks.get(service_network_identifier)
        if not service_network:
            raise ResourceNotFoundException(service_network_identifier)
        return service_network

    @paginate(pagination_model=PAGINATION_MODEL)
    def list_service_networks(self) -> list[VPCLatticeServiceNetwork]:
        return list(self.service_networks.values())

    def create_service_network_vpc_association(
        self,
        client_token: str,
        security_group_ids: list[str] | None,
        service_network_identifier: str,
        tags: dict[str, str] | None,
        vpc_identifier: str,
    ) -> VPCLatticeServiceNetworkVpcAssociation:
        sn = self.get_service_network(service_network_identifier)
        assoc = VPCLatticeServiceNetworkVpcAssociation(
            self.region_name,
            self.account_id,
            client_token,
            security_group_ids,
            sn,
            tags,
            vpc_identifier,
        )
        self.service_network_vpc_associations[assoc.id] = assoc
        self.tag_resource(assoc.arn, tags or {})
        return assoc

    def create_rule(
        self,
        action: dict[str, Any],
        client_token: str,
        listener_identifier: str,
        match: dict[str, Any],
        name: str,
        priority: int,
        service_identifier: str,
        tags: dict[str, str],
    ) -> VPCLatticeRule:
        rule = VPCLatticeRule(
            self.region_name,
            self.account_id,
            action,
            client_token,
            listener_identifier,
            match,
            name,
            priority,
            service_identifier,
            tags,
        )
        self.rules[rule.id] = rule
        self.tag_resource(rule.arn, tags or {})
        return rule

    def list_tags_for_resource(self, resource_arn: str) -> dict[str, str]:
        return self.tagger.get_tag_dict_for_resource(resource_arn)

    def create_access_log_subscription(
        self,
        resourceIdentifier: str,
        destinationArn: str,
        client_token: str | None,
        serviceNetworkLogType: str | None,
        tags: dict[str, str] | None,
    ) -> VPCLatticeAccessLogSubscription:
        resource: Any = None
        if resourceIdentifier.startswith("sn-"):
            resource = self.service_networks.get(resourceIdentifier)
        elif resourceIdentifier.startswith("svc-"):
            resource = self.services.get(resourceIdentifier)
        else:
            raise ValidationException(
                "Invalid parameter resourceIdentifier, must start with 'sn-' or 'svc-'"
            )

        if not resource:
            raise ResourceNotFoundException(f"Resource {resourceIdentifier} not found")

        sub = VPCLatticeAccessLogSubscription(
            self.region_name,
            self.account_id,
            destinationArn,
            resource.arn,
            resource.id,
            serviceNetworkLogType,
            tags,
        )

        self.access_log_subscriptions[sub.id] = sub
        self.tag_resource(sub.arn, tags or {})
        return sub

    def get_access_log_subscription(
        self, accessLogSubscriptionIdentifier: str
    ) -> VPCLatticeAccessLogSubscription:
        sub = self.access_log_subscriptions.get(accessLogSubscriptionIdentifier)
        if not sub:
            raise ResourceNotFoundException(
                f"Access Log Subscription {accessLogSubscriptionIdentifier} not found"
            )
        return sub

    def list_access_log_subscriptions(
        self,
        resourceIdentifier: str,
        maxResults: int | None = None,
        nextToken: str | None = None,
    ) -> list[VPCLatticeAccessLogSubscription]:
        return [
            sub
            for sub in self.access_log_subscriptions.values()
            if sub.resourceId == resourceIdentifier
            or sub.resourceArn == resourceIdentifier
        ][:maxResults]

    def update_access_log_subscription(
        self,
        accessLogSubscriptionIdentifier: str,
        destinationArn: str,
    ) -> VPCLatticeAccessLogSubscription:
        sub = self.access_log_subscriptions.get(accessLogSubscriptionIdentifier)
        if not sub:
            raise ResourceNotFoundException(
                f"Access Log Subscription {accessLogSubscriptionIdentifier} not found"
            )

        sub.destinationArn = destinationArn
        sub.last_updated_at = datetime.now(timezone.utc).isoformat()

        return sub

    def delete_access_log_subscription(
        self, accessLogSubscriptionIdentifier: str
    ) -> None:
        sub = self.access_log_subscriptions.get(accessLogSubscriptionIdentifier)
        if not sub:
            raise ResourceNotFoundException(
                f"Access Log Subscription {accessLogSubscriptionIdentifier} not found"
            )
        del self.access_log_subscriptions[accessLogSubscriptionIdentifier]

    def put_auth_policy(self, resourceIdentifier: str, policy: str) -> AuthPolicy:
        # ResourceConfig not supported yet, only handle Service and Service Network
        if resourceIdentifier.startswith("arn:"):
            resourceIdentifier = resourceIdentifier.rsplit("/", 1)[-1]

        resource: Any = None

        if resourceIdentifier.startswith("sn-"):
            resource = self.service_networks.get(resourceIdentifier)
        elif resourceIdentifier.startswith("svc-"):
            resource = self.services.get(resourceIdentifier)
        else:
            raise ValidationException(
                "Invalid parameter resourceIdentifier, must start with 'sn-' or 'svc-'"
            )

        if not resource:
            raise ResourceNotFoundException(f"Resource {resourceIdentifier} not found")

        # Handle state management
        state = "Inactive" if resource.auth_type == "NONE" else "Active"

        now_iso = datetime.now(timezone.utc).isoformat()
        if resourceIdentifier in self.auth_policies:
            auth_policy = self.auth_policies[resourceIdentifier]
            auth_policy.policy = policy
            auth_policy.last_updated_at = now_iso
            auth_policy.state = state

        else:
            auth_policy = AuthPolicy(
                policy_name=resourceIdentifier,
                policy=policy,
                created_at=now_iso,
                last_updated_at=now_iso,
                state=state,
            )

            self.auth_policies[resourceIdentifier] = auth_policy

        return auth_policy

    def get_auth_policy(self, resourceIdentifier: str) -> AuthPolicy:
        original_identifier = resourceIdentifier
        if resourceIdentifier.startswith("arn:"):
            resourceIdentifier = resourceIdentifier.rsplit("/", 1)[-1]

        auth_policy = self.auth_policies.get(resourceIdentifier)
        if not auth_policy:
            raise ResourceNotFoundException(f"Resource {original_identifier} not found")

        resource: Any
        if resourceIdentifier.startswith("sn-"):
            resource = self.service_networks.get(resourceIdentifier)
        else:
            resource = self.services.get(resourceIdentifier)

        auth_policy.state = (
            "Inactive" if not resource or resource.auth_type == "NONE" else "Active"
        )

        return auth_policy

    def delete_auth_policy(self, resourceIdentifier: str) -> None:
        original_identifier = resourceIdentifier

        if resourceIdentifier.startswith("arn:"):
            resourceIdentifier = resourceIdentifier.rsplit("/", 1)[-1]

        if resourceIdentifier not in self.auth_policies:
            raise ResourceNotFoundException(f"Resource {original_identifier} not found")

        del self.auth_policies[resourceIdentifier]

    def put_resource_policy(self, resourceArn: str, policy: str) -> None:
        self.resource_policies[resourceArn] = policy

    def get_resource_policy(self, resourceArn: str) -> str:
        if resourceArn not in self.resource_policies:
            raise ResourceNotFoundException(f"Resource {resourceArn} not found")

        return self.resource_policies[resourceArn]

    def delete_resource_policy(self, resourceArn: str) -> None:
        if resourceArn not in self.resource_policies:
            raise ResourceNotFoundException(f"Resource {resourceArn} not found")

        del self.resource_policies[resourceArn]

    @paginate(pagination_model=PAGINATION_MODEL)
    def list_service_network_vpc_associations(
        self,
        serviceNetworkIdentifier: str | None = None,
        vpcIdentifier: str | None = None,
    ) -> list[VPCLatticeServiceNetworkVpcAssociation]:
        associations = list(self.service_network_vpc_associations.values())
        if serviceNetworkIdentifier:
            associations = [
                assoc
                for assoc in associations
                if assoc.service_network_id == serviceNetworkIdentifier
                or assoc.service_network_arn == serviceNetworkIdentifier
            ]
        if vpcIdentifier:
            associations = [
                assoc for assoc in associations if assoc.vpc_id == vpcIdentifier
            ]
        return associations

    def create_listener(
        self,
        service_identifier: str,
        name: str,
        protocol: str,
        port: int | None,
        default_action: dict[str, Any],
        client_token: str,
        tags: dict[str, str] | None,
    ) -> VPCLatticeListener:
        service = self.get_service(service_identifier)

        listener = VPCLatticeListener(
            self.region_name,
            self.account_id,
            service.id,
            name,
            protocol,
            port,
            default_action,
            client_token,
            tags,
        )
        self.listeners[listener.id] = listener
        self.tag_resource(listener.arn, tags or {})
        return listener

    def get_listener(
        self, service_identifier: str, listener_identifier: str
    ) -> VPCLatticeListener:
        listener = self.listeners.get(listener_identifier)
        if not listener:
            raise ResourceNotFoundException(f"Listener {listener_identifier} not found")

        service = self.get_service(service_identifier)
        if listener.service_identifier != service.id:
            raise ResourceNotFoundException(
                f"Listener {listener_identifier} not found in service {service_identifier}"
            )

        return listener

    @paginate(pagination_model=PAGINATION_MODEL)
    def list_listeners(
        self,
        service_identifier: str,
    ) -> list[VPCLatticeListener]:
        service = self.get_service(service_identifier)
        return [
            listener
            for listener in self.listeners.values()
            if listener.service_identifier == service.id
        ]

    def create_target_group(
        self,
        name: str,
        target_type: str,
        config: dict[str, Any] | None,
        client_token: str,
        tags: dict[str, str] | None,
    ) -> VPCLatticeTargetGroup:
        target_group = VPCLatticeTargetGroup(
            self.region_name,
            self.account_id,
            name,
            target_type,
            config,
            client_token,
            tags,
        )
        self.target_groups[target_group.id] = target_group
        self.tag_resource(target_group.arn, tags or {})
        target_group.status = "ACTIVE"
        return target_group

    def get_target_group(self, target_group_identifier: str) -> VPCLatticeTargetGroup:
        target_group = self.target_groups.get(target_group_identifier)
        if not target_group:
            raise ResourceNotFoundException(
                f"Target group {target_group_identifier} not found"
            )
        return target_group

    @paginate(pagination_model=PAGINATION_MODEL)
    def list_target_groups(
        self,
        vpc_identifier: str | None = None,
        target_group_type: str | None = None,
    ) -> list[VPCLatticeTargetGroup]:
        target_groups = list(self.target_groups.values())

        if vpc_identifier:
            target_groups = [
                tg
                for tg in target_groups
                if tg.config.get("vpcIdentifier") == vpc_identifier
            ]

        if target_group_type:
            target_groups = [tg for tg in target_groups if tg.type == target_group_type]

        return target_groups

    def create_service_network_resource_association(
        self,
        resource_config_identifier: str,
        service_network_identifier: str,
        private_dns_enabled: bool,
        client_token: str,
        tags: dict[str, str] | None,
    ) -> VPCLatticeServiceNetworkResourceAssociation:
        service_network = self.get_service_network(service_network_identifier)
        # TODO: Once resource configs are implemented, do something with the resource_config_id

        assoc = VPCLatticeServiceNetworkResourceAssociation(
            self.region_name,
            self.account_id,
            service_network,
            # resource_config,
            private_dns_enabled,
            client_token,
            tags,
        )
        self.service_network_resource_associations[assoc.id] = assoc
        self.tag_resource(assoc.arn, tags or {})
        return assoc

    def get_service_network_resource_association(
        self, service_network_resource_association_identifier: str
    ) -> VPCLatticeServiceNetworkResourceAssociation:
        assoc = self.service_network_resource_associations.get(
            service_network_resource_association_identifier
        )
        if not assoc:
            raise ResourceNotFoundException(
                f"Service network resource association {service_network_resource_association_identifier} not found"
            )
        return assoc

    @paginate(pagination_model=PAGINATION_MODEL)
    def list_service_network_resource_associations(
        self,
        service_network_identifier: str | None = None,
        resource_configuration_identifier: str | None = None,
    ) -> list[VPCLatticeServiceNetworkResourceAssociation]:
        associations = list(self.service_network_resource_associations.values())
        if service_network_identifier:
            associations = [
                assoc
                for assoc in associations
                if assoc.service_network_id == service_network_identifier
            ]

        # TODO: Once resource configs are implemented, add this filter block
        # if resource_configuration_identifier:
        #     associations = [
        #         assoc
        #         for assoc in associations
        #         if assoc.resource_id == resource_configuration_identifier
        #     ]

        return associations

    def create_service_network_service_association(
        self,
        service_network_identifier: str,
        service_identifier: str,
        client_token: str,
        tags: dict[str, str] | None,
    ) -> VPCLatticeServiceNetworkServiceAssociation:
        service_network = self.get_service_network(service_network_identifier)
        service = self.get_service(service_identifier)

        assoc = VPCLatticeServiceNetworkServiceAssociation(
            self.region_name,
            self.account_id,
            service_network.id,
            service.id,
            client_token,
            tags,
        )
        self.service_network_service_associations[assoc.id] = assoc
        self.tag_resource(assoc.arn, tags or {})
        return assoc

    def get_service_network_service_association(
        self, service_network_service_association_identifier: str
    ) -> VPCLatticeServiceNetworkServiceAssociation:
        assoc = self.service_network_service_associations.get(
            service_network_service_association_identifier
        )
        if not assoc:
            raise ResourceNotFoundException(
                f"Service network service association {service_network_service_association_identifier} not found"
            )
        return assoc

    @paginate(pagination_model=PAGINATION_MODEL)
    def list_service_network_service_associations(
        self,
        service_network_identifier: str | None = None,
        service_identifier: str | None = None,
    ) -> list[VPCLatticeServiceNetworkServiceAssociation]:
        associations = list(self.service_network_service_associations.values())
        if service_network_identifier:
            associations = [
                assoc
                for assoc in associations
                if assoc.service_network_id == service_network_identifier
            ]

        if service_identifier:
            associations = [
                assoc
                for assoc in associations
                if assoc.service_id == service_identifier
            ]
        return associations

    # Resource Groups Tagging API (TaggableResourcesMixin method overrides)
    def iter_tagged_resources(self) -> Iterator[TaggedResource]:
        sources: dict[str, dict[str, Any]] = {
            "vpc-lattice:accesslogsubscription": self.access_log_subscriptions,
            "vpc-lattice:listener": self.listeners,
            "vpc-lattice:rule": self.rules,
            "vpc-lattice:service": self.services,
            "vpc-lattice:servicenetwork": self.service_networks,
            "vpc-lattice:servicenetworkresourceassociation": self.service_network_resource_associations,
            "vpc-lattice:servicenetworkserviceassociation": self.service_network_service_associations,
            "vpc-lattice:servicenetworkvpcassociation": self.service_network_vpc_associations,
            "vpc-lattice:targetgroup": self.target_groups,
        }
        for resource_type, items in sources.items():
            for item in items.values():
                yield TaggedResource(
                    arn=item.arn,
                    tags=self.tagger.get_tag_dict_for_resource(item.arn),
                    resource_type=resource_type,
                )

    def tag_resource(self, resource_arn: str, tags: dict[str, str]) -> None:
        tags_input = self.tagger.convert_dict_to_tags_input(tags or {})
        self.tagger.tag_resource(resource_arn, tags_input)

    def untag_resource(self, resource_arn: str, tag_keys: list[str]) -> None:
        self.tagger.untag_resource_using_names(resource_arn, tag_keys)


vpclattice_backends: BackendDict[VPCLatticeBackend] = BackendDict(
    VPCLatticeBackend, "vpc-lattice"
)
