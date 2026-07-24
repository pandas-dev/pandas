import json
from urllib.parse import unquote

from moto.core.responses import BaseResponse

from .models import VPCLatticeBackend, vpclattice_backends


class VPCLatticeResponse(BaseResponse):
    def __init__(self) -> None:
        super().__init__(service_name="vpc-lattice")
        self.automated_parameter_parsing = True

    @property
    def backend(self) -> VPCLatticeBackend:
        return vpclattice_backends[self.current_account][self.region]

    def create_service(self) -> str:
        service = self.backend.create_service(
            auth_type=self._get_param("authType"),
            certificate_arn=self._get_param("certificateArn"),
            client_token=self._get_param("clientToken"),
            custom_domain_name=self._get_param("customDomainName"),
            name=self._get_param("name"),
            tags=self._get_param("tags"),
        )
        return json.dumps(service.to_dict())

    def get_service(self) -> str:
        path = unquote(self.path)
        service = self.backend.get_service(service_identifier=path.split("/")[-1])
        return json.dumps(service.to_dict())

    def list_services(self) -> str:
        max_results = self._get_param("MaxResults")
        next_token = self._get_param("NextToken")
        services, next_token = self.backend.list_services(
            max_results=max_results, next_token=next_token
        )
        response = {
            "items": [service.to_dict() for service in services],
            "nextToken": next_token,
        }
        return json.dumps(response)

    def create_service_network(self) -> str:
        sn = self.backend.create_service_network(
            auth_type=self._get_param("authType"),
            client_token=self._get_param("clientToken"),
            name=self._get_param("name"),
            sharing_config=self._get_param("sharingConfig"),
            tags=self._get_param("tags"),
        )
        return json.dumps(sn.to_dict())

    def get_service_network(self) -> str:
        path = unquote(self.path)
        service_network = self.backend.get_service_network(
            service_network_identifier=path.split("/")[-1]
        )
        return json.dumps(service_network.to_dict())

    def list_service_networks(self) -> str:
        max_results = self._get_param("MaxResults")
        next_token = self._get_param("NextToken")
        service_networks, next_token = self.backend.list_service_networks(
            max_results=max_results, next_token=next_token
        )
        response = {
            "items": [
                service_network.to_dict() for service_network in service_networks
            ],
            "nextToken": next_token,
        }
        return json.dumps(response)

    def create_service_network_vpc_association(self) -> str:
        assoc = self.backend.create_service_network_vpc_association(
            client_token=self._get_param("clientToken"),
            security_group_ids=self._get_param("securityGroupIds"),
            service_network_identifier=self._get_param("serviceNetworkIdentifier"),
            tags=self._get_param("tags"),
            vpc_identifier=self._get_param("vpcIdentifier"),
        )
        return json.dumps(assoc.to_dict())

    def create_rule(self) -> str:
        rule = self.backend.create_rule(
            action=self._get_param("action"),
            client_token=self._get_param("clientToken"),
            listener_identifier=self._get_param("listenerIdentifier"),
            match=self._get_param("match"),
            name=self._get_param("name"),
            priority=self._get_param("priority"),
            service_identifier=self._get_param("serviceIdentifier"),
            tags=self._get_param("tags"),
        )
        return json.dumps(rule.to_dict())

    def list_tags_for_resource(self) -> str:
        resource_arn = unquote(self._get_param("resourceArn"))
        tags = self.backend.list_tags_for_resource(resource_arn)
        return json.dumps({"tags": tags})

    def tag_resource(self) -> str:
        resource_arn = unquote(self._get_param("resourceArn"))
        tags = self._get_param("tags")
        self.backend.tag_resource(resource_arn, tags)
        return json.dumps({})

    def untag_resource(self) -> str:
        resource_arn = unquote(self._get_param("resourceArn"))
        tag_keys = self._get_param("tagKeys", [])
        self.backend.untag_resource(resource_arn=resource_arn, tag_keys=tag_keys)
        return json.dumps({})

    def create_access_log_subscription(self) -> str:
        sub = self.backend.create_access_log_subscription(
            resourceIdentifier=self._get_param("resourceIdentifier"),
            destinationArn=self._get_param("destinationArn"),
            client_token=self._get_param("clientToken"),
            serviceNetworkLogType=self._get_param("serviceNetworkLogType"),
            tags=self._get_param("tags"),
        )

        return json.dumps(sub.to_dict())

    def get_access_log_subscription(self) -> str:
        path = unquote(self.path)

        sub = self.backend.get_access_log_subscription(
            accessLogSubscriptionIdentifier=path.split("/")[-1]
        )

        return json.dumps(sub.to_dict())

    def list_access_log_subscriptions(self) -> str:
        subs = self.backend.list_access_log_subscriptions(
            resourceIdentifier=self._get_param("resourceIdentifier"),
            maxResults=self._get_int_param("maxResults"),
            nextToken=self._get_param("nextToken"),
        )

        return json.dumps({"items": [s.to_dict() for s in subs], "nextToken": ""})

    def update_access_log_subscription(self) -> str:
        path = unquote(self.path)

        sub = self.backend.update_access_log_subscription(
            accessLogSubscriptionIdentifier=path.split("/")[-1],
            destinationArn=self._get_param("destinationArn"),
        )

        return json.dumps(sub.to_dict())

    def delete_access_log_subscription(self) -> str:
        path = unquote(self.path)

        self.backend.delete_access_log_subscription(
            accessLogSubscriptionIdentifier=path.split("/")[-1]
        )

        return json.dumps({})

    def put_auth_policy(self) -> str:
        resourceId = unquote(self._get_param("resourceIdentifier"))
        policy = self._get_param("policy")

        auth_policy = self.backend.put_auth_policy(
            resourceIdentifier=resourceId,
            policy=policy,
        )

        response = {
            "policy": auth_policy.policy,
            "state": auth_policy.state,
        }
        return json.dumps(response)

    def get_auth_policy(self) -> str:
        resourceId = unquote(self._get_param("resourceIdentifier"))
        auth_policy = self.backend.get_auth_policy(resourceIdentifier=resourceId)

        response = {
            "policy": auth_policy.policy,
            "state": auth_policy.state,
            "createdAt": auth_policy.created_at,
            "lastUpdatedAt": auth_policy.last_updated_at,
        }
        return json.dumps(response)

    def delete_auth_policy(self) -> str:
        resourceId = unquote(self._get_param("resourceIdentifier"))

        self.backend.delete_auth_policy(resourceIdentifier=resourceId)

        return "{}"

    def put_resource_policy(self) -> str:
        resource_arn = unquote(self._get_param("resourceArn"))
        policy = self._get_param("policy")
        self.backend.put_resource_policy(
            resourceArn=resource_arn,
            policy=policy,
        )

        return "{}"

    def get_resource_policy(self) -> str:
        resource_arn = unquote(self._get_param("resourceArn"))

        resource_policy = self.backend.get_resource_policy(resourceArn=resource_arn)

        return json.dumps({"policy": resource_policy})

    def delete_resource_policy(self) -> str:
        resource_arn = unquote(self._get_param("resourceArn"))

        self.backend.delete_resource_policy(resourceArn=resource_arn)

        return "{}"

    def list_service_network_vpc_associations(self) -> str:
        service_network_identifier = self._get_param("serviceNetworkIdentifier")
        if service_network_identifier:
            service_network_identifier = unquote(service_network_identifier)

        vpc_identifier = self._get_param("vpcIdentifier")
        if vpc_identifier:
            vpc_identifier = unquote(vpc_identifier)

        max_results = self._get_int_param("maxResults")
        next_token = self._get_param("nextToken")

        associations, next_token = self.backend.list_service_network_vpc_associations(
            serviceNetworkIdentifier=service_network_identifier,
            vpcIdentifier=vpc_identifier,
            max_results=max_results,
            next_token=next_token,
        )
        return json.dumps(
            {
                "items": [assoc.to_dict() for assoc in associations],
                "nextToken": next_token,
            }
        )

    def create_listener(self) -> str:
        service = self.backend.create_listener(
            service_identifier=self._get_param("serviceIdentifier"),
            name=self._get_param("name"),
            protocol=self._get_param("protocol"),
            port=self._get_param("port"),
            default_action=self._get_param("defaultAction"),
            client_token=self._get_param("clientToken"),
            tags=self._get_param("tags"),
        )
        return json.dumps(service.to_dict())

    def get_listener(self) -> str:
        path = unquote(self.path)
        service_listener = self.backend.get_listener(
            service_identifier=path.split("/")[-3],
            listener_identifier=path.split("/")[-1],
        )
        return json.dumps(service_listener.to_dict())

    def list_listeners(self) -> str:
        service_identifier = self._get_param("serviceIdentifier")
        max_results = self._get_param("maxResults")
        next_token = self._get_param("nextToken")
        service_listeners, next_token = self.backend.list_listeners(
            service_identifier=service_identifier,
            max_results=max_results,
            next_token=next_token,
        )
        response = {
            "items": [
                service_listener.to_dict() for service_listener in service_listeners
            ],
            "nextToken": next_token,
        }
        return json.dumps(response)

    def create_target_group(self) -> str:
        target_group = self.backend.create_target_group(
            name=self._get_param("name"),
            target_type=self._get_param("type"),
            config=self._get_param("config"),
            client_token=self._get_param("clientToken"),
            tags=self._get_param("tags"),
        )
        return json.dumps(target_group.to_dict())

    def get_target_group(self) -> str:
        path = unquote(self.path)
        target_group = self.backend.get_target_group(
            target_group_identifier=path.split("/")[-1]
        )
        return json.dumps(target_group.to_dict())

    def list_target_groups(self) -> str:
        vpc_identifier = self._get_param("vpcIdentifier")
        target_group_type = self._get_param("targetGroupType")
        max_results = self._get_param("maxResults")
        next_token = self._get_param("nextToken")
        target_groups, next_token = self.backend.list_target_groups(
            vpc_identifier=vpc_identifier,
            target_group_type=target_group_type,
            max_results=max_results,
            next_token=next_token,
        )
        response = {
            "items": [target_group.to_dict() for target_group in target_groups],
            "nextToken": next_token,
        }
        return json.dumps(response)

    def create_service_network_resource_association(self) -> str:
        assoc = self.backend.create_service_network_resource_association(
            resource_config_identifier=self._get_param(
                "resourceConfigurationIdentifier"
            ),
            service_network_identifier=self._get_param("serviceNetworkIdentifier"),
            private_dns_enabled=self._get_param("privateDnsEnabled"),
            client_token=self._get_param("clientToken"),
            tags=self._get_param("tags"),
        )
        return json.dumps(assoc.to_dict())

    def get_service_network_resource_association(self) -> str:
        path = unquote(self.path)
        assoc = self.backend.get_service_network_resource_association(
            service_network_resource_association_identifier=path.split("/")[-1]
        )
        return json.dumps(assoc.to_dict())

    def list_service_network_resource_associations(self) -> str:
        resource_configuration_identifier = self._get_param(
            "resourceConfigurationIdentifier"
        )
        service_network_identifier = self._get_param("serviceNetworkIdentifier")
        max_results = self._get_param("maxResults")
        next_token = self._get_param("nextToken")
        assocs, next_token = self.backend.list_service_network_resource_associations(
            resource_configuration_identifier=resource_configuration_identifier,
            service_network_identifier=service_network_identifier,
            max_results=max_results,
            next_token=next_token,
        )
        response = {
            "items": [assoc.to_dict() for assoc in assocs],
            "nextToken": next_token,
        }
        return json.dumps(response)

    def create_service_network_service_association(self) -> str:
        assoc = self.backend.create_service_network_service_association(
            service_identifier=self._get_param("serviceIdentifier"),
            service_network_identifier=self._get_param("serviceNetworkIdentifier"),
            client_token=self._get_param("clientToken"),
            tags=self._get_param("tags"),
        )
        return json.dumps(assoc.to_dict())

    def get_service_network_service_association(self) -> str:
        path = unquote(self.path)
        assoc = self.backend.get_service_network_service_association(
            service_network_service_association_identifier=path.split("/")[-1]
        )
        return json.dumps(assoc.to_dict())

    def list_service_network_service_associations(self) -> str:
        service_identifier = self._get_param("serviceIdentifier")
        service_network_identifier = self._get_param("serviceNetworkIdentifier")
        max_results = self._get_param("maxResults")
        next_token = self._get_param("nextToken")
        assocs, next_token = self.backend.list_service_network_service_associations(
            service_network_identifier=service_network_identifier,
            service_identifier=service_identifier,
            max_results=max_results,
            next_token=next_token,
        )
        response = {
            "items": [assoc.to_dict() for assoc in assocs],
            "nextToken": next_token,
        }
        return json.dumps(response)
