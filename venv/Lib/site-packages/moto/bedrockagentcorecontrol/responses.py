"""BedrockAgentCoreControl responses."""

from typing import Any

from moto.core.responses import ActionResult, BaseResponse, EmptyResult

from .models import BedrockAgentCoreControlBackend, bedrockagentcorecontrol_backends


class BedrockAgentCoreControlResponse(BaseResponse):
    def __init__(self) -> None:
        super().__init__(service_name="bedrock-agentcore-control")
        self.automated_parameter_parsing = True

    @property
    def backend(self) -> BedrockAgentCoreControlBackend:
        return bedrockagentcorecontrol_backends[self.current_account][self.region]

    def create_agent_runtime(self) -> ActionResult:
        params = self._get_params()
        runtime = self.backend.create_agent_runtime(
            agent_runtime_name=params["agentRuntimeName"],
            agent_runtime_artifact=params["agentRuntimeArtifact"],
            role_arn=params["roleArn"],
            network_configuration=params["networkConfiguration"],
            description=params.get("description"),
            authorizer_configuration=params.get("authorizerConfiguration"),
            request_header_configuration=params.get("requestHeaderConfiguration"),
            protocol_configuration=params.get("protocolConfiguration"),
            lifecycle_configuration=params.get("lifecycleConfiguration"),
            environment_variables=params.get("environmentVariables"),
            tags=params.get("tags"),
        )
        return ActionResult(
            {
                "agentRuntimeArn": runtime.agent_runtime_arn,
                "workloadIdentityDetails": runtime.workload_identity_details,
                "agentRuntimeId": runtime.agent_runtime_id,
                "agentRuntimeVersion": runtime.agent_runtime_version,
                "createdAt": runtime.created_at,
                "status": "CREATING",
            }
        )

    def get_agent_runtime(self) -> ActionResult:
        agent_runtime_id = self._get_param("agentRuntimeId")
        runtime = self.backend.get_agent_runtime(agent_runtime_id)
        result: dict[str, Any] = {
            "agentRuntimeArn": runtime.agent_runtime_arn,
            "agentRuntimeName": runtime.agent_runtime_name,
            "agentRuntimeId": runtime.agent_runtime_id,
            "agentRuntimeVersion": runtime.agent_runtime_version,
            "createdAt": runtime.created_at,
            "lastUpdatedAt": runtime.last_updated_at,
            "roleArn": runtime.role_arn,
            "networkConfiguration": runtime.network_configuration,
            "status": runtime.status,
            "lifecycleConfiguration": runtime.lifecycle_configuration,
            "description": runtime.description,
            "workloadIdentityDetails": runtime.workload_identity_details,
            "agentRuntimeArtifact": runtime.agent_runtime_artifact,
        }
        if runtime.protocol_configuration:
            result["protocolConfiguration"] = runtime.protocol_configuration
        if runtime.environment_variables:
            result["environmentVariables"] = runtime.environment_variables
        if runtime.authorizer_configuration:
            result["authorizerConfiguration"] = runtime.authorizer_configuration
        if runtime.request_header_configuration:
            result["requestHeaderConfiguration"] = runtime.request_header_configuration
        return ActionResult(result)

    def update_agent_runtime(self) -> ActionResult:
        agent_runtime_id = self._get_param("agentRuntimeId")
        params = self._get_params()
        runtime = self.backend.update_agent_runtime(
            agent_runtime_id=agent_runtime_id,
            agent_runtime_artifact=params["agentRuntimeArtifact"],
            role_arn=params["roleArn"],
            network_configuration=params["networkConfiguration"],
            description=params.get("description"),
            authorizer_configuration=params.get("authorizerConfiguration"),
            request_header_configuration=params.get("requestHeaderConfiguration"),
            protocol_configuration=params.get("protocolConfiguration"),
            lifecycle_configuration=params.get("lifecycleConfiguration"),
            environment_variables=params.get("environmentVariables"),
        )
        return ActionResult(
            {
                "agentRuntimeArn": runtime.agent_runtime_arn,
                "agentRuntimeId": runtime.agent_runtime_id,
                "workloadIdentityDetails": runtime.workload_identity_details,
                "agentRuntimeVersion": runtime.agent_runtime_version,
                "createdAt": runtime.created_at,
                "lastUpdatedAt": runtime.last_updated_at,
                "status": "UPDATING",
            }
        )

    def delete_agent_runtime(self) -> ActionResult:
        agent_runtime_id = self._get_param("agentRuntimeId")
        runtime = self.backend.delete_agent_runtime(agent_runtime_id)
        return ActionResult(
            {
                "status": "DELETING",
                "agentRuntimeId": runtime.agent_runtime_id,
            }
        )

    def list_agent_runtimes(self) -> ActionResult:
        runtimes = self.backend.list_agent_runtimes()
        return ActionResult(
            {
                "agentRuntimes": [r.to_summary() for r in runtimes],
            }
        )

    def list_agent_runtime_versions(self) -> ActionResult:
        agent_runtime_id = self._get_param("agentRuntimeId")
        versions = self.backend.list_agent_runtime_versions(agent_runtime_id)
        return ActionResult(
            {
                "agentRuntimes": versions,
            }
        )

    def create_agent_runtime_endpoint(self) -> ActionResult:
        agent_runtime_id = self._get_param("agentRuntimeId")
        params = self._get_params()
        endpoint = self.backend.create_agent_runtime_endpoint(
            agent_runtime_id=agent_runtime_id,
            name=params["name"],
            agent_runtime_version=params.get("agentRuntimeVersion"),
            description=params.get("description"),
            tags=params.get("tags"),
        )
        return ActionResult(
            {
                "targetVersion": endpoint.target_version,
                "agentRuntimeEndpointArn": endpoint.agent_runtime_endpoint_arn,
                "agentRuntimeArn": endpoint.agent_runtime_arn,
                "agentRuntimeId": endpoint.agent_runtime_id,
                "endpointName": endpoint.name,
                "status": "CREATING",
                "createdAt": endpoint.created_at,
            }
        )

    def get_agent_runtime_endpoint(self) -> ActionResult:
        agent_runtime_id = self._get_param("agentRuntimeId")
        endpoint_name = self._get_param("endpointName")
        endpoint = self.backend.get_agent_runtime_endpoint(
            agent_runtime_id, endpoint_name
        )
        result: dict[str, Any] = {
            "agentRuntimeEndpointArn": endpoint.agent_runtime_endpoint_arn,
            "agentRuntimeArn": endpoint.agent_runtime_arn,
            "status": endpoint.status,
            "createdAt": endpoint.created_at,
            "lastUpdatedAt": endpoint.last_updated_at,
            "name": endpoint.name,
            "id": endpoint.endpoint_id,
            "targetVersion": endpoint.target_version,
            "liveVersion": endpoint.live_version,
        }
        if endpoint.description:
            result["description"] = endpoint.description
        return ActionResult(result)

    def update_agent_runtime_endpoint(self) -> ActionResult:
        agent_runtime_id = self._get_param("agentRuntimeId")
        endpoint_name = self._get_param("endpointName")
        params = self._get_params()
        endpoint = self.backend.update_agent_runtime_endpoint(
            agent_runtime_id=agent_runtime_id,
            endpoint_name=endpoint_name,
            agent_runtime_version=params.get("agentRuntimeVersion"),
            description=params.get("description"),
        )
        return ActionResult(
            {
                "agentRuntimeEndpointArn": endpoint.agent_runtime_endpoint_arn,
                "agentRuntimeArn": endpoint.agent_runtime_arn,
                "status": "UPDATING",
                "createdAt": endpoint.created_at,
                "lastUpdatedAt": endpoint.last_updated_at,
                "liveVersion": endpoint.live_version,
                "targetVersion": endpoint.target_version,
            }
        )

    def delete_agent_runtime_endpoint(self) -> ActionResult:
        agent_runtime_id = self._get_param("agentRuntimeId")
        endpoint_name = self._get_param("endpointName")
        endpoint = self.backend.delete_agent_runtime_endpoint(
            agent_runtime_id, endpoint_name
        )
        return ActionResult(
            {
                "status": "DELETING",
                "agentRuntimeId": endpoint.agent_runtime_id,
                "endpointName": endpoint.name,
            }
        )

    def list_agent_runtime_endpoints(self) -> ActionResult:
        agent_runtime_id = self._get_param("agentRuntimeId")
        endpoints = self.backend.list_agent_runtime_endpoints(agent_runtime_id)
        return ActionResult(
            {
                "runtimeEndpoints": [ep.to_summary() for ep in endpoints],
            }
        )

    def create_gateway(self) -> ActionResult:
        params = self._get_params()
        gateway = self.backend.create_gateway(
            name=params["name"],
            role_arn=params["roleArn"],
            protocol_type=params["protocolType"],
            authorizer_type=params["authorizerType"],
            description=params.get("description"),
            protocol_configuration=params.get("protocolConfiguration"),
            authorizer_configuration=params.get("authorizerConfiguration"),
            kms_key_arn=params.get("kmsKeyArn"),
            interceptor_configurations=params.get("interceptorConfigurations"),
            policy_engine_configuration=params.get("policyEngineConfiguration"),
            exception_level=params.get("exceptionLevel"),
            tags=params.get("tags"),
        )
        result = gateway.to_dict()
        result["status"] = "CREATING"
        return ActionResult(result)

    def get_gateway(self) -> ActionResult:
        gateway_identifier = self._get_param("gatewayIdentifier")
        gateway = self.backend.get_gateway(gateway_identifier)
        return ActionResult(gateway.to_dict())

    def update_gateway(self) -> ActionResult:
        gateway_identifier = self._get_param("gatewayIdentifier")
        params = self._get_params()
        gateway = self.backend.update_gateway(
            gateway_identifier=gateway_identifier,
            name=params["name"],
            role_arn=params["roleArn"],
            protocol_type=params["protocolType"],
            authorizer_type=params["authorizerType"],
            description=params.get("description"),
            protocol_configuration=params.get("protocolConfiguration"),
            authorizer_configuration=params.get("authorizerConfiguration"),
            kms_key_arn=params.get("kmsKeyArn"),
            interceptor_configurations=params.get("interceptorConfigurations"),
            policy_engine_configuration=params.get("policyEngineConfiguration"),
            exception_level=params.get("exceptionLevel"),
        )
        result = gateway.to_dict()
        result["status"] = "UPDATING"
        return ActionResult(result)

    def delete_gateway(self) -> ActionResult:
        gateway_identifier = self._get_param("gatewayIdentifier")
        gateway = self.backend.delete_gateway(gateway_identifier)
        return ActionResult(
            {
                "gatewayId": gateway.gateway_id,
                "status": "DELETING",
            }
        )

    def list_gateways(self) -> ActionResult:
        gateways = self.backend.list_gateways()
        return ActionResult(
            {
                "items": [g.to_summary() for g in gateways],
            }
        )

    def create_gateway_target(self) -> ActionResult:
        gateway_identifier = self._get_param("gatewayIdentifier")
        params = self._get_params()
        target = self.backend.create_gateway_target(
            gateway_identifier=gateway_identifier,
            name=params["name"],
            target_configuration=params["targetConfiguration"],
            description=params.get("description"),
            credential_provider_configurations=params.get(
                "credentialProviderConfigurations"
            ),
            metadata_configuration=params.get("metadataConfiguration"),
        )
        result = target.to_dict()
        result["status"] = "CREATING"
        return ActionResult(result)

    def get_gateway_target(self) -> ActionResult:
        gateway_identifier = self._get_param("gatewayIdentifier")
        target_id = self._get_param("targetId")
        target = self.backend.get_gateway_target(gateway_identifier, target_id)
        return ActionResult(target.to_dict())

    def update_gateway_target(self) -> ActionResult:
        gateway_identifier = self._get_param("gatewayIdentifier")
        target_id = self._get_param("targetId")
        params = self._get_params()
        target = self.backend.update_gateway_target(
            gateway_identifier=gateway_identifier,
            target_id=target_id,
            name=params["name"],
            target_configuration=params["targetConfiguration"],
            description=params.get("description"),
            credential_provider_configurations=params.get(
                "credentialProviderConfigurations"
            ),
            metadata_configuration=params.get("metadataConfiguration"),
        )
        result = target.to_dict()
        result["status"] = "UPDATING"
        return ActionResult(result)

    def delete_gateway_target(self) -> ActionResult:
        gateway_identifier = self._get_param("gatewayIdentifier")
        target_id = self._get_param("targetId")
        target = self.backend.delete_gateway_target(gateway_identifier, target_id)
        return ActionResult(
            {
                "gatewayArn": target.gateway_arn,
                "targetId": target.target_id,
                "status": "DELETING",
            }
        )

    def list_gateway_targets(self) -> ActionResult:
        gateway_identifier = self._get_param("gatewayIdentifier")
        targets = self.backend.list_gateway_targets(gateway_identifier)
        return ActionResult(
            {
                "items": [t.to_summary() for t in targets],
            }
        )

    def create_memory(self) -> ActionResult:
        params = self._get_params()
        memory = self.backend.create_memory(
            name=params["name"],
            event_expiry_duration=params["eventExpiryDuration"],
            description=params.get("description"),
            encryption_key_arn=params.get("encryptionKeyArn"),
            memory_execution_role_arn=params.get("memoryExecutionRoleArn"),
            memory_strategies=params.get("memoryStrategies"),
            tags=params.get("tags"),
        )
        result = memory.to_dict()
        result["status"] = "CREATING"
        return ActionResult({"memory": result})

    def get_memory(self) -> ActionResult:
        memory_id = self._get_param("memoryId")
        memory = self.backend.get_memory(memory_id)
        return ActionResult({"memory": memory.to_dict()})

    def update_memory(self) -> ActionResult:
        memory_id = self._get_param("memoryId")
        params = self._get_params()
        memory = self.backend.update_memory(
            memory_id=memory_id,
            description=params.get("description"),
            event_expiry_duration=params.get("eventExpiryDuration"),
            memory_execution_role_arn=params.get("memoryExecutionRoleArn"),
            memory_strategies=params.get("memoryStrategies"),
        )
        return ActionResult({"memory": memory.to_dict()})

    def delete_memory(self) -> ActionResult:
        memory_id = self._get_param("memoryId")
        memory = self.backend.delete_memory(memory_id)
        return ActionResult(
            {
                "memoryId": memory.memory_id,
                "status": "DELETING",
            }
        )

    def list_memories(self) -> ActionResult:
        memories = self.backend.list_memories()
        return ActionResult(
            {
                "memories": [m.to_summary() for m in memories],
            }
        )

    def tag_resource(self) -> ActionResult:
        resource_arn = self._get_param("resourceArn")
        params = self._get_params()
        self.backend.tag_resource(resource_arn, params["tags"])
        return EmptyResult()

    def untag_resource(self) -> ActionResult:
        resource_arn = self._get_param("resourceArn")
        tag_keys = self._get_param("tagKeys", [])
        self.backend.untag_resource(resource_arn, tag_keys)
        return EmptyResult()

    def list_tags_for_resource(self) -> ActionResult:
        resource_arn = self._get_param("resourceArn")
        tags = self.backend.list_tags_for_resource(resource_arn)
        return ActionResult({"tags": tags})
