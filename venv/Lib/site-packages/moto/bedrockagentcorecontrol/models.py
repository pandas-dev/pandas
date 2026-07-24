"""BedrockAgentCoreControl models."""

from collections import OrderedDict
from typing import Any

from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.common_models import BaseModel
from moto.core.utils import utcnow
from moto.moto_api._internal import mock_random
from moto.moto_api._internal.managed_state_model import ManagedState
from moto.utilities.tagging_service import TaggingService
from moto.utilities.utils import get_partition

from .exceptions import ConflictException, ResourceNotFoundException


class AgentRuntime(BaseModel, ManagedState):
    def __init__(
        self,
        region_name: str,
        account_id: str,
        agent_runtime_name: str,
        agent_runtime_artifact: dict[str, Any],
        role_arn: str,
        network_configuration: dict[str, Any],
        description: str | None,
        authorizer_configuration: dict[str, Any] | None,
        request_header_configuration: dict[str, Any] | None,
        protocol_configuration: dict[str, Any] | None,
        lifecycle_configuration: dict[str, Any] | None,
        environment_variables: dict[str, str] | None,
    ):
        ManagedState.__init__(
            self,
            "bedrock-agentcore-control::agent_runtime",
            transitions=[("CREATING", "READY")],
        )
        self.region_name = region_name
        self.account_id = account_id
        self.agent_runtime_id = (
            f"a{mock_random.get_random_hex(9)}-{mock_random.get_random_hex(10)}"
        )
        self.agent_runtime_version = "1"
        runtime_uuid = str(mock_random.uuid4())
        self.agent_runtime_arn = f"arn:{get_partition(region_name)}:bedrock-agentcore:{region_name}:{account_id}:agent/{runtime_uuid}:{self.agent_runtime_version}"
        self.agent_runtime_name = agent_runtime_name
        self.agent_runtime_artifact = agent_runtime_artifact
        self.role_arn = role_arn
        self.network_configuration = network_configuration
        self.description = description or ""
        self.authorizer_configuration = authorizer_configuration
        self.request_header_configuration = request_header_configuration
        self.protocol_configuration = protocol_configuration
        self.lifecycle_configuration = lifecycle_configuration or {}
        self.environment_variables = environment_variables
        now = utcnow()
        self.created_at = now
        self.last_updated_at = now
        self.workload_identity_details = {
            "workloadIdentityArn": f"arn:{get_partition(region_name)}:bedrock-agentcore:{region_name}:{account_id}:workload-identity-directory/default/workload-identity/{self.agent_runtime_id}"
        }
        # Store version snapshots for ListAgentRuntimeVersions
        self.versions: list[dict[str, Any]] = []
        self._snapshot_version()

    def _snapshot_version(self) -> None:
        self.versions.append(
            {
                "agentRuntimeArn": self.agent_runtime_arn,
                "agentRuntimeId": self.agent_runtime_id,
                "agentRuntimeVersion": self.agent_runtime_version,
                "agentRuntimeName": self.agent_runtime_name,
                "description": self.description,
                "lastUpdatedAt": self.last_updated_at,
                "status": self.status,
            }
        )

    def update(
        self,
        agent_runtime_artifact: dict[str, Any],
        role_arn: str,
        network_configuration: dict[str, Any],
        description: str | None,
        authorizer_configuration: dict[str, Any] | None,
        request_header_configuration: dict[str, Any] | None,
        protocol_configuration: dict[str, Any] | None,
        lifecycle_configuration: dict[str, Any] | None,
        environment_variables: dict[str, str] | None,
    ) -> None:
        self.agent_runtime_artifact = agent_runtime_artifact
        self.role_arn = role_arn
        self.network_configuration = network_configuration
        if description is not None:
            self.description = description
        if authorizer_configuration is not None:
            self.authorizer_configuration = authorizer_configuration
        if request_header_configuration is not None:
            self.request_header_configuration = request_header_configuration
        if protocol_configuration is not None:
            self.protocol_configuration = protocol_configuration
        if lifecycle_configuration is not None:
            self.lifecycle_configuration = lifecycle_configuration
        if environment_variables is not None:
            self.environment_variables = environment_variables
        new_version = str(int(self.agent_runtime_version) + 1)
        self.agent_runtime_version = new_version
        runtime_uuid = self.agent_runtime_arn.split("agent/")[1].split(":")[0]
        self.agent_runtime_arn = f"arn:{get_partition(self.region_name)}:bedrock-agentcore:{self.region_name}:{self.account_id}:agent/{runtime_uuid}:{new_version}"
        self.last_updated_at = utcnow()
        self.status = "UPDATING"
        self._snapshot_version()

    def to_summary(self) -> dict[str, Any]:
        return {
            "agentRuntimeArn": self.agent_runtime_arn,
            "agentRuntimeId": self.agent_runtime_id,
            "agentRuntimeVersion": self.agent_runtime_version,
            "agentRuntimeName": self.agent_runtime_name,
            "description": self.description,
            "lastUpdatedAt": self.last_updated_at,
            "status": self.status,
        }


class AgentRuntimeEndpoint(BaseModel, ManagedState):
    def __init__(
        self,
        region_name: str,
        account_id: str,
        agent_runtime: "AgentRuntime",
        name: str,
        agent_runtime_version: str | None,
        description: str | None,
    ):
        ManagedState.__init__(
            self,
            "bedrock-agentcore-control::agent_runtime_endpoint",
            transitions=[("CREATING", "READY")],
        )
        self.region_name = region_name
        self.account_id = account_id
        self.name = name
        self.endpoint_id = (
            f"e{mock_random.get_random_hex(9)}-{mock_random.get_random_hex(10)}"
        )
        endpoint_uuid = str(mock_random.uuid4())
        self.agent_runtime_endpoint_arn = f"arn:{get_partition(region_name)}:bedrock-agentcore:{region_name}:{account_id}:agentEndpoint/{endpoint_uuid}"
        self.agent_runtime_arn = agent_runtime.agent_runtime_arn
        self.agent_runtime_id = agent_runtime.agent_runtime_id
        self.target_version = (
            agent_runtime_version or agent_runtime.agent_runtime_version
        )
        self.live_version = self.target_version
        self.description = description or ""
        now = utcnow()
        self.created_at = now
        self.last_updated_at = now

    def update(
        self,
        agent_runtime_version: str | None,
        description: str | None,
    ) -> None:
        if agent_runtime_version is not None:
            self.target_version = agent_runtime_version
        if description is not None:
            self.description = description
        self.last_updated_at = utcnow()
        self.status = "UPDATING"

    def to_summary(self) -> dict[str, Any]:
        return {
            "name": self.name,
            "agentRuntimeEndpointArn": self.agent_runtime_endpoint_arn,
            "agentRuntimeArn": self.agent_runtime_arn,
            "status": self.status,
            "id": self.endpoint_id,
            "description": self.description,
            "createdAt": self.created_at,
            "lastUpdatedAt": self.last_updated_at,
        }


class Gateway(BaseModel, ManagedState):
    def __init__(
        self,
        region_name: str,
        account_id: str,
        name: str,
        role_arn: str,
        protocol_type: str,
        authorizer_type: str,
        description: str | None,
        protocol_configuration: dict[str, Any] | None,
        authorizer_configuration: dict[str, Any] | None,
        kms_key_arn: str | None,
        interceptor_configurations: list[dict[str, Any]] | None,
        policy_engine_configuration: dict[str, Any] | None,
        exception_level: str | None,
    ):
        ManagedState.__init__(
            self,
            "bedrock-agentcore-control::gateway",
            transitions=[("CREATING", "READY")],
        )
        self.region_name = region_name
        self.account_id = account_id
        self.name = name
        self.gateway_id = mock_random.get_random_hex(10)
        self.gateway_arn = f"arn:{get_partition(region_name)}:bedrock-agentcore:{region_name}:{account_id}:gateway/{self.gateway_id}"
        self.gateway_url = f"https://{self.gateway_id}.gateway.bedrock-agentcore.{region_name}.amazonaws.com"
        self.role_arn = role_arn
        self.protocol_type = protocol_type
        self.authorizer_type = authorizer_type
        self.description = description or ""
        self.protocol_configuration = protocol_configuration
        self.authorizer_configuration = authorizer_configuration
        self.kms_key_arn = kms_key_arn
        self.interceptor_configurations = interceptor_configurations
        self.policy_engine_configuration = policy_engine_configuration
        self.exception_level = exception_level
        now = utcnow()
        self.created_at = now
        self.updated_at = now
        self.workload_identity_details = {
            "workloadIdentityArn": f"arn:{get_partition(region_name)}:bedrock-agentcore:{region_name}:{account_id}:workload-identity-directory/default/workload-identity/{self.gateway_id}"
        }

    def update(
        self,
        name: str,
        role_arn: str,
        protocol_type: str,
        authorizer_type: str,
        description: str | None,
        protocol_configuration: dict[str, Any] | None,
        authorizer_configuration: dict[str, Any] | None,
        kms_key_arn: str | None,
        interceptor_configurations: list[dict[str, Any]] | None,
        policy_engine_configuration: dict[str, Any] | None,
        exception_level: str | None,
    ) -> None:
        self.name = name
        self.role_arn = role_arn
        self.protocol_type = protocol_type
        self.authorizer_type = authorizer_type
        if description is not None:
            self.description = description
        if protocol_configuration is not None:
            self.protocol_configuration = protocol_configuration
        if authorizer_configuration is not None:
            self.authorizer_configuration = authorizer_configuration
        if kms_key_arn is not None:
            self.kms_key_arn = kms_key_arn
        if interceptor_configurations is not None:
            self.interceptor_configurations = interceptor_configurations
        if policy_engine_configuration is not None:
            self.policy_engine_configuration = policy_engine_configuration
        if exception_level is not None:
            self.exception_level = exception_level
        self.updated_at = utcnow()
        self.status = "UPDATING"

    def to_dict(self) -> dict[str, Any]:
        result: dict[str, Any] = {
            "gatewayArn": self.gateway_arn,
            "gatewayId": self.gateway_id,
            "gatewayUrl": self.gateway_url,
            "createdAt": self.created_at,
            "updatedAt": self.updated_at,
            "status": self.status,
            "name": self.name,
            "roleArn": self.role_arn,
            "protocolType": self.protocol_type,
            "authorizerType": self.authorizer_type,
            "description": self.description,
            "workloadIdentityDetails": self.workload_identity_details,
        }
        if self.protocol_configuration:
            result["protocolConfiguration"] = self.protocol_configuration
        if self.authorizer_configuration:
            result["authorizerConfiguration"] = self.authorizer_configuration
        if self.kms_key_arn:
            result["kmsKeyArn"] = self.kms_key_arn
        if self.interceptor_configurations:
            result["interceptorConfigurations"] = self.interceptor_configurations
        if self.policy_engine_configuration:
            result["policyEngineConfiguration"] = self.policy_engine_configuration
        if self.exception_level:
            result["exceptionLevel"] = self.exception_level
        return result

    def to_summary(self) -> dict[str, Any]:
        result: dict[str, Any] = {
            "gatewayId": self.gateway_id,
            "name": self.name,
            "status": self.status,
            "createdAt": self.created_at,
            "updatedAt": self.updated_at,
            "authorizerType": self.authorizer_type,
            "protocolType": self.protocol_type,
        }
        if self.description:
            result["description"] = self.description
        return result


class GatewayTarget(BaseModel, ManagedState):
    def __init__(
        self,
        region_name: str,
        account_id: str,
        gateway: Gateway,
        name: str,
        target_configuration: dict[str, Any],
        description: str | None,
        credential_provider_configurations: list[dict[str, Any]] | None,
        metadata_configuration: dict[str, Any] | None,
    ):
        ManagedState.__init__(
            self,
            "bedrock-agentcore-control::gateway_target",
            transitions=[("CREATING", "READY")],
        )
        self.region_name = region_name
        self.account_id = account_id
        self.name = name
        self.target_id = mock_random.get_random_hex(10)
        self.gateway_arn = gateway.gateway_arn
        self.gateway_id = gateway.gateway_id
        self.target_configuration = target_configuration
        self.description = description or ""
        self.credential_provider_configurations = (
            credential_provider_configurations or []
        )
        self.metadata_configuration = metadata_configuration
        now = utcnow()
        self.created_at = now
        self.updated_at = now

    def update(
        self,
        name: str,
        target_configuration: dict[str, Any],
        description: str | None,
        credential_provider_configurations: list[dict[str, Any]] | None,
        metadata_configuration: dict[str, Any] | None,
    ) -> None:
        self.name = name
        self.target_configuration = target_configuration
        if description is not None:
            self.description = description
        if credential_provider_configurations is not None:
            self.credential_provider_configurations = credential_provider_configurations
        if metadata_configuration is not None:
            self.metadata_configuration = metadata_configuration
        self.updated_at = utcnow()
        self.status = "UPDATING"

    def to_dict(self) -> dict[str, Any]:
        result: dict[str, Any] = {
            "gatewayArn": self.gateway_arn,
            "targetId": self.target_id,
            "createdAt": self.created_at,
            "updatedAt": self.updated_at,
            "status": self.status,
            "name": self.name,
            "targetConfiguration": self.target_configuration,
            "credentialProviderConfigurations": self.credential_provider_configurations,
        }
        if self.description:
            result["description"] = self.description
        if self.metadata_configuration:
            result["metadataConfiguration"] = self.metadata_configuration
        return result

    def to_summary(self) -> dict[str, Any]:
        result: dict[str, Any] = {
            "targetId": self.target_id,
            "name": self.name,
            "status": self.status,
            "createdAt": self.created_at,
            "updatedAt": self.updated_at,
        }
        if self.description:
            result["description"] = self.description
        return result


class Memory(BaseModel, ManagedState):
    def __init__(
        self,
        region_name: str,
        account_id: str,
        name: str,
        event_expiry_duration: int,
        description: str | None,
        encryption_key_arn: str | None,
        memory_execution_role_arn: str | None,
        memory_strategies: list[dict[str, Any]] | None,
    ):
        ManagedState.__init__(
            self,
            "bedrock-agentcore-control::memory",
            transitions=[("CREATING", "ACTIVE")],
        )
        self.region_name = region_name
        self.account_id = account_id
        self.name = name
        self.memory_id = (
            f"m{mock_random.get_random_hex(9)}-{mock_random.get_random_hex(10)}"
        )
        self.memory_arn = f"arn:{get_partition(region_name)}:bedrock-agentcore:{region_name}:{account_id}:memory/{self.memory_id}"
        self.event_expiry_duration = event_expiry_duration
        self.description = description or ""
        self.encryption_key_arn = encryption_key_arn
        self.memory_execution_role_arn = memory_execution_role_arn
        self.strategies: list[dict[str, Any]] = []
        if memory_strategies:
            for strategy_input in memory_strategies:
                self._add_strategy(strategy_input)
        now = utcnow()
        self.created_at = now
        self.updated_at = now

    def _add_strategy(self, strategy_input: dict[str, Any]) -> None:
        strategy_id = (
            f"s{mock_random.get_random_hex(9)}-{mock_random.get_random_hex(10)}"
        )
        strategy_type_key = next(iter(strategy_input))
        strategy_data = strategy_input[strategy_type_key]
        type_map = {
            "semanticMemoryStrategy": "SEMANTIC",
            "summaryMemoryStrategy": "SUMMARIZATION",
            "userPreferenceMemoryStrategy": "USER_PREFERENCE",
            "customMemoryStrategy": "CUSTOM",
            "episodicMemoryStrategy": "EPISODIC",
        }
        now = utcnow()
        self.strategies.append(
            {
                "strategyId": strategy_id,
                "name": strategy_data.get("name", ""),
                "description": strategy_data.get("description", ""),
                "type": type_map.get(strategy_type_key, "CUSTOM"),
                "namespaces": strategy_data.get("namespaces", []),
                "status": "ACTIVE",
                "createdAt": now,
                "updatedAt": now,
            }
        )

    def update(
        self,
        description: str | None,
        event_expiry_duration: int | None,
        memory_execution_role_arn: str | None,
        memory_strategies: dict[str, Any] | None,
    ) -> None:
        if description is not None:
            self.description = description
        if event_expiry_duration is not None:
            self.event_expiry_duration = event_expiry_duration
        if memory_execution_role_arn is not None:
            self.memory_execution_role_arn = memory_execution_role_arn
        if memory_strategies:
            for strategy_input in memory_strategies.get("addMemoryStrategies", []):
                self._add_strategy(strategy_input)
            for modify in memory_strategies.get("modifyMemoryStrategies", []):
                sid = modify["memoryStrategyId"]
                for s in self.strategies:
                    if s["strategyId"] == sid:
                        if "description" in modify:
                            s["description"] = modify["description"]
                        if "namespaces" in modify:
                            s["namespaces"] = modify["namespaces"]
                        s["updatedAt"] = utcnow()
                        break
            for delete in memory_strategies.get("deleteMemoryStrategies", []):
                sid = delete["memoryStrategyId"]
                self.strategies = [s for s in self.strategies if s["strategyId"] != sid]
        self.updated_at = utcnow()
        self.status = "ACTIVE"

    def to_dict(self) -> dict[str, Any]:
        result: dict[str, Any] = {
            "arn": self.memory_arn,
            "id": self.memory_id,
            "name": self.name,
            "eventExpiryDuration": self.event_expiry_duration,
            "status": self.status,
            "createdAt": self.created_at,
            "updatedAt": self.updated_at,
        }
        if self.description:
            result["description"] = self.description
        if self.encryption_key_arn:
            result["encryptionKeyArn"] = self.encryption_key_arn
        if self.memory_execution_role_arn:
            result["memoryExecutionRoleArn"] = self.memory_execution_role_arn
        if self.strategies:
            result["strategies"] = self.strategies
        return result

    def to_summary(self) -> dict[str, Any]:
        return {
            "arn": self.memory_arn,
            "id": self.memory_id,
            "status": self.status,
            "createdAt": self.created_at,
            "updatedAt": self.updated_at,
        }


class BedrockAgentCoreControlBackend(BaseBackend):
    def __init__(self, region_name: str, account_id: str):
        super().__init__(region_name, account_id)
        self.agent_runtimes: dict[str, AgentRuntime] = OrderedDict()
        # endpoints keyed by (agent_runtime_id, endpoint_name)
        self.agent_runtime_endpoints: dict[tuple[str, str], AgentRuntimeEndpoint] = (
            OrderedDict()
        )
        self.gateways: dict[str, Gateway] = OrderedDict()
        self.gateway_targets: dict[tuple[str, str], GatewayTarget] = OrderedDict()
        self.memories: dict[str, Memory] = OrderedDict()
        self.tagger = TaggingService()

    def _get_runtime(self, agent_runtime_id: str) -> AgentRuntime:
        if agent_runtime_id not in self.agent_runtimes:
            raise ResourceNotFoundException(
                f"Could not find Agent Runtime with ID {agent_runtime_id}"
            )
        return self.agent_runtimes[agent_runtime_id]

    def create_agent_runtime(
        self,
        agent_runtime_name: str,
        agent_runtime_artifact: dict[str, Any],
        role_arn: str,
        network_configuration: dict[str, Any],
        description: str | None,
        authorizer_configuration: dict[str, Any] | None,
        request_header_configuration: dict[str, Any] | None,
        protocol_configuration: dict[str, Any] | None,
        lifecycle_configuration: dict[str, Any] | None,
        environment_variables: dict[str, str] | None,
        tags: dict[str, str] | None,
    ) -> AgentRuntime:
        runtime = AgentRuntime(
            region_name=self.region_name,
            account_id=self.account_id,
            agent_runtime_name=agent_runtime_name,
            agent_runtime_artifact=agent_runtime_artifact,
            role_arn=role_arn,
            network_configuration=network_configuration,
            description=description,
            authorizer_configuration=authorizer_configuration,
            request_header_configuration=request_header_configuration,
            protocol_configuration=protocol_configuration,
            lifecycle_configuration=lifecycle_configuration,
            environment_variables=environment_variables,
        )
        self.agent_runtimes[runtime.agent_runtime_id] = runtime
        if tags:
            self.tagger.tag_resource(
                runtime.agent_runtime_arn,
                [{"Key": k, "Value": v} for k, v in tags.items()],
            )
        return runtime

    def get_agent_runtime(self, agent_runtime_id: str) -> AgentRuntime:
        runtime = self._get_runtime(agent_runtime_id)
        runtime.advance()
        return runtime

    def update_agent_runtime(
        self,
        agent_runtime_id: str,
        agent_runtime_artifact: dict[str, Any],
        role_arn: str,
        network_configuration: dict[str, Any],
        description: str | None,
        authorizer_configuration: dict[str, Any] | None,
        request_header_configuration: dict[str, Any] | None,
        protocol_configuration: dict[str, Any] | None,
        lifecycle_configuration: dict[str, Any] | None,
        environment_variables: dict[str, str] | None,
    ) -> AgentRuntime:
        runtime = self._get_runtime(agent_runtime_id)
        runtime.update(
            agent_runtime_artifact=agent_runtime_artifact,
            role_arn=role_arn,
            network_configuration=network_configuration,
            description=description,
            authorizer_configuration=authorizer_configuration,
            request_header_configuration=request_header_configuration,
            protocol_configuration=protocol_configuration,
            lifecycle_configuration=lifecycle_configuration,
            environment_variables=environment_variables,
        )
        return runtime

    def delete_agent_runtime(self, agent_runtime_id: str) -> AgentRuntime:
        runtime = self._get_runtime(agent_runtime_id)
        # Remove all endpoints for this runtime
        keys_to_remove = [
            k for k in self.agent_runtime_endpoints if k[0] == agent_runtime_id
        ]
        for key in keys_to_remove:
            self.agent_runtime_endpoints.pop(key)
        self.agent_runtimes.pop(agent_runtime_id)
        return runtime

    def list_agent_runtimes(self) -> list[AgentRuntime]:
        return list(self.agent_runtimes.values())

    def list_agent_runtime_versions(
        self, agent_runtime_id: str
    ) -> list[dict[str, Any]]:
        runtime = self._get_runtime(agent_runtime_id)
        return list(reversed(runtime.versions))

    def create_agent_runtime_endpoint(
        self,
        agent_runtime_id: str,
        name: str,
        agent_runtime_version: str | None,
        description: str | None,
        tags: dict[str, str] | None,
    ) -> AgentRuntimeEndpoint:
        runtime = self._get_runtime(agent_runtime_id)
        key = (agent_runtime_id, name)
        if key in self.agent_runtime_endpoints:
            raise ConflictException(
                f"Endpoint {name} already exists for Agent Runtime {agent_runtime_id}"
            )
        endpoint = AgentRuntimeEndpoint(
            region_name=self.region_name,
            account_id=self.account_id,
            agent_runtime=runtime,
            name=name,
            agent_runtime_version=agent_runtime_version,
            description=description,
        )
        self.agent_runtime_endpoints[key] = endpoint
        if tags:
            self.tagger.tag_resource(
                endpoint.agent_runtime_endpoint_arn,
                [{"Key": k, "Value": v} for k, v in tags.items()],
            )
        return endpoint

    def get_agent_runtime_endpoint(
        self, agent_runtime_id: str, endpoint_name: str
    ) -> AgentRuntimeEndpoint:
        key = (agent_runtime_id, endpoint_name)
        if key not in self.agent_runtime_endpoints:
            raise ResourceNotFoundException(
                f"Could not find endpoint {endpoint_name} for Agent Runtime {agent_runtime_id}"
            )
        endpoint = self.agent_runtime_endpoints[key]
        endpoint.advance()
        return endpoint

    def update_agent_runtime_endpoint(
        self,
        agent_runtime_id: str,
        endpoint_name: str,
        agent_runtime_version: str | None,
        description: str | None,
    ) -> AgentRuntimeEndpoint:
        endpoint = self.get_agent_runtime_endpoint(agent_runtime_id, endpoint_name)
        endpoint.update(
            agent_runtime_version=agent_runtime_version,
            description=description,
        )
        return endpoint

    def delete_agent_runtime_endpoint(
        self, agent_runtime_id: str, endpoint_name: str
    ) -> AgentRuntimeEndpoint:
        endpoint = self.get_agent_runtime_endpoint(agent_runtime_id, endpoint_name)
        key = (agent_runtime_id, endpoint_name)
        self.agent_runtime_endpoints.pop(key)
        return endpoint

    def list_agent_runtime_endpoints(
        self, agent_runtime_id: str
    ) -> list[AgentRuntimeEndpoint]:
        self._get_runtime(agent_runtime_id)
        return [
            ep
            for (rid, _), ep in self.agent_runtime_endpoints.items()
            if rid == agent_runtime_id
        ]

    def _get_gateway(self, gateway_identifier: str) -> Gateway:
        if gateway_identifier not in self.gateways:
            raise ResourceNotFoundException(
                f"Could not find Gateway with ID {gateway_identifier}"
            )
        return self.gateways[gateway_identifier]

    def create_gateway(
        self,
        name: str,
        role_arn: str,
        protocol_type: str,
        authorizer_type: str,
        description: str | None,
        protocol_configuration: dict[str, Any] | None,
        authorizer_configuration: dict[str, Any] | None,
        kms_key_arn: str | None,
        interceptor_configurations: list[dict[str, Any]] | None,
        policy_engine_configuration: dict[str, Any] | None,
        exception_level: str | None,
        tags: dict[str, str] | None,
    ) -> Gateway:
        gateway = Gateway(
            region_name=self.region_name,
            account_id=self.account_id,
            name=name,
            role_arn=role_arn,
            protocol_type=protocol_type,
            authorizer_type=authorizer_type,
            description=description,
            protocol_configuration=protocol_configuration,
            authorizer_configuration=authorizer_configuration,
            kms_key_arn=kms_key_arn,
            interceptor_configurations=interceptor_configurations,
            policy_engine_configuration=policy_engine_configuration,
            exception_level=exception_level,
        )
        self.gateways[gateway.gateway_id] = gateway
        if tags:
            self.tagger.tag_resource(
                gateway.gateway_arn,
                [{"Key": k, "Value": v} for k, v in tags.items()],
            )
        return gateway

    def get_gateway(self, gateway_identifier: str) -> Gateway:
        gateway = self._get_gateway(gateway_identifier)
        gateway.advance()
        return gateway

    def update_gateway(
        self,
        gateway_identifier: str,
        name: str,
        role_arn: str,
        protocol_type: str,
        authorizer_type: str,
        description: str | None,
        protocol_configuration: dict[str, Any] | None,
        authorizer_configuration: dict[str, Any] | None,
        kms_key_arn: str | None,
        interceptor_configurations: list[dict[str, Any]] | None,
        policy_engine_configuration: dict[str, Any] | None,
        exception_level: str | None,
    ) -> Gateway:
        gateway = self._get_gateway(gateway_identifier)
        gateway.update(
            name=name,
            role_arn=role_arn,
            protocol_type=protocol_type,
            authorizer_type=authorizer_type,
            description=description,
            protocol_configuration=protocol_configuration,
            authorizer_configuration=authorizer_configuration,
            kms_key_arn=kms_key_arn,
            interceptor_configurations=interceptor_configurations,
            policy_engine_configuration=policy_engine_configuration,
            exception_level=exception_level,
        )
        return gateway

    def delete_gateway(self, gateway_identifier: str) -> Gateway:
        gateway = self._get_gateway(gateway_identifier)
        keys_to_remove = [k for k in self.gateway_targets if k[0] == gateway.gateway_id]
        for key in keys_to_remove:
            self.gateway_targets.pop(key)
        self.gateways.pop(gateway.gateway_id)
        return gateway

    def list_gateways(self) -> list[Gateway]:
        return list(self.gateways.values())

    def _get_gateway_target(self, gateway_id: str, target_id: str) -> GatewayTarget:
        key = (gateway_id, target_id)
        if key not in self.gateway_targets:
            raise ResourceNotFoundException(
                f"Could not find Target {target_id} for Gateway {gateway_id}"
            )
        return self.gateway_targets[key]

    def create_gateway_target(
        self,
        gateway_identifier: str,
        name: str,
        target_configuration: dict[str, Any],
        description: str | None,
        credential_provider_configurations: list[dict[str, Any]] | None,
        metadata_configuration: dict[str, Any] | None,
    ) -> GatewayTarget:
        gateway = self._get_gateway(gateway_identifier)
        target = GatewayTarget(
            region_name=self.region_name,
            account_id=self.account_id,
            gateway=gateway,
            name=name,
            target_configuration=target_configuration,
            description=description,
            credential_provider_configurations=credential_provider_configurations,
            metadata_configuration=metadata_configuration,
        )
        self.gateway_targets[(gateway.gateway_id, target.target_id)] = target
        return target

    def get_gateway_target(
        self, gateway_identifier: str, target_id: str
    ) -> GatewayTarget:
        gateway = self._get_gateway(gateway_identifier)
        target = self._get_gateway_target(gateway.gateway_id, target_id)
        target.advance()
        return target

    def update_gateway_target(
        self,
        gateway_identifier: str,
        target_id: str,
        name: str,
        target_configuration: dict[str, Any],
        description: str | None,
        credential_provider_configurations: list[dict[str, Any]] | None,
        metadata_configuration: dict[str, Any] | None,
    ) -> GatewayTarget:
        gateway = self._get_gateway(gateway_identifier)
        target = self._get_gateway_target(gateway.gateway_id, target_id)
        target.update(
            name=name,
            target_configuration=target_configuration,
            description=description,
            credential_provider_configurations=credential_provider_configurations,
            metadata_configuration=metadata_configuration,
        )
        return target

    def delete_gateway_target(
        self, gateway_identifier: str, target_id: str
    ) -> GatewayTarget:
        gateway = self._get_gateway(gateway_identifier)
        target = self._get_gateway_target(gateway.gateway_id, target_id)
        self.gateway_targets.pop((gateway.gateway_id, target_id))
        return target

    def list_gateway_targets(self, gateway_identifier: str) -> list[GatewayTarget]:
        gateway = self._get_gateway(gateway_identifier)
        return [
            t
            for (gid, _), t in self.gateway_targets.items()
            if gid == gateway.gateway_id
        ]

    def _get_memory(self, memory_id: str) -> Memory:
        if memory_id not in self.memories:
            raise ResourceNotFoundException(
                f"Could not find Memory with ID {memory_id}"
            )
        return self.memories[memory_id]

    def create_memory(
        self,
        name: str,
        event_expiry_duration: int,
        description: str | None,
        encryption_key_arn: str | None,
        memory_execution_role_arn: str | None,
        memory_strategies: list[dict[str, Any]] | None,
        tags: dict[str, str] | None,
    ) -> Memory:
        memory = Memory(
            region_name=self.region_name,
            account_id=self.account_id,
            name=name,
            event_expiry_duration=event_expiry_duration,
            description=description,
            encryption_key_arn=encryption_key_arn,
            memory_execution_role_arn=memory_execution_role_arn,
            memory_strategies=memory_strategies,
        )
        self.memories[memory.memory_id] = memory
        if tags:
            self.tagger.tag_resource(
                memory.memory_arn,
                [{"Key": k, "Value": v} for k, v in tags.items()],
            )
        return memory

    def get_memory(self, memory_id: str) -> Memory:
        memory = self._get_memory(memory_id)
        memory.advance()
        return memory

    def update_memory(
        self,
        memory_id: str,
        description: str | None,
        event_expiry_duration: int | None,
        memory_execution_role_arn: str | None,
        memory_strategies: dict[str, Any] | None,
    ) -> Memory:
        memory = self._get_memory(memory_id)
        memory.update(
            description=description,
            event_expiry_duration=event_expiry_duration,
            memory_execution_role_arn=memory_execution_role_arn,
            memory_strategies=memory_strategies,
        )
        return memory

    def delete_memory(self, memory_id: str) -> Memory:
        memory = self._get_memory(memory_id)
        self.memories.pop(memory_id)
        return memory

    def list_memories(self) -> list[Memory]:
        return list(self.memories.values())

    def tag_resource(self, resource_arn: str, tags: dict[str, str]) -> None:
        self.tagger.tag_resource(
            resource_arn, [{"Key": k, "Value": v} for k, v in tags.items()]
        )

    def untag_resource(self, resource_arn: str, tag_keys: list[str]) -> None:
        self.tagger.untag_resource_using_names(resource_arn, tag_keys)

    def list_tags_for_resource(self, resource_arn: str) -> dict[str, str]:
        return self.tagger.get_tag_dict_for_resource(resource_arn)


bedrockagentcorecontrol_backends = BackendDict(
    BedrockAgentCoreControlBackend,
    "bedrock-agentcore-control",
    # botocore does not yet include endpoint data for this service
    use_boto3_regions=False,
    additional_regions=[
        "us-east-1",
        "us-east-2",
        "us-west-1",
        "us-west-2",
        "ap-south-1",
        "ap-northeast-1",
        "ap-northeast-2",
        "ap-southeast-1",
        "ap-southeast-2",
        "ca-central-1",
        "eu-central-1",
        "eu-west-1",
        "eu-west-2",
        "eu-west-3",
        "eu-north-1",
        "sa-east-1",
    ],
)
