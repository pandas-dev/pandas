"""Handles incoming bedrockagent requests, invokes methods, returns responses."""

import json
from urllib.parse import unquote

from moto.core.responses import BaseResponse

from .models import AgentsforBedrockBackend, bedrockagent_backends


class AgentsforBedrockResponse(BaseResponse):
    """Handler for AgentsforBedrock requests and responses."""

    def __init__(self) -> None:
        super().__init__(service_name="bedrock-agent")

    @property
    def bedrockagent_backend(self) -> AgentsforBedrockBackend:
        """Return backend instance specific for this region."""
        return bedrockagent_backends[self.current_account][self.region]

    def create_agent(self) -> str:
        params = json.loads(self.body)
        agent_name = params.get("agentName")
        client_token = params.get("clientToken")
        instruction = params.get("instruction")
        foundation_model = params.get("foundationModel")
        description = params.get("description")
        idle_session_ttl_in_seconds = params.get("idleSessionTTLInSeconds")
        agent_resource_role_arn = params.get("agentResourceRoleArn")
        customer_encryption_key_arn = params.get("customerEncryptionKeyArn")
        tags = params.get("tags")
        prompt_override_configuration = params.get("promptOverrideConfiguration")
        agent = self.bedrockagent_backend.create_agent(
            agent_name=agent_name,
            client_token=client_token,
            instruction=instruction,
            foundation_model=foundation_model,
            description=description,
            idle_session_ttl_in_seconds=idle_session_ttl_in_seconds,
            agent_resource_role_arn=agent_resource_role_arn,
            customer_encryption_key_arn=customer_encryption_key_arn,
            tags=tags,
            prompt_override_configuration=prompt_override_configuration,
        )
        return json.dumps({"agent": dict(agent.to_dict())})

    def get_agent(self) -> str:
        agent_id = self.path.split("/")[-2]
        agent = self.bedrockagent_backend.get_agent(agent_id=agent_id)
        return json.dumps({"agent": dict(agent.to_dict())})

    def list_agents(self) -> str:
        params = json.loads(self.body)
        max_results = params.get("maxResults")
        next_token = params.get("nextToken")
        max_results = int(max_results) if max_results else None
        agents, next_token = self.bedrockagent_backend.list_agents(
            max_results=max_results,
            next_token=next_token,
        )
        return json.dumps(
            {
                "agentSummaries": agents,
                "nextToken": next_token,
            }
        )

    def delete_agent(self) -> str:
        params = self._get_params()
        skip_resource_in_use_check = params.get("skipResourceInUseCheck")
        agent_id = self.path.split("/")[-2]
        agent_id, agent_status = self.bedrockagent_backend.delete_agent(
            agent_id=agent_id, skip_resource_in_use_check=skip_resource_in_use_check
        )
        return json.dumps({"agentId": agent_id, "agentStatus": agent_status})

    def create_knowledge_base(self) -> str:
        params = json.loads(self.body)
        client_token = params.get("clientToken")
        name = params.get("name")
        description = params.get("description")
        role_arn = params.get("roleArn")
        knowledge_base_configuration = params.get("knowledgeBaseConfiguration")
        storage_configuration = params.get("storageConfiguration")
        tags = params.get("tags")
        knowledge_base = self.bedrockagent_backend.create_knowledge_base(
            client_token=client_token,
            name=name,
            description=description,
            role_arn=role_arn,
            knowledge_base_configuration=knowledge_base_configuration,
            storage_configuration=storage_configuration,
            tags=tags,
        )
        return json.dumps({"knowledgeBase": dict(knowledge_base.to_dict())})

    def list_knowledge_bases(self) -> str:
        params = json.loads(self.body)
        max_results = params.get("maxResults")
        next_token = params.get("nextToken")
        max_results = int(max_results) if max_results else None
        knowledge_bases, next_token = self.bedrockagent_backend.list_knowledge_bases(
            max_results=max_results,
            next_token=next_token,
        )
        return json.dumps(
            {
                "knowledgeBaseSummaries": knowledge_bases,
                "nextToken": next_token,
            }
        )

    def delete_knowledge_base(self) -> str:
        knowledge_base_id = self.path.split("/")[-1]
        (
            knowledge_base_id,
            knowledge_base_status,
        ) = self.bedrockagent_backend.delete_knowledge_base(
            knowledge_base_id=knowledge_base_id
        )
        return json.dumps(
            {"knowledgeBaseId": knowledge_base_id, "status": knowledge_base_status}
        )

    def get_knowledge_base(self) -> str:
        knowledge_base_id = self.path.split("/")[-1]
        knowledge_base = self.bedrockagent_backend.get_knowledge_base(
            knowledge_base_id=knowledge_base_id
        )
        return json.dumps({"knowledgeBase": knowledge_base.to_dict()})

    def tag_resource(self) -> str:
        params = json.loads(self.body)
        resource_arn = unquote(self.path.split("/tags/")[-1])
        tags = params.get("tags")
        self.bedrockagent_backend.tag_resource(resource_arn=resource_arn, tags=tags)
        return json.dumps(dict())

    def untag_resource(self) -> str:
        resource_arn = unquote(self.path.split("/tags/")[-1])
        tag_keys = self.querystring.get("tagKeys", [])
        self.bedrockagent_backend.untag_resource(
            resource_arn=resource_arn, tag_keys=tag_keys
        )
        return json.dumps(dict())

    def list_tags_for_resource(self) -> str:
        resource_arn = unquote(self.path.split("/tags/")[-1])
        tags = self.bedrockagent_backend.list_tags_for_resource(
            resource_arn=resource_arn
        )
        return json.dumps(dict(tags=tags))
