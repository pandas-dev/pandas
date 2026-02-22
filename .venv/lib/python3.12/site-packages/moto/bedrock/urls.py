"""bedrock base URL and path."""

from ..bedrockagent.responses import AgentsforBedrockResponse
from .responses import BedrockResponse

url_bases = [
    r"https?://bedrock\.(.+)\.amazonaws\.com",
]

url_paths = {
    "{0}/.*$": BedrockResponse.dispatch,
    "{0}/agents/?$": AgentsforBedrockResponse.dispatch,
    "{0}/agents/(?P<agent_name>[^/]+)/$": AgentsforBedrockResponse.dispatch,
    "{0}/custom-models$": BedrockResponse.dispatch,
    "{0}/custom-models/(?P<modelIdentifier>[^/]+)$": BedrockResponse.dispatch,
    "{0}/custom-models/(?P<arn_prefix>[^/]+)/(?P<jobIdentifier>[^/]+)$": BedrockResponse.dispatch,
    "{0}/knowledgebases$": AgentsforBedrockResponse.dispatch,
    "{0}/knowledgebases/(?P<kb_name>[^/]+)$": AgentsforBedrockResponse.dispatch,
    "{0}/knowledgebases/(?P<kb_name>[^/]+)/$": AgentsforBedrockResponse.dispatch,
    "{0}/listTagsForResource$": BedrockResponse.dispatch,
    "{0}/logging/modelinvocations$": BedrockResponse.dispatch,
    "{0}/model-customization-jobs$": BedrockResponse.dispatch,
    "{0}/model-customization-jobs/(?P<jobIdentifier>[^/]+)$": BedrockResponse.dispatch,
    "{0}/model-customization-jobs/(?P<jobIdentifier>[^/]+)/stop$": BedrockResponse.dispatch,
    "{0}/tags/(?P<resource_arn>[^/]+)$": AgentsforBedrockResponse.dispatch,
    "{0}/tags/(?P<arn_prefix>[^/]+)/(?P<name>[^/]+)$": AgentsforBedrockResponse.dispatch,
    "{0}/tagResource$": BedrockResponse.dispatch,
    "{0}/untagResource$": BedrockResponse.dispatch,
}
