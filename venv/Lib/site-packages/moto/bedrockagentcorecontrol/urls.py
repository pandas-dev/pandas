"""BedrockAgentCoreControl URLs."""

from .responses import BedrockAgentCoreControlResponse

url_bases = [
    r"https?://bedrock-agentcore-control\.(.+)\.amazonaws\.com",
]

url_paths = {
    "{0}/.*$": BedrockAgentCoreControlResponse.dispatch,
}
