"""bedrockagent base URL and path."""

from .responses import AgentsforBedrockResponse

url_bases = [
    r"https?://bedrock-agent\.(.+)\.amazonaws\.com",
]

url_paths = {
    "{0}/.*$": AgentsforBedrockResponse.dispatch,
}
