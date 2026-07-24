"""bedrockruntime base URL and path."""

from .responses import BedrockRuntimeResponse

url_bases = [
    r"https?://bedrock-runtime\.(.+)\.amazonaws\.com",
]

url_paths = {
    "{0}/.*$": BedrockRuntimeResponse.dispatch,
}
