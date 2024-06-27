"""signer base URL and path."""

from .responses import SignerResponse

url_bases = [
    r"https?://signer\.(.+)\.amazonaws\.com",
]


url_paths = {
    "{0}/tags/(?P<profile_arn>[^/]+)$": SignerResponse.dispatch,
    "{0}/tags/(?P<arn_prefix>[^/]+)/signing-profiles/(?P<profile_name>[^/]+)$": SignerResponse.method_dispatch(
        SignerResponse.tags  # type: ignore
    ),
    "{0}/signing-profiles/(?P<profile_name>[^/]+)$": SignerResponse.dispatch,
    "{0}/signing-platforms$": SignerResponse.dispatch,
}
