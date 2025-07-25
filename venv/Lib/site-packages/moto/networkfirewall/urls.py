"""networkfirewall base URL and path."""

from .responses import NetworkFirewallResponse

url_bases = [
    r"https?://network-firewall\.(.+)\.amazonaws\.com",
]

url_paths = {
    "{0}/$": NetworkFirewallResponse.dispatch,
}
