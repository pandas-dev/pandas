"""ebs base URL and path."""

from .responses import EBSResponse

url_bases = [r"https?://ebs\.(.+)\.amazonaws\.com"]


url_paths = {
    "{0}/snapshots$": EBSResponse.dispatch,
    "{0}/snapshots/completion/(?P<snapshot_id>[^/]+)$": EBSResponse.dispatch,
    "{0}/snapshots/(?P<snapshot_id>[^/]+)/changedblocks$": EBSResponse.dispatch,
    "{0}/snapshots/(?P<snapshot_id>[^/]+)/blocks$": EBSResponse.dispatch,
    "{0}/snapshots/(?P<snapshot_id>[^/]+)/blocks/(?P<block_idx>[^/]+)$": EBSResponse.dispatch,
}
