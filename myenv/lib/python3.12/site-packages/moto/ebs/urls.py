"""ebs base URL and path."""

from .responses import EBSResponse

url_bases = [r"https?://ebs\.(.+)\.amazonaws\.com"]


url_paths = {
    "{0}/snapshots$": EBSResponse.method_dispatch(EBSResponse.snapshots),
    "{0}/snapshots/completion/(?P<snapshot_id>[^/]+)$": EBSResponse.method_dispatch(
        EBSResponse.complete_snapshot
    ),
    "{0}/snapshots/(?P<snapshot_id>[^/]+)/changedblocks$": EBSResponse.method_dispatch(
        EBSResponse.snapshot_changed_blocks
    ),
    "{0}/snapshots/(?P<snapshot_id>[^/]+)/blocks$": EBSResponse.method_dispatch(
        EBSResponse.snapshot_blocks
    ),
    "{0}/snapshots/(?P<snapshot_id>[^/]+)/blocks/(?P<block_idx>[^/]+)$": EBSResponse.method_dispatch(
        EBSResponse.snapshot_block
    ),
}
