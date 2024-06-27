from .responses import ManagedBlockchainResponse

url_bases = [r"https?://managedblockchain\.(.+)\.amazonaws.com"]

url_paths = {
    "{0}/networks$": ManagedBlockchainResponse.dispatch,
    "{0}/networks/(?P<networkid>[^/.]+)$": ManagedBlockchainResponse.dispatch,
    "{0}/networks/(?P<networkid>[^/.]+)/proposals$": ManagedBlockchainResponse.dispatch,
    "{0}/networks/(?P<networkid>[^/.]+)/proposals/(?P<proposalid>[^/.]+)$": ManagedBlockchainResponse.dispatch,
    "{0}/networks/(?P<networkid>[^/.]+)/proposals/(?P<proposalid>[^/.]+)/votes$": ManagedBlockchainResponse.dispatch,
    "{0}/invitations$": ManagedBlockchainResponse.dispatch,
    "{0}/invitations/(?P<invitationid>[^/.]+)$": ManagedBlockchainResponse.dispatch,
    "{0}/networks/(?P<networkid>[^/.]+)/members$": ManagedBlockchainResponse.dispatch,
    "{0}/networks/(?P<networkid>[^/.]+)/members/(?P<memberid>[^/.]+)$": ManagedBlockchainResponse.dispatch,
    "{0}/networks/(?P<networkid>[^/.]+)/members/(?P<memberid>[^/.]+)/nodes$": ManagedBlockchainResponse.dispatch,
    "{0}/networks/(?P<networkid>[^/.]+)/members/(?P<memberid>[^/.]+)/nodes?(?P<querys>[^/.]+)$": ManagedBlockchainResponse.dispatch,
    "{0}/networks/(?P<networkid>[^/.]+)/members/(?P<memberid>[^/.]+)/nodes/(?P<nodeid>[^/.]+)$": ManagedBlockchainResponse.dispatch,
    # >= botocore 1.19.41 (API change - memberId is now part of query-string or body)
    "{0}/networks/(?P<networkid>[^/.]+)/nodes$": ManagedBlockchainResponse.dispatch,
    "{0}/networks/(?P<networkid>[^/.]+)/nodes/(?P<nodeid>[^/.]+)$": ManagedBlockchainResponse.dispatch,
}
