from .responses import EKSResponse

url_bases = [
    r"https?://eks\.(.+)\.amazonaws.com",
]


response = EKSResponse()


url_paths = {
    "{0}/clusters$": response.dispatch,
    "{0}/clusters/(?P<name>[^/]+)$": response.dispatch,
    "{0}/clusters/(?P<name>[^/]+)/node-groups$": response.dispatch,
    "{0}/clusters/(?P<name>[^/]+)/node-groups/(?P<nodegroupName>[^/]+)$": response.dispatch,
    "{0}/clusters/(?P<name>[^/]+)/fargate-profiles$": response.dispatch,
    "{0}/clusters/(?P<name>[^/]+)/fargate-profiles/(?P<fargateProfileName>[^/]+)$": response.dispatch,
    "{0}/tags/(?P<resourceArn>.+)$": response.tags,
}
