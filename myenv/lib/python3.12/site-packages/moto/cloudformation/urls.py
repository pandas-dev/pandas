from .responses import CloudFormationResponse

url_bases = [r"https?://cloudformation\.(.+)\.amazonaws\.com"]

url_paths = {
    "{0}/$": CloudFormationResponse.dispatch,
    "{0}/cloudformation_(?P<region>[^/]+)/cfnresponse$": CloudFormationResponse.cfnresponse,
}
