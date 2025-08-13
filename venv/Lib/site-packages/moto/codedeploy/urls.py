"""codedeploy base URL and path."""

from .responses import CodeDeployResponse

url_bases = [
    r"https?://codedeploy\.(.+)\.amazonaws\.com",
]

url_paths = {
    "{0}/$": CodeDeployResponse.dispatch,
    "{0}/list-tags-for-resource$": CodeDeployResponse.list_tags_for_resource,
    "{0}/tag-resource$": CodeDeployResponse.tag_resource,
    "{0}/untag-resource$": CodeDeployResponse.untag_resource,
}
