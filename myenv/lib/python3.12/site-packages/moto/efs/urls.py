from .responses import EFSResponse

url_bases = [
    r"https?://elasticfilesystem\.(.+)\.amazonaws.com",
    r"https?://elasticfilesystem\.amazonaws.com",
]


response = EFSResponse()


url_paths = {
    "{0}/.*?$": response.dispatch,
    "/2015-02-01/access-points": response.dispatch,
    "/2015-02-01/access-points/<access_point_id>": response.dispatch,
    "/2015-02-01/file-systems": response.dispatch,
    "/2015-02-01/file-systems/<file_system_id>": response.dispatch,
    "/2015-02-01/file-systems/<file_system_id>/backup-policy": response.dispatch,
    "/2015-02-01/file-systems/<file_system_id>/lifecycle-configuration": response.dispatch,
    "/2015-02-01/mount-targets": response.dispatch,
    "/2015-02-01/mount-targets/<mount_target_id>": response.dispatch,
    "/2015-02-01/mount-targets/<mount_target_id>/security-groups": response.dispatch,
    "/2015-02-01/resource-tags/<resource_id>": response.dispatch,
}
