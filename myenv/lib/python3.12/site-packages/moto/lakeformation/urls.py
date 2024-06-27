"""lakeformation base URL and path."""

from .responses import LakeFormationResponse

url_bases = [
    r"https?://lakeformation\.(.+)\.amazonaws\.com",
]


response = LakeFormationResponse()


url_paths = {
    "{0}/DescribeResource$": response.dispatch,
    "{0}/DeregisterResource$": response.dispatch,
    "{0}/RegisterResource$": response.dispatch,
    "{0}/ListResources$": response.dispatch,
    "{0}/GetDataLakeSettings$": response.dispatch,
    "{0}/PutDataLakeSettings$": response.dispatch,
    "{0}/GrantPermissions$": response.dispatch,
    "{0}/ListPermissions$": response.dispatch,
    "{0}/RevokePermissions$": response.dispatch,
    "{0}/CreateLFTag$": response.dispatch,
    "{0}/GetLFTag$": response.dispatch,
    "{0}/DeleteLFTag$": response.dispatch,
    "{0}/UpdateLFTag": response.dispatch,
    "{0}/ListLFTags$": response.dispatch,
    "{0}/AddLFTagsToResource": response.dispatch,
    "{0}/RemoveLFTagsFromResource": response.dispatch,
    "{0}/GetResourceLFTags": response.dispatch,
    "{0}/ListDataCellsFilter$": response.dispatch,
    "{0}/BatchGrantPermissions$": response.dispatch,
    "{0}/BatchRevokePermissions$": response.dispatch,
}
