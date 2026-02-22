from .responses import DatabaseMigrationServiceResponse

url_bases = [r"https?://dms\.(.+)\.amazonaws\.com"]


url_paths = {
    "{0}/$": DatabaseMigrationServiceResponse.dispatch,
}
