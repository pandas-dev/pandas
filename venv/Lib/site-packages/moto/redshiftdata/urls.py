from .responses import RedshiftDataAPIServiceResponse

url_bases = [
    r"https?://redshift-data\.(.+)\.amazonaws\.com",
]


url_paths = {
    "{0}/$": RedshiftDataAPIServiceResponse.dispatch,
}
