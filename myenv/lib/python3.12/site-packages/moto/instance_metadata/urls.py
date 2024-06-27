from .responses import InstanceMetadataResponse

url_bases = ["http://169.254.169.254"]

url_paths = {
    "{0}/(?P<path>.+)": InstanceMetadataResponse.method_dispatch(
        InstanceMetadataResponse.metadata_response
    )
}
