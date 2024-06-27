from .responses import ElasticTranscoderResponse

url_bases = [
    r"https?://elastictranscoder\.(.+)\.amazonaws.com",
]


url_paths = {
    r"{0}/(?P<api_version>[^/]+)/pipelines/?$": ElasticTranscoderResponse.method_dispatch(
        ElasticTranscoderResponse.pipelines
    ),
    r"{0}/(?P<api_version>[^/]+)/pipelines/(?P<pipeline_id>[^/]+)/?$": ElasticTranscoderResponse.method_dispatch(
        ElasticTranscoderResponse.individual_pipeline
    ),
}
