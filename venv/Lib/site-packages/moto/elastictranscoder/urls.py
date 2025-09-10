from .responses import ElasticTranscoderResponse

url_bases = [
    r"https?://elastictranscoder\.(.+)\.amazonaws.com",
]


url_paths = {
    r"{0}/(?P<api_version>[^/]+)/pipelines/?$": ElasticTranscoderResponse.dispatch,
    r"{0}/(?P<api_version>[^/]+)/pipelines/(?P<pipeline_id>[^/]+)/?$": ElasticTranscoderResponse.dispatch,
}
