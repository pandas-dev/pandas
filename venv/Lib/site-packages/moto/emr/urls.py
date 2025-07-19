from .responses import ElasticMapReduceResponse

url_bases = [
    r"https?://(.+)\.elasticmapreduce\.amazonaws.com",
    r"https?://elasticmapreduce\.(.+)\.amazonaws.com",
]

url_paths = {"{0}/$": ElasticMapReduceResponse.dispatch}
