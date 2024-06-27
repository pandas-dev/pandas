from moto.opensearch.responses import OpenSearchServiceResponse

from .responses import ElasticsearchServiceResponse

url_bases = [
    r"https?://es\.(.+)\.amazonaws\.com",
]


url_paths = {
    "{0}/2015-01-01/domain$": ElasticsearchServiceResponse.list_domains,
    "{0}/2015-01-01/es/domain$": ElasticsearchServiceResponse.domains,
    "{0}/2015-01-01/es/domain/(?P<domainname>[^/]+)": ElasticsearchServiceResponse.domain,
    "{0}/2021-01-01/domain$": OpenSearchServiceResponse.dispatch,
    "{0}/2021-01-01/opensearch/compatibleVersions": OpenSearchServiceResponse.dispatch,
    "{0}/2021-01-01/opensearch/domain": OpenSearchServiceResponse.dispatch,
    "{0}/2021-01-01/opensearch/domain/(?P<domainname>[^/]+)": OpenSearchServiceResponse.dispatch,
    "{0}/2021-01-01/opensearch/domain/(?P<domainname>[^/]+)/config": OpenSearchServiceResponse.dispatch,
    "{0}/2021-01-01/tags/?": OpenSearchServiceResponse.dispatch,
    "{0}/2021-01-01/tags-removal/": OpenSearchServiceResponse.dispatch,
}
