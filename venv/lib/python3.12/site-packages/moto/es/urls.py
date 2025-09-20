from moto.opensearch.responses import OpenSearchServiceResponse

url_bases = [
    r"https?://es\.(.+)\.amazonaws\.com",
]


url_paths = {
    "{0}/2015-01-01/es/domain/(?P<domainname>[^/]+)/config": OpenSearchServiceResponse.describe_es_domain_config,
    "{0}/2015-01-01/domain$": OpenSearchServiceResponse.list_domains,
    "{0}/2015-01-01/es/domain$": OpenSearchServiceResponse.domains,
    "{0}/2015-01-01/es/domain/(?P<domainname>[^/]+)": OpenSearchServiceResponse.domain,
    "{0}/2015-01-01/es/domain-info$": OpenSearchServiceResponse.list_domains,
    "{0}/2015-01-01/tags/?$": OpenSearchServiceResponse.tags,
    "{0}/2015-01-01/tags-removal/?": OpenSearchServiceResponse.tag_removal,
    "{0}/2021-01-01/domain$": OpenSearchServiceResponse.dispatch,
    "{0}/2021-01-01/opensearch/compatibleVersions": OpenSearchServiceResponse.dispatch,
    "{0}/2021-01-01/opensearch/domain": OpenSearchServiceResponse.dispatch,
    "{0}/2021-01-01/opensearch/domain/(?P<domainname>[^/]+)": OpenSearchServiceResponse.dispatch,
    "{0}/2021-01-01/opensearch/domain/(?P<domainname>[^/]+)/config": OpenSearchServiceResponse.dispatch,
    "{0}/2021-01-01/opensearch/domain-info$": OpenSearchServiceResponse.dispatch,
    "{0}/2021-01-01/tags/?$": OpenSearchServiceResponse.dispatch,
    "{0}/2021-01-01/tags-removal/?": OpenSearchServiceResponse.dispatch,
}
