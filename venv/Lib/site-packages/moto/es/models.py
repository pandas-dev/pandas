from moto.core.base_backend import BackendDict
from moto.opensearch.models import OpenSearchServiceBackend


class ElasticsearchServiceBackend(OpenSearchServiceBackend):
    def create_elasticsearch_domain(self) -> None:
        # Functionality is part of OpenSearch, as that includes all of ES functionality + more
        # Method is kept here to make sure our documentation register this as supported
        pass

    def delete_elasticsearch_domain(self) -> None:
        # Functionality is part of OpenSearch, as that includes all of ES functionality + more
        # Method is kept here to make sure our documentation register this as supported
        pass

    def describe_elasticsearch_domain(self) -> None:
        # Functionality is part of OpenSearch, as that includes all of ES functionality + more
        # Method is kept here to make sure our documentation register this as supported
        pass

    def describe_elasticsearch_domains(self) -> None:
        # Functionality is part of OpenSearch, as that includes all of ES functionality + more
        # Method is kept here to make sure our documentation register this as supported
        pass


es_backends = BackendDict(ElasticsearchServiceBackend, "es")
