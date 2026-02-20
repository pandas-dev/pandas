"""osis base URL and path."""

from .responses import OpenSearchIngestionResponse

url_bases = [
    r"https?://osis\.(.+)\.amazonaws\.com",
]

url_paths = {
    "{0}/2022-01-01/osis/createPipeline$": OpenSearchIngestionResponse.dispatch,
    "{0}/2022-01-01/osis/deletePipeline/(?P<PipelineName>[^/]+)$": OpenSearchIngestionResponse.dispatch,
    "{0}/2022-01-01/osis/getPipeline/(?P<PipelineName>[^/]+)$": OpenSearchIngestionResponse.dispatch,
    "{0}/2022-01-01/osis/startPipeline/(?P<PipelineName>[^/]+)$": OpenSearchIngestionResponse.dispatch,
    "{0}/2022-01-01/osis/stopPipeline/(?P<PipelineName>[^/]+)$": OpenSearchIngestionResponse.dispatch,
    "{0}/2022-01-01/osis/listPipelines$": OpenSearchIngestionResponse.dispatch,
    "{0}/2022-01-01/osis/listTagsForResource/$": OpenSearchIngestionResponse.dispatch,
    "{0}/2022-01-01/osis/updatePipeline/(?P<PipelineName>[^/]+)$": OpenSearchIngestionResponse.dispatch,
    "{0}/2022-01-01/osis/tagResource/$": OpenSearchIngestionResponse.dispatch,
    "{0}/2022-01-01/osis/untagResource/$": OpenSearchIngestionResponse.dispatch,
}
