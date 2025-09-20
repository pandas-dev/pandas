from .responses import SQSResponse

url_bases = [r"https?://(.*\.)?(queue|sqs)\.(.*\.)?amazonaws\.com"]

dispatch = SQSResponse().dispatch

url_paths = {
    "{0}/$": dispatch,
    r"{0}/(?P<account_id>\d+)/(?P<queue_name>[a-zA-Z0-9\-_\.]+)": dispatch,
}
