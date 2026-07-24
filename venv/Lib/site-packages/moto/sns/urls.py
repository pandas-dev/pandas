from .responses import SNSResponse

url_bases = [r"https?://sns\.(.+)\.amazonaws\.com"]

url_paths = {
    "{0}/$": SNSResponse.dispatch,
    "{0}/SimpleNotificationService-(?P<pem_id>[^/]+).pem": SNSResponse.serve_pem,
}
