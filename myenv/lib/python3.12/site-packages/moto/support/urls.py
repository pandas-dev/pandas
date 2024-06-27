from .responses import SupportResponse

url_bases = [r"https?://support\.(.+)\.amazonaws\.com"]


url_paths = {
    "{0}/$": SupportResponse.dispatch,
}
