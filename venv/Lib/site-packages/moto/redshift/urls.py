from .responses import RedshiftResponse

url_bases = [r"https?://redshift\.(.+)\.amazonaws\.com"]

url_paths = {"{0}/$": RedshiftResponse.dispatch}
