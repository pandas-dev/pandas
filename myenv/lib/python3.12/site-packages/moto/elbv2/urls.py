from ..elb.urls import api_version_elb_backend

url_bases = [r"https?://elasticloadbalancing\.(.+)\.amazonaws.com"]

url_paths = {"{0}/$": api_version_elb_backend}
