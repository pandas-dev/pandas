from urllib.parse import urlparse


def region_from_applicationautoscaling_url(url: str) -> str:
    domain = urlparse(url).netloc

    if "." in domain:
        return domain.split(".")[1]
    else:
        return "us-east-1"
