"""kinesisanalyticsv2 base URL and path."""

from .responses import KinesisAnalyticsV2Response

url_bases = [
    r"https?://kinesisanalytics\.(.+)\.amazonaws\.com",
]

url_paths = {
    "{0}/$": KinesisAnalyticsV2Response.dispatch,
}
