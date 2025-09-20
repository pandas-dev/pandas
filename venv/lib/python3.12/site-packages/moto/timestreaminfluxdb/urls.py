"""timestreaminfluxdb base URL and path."""

from .responses import TimestreamInfluxDBResponse

url_bases = [
    r"https?://timestream-influxdb\.(.+)\.amazonaws\.com",
]

url_paths = {
    "{0}/$": TimestreamInfluxDBResponse.dispatch,
}
