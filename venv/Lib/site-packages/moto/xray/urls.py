from .responses import XRayResponse

url_bases = [r"https?://xray\.(.+)\.amazonaws.com"]

url_paths = {
    "{0}/TelemetryRecords$": XRayResponse.dispatch,
    "{0}/TraceSegments$": XRayResponse.dispatch,
    "{0}/Traces$": XRayResponse.dispatch,
    "{0}/ServiceGraph$": XRayResponse.dispatch,
    "{0}/TraceGraph$": XRayResponse.dispatch,
    "{0}/TraceSummaries$": XRayResponse.dispatch,
}
