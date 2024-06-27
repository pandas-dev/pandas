"""apigatewaymanagementapi base URL and path."""

from .responses import ApiGatewayManagementApiResponse

url_bases = [r"https?://execute-api\.(.+)\.amazonaws\.com"]


response = ApiGatewayManagementApiResponse()


url_paths = {
    "{0}/@connections/(?P<connection_id>[^/]+)$": response.dispatch,
}
