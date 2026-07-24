"""apigatewaymanagementapi base URL and path."""

from .responses import ApiGatewayManagementApiResponse

# execute-api.us-east-1.amazonaws.com
# api_id.execute-api.us-east-1.amazonaws.com
url_bases = [r"https?://([^.]+\.)*execute-api\.[^.]+\.amazonaws\.com"]


response = ApiGatewayManagementApiResponse()


url_paths = {
    "{0}/@connections/(?P<connection_id>[^/]+)$": response.dispatch,
    "{0}/(?P<stage_name>.+/)+@connections/(?P<connection_id>[^/]+)$": response.connect_to_apigateway,
}
