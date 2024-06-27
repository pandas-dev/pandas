"""apigatewayv2 base URL and path."""

from .responses import ApiGatewayV2Response

url_bases = [
    r"https?://apigateway\.(.+)\.amazonaws\.com",
]


url_paths = {
    "{0}/v2/apis$": ApiGatewayV2Response.dispatch,
    "{0}/v2/apis/(?P<api_id>[^/]+)$": ApiGatewayV2Response.dispatch,
    "{0}/v2/apis/(?P<api_id>[^/]+)/authorizers$": ApiGatewayV2Response.dispatch,
    "{0}/v2/apis/(?P<api_id>[^/]+)/authorizers/(?P<authorizer_id>[^/]+)$": ApiGatewayV2Response.dispatch,
    "{0}/v2/apis/(?P<api_id>[^/]+)/cors$": ApiGatewayV2Response.dispatch,
    "{0}/v2/apis/(?P<api_id>[^/]+)/integrations$": ApiGatewayV2Response.dispatch,
    "{0}/v2/apis/(?P<api_id>[^/]+)/integrations/(?P<integration_id>[^/]+)$": ApiGatewayV2Response.dispatch,
    "{0}/v2/apis/(?P<api_id>[^/]+)/integrations/(?P<integration_id>[^/]+)/integrationresponses$": ApiGatewayV2Response.dispatch,
    "{0}/v2/apis/(?P<api_id>[^/]+)/integrations/(?P<integration_id>[^/]+)/integrationresponses/(?P<integration_response_id>[^/]+)$": ApiGatewayV2Response.dispatch,
    "{0}/v2/apis/(?P<api_id>[^/]+)/models$": ApiGatewayV2Response.dispatch,
    "{0}/v2/apis/(?P<api_id>[^/]+)/models/(?P<model_id>[^/]+)$": ApiGatewayV2Response.dispatch,
    "{0}/v2/apis/(?P<api_id>[^/]+)/routes$": ApiGatewayV2Response.dispatch,
    "{0}/v2/apis/(?P<api_id>[^/]+)/routes/(?P<route_id>[^/]+)$": ApiGatewayV2Response.dispatch,
    "{0}/v2/apis/(?P<api_id>[^/]+)/routes/(?P<route_id>[^/]+)/routeresponses$": ApiGatewayV2Response.dispatch,
    "{0}/v2/apis/(?P<api_id>[^/]+)/routes/(?P<route_id>[^/]+)/routeresponses/(?P<route_response_id>[^/]+)$": ApiGatewayV2Response.dispatch,
    "{0}/v2/apis/(?P<api_id>[^/]+)/routes/(?P<route_id>[^/]+)/requestparameters/(?P<request_parameter>[^/]+)$": ApiGatewayV2Response.dispatch,
    "{0}/v2/apis/(?P<api_id>[^/]+)/stages$": ApiGatewayV2Response.dispatch,
    "{0}/v2/apis/(?P<api_id>[^/]+)/stages/(?P<stage_name>[^/]+)$": ApiGatewayV2Response.dispatch,
    "{0}/v2/tags/(?P<resource_arn>[^/]+)$": ApiGatewayV2Response.dispatch,
    "{0}/v2/tags/(?P<resource_arn_pt1>[^/]+)/apis/(?P<resource_arn_pt2>[^/]+)$": ApiGatewayV2Response.dispatch,
    "{0}/v2/tags/(?P<resource_arn_pt1>[^/]+)/vpclinks/(?P<resource_arn_pt2>[^/]+)$": ApiGatewayV2Response.dispatch,
    "{0}/v2/vpclinks$": ApiGatewayV2Response.dispatch,
    "{0}/v2/vpclinks/(?P<vpc_link_id>[^/]+)$": ApiGatewayV2Response.dispatch,
    "{0}/v2/domainnames$": ApiGatewayV2Response.dispatch,
    "{0}/v2/domainnames/(?P<domain_name>[^/]+)$": ApiGatewayV2Response.dispatch,
    "{0}/v2/domainnames/(?P<domain_name>[^/]+)/apimappings$": ApiGatewayV2Response.dispatch,
    "{0}/v2/domainnames/(?P<domain_name>[^/]+)/apimappings/(?P<api_mapping_id>[^/]+)$": ApiGatewayV2Response.dispatch,
}
