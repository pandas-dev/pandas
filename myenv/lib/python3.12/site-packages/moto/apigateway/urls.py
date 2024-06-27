from ..apigatewayv2.urls import url_paths as url_paths_v2
from .responses import APIGatewayResponse

url_bases = [r"https?://apigateway\.(.+)\.amazonaws.com"]

url_paths = {
    "{0}/restapis$": APIGatewayResponse.dispatch,
    "{0}/restapis/(?P<function_id>[^/]+)/?$": APIGatewayResponse.dispatch,
    "{0}/restapis/(?P<function_id>[^/]+)/resources$": APIGatewayResponse.dispatch,
    "{0}/restapis/(?P<function_id>[^/]+)/authorizers$": APIGatewayResponse.dispatch,
    "{0}/restapis/(?P<function_id>[^/]+)/authorizers/(?P<authorizer_id>[^/]+)/?$": APIGatewayResponse.dispatch,
    "{0}/restapis/(?P<function_id>[^/]+)/stages$": APIGatewayResponse.dispatch,
    "{0}/tags/(?P<resourceArn>[^/]+)$": APIGatewayResponse.dispatch,
    "{0}/tags/arn:(?P<partition>[^/]+):apigateway:(?P<region_name>[^/]+)::/restapis/(?P<function_id>[^/]+)/stages/(?P<stage_name>[^/]+)$": APIGatewayResponse.dispatch,
    "{0}/restapis/(?P<function_id>[^/]+)/stages/(?P<stage_name>[^/]+)/?$": APIGatewayResponse.dispatch,
    "{0}/restapis/(?P<function_id>[^/]+)/stages/(?P<stage_name>[^/]+)/exports/(?P<export_type>[^/]+)/?$": APIGatewayResponse.dispatch,
    "{0}/restapis/(?P<function_id>[^/]+)/deployments$": APIGatewayResponse.dispatch,
    "{0}/restapis/(?P<function_id>[^/]+)/deployments/(?P<deployment_id>[^/]+)/?$": APIGatewayResponse.dispatch,
    "{0}/restapis/(?P<function_id>[^/]+)/resources/(?P<resource_id>[^/]+)/?$": APIGatewayResponse.dispatch,
    "{0}/restapis/(?P<function_id>[^/]+)/resources/(?P<resource_id>[^/]+)/methods/(?P<method_name>[^/]+)/?$": APIGatewayResponse.dispatch,
    r"{0}/restapis/(?P<function_id>[^/]+)/resources/(?P<resource_id>[^/]+)/methods/(?P<method_name>[^/]+)/responses/(?P<status_code>\d+)$": APIGatewayResponse.dispatch,
    r"{0}/restapis/(?P<function_id>[^/]+)/resources/(?P<resource_id>[^/]+)/methods/(?P<method_name>[^/]+)/integration$": APIGatewayResponse.dispatch,
    r"{0}/restapis/(?P<function_id>[^/]+)/resources/(?P<resource_id>[^/]+)/methods/(?P<method_name>[^/]+)/integration/responses/(?P<status_code>\d+)$": APIGatewayResponse.dispatch,
    r"{0}/restapis/(?P<function_id>[^/]+)/resources/(?P<resource_id>[^/]+)/methods/(?P<method_name>[^/]+)/integration/responses/(?P<status_code>\d+)/$": APIGatewayResponse.dispatch,
    "{0}/apikeys$": APIGatewayResponse.dispatch,
    "{0}/apikeys/(?P<apikey>[^/]+)": APIGatewayResponse.dispatch,
    "{0}/usageplans$": APIGatewayResponse.dispatch,
    "{0}/domainnames$": APIGatewayResponse.dispatch,
    "{0}/restapis/(?P<function_id>[^/]+)/models$": APIGatewayResponse.dispatch,
    "{0}/restapis/(?P<function_id>[^/]+)/models/(?P<model_name>[^/]+)/?$": APIGatewayResponse.dispatch,
    "{0}/domainnames/(?P<domain_name>[^/]+)/?$": APIGatewayResponse.dispatch,
    "{0}/domainnames/(?P<domain_name>[^/]+)/basepathmappings$": APIGatewayResponse.dispatch,
    "{0}/domainnames/(?P<domain_name>[^/]+)/basepathmappings/(?P<base_path_mapping>[^/]+)$": APIGatewayResponse.dispatch,
    "{0}/usageplans/(?P<usage_plan_id>[^/]+)/?$": APIGatewayResponse.dispatch,
    "{0}/usageplans/(?P<usage_plan_id>[^/]+)/keys$": APIGatewayResponse.dispatch,
    "{0}/usageplans/(?P<usage_plan_id>[^/]+)/keys/(?P<api_key_id>[^/]+)/?$": APIGatewayResponse.dispatch,
    "{0}/restapis/(?P<function_id>[^/]+)/requestvalidators$": APIGatewayResponse.dispatch,
    "{0}/restapis/(?P<api_id>[^/]+)/requestvalidators/(?P<validator_id>[^/]+)/?$": APIGatewayResponse.dispatch,
    "{0}/restapis/(?P<api_id>[^/]+)/gatewayresponses/?$": APIGatewayResponse.dispatch,
    "{0}/restapis/(?P<api_id>[^/]+)/gatewayresponses/(?P<response_type>[^/]+)/?$": APIGatewayResponse.dispatch,
    "{0}/vpclinks$": APIGatewayResponse.dispatch,
    "{0}/vpclinks/(?P<vpclink_id>[^/]+)": APIGatewayResponse.dispatch,
}

# Also manages the APIGatewayV2
url_paths.update(url_paths_v2)  # type: ignore[arg-type]
