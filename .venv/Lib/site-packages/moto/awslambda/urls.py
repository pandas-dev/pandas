from .responses import LambdaResponse

url_bases = [r"https?://lambda\.(.+)\.amazonaws\.com"]


url_paths = {
    r"{0}/(?P<api_version>[^/]+)/functions$": LambdaResponse.dispatch,
    r"{0}/(?P<api_version>[^/]+)/functions/$": LambdaResponse.dispatch,
    r"{0}/(?P<api_version>[^/]+)/functions/(?P<function_name>[\w_:%-]+)/?$": LambdaResponse.dispatch,
    r"{0}/(?P<api_version>[^/]+)/functions/(?P<function_name>[\w_:%-]+)/aliases$": LambdaResponse.dispatch,
    r"{0}/(?P<api_version>[^/]+)/functions/(?P<function_name>[\w_:%-]+)/aliases/(?P<alias_name>[\w_-]+)$": LambdaResponse.dispatch,
    r"{0}/(?P<api_version>[^/]+)/functions/(?P<function_name>[\w_:%-]+)/versions/?$": LambdaResponse.dispatch,
    r"{0}/(?P<api_version>[^/]+)/event-source-mappings/$": LambdaResponse.dispatch,
    r"{0}/(?P<api_version>[^/]+)/event-source-mappings/(?P<UUID>[\w_-]+)/?$": LambdaResponse.dispatch,
    r"{0}/(?P<api_version>[^/]+)/functions/(?P<function_name>[\w_-]+)/invocations/?$": LambdaResponse.dispatch,
    r"{0}/(?P<api_version>[^/]+)/functions/(?P<resource_arn>.+)/invocations/?$": LambdaResponse.dispatch,
    r"{0}/(?P<api_version>[^/]+)/functions/(?P<function_name>[\w_:%-]+)/invoke-async$": LambdaResponse.dispatch,
    r"{0}/(?P<api_version>[^/]+)/functions/(?P<function_name>[\w_:%-]+)/invoke-async/$": LambdaResponse.dispatch,
    r"{0}/(?P<api_version>[^/]+)/tags/(?P<resource_arn>.+)": LambdaResponse.dispatch,
    r"{0}/(?P<api_version>[^/]+)/functions/(?P<function_name>[\w_:%-]+)/policy/(?P<statement_id>[\w_-]+)$": LambdaResponse.dispatch,
    r"{0}/(?P<api_version>[^/]+)/functions/(?P<function_name>[\w_:%-]+)/policy/?$": LambdaResponse.dispatch,
    r"{0}/(?P<api_version>[^/]+)/functions/(?P<function_name>[\w_:%-]+)/configuration/?$": LambdaResponse.dispatch,
    r"{0}/(?P<api_version>[^/]+)/functions/(?P<function_name>[\w_:%-]+)/code/?$": LambdaResponse.dispatch,
    r"{0}/(?P<api_version>[^/]+)/functions/(?P<function_name>[\w_:%-]+)/code-signing-config$": LambdaResponse.dispatch,
    r"{0}/(?P<api_version>[^/]+)/functions/(?P<function_name>[\w_:%-]+)/concurrency/?$": LambdaResponse.dispatch,
    r"{0}/(?P<api_version>[^/]+)/functions/(?P<function_name>[\w_:%-]+)/url/?$": LambdaResponse.dispatch,
    r"{0}/(?P<api_version>[^/]+)/layers$": LambdaResponse.dispatch,
    r"{0}/(?P<api_version>[^/]+)/layers/$": LambdaResponse.dispatch,
    r"{0}/(?P<api_version>[^/]+)/layers/(?P<layer_name>.+)/versions$": LambdaResponse.dispatch,
    r"{0}/(?P<api_version>[^/]+)/layers/(?P<layer_name>.+)/versions/$": LambdaResponse.dispatch,
    r"{0}/(?P<api_version>[^/]+)/layers/(?P<layer_name>.+)/versions/(?P<layer_version>[\w_-]+)$": LambdaResponse.dispatch,
    r"{0}/(?P<api_version>[^/]+)/functions/(?P<function_name>[\w_:%-]+)/event-invoke-config/?$": LambdaResponse.dispatch,
    r"{0}/(?P<api_version>[^/]+)/functions/(?P<function_name>[\w_:%-]+)/event-invoke-config/list$": LambdaResponse.dispatch,
}
