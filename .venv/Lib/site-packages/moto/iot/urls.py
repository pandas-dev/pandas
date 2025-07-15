from .responses import IoTResponse

url_bases = [r"https?://iot\.(.+)\.amazonaws\.com"]


url_paths = {
    #
    # Paths for :class:`moto.core.models.MockAWS`
    #
    # This route requires special handling.
    "{0}/attached-policies/(?P<target>.*)$": IoTResponse.dispatch_attached_policies,
    # The remaining routes can be handled by the default dispatcher.
    "{0}/.*$": IoTResponse.dispatch,
    #
    # (Flask) Paths for :class:`moto.core.models.ServerModeMockAWS`
    #
    # This route requires special handling.
    "{0}/attached-policies/<path:target>$": IoTResponse.dispatch_attached_policies,
    # The remaining routes can be handled by the default dispatcher.
    "{0}/<path:route>$": IoTResponse.dispatch,
}
