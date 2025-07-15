"""inspector2 base URL and path."""

from .responses import Inspector2Response

url_bases = [
    r"https?://inspector2\.(.+)\.amazonaws\.com",
]


url_paths = {
    "{0}/delegatedadminaccounts/disable$": Inspector2Response.dispatch,
    "{0}/delegatedadminaccounts/enable$": Inspector2Response.dispatch,
    "{0}/delegatedadminaccounts/list$": Inspector2Response.dispatch,
    "{0}/disable$": Inspector2Response.dispatch,
    "{0}/enable$": Inspector2Response.dispatch,
    "{0}/filters/create$": Inspector2Response.dispatch,
    "{0}/filters/delete$": Inspector2Response.dispatch,
    "{0}/filters/list$": Inspector2Response.dispatch,
    "{0}/findings/list$": Inspector2Response.dispatch,
    "{0}/members/associate$": Inspector2Response.dispatch,
    "{0}/members/get": Inspector2Response.dispatch,
    "{0}/members/list": Inspector2Response.dispatch,
    "{0}/members/disassociate$": Inspector2Response.dispatch,
    "{0}/status/batch/get$": Inspector2Response.dispatch,
    "{0}/organizationconfiguration/describe$": Inspector2Response.dispatch,
    "{0}/organizationconfiguration/update$": Inspector2Response.dispatch,
    "{0}/tags/(?P<resource_arn>.+)$": Inspector2Response.dispatch,
}
