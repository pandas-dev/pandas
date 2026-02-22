"""lexv2models base URL and path."""

from .responses import LexModelsV2Response

url_bases = [
    r"https?://lex\.(.+)\.amazonaws\.com",
    r"https?://models-v2-lex\.(.+)\.amazonaws\.com",
]

url_paths = {
    "{0}/bots/$": LexModelsV2Response.dispatch,
    "{0}/bots/(?P<botId>[^/]+)/$": LexModelsV2Response.dispatch,
    "{0}/bots/(?P<botId>[^/]+)/botaliases/$": LexModelsV2Response.dispatch,
    "{0}/bots/(?P<botId>[^/]+)/botaliases/(?P<botAliasId>[^/]+)/$": LexModelsV2Response.dispatch,
    "{0}/policy/(?P<resourceArn>.+)/$": LexModelsV2Response.dispatch,
    "{0}/tags/(?P<resourceARN>[^/]+)$": LexModelsV2Response.dispatch,
    "{0}/tags/(?P<arn_prefix>[^/]+)/(?P<bot_id>[^/]+)$": LexModelsV2Response.dispatch,
}
