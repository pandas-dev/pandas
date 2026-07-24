from .responses import DataBrewResponse

url_bases = [r"https?://databrew\.(.+)\.amazonaws.com"]

url_paths = {
    "{0}/recipeVersions$": DataBrewResponse.dispatch,
    "{0}/recipes$": DataBrewResponse.dispatch,
    "{0}/recipes/(?P<recipe_name>[^/]+)$": DataBrewResponse.dispatch,
    "{0}/recipes/(?P<recipe_name>[^/]+)/recipeVersion/(?P<recipe_version>[^/]+)": DataBrewResponse.dispatch,
    "{0}/recipes/(?P<recipe_name>[^/]+)/publishRecipe$": DataBrewResponse.dispatch,
    "{0}/rulesets$": DataBrewResponse.dispatch,
    "{0}/rulesets/(?P<ruleset_name>[^/]+)$": DataBrewResponse.dispatch,
    "{0}/datasets$": DataBrewResponse.dispatch,
    "{0}/datasets/(?P<dataset_name>[^/]+)$": DataBrewResponse.dispatch,
    "{0}/jobs$": DataBrewResponse.dispatch,
    "{0}/jobs/(?P<job_name>[^/]+)$": DataBrewResponse.dispatch,
    "{0}/profileJobs$": DataBrewResponse.dispatch,
    "{0}/recipeJobs$": DataBrewResponse.dispatch,
    "{0}/profileJobs/(?P<job_name>[^/]+)$": DataBrewResponse.dispatch,
    "{0}/recipeJobs/(?P<job_name>[^/]+)$": DataBrewResponse.dispatch,
}
