"""appmesh base URL and path."""

from .responses import AppMeshResponse

url_bases = [
    r"https?://appmesh\.(.+)\.amazonaws\.com",
]

url_paths = {
    "{0}/v20190125/meshes$": AppMeshResponse.dispatch,
    "{0}/v20190125/meshes/(?P<meshName>[^/]+)$": AppMeshResponse.dispatch,
    "{0}/v20190125/tags$": AppMeshResponse.dispatch,
    "{0}/v20190125/tag$": AppMeshResponse.dispatch,
    "{0}/v20190125/meshes/(?P<meshName>.*)/virtualRouters/(?P<virtualRouterName>[^/]+)$": AppMeshResponse.dispatch,
    "{0}/v20190125/meshes/(?P<meshName>.*)/virtualRouters$": AppMeshResponse.dispatch,
    "{0}/v20190125/meshes/(?P<meshName>.*)/virtualRouter/(?P<virtualRouterName>.*)/routes$": AppMeshResponse.dispatch,
    "{0}/v20190125/meshes/(?P<meshName>.*)/virtualRouter/(?P<virtualRouterName>.*)/routes/(?P<routeName>[^/]+)$": AppMeshResponse.dispatch,
    "{0}/v20190125/meshes/(?P<meshName>.*)/virtualNodes/(?P<virtualNodeName>[^/]+)$": AppMeshResponse.dispatch,
    "{0}/v20190125/meshes/(?P<meshName>.*)/virtualNodes$": AppMeshResponse.dispatch,
}
