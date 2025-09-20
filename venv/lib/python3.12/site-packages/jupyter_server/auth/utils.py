"""A module with various utility methods for authorization in Jupyter Server."""

# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.
import importlib
import random
import re
import warnings


def warn_disabled_authorization():
    """DEPRECATED, does nothing"""
    warnings.warn(
        "jupyter_server.auth.utils.warn_disabled_authorization is deprecated",
        DeprecationWarning,
        stacklevel=2,
    )


HTTP_METHOD_TO_AUTH_ACTION = {
    "GET": "read",
    "HEAD": "read",
    "OPTIONS": "read",
    "POST": "write",
    "PUT": "write",
    "PATCH": "write",
    "DELETE": "write",
    "WEBSOCKET": "execute",
}


def get_regex_to_resource_map():
    """Returns a dictionary with all of Jupyter Server's
    request handler URL regex patterns mapped to
    their resource name.

    e.g.
    { "/api/contents/<regex_pattern>": "contents", ...}
    """
    from jupyter_server.serverapp import JUPYTER_SERVICE_HANDLERS

    modules = []
    for mod_name in JUPYTER_SERVICE_HANDLERS.values():
        if mod_name:
            modules.extend(mod_name)
    resource_map = {}
    for handler_module in modules:
        mod = importlib.import_module(handler_module)
        name = mod.AUTH_RESOURCE
        for handler in mod.default_handlers:
            url_regex = handler[0]
            resource_map[url_regex] = name
    # terminal plugin doesn't have importable url patterns
    # get these from terminal/__init__.py
    for url_regex in [
        r"/terminals/websocket/(\w+)",
        "/api/terminals",
        r"/api/terminals/(\w+)",
    ]:
        resource_map[url_regex] = "terminals"
    return resource_map


def match_url_to_resource(url, regex_mapping=None):
    """Finds the JupyterHandler regex pattern that would
    match the given URL and returns the resource name (str)
    of that handler.

    e.g.
    /api/contents/... returns "contents"
    """
    if not regex_mapping:
        regex_mapping = get_regex_to_resource_map()
    for regex, auth_resource in regex_mapping.items():
        pattern = re.compile(regex)
        if pattern.fullmatch(url):
            return auth_resource


# From https://en.wikipedia.org/wiki/Moons_of_Jupiter
moons_of_jupyter = [
    "Metis",
    "Adrastea",
    "Amalthea",
    "Thebe",
    "Io",
    "Europa",
    "Ganymede",
    "Callisto",
    "Themisto",
    "Leda",
    "Ersa",
    "Pandia",
    "Himalia",
    "Lysithea",
    "Elara",
    "Dia",
    "Carpo",
    "Valetudo",
    "Euporie",
    "Eupheme",
    # 'S/2003 J 18',
    # 'S/2010 J 2',
    "Helike",
    # 'S/2003 J 16',
    # 'S/2003 J 2',
    "Euanthe",
    # 'S/2017 J 7',
    "Hermippe",
    "Praxidike",
    "Thyone",
    "Thelxinoe",
    # 'S/2017 J 3',
    "Ananke",
    "Mneme",
    # 'S/2016 J 1',
    "Orthosie",
    "Harpalyke",
    "Iocaste",
    # 'S/2017 J 9',
    # 'S/2003 J 12',
    # 'S/2003 J 4',
    "Erinome",
    "Aitne",
    "Herse",
    "Taygete",
    # 'S/2017 J 2',
    # 'S/2017 J 6',
    "Eukelade",
    "Carme",
    # 'S/2003 J 19',
    "Isonoe",
    # 'S/2003 J 10',
    "Autonoe",
    "Philophrosyne",
    "Cyllene",
    "Pasithee",
    # 'S/2010 J 1',
    "Pasiphae",
    "Sponde",
    # 'S/2017 J 8',
    "Eurydome",
    # 'S/2017 J 5',
    "Kalyke",
    "Hegemone",
    "Kale",
    "Kallichore",
    # 'S/2011 J 1',
    # 'S/2017 J 1',
    "Chaldene",
    "Arche",
    "Eirene",
    "Kore",
    # 'S/2011 J 2',
    # 'S/2003 J 9',
    "Megaclite",
    "Aoede",
    # 'S/2003 J 23',
    "Callirrhoe",
    "Sinope",
]


def get_anonymous_username() -> str:
    """
    Get a random user-name based on the moons of Jupyter.
    This function returns names like "Anonymous Io" or "Anonymous Metis".
    """
    return moons_of_jupyter[random.randint(0, len(moons_of_jupyter) - 1)]  # noqa: S311
