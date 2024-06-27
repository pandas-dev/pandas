"""Tornado handlers for listing extensions."""

# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.
from __future__ import annotations

import json
from logging import Logger

import requests
import tornado
from jupyter_server.base.handlers import APIHandler

LISTINGS_URL_SUFFIX = "@jupyterlab/extensionmanager-extension/listings.json"


def fetch_listings(logger: Logger | None) -> None:
    """Fetch the listings for the extension manager."""
    if not logger:
        from traitlets import log

        logger = log.get_logger()  # type:ignore[assignment]
    assert logger is not None
    if len(ListingsHandler.blocked_extensions_uris) > 0:
        blocked_extensions = []
        for blocked_extensions_uri in ListingsHandler.blocked_extensions_uris:
            logger.info(
                "Fetching blocked_extensions from %s", ListingsHandler.blocked_extensions_uris
            )
            r = requests.request(
                "GET", blocked_extensions_uri, **ListingsHandler.listings_request_opts
            )
            j = json.loads(r.text)
            for b in j["blocked_extensions"]:
                blocked_extensions.append(b)
            ListingsHandler.blocked_extensions = blocked_extensions
    if len(ListingsHandler.allowed_extensions_uris) > 0:
        allowed_extensions = []
        for allowed_extensions_uri in ListingsHandler.allowed_extensions_uris:
            logger.info(
                "Fetching allowed_extensions from %s", ListingsHandler.allowed_extensions_uris
            )
            r = requests.request(
                "GET", allowed_extensions_uri, **ListingsHandler.listings_request_opts
            )
            j = json.loads(r.text)
            for w in j["allowed_extensions"]:
                allowed_extensions.append(w)
        ListingsHandler.allowed_extensions = allowed_extensions
    ListingsHandler.listings = json.dumps(  # type:ignore[attr-defined]
        {
            "blocked_extensions_uris": list(ListingsHandler.blocked_extensions_uris),
            "allowed_extensions_uris": list(ListingsHandler.allowed_extensions_uris),
            "blocked_extensions": ListingsHandler.blocked_extensions,
            "allowed_extensions": ListingsHandler.allowed_extensions,
        }
    )


class ListingsHandler(APIHandler):
    """An handler that returns the listings specs."""

    """Below fields are class level fields that are accessed and populated
    by the initialization and the fetch_listings methods.
    Some fields are initialized before the handler creation in the
    handlers.py#add_handlers method.
    Having those fields predefined reduces the guards in the methods using
    them.
    """
    # The list of blocked_extensions URIS.
    blocked_extensions_uris: set = set()
    # The list of allowed_extensions URIS.
    allowed_extensions_uris: set = set()
    # The blocked extensions extensions.
    blocked_extensions: list = []
    # The allowed extensions extensions.
    allowed_extensions: list = []
    # The provider request options to be used for the request library.
    listings_request_opts: dict = {}
    # The callback time for the periodic callback in seconds.
    listings_refresh_seconds: int
    # The PeriodicCallback that schedule the call to fetch_listings method.
    pc = None

    def get(self, path: str) -> None:
        """Get the listings for the extension manager."""
        self.set_header("Content-Type", "application/json")
        if path == LISTINGS_URL_SUFFIX:
            self.write(ListingsHandler.listings)  # type:ignore[attr-defined]
        else:
            raise tornado.web.HTTPError(400)
