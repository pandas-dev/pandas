"""Announcements handler for JupyterLab."""

# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.

import abc
import hashlib
import json
import xml.etree.ElementTree as ET
from collections.abc import Awaitable
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from typing import Optional, Union

from jupyter_server.base.handlers import APIHandler
from jupyterlab_server.translation_utils import translator
from packaging.version import parse
from tornado import httpclient, web

from jupyterlab._version import __version__

ISO8601_FORMAT = "%Y-%m-%dT%H:%M:%S%z"
JUPYTERLAB_LAST_RELEASE_URL = "https://pypi.org/pypi/jupyterlab/json"
JUPYTERLAB_RELEASE_URL = "https://github.com/jupyterlab/jupyterlab/releases/tag/v"


def format_datetime(dt_str: str):
    return datetime.fromisoformat(dt_str).timestamp() * 1000


@dataclass(frozen=True)
class Notification:
    """Notification

    Attributes:
        createdAt: Creation date
        message: Notification message
        modifiedAt: Modification date
        type: Notification type — ["default", "error", "info", "success", "warning"]
        link: Notification link button as a tuple (label, URL)
        options: Notification options
    """

    createdAt: float  # noqa
    message: str
    modifiedAt: float  # noqa
    type: str = "default"
    link: tuple[str, str] = field(default_factory=tuple)
    options: dict = field(default_factory=dict)


class CheckForUpdateABC(abc.ABC):
    """Abstract class to check for update.

    Args:
        version: Current JupyterLab version

    Attributes:
        version - str: Current JupyterLab version
        logger - logging.Logger: Server logger
    """

    def __init__(self, version: str) -> None:
        self.version = version

    @abc.abstractmethod
    async def __call__(self) -> Awaitable[Union[None, str, tuple[str, tuple[str, str]]]]:
        """Get the notification message if a new version is available.

        Returns:
            None if there is not update.
            or the notification message
            or the notification message and a tuple(label, URL link) for the user to get more information
        """
        msg = "CheckForUpdateABC.__call__ is not implemented"
        raise NotImplementedError(msg)


class CheckForUpdate(CheckForUpdateABC):
    """Default class to check for update.

    Args:
        version: Current JupyterLab version

    Attributes:
        version - str: Current JupyterLab version
        logger - logging.Logger: Server logger
    """

    async def __call__(self) -> Awaitable[tuple[str, tuple[str, str]]]:
        """Get the notification message if a new version is available.

        Returns:
            None if there is no update.
            or the notification message
            or the notification message and a tuple(label, URL link) for the user to get more information
        """
        http_client = httpclient.AsyncHTTPClient()
        try:
            response = await http_client.fetch(
                JUPYTERLAB_LAST_RELEASE_URL,
                headers={"Content-Type": "application/json"},
            )
            data = json.loads(response.body).get("info")
            last_version = data["version"]
        except Exception as e:
            self.logger.debug("Failed to get latest version", exc_info=e)
            return None
        else:
            if parse(self.version) < parse(last_version):
                trans = translator.load("jupyterlab")
                return (
                    trans.__(f"A newer version ({last_version}) of JupyterLab is available."),
                    (trans.__("Read more…"), f"{JUPYTERLAB_RELEASE_URL}{last_version}"),
                )
            else:
                return None


class NeverCheckForUpdate(CheckForUpdateABC):
    """Check update version that does nothing.

    This is provided for administrators that want to
    turn off requesting external resources.

    Args:
        version: Current JupyterLab version

    Attributes:
        version - str: Current JupyterLab version
        logger - logging.Logger: Server logger
    """

    async def __call__(self) -> Awaitable[None]:
        """Get the notification message if a new version is available.

        Returns:
            None if there is no update.
            or the notification message
            or the notification message and a tuple(label, URL link) for the user to get more information
        """
        return None


class CheckForUpdateHandler(APIHandler):
    """Check for Updates API handler.

    Args:
        update_check: The class checking for a new version
    """

    def initialize(
        self,
        update_checker: Optional[CheckForUpdate] = None,
    ) -> None:
        super().initialize()
        self.update_checker = (
            NeverCheckForUpdate(__version__) if update_checker is None else update_checker
        )
        self.update_checker.logger = self.log

    @web.authenticated
    async def get(self):
        """Check for updates.
        Response:
            {
                "notification": Optional[Notification]
            }
        """
        notification = None
        out = await self.update_checker()
        if out:
            message, link = (out, ()) if isinstance(out, str) else out
            now = datetime.now(tz=timezone.utc).timestamp() * 1000.0
            hash_ = hashlib.sha1(message.encode()).hexdigest()  # noqa: S324
            notification = Notification(
                message=message,
                createdAt=now,
                modifiedAt=now,
                type="info",
                link=link,
                options={"data": {"id": hash_, "tags": ["update"]}},
            )

        self.set_status(200)
        self.finish(
            json.dumps({"notification": None if notification is None else asdict(notification)})
        )


class NewsHandler(APIHandler):
    """News API handler.

    Args:
        news_url: The Atom feed to fetch for news
    """

    def initialize(
        self,
        news_url: Optional[str] = None,
    ) -> None:
        super().initialize()
        self.news_url = news_url

    @web.authenticated
    async def get(self):
        """Get the news.

        Response:
            {
                "news": List[Notification]
            }
        """
        news = []

        http_client = httpclient.AsyncHTTPClient()

        if self.news_url is not None:
            trans = translator.load("jupyterlab")

            # Those registrations are global, naming them to reduce chance of clashes
            xml_namespaces = {"atom": "http://www.w3.org/2005/Atom"}
            for key, spec in xml_namespaces.items():
                ET.register_namespace(key, spec)

            try:
                response = await http_client.fetch(
                    self.news_url,
                    headers={"Content-Type": "application/atom+xml"},
                )
                tree = ET.fromstring(response.body)  # noqa S314

                def build_entry(node):
                    def get_xml_text(attr: str, default: Optional[str] = None) -> str:
                        node_item = node.find(f"atom:{attr}", xml_namespaces)
                        if node_item is not None:
                            return node_item.text
                        elif default is not None:
                            return default
                        else:
                            error_m = (
                                f"atom feed entry does not contain a required attribute: {attr}"
                            )
                            raise KeyError(error_m)

                    entry_title = get_xml_text("title")
                    entry_id = get_xml_text("id")
                    entry_updated = get_xml_text("updated")
                    entry_published = get_xml_text("published", entry_updated)
                    entry_summary = get_xml_text("summary", default="")
                    links = node.findall("atom:link", xml_namespaces)
                    if len(links) > 1:
                        alternate = list(filter(lambda elem: elem.get("rel") == "alternate", links))
                        link_node = alternate[0] if alternate else links[0]
                    else:
                        link_node = links[0] if len(links) == 1 else None
                    entry_link = link_node.get("href") if link_node is not None else None

                    message = (
                        "\n".join([entry_title, entry_summary]) if entry_summary else entry_title
                    )
                    modified_at = format_datetime(entry_updated)
                    created_at = format_datetime(entry_published)
                    notification = Notification(
                        message=message,
                        createdAt=created_at,
                        modifiedAt=modified_at,
                        type="info",
                        link=None
                        if entry_link is None
                        else (
                            trans.__("Open full post"),
                            entry_link,
                        ),
                        options={
                            "data": {
                                "id": entry_id,
                                "tags": ["news"],
                            }
                        },
                    )
                    return notification

                entries = map(build_entry, tree.findall("atom:entry", xml_namespaces))
                news.extend(entries)
            except Exception as e:
                self.log.debug(
                    f"Failed to get announcements from Atom feed: {self.news_url}",
                    exc_info=e,
                )

        self.set_status(200)
        self.finish(json.dumps({"news": list(map(asdict, news))}))


news_handler_path = r"/lab/api/news"
check_update_handler_path = r"/lab/api/update"
