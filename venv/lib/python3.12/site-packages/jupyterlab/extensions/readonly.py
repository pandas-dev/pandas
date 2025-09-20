"""Extension manager without installation capabilities."""

# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.

import sys
from typing import Optional

from jupyterlab_server.translation_utils import translator

from .manager import ActionResult, ExtensionManager, ExtensionManagerMetadata, ExtensionPackage


class ReadOnlyExtensionManager(ExtensionManager):
    """Extension manager without installation capabilities."""

    @property
    def metadata(self) -> ExtensionManagerMetadata:
        """Extension manager metadata."""
        return ExtensionManagerMetadata("read-only", install_path=sys.prefix)

    async def get_latest_version(self, pkg: str) -> Optional[str]:
        """Return the latest available version for a given extension.

        Args:
            pkg: The extension to search for
        Returns:
            The latest available version
        """
        return None

    async def list_packages(
        self, query: str, page: int, per_page: int
    ) -> tuple[dict[str, ExtensionPackage], Optional[int]]:
        """List the available extensions.

        Args:
            query: The search extension query
            page: The result page
            per_page: The number of results per page
        Returns:
            The available extensions in a mapping {name: metadata}
            The results last page; None if the manager does not support pagination
        """
        return {}, None

    async def install(self, extension: str, version: Optional[str] = None) -> ActionResult:
        """Install the required extension.

        Note:
            If the user must be notified with a message (like asking to restart the
            server), the result should be
            {"status": "warning", "message": "<explanation for the user>"}

        Args:
            extension: The extension name
            version: The version to install; default None (i.e. the latest possible)
        Returns:
            The action result
        """
        trans = translator.load("jupyterlab")
        return ActionResult(
            status="error", message=trans.gettext("Extension installation not supported.")
        )

    async def uninstall(self, extension: str) -> ActionResult:
        """Uninstall the required extension.

        Note:
            If the user must be notified with a message (like asking to restart the
            server), the result should be
            {"status": "warning", "message": "<explanation for the user>"}

        Args:
            extension: The extension name
        Returns:
            The action result
        """
        trans = translator.load("jupyterlab")
        return ActionResult(
            status="error", message=trans.gettext("Extension removal not supported.")
        )
