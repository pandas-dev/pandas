###############################################################################
#
# Url - A class to represent URLs in Excel.
#
# SPDX-License-Identifier: BSD-2-Clause
#
# Copyright (c) 2013-2025, John McNamara, jmcnamara@cpan.org
#


import re
from enum import Enum
from typing import Optional


class UrlTypes(Enum):
    """
    Enum to represent different types of URLS.
    """

    UNKNOWN = 0
    URL = 1
    INTERNAL = 2
    EXTERNAL = 3


class Url:
    """
    A class to represent URLs in Excel.

    """

    MAX_URL_LEN = 2080
    MAX_PARAMETER_LEN = 255

    def __init__(self, link: str):
        self._link_type: UrlTypes = UrlTypes.UNKNOWN
        self._original_url: str = link
        self._link: str = link
        self._relationship_link: str = link
        self._text: str = ""
        self._tip: str = ""
        self._anchor: str = ""
        self._is_object_link: bool = False
        self._rel_index: int = 0

        self._parse_url()

        if len(self._link) > self.MAX_URL_LEN:
            raise ValueError("URL exceeds Excel's maximum length.")

        if len(self._anchor) > self.MAX_URL_LEN:
            raise ValueError("Anchor segment or url exceeds Excel's maximum length.")

        if len(self._tip) > self.MAX_PARAMETER_LEN:
            raise ValueError("Hyperlink tool tip exceeds Excel's maximum length.")

        self._escape_strings()

    def __repr__(self):
        """
        Return a string representation of the Url instance.

        """
        return (
            "\n"
            f"Url:\n"
            f"  _link_type         = {self._link_type.name}\n"
            f"  _original_url      = {self._original_url}\n"
            f"  _link              = {self._link}\n"
            f"  _relationship_link = {self._relationship_link}\n"
            f"  _text              = {self._text}\n"
            f"  _tip               = {self._tip}\n"
            f"  _anchor            = {self._anchor}\n"
            f"  _is_object_link    = {self._is_object_link}\n"
            f"  _rel_index         = {self._rel_index}\n"
        )

    @classmethod
    def from_options(cls, options: dict) -> Optional["Url"]:
        """
        For backward compatibility, convert the 'url' key and 'tip' keys in an
        options dictionary to a Url object, or return the Url object if already
        an instance.

        Args:
            options (dict): A dictionary that may contain a 'url' key.

        Returns:
            url: A Url object or None.

        """
        if not isinstance(options, dict):
            raise TypeError("The 'options' parameter must be a dictionary.")

        url = options.get("url")

        if isinstance(url, str):
            url = cls(options["url"])
            if options.get("tip"):
                url._tip = options["tip"]

        return url

    @property
    def text(self) -> str:
        """Get the alternative, user-friendly, text for the URL."""
        return self._text

    @text.setter
    def text(self, value: str):
        """Set the alternative, user-friendly, text for the URL."""
        self._text = value

    @property
    def tip(self) -> str:
        """Get the screen tip for the URL."""
        return self._tip

    @tip.setter
    def tip(self, value: str):
        """Set the screen tip for the URL."""
        self._tip = value

    def _parse_url(self):
        """Parse the URL and determine its type."""

        # Handle mail address links.
        if self._link.startswith("mailto:"):
            self._link_type = UrlTypes.URL

            if not self._text:
                self._text = self._link.replace("mailto:", "", 1)

        # Handle links to cells within the workbook.
        elif self._link.startswith("internal:"):
            self._link_type = UrlTypes.INTERNAL
            self._relationship_link = self._link.replace("internal:", "#", 1)
            self._link = self._link.replace("internal:", "", 1)
            self._anchor = self._link

            if not self._text:
                self._text = self._anchor

        # Handle links to other files or cells in other Excel files.
        elif self._link.startswith("file://") or self._link.startswith("external:"):
            self._link_type = UrlTypes.EXTERNAL

            # Handle backward compatibility with external: links.
            file_url = self._original_url.replace("external:", "file:///", 1)

            link_path = file_url
            link_path = link_path.replace("file:///", "", 1)
            link_path = link_path.replace("file://", "", 1)
            link_path = link_path.replace("/", "\\")

            if self._is_relative_path(link_path):
                self._link = link_path
            else:
                self._link = "file:///" + link_path

            if not self._text:
                self._text = link_path

            if "#" in self._link:
                self._link, self._anchor = self._link.split("#", 1)

            # Set up the relationship link. This doesn't usually contain the
            # anchor unless it is a link from an object like an image.
            if self._is_object_link:
                if self._is_relative_path(link_path):
                    self._relationship_link = self._link.replace("\\", "/")
                else:
                    self._relationship_link = file_url

            else:
                self._relationship_link = self._link

            # Convert a .\dir\file.xlsx link to dir\file.xlsx.
            if self._relationship_link.startswith(".\\"):
                self._relationship_link = self._relationship_link.replace(".\\", "", 1)

        # Handle standard Excel links like http://, https://, ftp://, ftps://
        # but also allow custom "foo://bar" URLs.
        elif "://" in self._link:
            self._link_type = UrlTypes.URL

            if not self._text:
                self._text = self._link

            if "#" in self._link:
                self._link, self._anchor = self._link.split("#", 1)

            # Set up the relationship link. This doesn't usually contain the
            # anchor unless it is a link from an object like an image.
            if self._is_object_link:
                self._relationship_link = self._original_url
            else:
                self._relationship_link = self._link

        else:
            raise ValueError(f"Unknown URL type: {self._original_url}")

    def _set_object_link(self):
        """
        Set the _is_object_link flag and re-parse the URL since the relationship
        link is different for object links.

        """
        self._is_object_link = True
        self._link = self._original_url
        self._parse_url()
        self._escape_strings()

    def _escape_strings(self):
        """Escape special characters in the URL strings."""

        if self._link_type != UrlTypes.INTERNAL:
            self._link = self._escape_url(self._link)
            self._relationship_link = self._escape_url(self._relationship_link)

        # Excel additionally escapes # to %23 in file paths.
        if self._link_type == UrlTypes.EXTERNAL:
            self._relationship_link = self._relationship_link.replace("#", "%23")

    def _target(self) -> str:
        """Get the target for relationship IDs."""
        return self._relationship_link

    def _target_mode(self) -> str:
        """Get the target mode for relationship IDs."""
        if self._link_type == UrlTypes.INTERNAL:
            return ""

        return "External"

    @staticmethod
    def _is_relative_path(url: str) -> bool:
        """Check if a URL is a relative path."""
        if url.startswith(r"\\"):
            return False

        if url[0].isalpha() and url[1] == ":":
            return False

        return True

    @staticmethod
    def _escape_url(url: str) -> str:
        """Escape special characters in a URL."""
        # Don't escape URL if it looks already escaped.
        if re.search("%[0-9a-fA-F]{2}", url):
            return url

        # Can't use url.quote() here because it doesn't match Excel.
        return (
            url.replace("%", "%25")
            .replace('"', "%22")
            .replace(" ", "%20")
            .replace("<", "%3c")
            .replace(">", "%3e")
            .replace("[", "%5b")
            .replace("]", "%5d")
            .replace("^", "%5e")
            .replace("`", "%60")
            .replace("{", "%7b")
            .replace("}", "%7d")
        )
