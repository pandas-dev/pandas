from __future__ import annotations

import os
from pathlib import Path
from typing import TYPE_CHECKING, Any, List, Optional

from rattler.rattler import PyAboutJson

if TYPE_CHECKING:
    from rattler.networking.client import Client


class AboutJson:
    """
    The `about.json` file contains metadata about the package.
    """

    _inner: PyAboutJson

    @staticmethod
    def from_path(path: os.PathLike[str]) -> AboutJson:
        """
        Parses the object from a file specified by a `path`, using a format
        appropriate for the file type.

        For example, if the file is in JSON format, this function reads the data
        from the file at the specified path, parse the JSON string and return the
        resulting object. If the file is not in a parsable format or if the file
        could not read, this function returns an error.

        Examples
        --------
        ```python
        >>> about = AboutJson.from_path("../test-data/dummy-about.json")
        >>> about
        AboutJson()
        >>>
        ```
        """
        return AboutJson._from_py_about_json(PyAboutJson.from_path(Path(path)))

    @staticmethod
    def from_package_directory(path: os.PathLike[str]) -> AboutJson:
        """
        Parses the object by looking up the appropriate file from the root of the
        specified Conda archive directory, using a format appropriate for the file
        type.

        For example, if the file is in JSON format, this function reads the
        appropriate file from the archive, parse the JSON string and return the
        resulting object. If the file is not in a parsable format or if the file
        could not be read, this function returns an error.
        """
        return AboutJson._from_py_about_json(PyAboutJson.from_package_directory(Path(path)))

    @staticmethod
    def from_str(string: str) -> AboutJson:
        """
        Parses the object from a string, using a format appropriate for the file
        type.

        For example, if the file is in JSON format, this function parses the JSON
        string and returns the resulting object. If the file is not in a parsable
        format, this function returns an error.

        Examples
        --------
        ```python
        >>> import json
        >>> with open("../test-data/dummy-about.json", 'r') as file:
        ...     json_str = json.dumps(json.load(file))
        >>> about = AboutJson.from_str(json_str)
        >>> about
        AboutJson()
        >>>
        ```
        """
        return AboutJson._from_py_about_json(PyAboutJson.from_str(string))

    @classmethod
    async def from_remote_url(cls, client: Client, url: str) -> Optional[AboutJson]:
        """
        Fetches `info/about.json` from a remote package archive URL.
        """
        py_about_json = await PyAboutJson.from_remote_url(client._client, url)
        if py_about_json is None:
            return None
        return cls._from_py_about_json(py_about_json)

    @staticmethod
    def package_path() -> Path:
        """
        Returns the path to the file within the Conda archive.

        The path is relative to the root of the archive and includes any necessary
        directories.

        Examples
        --------
        ```python
        >>> str(AboutJson.package_path())
        'info/about.json'
        >>>
        ```
        """
        return PyAboutJson.package_path()

    @property
    def channels(self) -> List[str]:
        """
        A list of channels that where used during the build.

        Examples
        --------
        ```python
        >>> about = AboutJson.from_path("../test-data/dummy-about.json")
        >>> about.channels
        ['https://conda.anaconda.org/conda-forge']
        >>> about.channels = ['https://test.channel']
        >>> about.channels
        ['https://test.channel']
        >>>
        ```
        """
        return self._inner.channels

    @channels.setter
    def channels(self, value: List[str]) -> None:
        self._inner.channels = value

    @property
    def description(self) -> Optional[str]:
        """
        Description of the package.

        Examples
        --------
        ```python
        >>> about = AboutJson.from_path("../test-data/dummy-about.json")
        >>> about.description
        'A dummy description.'
        >>> about.description = 'New description'
        >>> about.description
        'New description'
        >>>
        ```
        """
        if description := self._inner.description:
            return description

        return None

    @description.setter
    def description(self, value: Optional[str]) -> None:
        self._inner.description = value

    @property
    def dev_url(self) -> List[str]:
        """
        A list of URLs to the development page of the package.

        Examples
        --------
        ```python
        >>> about = AboutJson.from_path("../test-data/dummy-about.json")
        >>> about.dev_url
        ['https://github.com/conda/rattler']
        >>> about.dev_url = ['https://test.dev']
        >>> about.dev_url
        ['https://test.dev/']
        >>>
        ```
        """
        return self._inner.dev_url

    @dev_url.setter
    def dev_url(self, value: List[str]) -> None:
        self._inner.dev_url = value

    @property
    def doc_url(self) -> List[str]:
        """
        A list of URLs to the documentation of the package.

        Examples
        --------
        ```python
        >>> about = AboutJson.from_path("../test-data/dummy-about.json")
        >>> about.doc_url
        ['https://conda.github.io/rattler/py-rattler/']
        >>> about.doc_url = ['https://test.docs']
        >>> about.doc_url
        ['https://test.docs/']
        >>>
        ```
        """
        return self._inner.doc_url

    @doc_url.setter
    def doc_url(self, value: List[str]) -> None:
        self._inner.doc_url = value

    @property
    def home(self) -> List[str]:
        """
        A list URL to the homepage of the package.

        Examples
        --------
        ```python
        >>> about = AboutJson.from_path("../test-data/dummy-about.json")
        >>> about.home
        ['http://github.com/conda/rattler']
        >>> about.home = ['https://test.home']
        >>> about.home
        ['https://test.home/']
        >>>
        ```
        """
        return self._inner.home

    @home.setter
    def home(self, value: List[str]) -> None:
        self._inner.home = value

    @property
    def extra(self) -> dict[str, Any]:
        """
        The JSON-serializable `extra` metadata attached to `about.json`.

        Examples
        --------
        ```python
        >>> about = AboutJson.from_str('{"extra": {"flow_id": "2024.08.13"}}')
        >>> about.extra
        {'flow_id': '2024.08.13'}
        >>> about.extra = {"nested": {"value": [1, 2, None]}}
        >>> about.extra
        {'nested': {'value': [1, 2, None]}}
        >>>
        ```
        """
        return self._inner.extra

    @extra.setter
    def extra(self, value: dict[str, Any]) -> None:
        self._inner.extra = value

    @property
    def license(self) -> Optional[str]:
        """
        The license of the package.

        Examples
        --------
        ```python
        >>> about = AboutJson.from_path("../test-data/dummy-about.json")
        >>> about.license
        'BSD-3-Clause'
        >>> about.license = 'MIT'
        >>> about.license
        'MIT'
        >>>
        ```
        """
        if license := self._inner.license:
            return license

        return None

    @license.setter
    def license(self, value: Optional[str]) -> None:
        self._inner.license = value

    @property
    def license_family(self) -> Optional[str]:
        """
        The license family of the package.

        Examples
        --------
        ```python
        >>> about = AboutJson.from_path("../test-data/dummy-about.json")
        >>> about.license_family
        >>> type(about.license_family)
        <class 'NoneType'>
        >>> about.license_family = 'BSD'
        >>> about.license_family
        'BSD'
        >>>
        ```
        """
        if license_family := self._inner.license_family:
            return license_family

        return None

    @license_family.setter
    def license_family(self, value: Optional[str]) -> None:
        self._inner.license_family = value

    @property
    def source_url(self) -> Optional[str]:
        """
        The URL to the latest source code of the package.

        Examples
        --------
        ```python
        >>> about = AboutJson.from_path("../test-data/dummy-about.json")
        >>> about.source_url
        'https://github.com/conda/rattler'
        >>> about.source_url = 'https://test.source'
        >>> about.source_url
        'https://test.source/'
        >>>
        ```
        """
        if source_url := self._inner.source_url:
            return source_url

        return None

    @source_url.setter
    def source_url(self, value: Optional[str]) -> None:
        self._inner.source_url = value

    @property
    def summary(self) -> Optional[str]:
        """
        A Short summary description.

        Examples
        --------
        ```python
        >>> about = AboutJson.from_path("../test-data/dummy-about.json")
        >>> about.summary
        'A dummy summary.'
        >>> about.summary = 'New summary'
        >>> about.summary
        'New summary'
        >>>
        ```
        """
        if summary := self._inner.summary:
            return summary

        return None

    @summary.setter
    def summary(self, value: Optional[str]) -> None:
        self._inner.summary = value

    @classmethod
    def _from_py_about_json(cls, py_about_json: PyAboutJson) -> AboutJson:
        about_json = cls.__new__(cls)
        about_json._inner = py_about_json

        return about_json

    def __repr__(self) -> str:
        """
        Returns a representation of the AboutJson.
        """
        return "AboutJson()"
