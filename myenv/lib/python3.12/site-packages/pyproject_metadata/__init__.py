# SPDX-License-Identifier: MIT

from __future__ import annotations

import collections
import copy
import dataclasses
import email.utils
import os
import os.path
import pathlib
import sys
import typing


if typing.TYPE_CHECKING:
    from collections.abc import Mapping
    from typing import Any

    from packaging.requirements import Requirement

    if sys.version_info < (3, 11):
        from typing_extensions import Self
    else:
        from typing import Self

import packaging.markers
import packaging.requirements
import packaging.specifiers
import packaging.utils
import packaging.version


__version__ = '0.8.0'

KNOWN_METADATA_VERSIONS = {'2.1', '2.2', '2.3'}


class ConfigurationError(Exception):
    '''Error in the backend metadata.'''
    def __init__(self, msg: str, *, key: str | None = None):
        super().__init__(msg)
        self._key = key

    @property
    def key(self) -> str | None:  # pragma: no cover
        return self._key


class RFC822Message:
    '''Python-flavored RFC 822 message implementation.'''

    def __init__(self) -> None:
        self.headers: collections.OrderedDict[str, list[str]] = collections.OrderedDict()
        self.body: str | None = None

    def __setitem__(self, name: str, value: str | None) -> None:
        if not value:
            return
        if name not in self.headers:
            self.headers[name] = []
        self.headers[name].append(value)

    def __str__(self) -> str:
        text = ''
        for name, entries in self.headers.items():
            for entry in entries:
                lines = entry.strip('\n').split('\n')
                text += f'{name}: {lines[0]}\n'
                for line in lines[1:]:
                    text += ' ' * 8 + line + '\n'
        if self.body:
            text += '\n' + self.body
        return text

    def __bytes__(self) -> bytes:
        return str(self).encode()


class DataFetcher:
    def __init__(self, data: Mapping[str, Any]) -> None:
        self._data = data

    def __contains__(self, key: Any) -> bool:
        if not isinstance(key, str):
            return False
        val = self._data
        try:
            for part in key.split('.'):
                val = val[part]
        except KeyError:
            return False
        return True

    def get(self, key: str) -> Any:
        val = self._data
        for part in key.split('.'):
            val = val[part]
        return val

    def get_str(self, key: str) -> str | None:
        try:
            val = self.get(key)
            if not isinstance(val, str):
                msg = f'Field "{key}" has an invalid type, expecting a string (got "{val}")'
                raise ConfigurationError(msg, key=key)
            return val
        except KeyError:
            return None

    def get_list(self, key: str) -> list[str]:
        try:
            val = self.get(key)
            if not isinstance(val, list):
                msg = f'Field "{key}" has an invalid type, expecting a list of strings (got "{val}")'
                raise ConfigurationError(msg, key=val)
            for item in val:
                if not isinstance(item, str):
                    msg = f'Field "{key}" contains item with invalid type, expecting a string (got "{item}")'
                    raise ConfigurationError(msg, key=key)
            return val
        except KeyError:
            return []

    def get_dict(self, key: str) -> dict[str, str]:
        try:
            val = self.get(key)
            if not isinstance(val, dict):
                msg = f'Field "{key}" has an invalid type, expecting a dictionary of strings (got "{val}")'
                raise ConfigurationError(msg, key=key)
            for subkey, item in val.items():
                if not isinstance(item, str):
                    msg = f'Field "{key}.{subkey}" has an invalid type, expecting a string (got "{item}")'
                    raise ConfigurationError(msg, key=f'{key}.{subkey}')
            return val
        except KeyError:
            return {}

    def get_people(self, key: str) -> list[tuple[str, str]]:
        try:
            val = self.get(key)
            if not (
                isinstance(val, list)
                and all(isinstance(x, dict) for x in val)
                and all(
                    isinstance(item, str)
                    for items in [_dict.values() for _dict in val]
                    for item in items
                )
            ):
                msg = (
                    f'Field "{key}" has an invalid type, expecting a list of '
                    f'dictionaries containing the "name" and/or "email" keys (got "{val}")'
                )
                raise ConfigurationError(msg, key=key)
            return [
                (entry.get('name', 'Unknown'), entry.get('email'))
                for entry in val
            ]
        except KeyError:
            return []


class License(typing.NamedTuple):
    text: str
    file: pathlib.Path | None


class Readme(typing.NamedTuple):
    text: str
    file: pathlib.Path | None
    content_type: str


@dataclasses.dataclass
class StandardMetadata:
    name: str
    version: packaging.version.Version | None = None
    description: str | None = None
    license: License | None = None
    readme: Readme | None = None
    requires_python: packaging.specifiers.SpecifierSet | None = None
    dependencies: list[Requirement] = dataclasses.field(default_factory=list)
    optional_dependencies: dict[str, list[Requirement]] = dataclasses.field(default_factory=dict)
    entrypoints: dict[str, dict[str, str]] = dataclasses.field(default_factory=dict)
    authors: list[tuple[str, str]] = dataclasses.field(default_factory=list)
    maintainers: list[tuple[str, str]] = dataclasses.field(default_factory=list)
    urls: dict[str, str] = dataclasses.field(default_factory=dict)
    classifiers: list[str] = dataclasses.field(default_factory=list)
    keywords: list[str] = dataclasses.field(default_factory=list)
    scripts: dict[str, str] = dataclasses.field(default_factory=dict)
    gui_scripts: dict[str, str] = dataclasses.field(default_factory=dict)
    dynamic: list[str] = dataclasses.field(default_factory=list)

    _metadata_version: str | None = None

    @property
    def metadata_version(self) -> str:
        if self._metadata_version is None:
            return '2.2' if self.dynamic else '2.1'
        return self._metadata_version

    @property
    def canonical_name(self) -> str:
        return packaging.utils.canonicalize_name(self.name)

    @classmethod
    def from_pyproject(
        cls,
        data: Mapping[str, Any],
        project_dir: str | os.PathLike[str] = os.path.curdir,
        metadata_version: str | None = None,
    ) -> Self:
        fetcher = DataFetcher(data)
        project_dir = pathlib.Path(project_dir)

        if 'project' not in fetcher:
            msg = 'Section "project" missing in pyproject.toml'
            raise ConfigurationError(msg)

        dynamic = fetcher.get_list('project.dynamic')
        if 'name' in dynamic:
            msg = 'Unsupported field "name" in "project.dynamic"'
            raise ConfigurationError(msg)

        for field in dynamic:
            if field in data['project']:
                msg = f'Field "project.{field}" declared as dynamic in "project.dynamic" but is defined'
                raise ConfigurationError(msg)

        name = fetcher.get_str('project.name')
        if not name:
            msg = 'Field "project.name" missing'
            raise ConfigurationError(msg)

        version_string = fetcher.get_str('project.version')
        requires_python_string = fetcher.get_str('project.requires-python')
        version = packaging.version.Version(version_string) if version_string else None

        if version is None and 'version' not in dynamic:
            msg = 'Field "project.version" missing and "version" not specified in "project.dynamic"'
            raise ConfigurationError(msg)

        # Description fills Summary, which cannot be multiline
        # However, throwing an error isn't backward compatible,
        # so leave it up to the users for now.
        description = fetcher.get_str('project.description')

        if metadata_version and metadata_version not in KNOWN_METADATA_VERSIONS:
            msg = f'The metadata_version must be one of {KNOWN_METADATA_VERSIONS} or None (default)'
            raise ConfigurationError(msg)

        return cls(
            name,
            version,
            description,
            cls._get_license(fetcher, project_dir),
            cls._get_readme(fetcher, project_dir),
            packaging.specifiers.SpecifierSet(requires_python_string) if requires_python_string else None,
            cls._get_dependencies(fetcher),
            cls._get_optional_dependencies(fetcher),
            cls._get_entrypoints(fetcher),
            fetcher.get_people('project.authors'),
            fetcher.get_people('project.maintainers'),
            fetcher.get_dict('project.urls'),
            fetcher.get_list('project.classifiers'),
            fetcher.get_list('project.keywords'),
            fetcher.get_dict('project.scripts'),
            fetcher.get_dict('project.gui-scripts'),
            dynamic,
            metadata_version,
        )

    def _update_dynamic(self, value: Any) -> None:
        if value and 'version' in self.dynamic:
            self.dynamic.remove('version')

    def __setattr__(self, name: str, value: Any) -> None:
        # update dynamic when version is set
        if name == 'version' and hasattr(self, 'dynamic'):
            self._update_dynamic(value)
        super().__setattr__(name, value)

    def as_rfc822(self) -> RFC822Message:
        message = RFC822Message()
        self.write_to_rfc822(message)
        return message

    def write_to_rfc822(self, message: RFC822Message) -> None:  # noqa: C901
        message['Metadata-Version'] = self.metadata_version
        message['Name'] = self.name
        if not self.version:
            msg = 'Missing version field'
            raise ConfigurationError(msg)
        message['Version'] = str(self.version)
        # skip 'Platform'
        # skip 'Supported-Platform'
        if self.description:
            message['Summary'] = self.description
        message['Keywords'] = ','.join(self.keywords)
        if 'homepage' in self.urls:
            message['Home-page'] = self.urls['homepage']
        # skip 'Download-URL'
        message['Author'] = self._name_list(self.authors)
        message['Author-Email'] = self._email_list(self.authors)
        message['Maintainer'] = self._name_list(self.maintainers)
        message['Maintainer-Email'] = self._email_list(self.maintainers)
        if self.license:
            message['License'] = self.license.text
        for classifier in self.classifiers:
            message['Classifier'] = classifier
        # skip 'Provides-Dist'
        # skip 'Obsoletes-Dist'
        # skip 'Requires-External'
        for name, url in self.urls.items():
            message['Project-URL'] = f'{name.capitalize()}, {url}'
        if self.requires_python:
            message['Requires-Python'] = str(self.requires_python)
        for dep in self.dependencies:
            message['Requires-Dist'] = str(dep)
        for extra, requirements in self.optional_dependencies.items():
            norm_extra = extra.replace('.', '-').replace('_', '-').lower()
            message['Provides-Extra'] = norm_extra
            for requirement in requirements:
                message['Requires-Dist'] = str(self._build_extra_req(norm_extra, requirement))
        if self.readme:
            if self.readme.content_type:
                message['Description-Content-Type'] = self.readme.content_type
            message.body = self.readme.text
        # Core Metadata 2.2
        if self.metadata_version != '2.1':
            for field in self.dynamic:
                if field in ('name', 'version'):
                    msg = f'Field cannot be dynamic: {field}'
                    raise ConfigurationError(msg)
                message['Dynamic'] = field

    def _name_list(self, people: list[tuple[str, str]]) -> str:
        return ', '.join(
            name
            for name, email_ in people
            if not email_
        )

    def _email_list(self, people: list[tuple[str, str]]) -> str:
        return ', '.join(
            email.utils.formataddr((name, _email))
            for name, _email in people
            if _email
        )

    def _build_extra_req(
        self,
        extra: str,
        requirement: Requirement,
    ) -> Requirement:
        # append or add our extra marker
        requirement = copy.copy(requirement)
        if requirement.marker:
            if 'or' in requirement.marker._markers:
                requirement.marker = packaging.markers.Marker(
                    f'({requirement.marker}) and extra == "{extra}"'
                )
            else:
                requirement.marker = packaging.markers.Marker(
                    f'{requirement.marker} and extra == "{extra}"'
                )
        else:
            requirement.marker = packaging.markers.Marker(f'extra == "{extra}"')
        return requirement

    @staticmethod
    def _get_license(fetcher: DataFetcher, project_dir: pathlib.Path) -> License | None:
        if 'project.license' not in fetcher:
            return None

        _license = fetcher.get_dict('project.license')
        for field in _license:
            if field not in ('file', 'text'):
                msg = f'Unexpected field "project.license.{field}"'
                raise ConfigurationError(msg, key=f'project.license.{field}')

        file: pathlib.Path | None = None
        filename = fetcher.get_str('project.license.file')
        text = fetcher.get_str('project.license.text')

        if (filename and text) or (not filename and not text):
            msg = f'Invalid "project.license" value, expecting either "file" or "text" (got "{_license}")'
            raise ConfigurationError(msg, key='project.license')

        if filename:
            file = project_dir.joinpath(filename)
            if not file.is_file():
                msg = f'License file not found ("{filename}")'
                raise ConfigurationError(msg, key='project.license.file')
            text = file.read_text(encoding='utf-8')

        assert text is not None
        return License(text, file)

    @staticmethod
    def _get_readme(fetcher: DataFetcher, project_dir: pathlib.Path) -> Readme | None:  # noqa: C901
        if 'project.readme' not in fetcher:
            return None

        filename: str | None
        file: pathlib.Path | None = None
        text: str | None
        content_type: str | None

        readme = fetcher.get('project.readme')
        if isinstance(readme, str):
            # readme is a file
            text = None
            filename = readme
            if filename.endswith('.md'):
                content_type = 'text/markdown'
            elif filename.endswith('.rst'):
                content_type = 'text/x-rst'
            else:
                msg = f'Could not infer content type for readme file "{filename}"'
                raise ConfigurationError(msg, key='project.readme')
        elif isinstance(readme, dict):
            # readme is a dict containing either 'file' or 'text', and content-type
            for field in readme:
                if field not in ('content-type', 'file', 'text'):
                    msg = f'Unexpected field "project.readme.{field}"'
                    raise ConfigurationError(msg, key=f'project.readme.{field}')
            content_type = fetcher.get_str('project.readme.content-type')
            filename = fetcher.get_str('project.readme.file')
            text = fetcher.get_str('project.readme.text')
            if (filename and text) or (not filename and not text):
                msg = f'Invalid "project.readme" value, expecting either "file" or "text" (got "{readme}")'
                raise ConfigurationError(msg, key='project.readme')
            if not content_type:
                msg = 'Field "project.readme.content-type" missing'
                raise ConfigurationError(msg, key='project.readme.content-type')
        else:
            msg = (
                f'Field "project.readme" has an invalid type, expecting either, '
                f'a string or dictionary of strings (got "{readme}")'
            )
            raise ConfigurationError(msg, key='project.readme')

        if filename:
            file = project_dir.joinpath(filename)
            if not file.is_file():
                msg = f'Readme file not found ("{filename}")'
                raise ConfigurationError(msg, key='project.readme.file')
            text = file.read_text(encoding='utf-8')

        assert text is not None
        return Readme(text, file, content_type)

    @staticmethod
    def _get_dependencies(fetcher: DataFetcher) -> list[Requirement]:
        try:
            requirement_strings = fetcher.get_list('project.dependencies')
        except KeyError:
            return []

        requirements: list[Requirement] = []
        for req in requirement_strings:
            try:
                requirements.append(packaging.requirements.Requirement(req))
            except packaging.requirements.InvalidRequirement as e:
                msg = (
                    'Field "project.dependencies" contains an invalid PEP 508 '
                    f'requirement string "{req}" ("{e}")'
                )
                raise ConfigurationError(msg) from None
        return requirements

    @staticmethod
    def _get_optional_dependencies(fetcher: DataFetcher) -> dict[str, list[Requirement]]:
        try:
            val = fetcher.get('project.optional-dependencies')
        except KeyError:
            return {}

        requirements_dict: dict[str, list[Requirement]] = {}
        if not isinstance(val, dict):
            msg = (
                'Field "project.optional-dependencies" has an invalid type, expecting a '
                f'dictionary of PEP 508 requirement strings (got "{val}")'
            )
            raise ConfigurationError(msg)
        for extra, requirements in val.copy().items():
            assert isinstance(extra, str)
            if not isinstance(requirements, list):
                msg = (
                    f'Field "project.optional-dependencies.{extra}" has an invalid type, expecting a '
                    f'dictionary PEP 508 requirement strings (got "{requirements}")'
                )
                raise ConfigurationError(msg)
            requirements_dict[extra] = []
            for req in requirements:
                if not isinstance(req, str):
                    msg = (
                        f'Field "project.optional-dependencies.{extra}" has an invalid type, '
                        f'expecting a PEP 508 requirement string (got "{req}")'
                    )
                    raise ConfigurationError(msg)
                try:
                    requirements_dict[extra].append(packaging.requirements.Requirement(req))
                except packaging.requirements.InvalidRequirement as e:
                    msg = (
                        f'Field "project.optional-dependencies.{extra}" contains '
                        f'an invalid PEP 508 requirement string "{req}" ("{e}")'
                    )
                    raise ConfigurationError(msg) from None
        return dict(requirements_dict)

    @staticmethod
    def _get_entrypoints(fetcher: DataFetcher) -> dict[str, dict[str, str]]:
        try:
            val = fetcher.get('project.entry-points')
        except KeyError:
            return {}
        if not isinstance(val, dict):
            msg = (
                'Field "project.entry-points" has an invalid type, expecting a '
                f'dictionary of entrypoint sections (got "{val}")'
            )
            raise ConfigurationError(msg)
        for section, entrypoints in val.items():
            assert isinstance(section, str)
            if not isinstance(entrypoints, dict):
                msg = (
                    f'Field "project.entry-points.{section}" has an invalid type, expecting a '
                    f'dictionary of entrypoints (got "{entrypoints}")'
                )
                raise ConfigurationError(msg)
            for name, entrypoint in entrypoints.items():
                assert isinstance(name, str)
                if not isinstance(entrypoint, str):
                    msg = (
                        f'Field "project.entry-points.{section}.{name}" has an invalid type, '
                        f'expecting a string (got "{entrypoint}")'
                    )
                    raise ConfigurationError(msg)
        return val
