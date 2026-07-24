"""Manager to read and modify config data in JSON files."""

# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.
from __future__ import annotations

import copy
import errno
import glob
import json
import os
import typing as t

from traitlets.config import LoggingConfigurable
from traitlets.traitlets import Bool, Unicode

StrDict = dict[str, t.Any]


def recursive_update(target: StrDict, new: StrDict) -> None:
    """Recursively update one dictionary using another.

    None values will delete their keys.
    """
    for k, v in new.items():
        if isinstance(v, dict):
            if k not in target:
                target[k] = {}
            recursive_update(target[k], v)
            if not target[k]:
                # Prune empty subdicts
                del target[k]

        elif v is None:
            target.pop(k, None)

        else:
            target[k] = v


def remove_defaults(data: StrDict, defaults: StrDict) -> None:
    """Recursively remove items from dict that are already in defaults"""
    # copy the iterator, since data will be modified
    for key, value in list(data.items()):
        if key in defaults:
            if isinstance(value, dict):
                remove_defaults(data[key], defaults[key])
                if not data[key]:  # prune empty subdicts
                    del data[key]
            elif value == defaults[key]:
                del data[key]


class BaseJSONConfigManager(LoggingConfigurable):
    """General JSON config manager

    Deals with persisting/storing config in a json file with optionally
    default values in a {section_name}.d directory.
    """

    config_dir = Unicode(".")
    read_directory = Bool(True)

    def ensure_config_dir_exists(self) -> None:
        """Will try to create the config_dir directory."""
        try:
            os.makedirs(self.config_dir, 0o755)
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise

    def file_name(self, section_name: str) -> str:
        """Returns the json filename for the section_name: {config_dir}/{section_name}.json"""
        return os.path.join(self.config_dir, section_name + ".json")

    def directory(self, section_name: str) -> str:
        """Returns the directory name for the section name: {config_dir}/{section_name}.d"""
        return os.path.join(self.config_dir, section_name + ".d")

    def get(self, section_name: str, include_root: bool = True) -> dict[str, t.Any]:
        """Retrieve the config data for the specified section.

        Returns the data as a dictionary, or an empty dictionary if the file
        doesn't exist.

        When include_root is False, it will not read the root .json file,
        effectively returning the default values.
        """
        paths = [self.file_name(section_name)] if include_root else []
        if self.read_directory:
            pattern = os.path.join(self.directory(section_name), "*.json")
            # These json files should be processed first so that the
            # {section_name}.json take precedence.
            # The idea behind this is that installing a Python package may
            # put a json file somewhere in the a .d directory, while the
            # .json file is probably a user configuration.
            paths = sorted(glob.glob(pattern)) + paths
        self.log.debug(
            "Paths used for configuration of %s: \n\t%s",
            section_name,
            "\n\t".join(paths),
        )
        data: dict[str, t.Any] = {}
        for path in paths:
            if os.path.isfile(path) and os.path.getsize(path):
                with open(path, encoding="utf-8") as f:
                    try:
                        recursive_update(data, json.load(f))
                    except json.decoder.JSONDecodeError:
                        self.log.warning("Invalid JSON in %s, skipping", path)
        return data

    def set(self, section_name: str, data: t.Any) -> None:
        """Store the given config data."""
        filename = self.file_name(section_name)
        self.ensure_config_dir_exists()

        if self.read_directory:
            # we will modify data in place, so make a copy
            data = copy.deepcopy(data)
            defaults = self.get(section_name, include_root=False)
            remove_defaults(data, defaults)

        # Generate the JSON up front, since it could raise an exception,
        # in order to avoid writing half-finished corrupted data to disk.
        json_content = json.dumps(data, indent=2)
        with open(filename, "w", encoding="utf-8") as f:
            f.write(json_content)

    def update(self, section_name: str, new_data: t.Any) -> dict[str, t.Any]:
        """Modify the config section by recursively updating it with new_data.

        Returns the modified config data as a dictionary.
        """
        data = self.get(section_name)
        recursive_update(data, new_data)
        self.set(section_name, data)
        return data
