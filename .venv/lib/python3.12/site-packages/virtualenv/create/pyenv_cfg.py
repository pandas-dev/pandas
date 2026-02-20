from __future__ import annotations

import logging
import os
from collections import OrderedDict

LOGGER = logging.getLogger(__name__)


class PyEnvCfg:
    def __init__(self, content, path) -> None:
        self.content = content
        self.path = path

    @classmethod
    def from_folder(cls, folder):
        return cls.from_file(folder / "pyvenv.cfg")

    @classmethod
    def from_file(cls, path):
        content = cls._read_values(path) if path.exists() else OrderedDict()
        return PyEnvCfg(content, path)

    @staticmethod
    def _read_values(path):
        content = OrderedDict()
        for line in path.read_text(encoding="utf-8").splitlines():
            equals_at = line.index("=")
            key = line[:equals_at].strip()
            value = line[equals_at + 1 :].strip()
            if len(value) > 1 and value[0] in {"'", '"'} and value[0] == value[-1]:
                value = value[1:-1]
            content[key] = value
        return content

    def write(self):
        LOGGER.debug("write %s", self.path)
        text = ""
        for key, value in self.content.items():
            # Use abspath to normalize relative paths but preserve symlinks (match venv behavior)
            # See issue #2770 - realpath resolves symlinks which breaks prefix symlinks
            if key == "prompt" and value:
                normalized_value = f'"{value}"'
            else:
                normalized_value = os.path.abspath(value) if value and os.path.exists(value) else value
            line = f"{key} = {normalized_value}"
            LOGGER.debug("\t%s", line)
            text += line
            text += "\n"
        self.path.write_text(text, encoding="utf-8")

    def refresh(self):
        self.content = self._read_values(self.path)
        return self.content

    def __setitem__(self, key, value) -> None:
        self.content[key] = value

    def __getitem__(self, key):
        return self.content[key]

    def __contains__(self, item) -> bool:
        return item in self.content

    def update(self, other):
        self.content.update(other)
        return self

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}(path={self.path})"


__all__ = [
    "PyEnvCfg",
]
