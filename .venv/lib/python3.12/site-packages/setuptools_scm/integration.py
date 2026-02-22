from __future__ import annotations

import logging
import textwrap

from pathlib import Path

from . import _types as _t

log = logging.getLogger(__name__)


def data_from_mime(path: _t.PathT, content: str | None = None) -> dict[str, str]:
    """return a mapping from mime/pseudo-mime content
    :param path: path to the mime file
    :param content: content of the mime file, if None, read from path
    :rtype: dict[str, str]

    """

    if content is None:
        content = Path(path).read_text(encoding="utf-8")
    log.debug("mime %s content:\n%s", path, textwrap.indent(content, "    "))

    from email.parser import HeaderParser

    parser = HeaderParser()
    message = parser.parsestr(content)
    data = dict(message.items())
    log.debug("mime %s data:\n%s", path, data)
    return data
