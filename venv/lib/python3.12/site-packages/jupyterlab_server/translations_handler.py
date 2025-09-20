"""
Translation handler.
"""
# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.
from __future__ import annotations

import json
import traceback
from functools import partial

import tornado

from .settings_utils import SchemaHandler
from .translation_utils import (
    SYS_LOCALE,
    get_language_pack,
    get_language_packs,
    is_valid_locale,
    translator,
)


class TranslationsHandler(SchemaHandler):
    """An API handler for translations."""

    @tornado.web.authenticated
    async def get(self, locale: str | None = None) -> None:
        """
        Get installed language packs.

        If `locale` is equals to "default", the default locale will be used.

        Parameters
        ----------
        locale: str, optional
            If no locale is provided, it will list all the installed language packs.
            Default is `None`.
        """
        data: dict
        data, message = {}, ""
        try:
            current_loop = tornado.ioloop.IOLoop.current()
            if locale is None:
                data, message = await current_loop.run_in_executor(
                    None,
                    partial(get_language_packs, display_locale=self.get_current_locale()),
                )
            else:
                locale = locale or SYS_LOCALE
                if locale == "default":
                    locale = SYS_LOCALE
                data, message = await current_loop.run_in_executor(
                    None, partial(get_language_pack, locale)
                )
                if data == {} and not message:
                    if is_valid_locale(locale):
                        message = f"Language pack '{locale}' not installed!"
                    else:
                        message = f"Language pack '{locale}' not valid!"
                elif is_valid_locale(locale):
                    # only change locale if the language pack is installed and valid
                    translator.set_locale(locale)
        except Exception:
            message = traceback.format_exc()

        self.set_status(200)
        self.finish(json.dumps({"data": data, "message": message}))
