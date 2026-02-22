"""Application data stored by virtualenv."""

from __future__ import annotations

import logging
import os
import shutil

from platformdirs import user_cache_dir, user_data_dir

from .na import AppDataDisabled
from .read_only import ReadOnlyAppData
from .via_disk_folder import AppDataDiskFolder
from .via_tempdir import TempAppData

LOGGER = logging.getLogger(__name__)


def _default_app_data_dir(env):
    key = "VIRTUALENV_OVERRIDE_APP_DATA"
    if key in env:
        return env[key]
    return _cache_dir_with_migration()


def _cache_dir_with_migration():
    new_dir = user_cache_dir(appname="virtualenv", appauthor="pypa")
    old_dir = user_data_dir(appname="virtualenv", appauthor="pypa")
    if new_dir == old_dir:
        return new_dir
    if os.path.isdir(old_dir) and not os.path.isdir(new_dir):
        LOGGER.info("migrating app data from %s to %s", old_dir, new_dir)
        try:
            shutil.move(old_dir, new_dir)
        except OSError as exception:
            LOGGER.warning(
                "could not migrate app data from %s to %s: %r, using old location", old_dir, new_dir, exception
            )
            return old_dir
    return new_dir


def make_app_data(folder, **kwargs):
    is_read_only = kwargs.pop("read_only")
    env = kwargs.pop("env")
    if kwargs:  # py3+ kwonly
        msg = "unexpected keywords: {}"
        raise TypeError(msg)

    if folder is None:
        folder = _default_app_data_dir(env)
    folder = os.path.abspath(folder)

    if is_read_only:
        return ReadOnlyAppData(folder)

    try:
        os.makedirs(folder, exist_ok=True)
        LOGGER.debug("created app data folder %s", folder)
    except OSError as exception:
        LOGGER.info("could not create app data folder %s due to %r", folder, exception)

    if os.access(folder, os.W_OK):
        return AppDataDiskFolder(folder)
    LOGGER.debug("app data folder %s has no write access", folder)
    return TempAppData()


__all__ = (
    "AppDataDisabled",
    "AppDataDiskFolder",
    "ReadOnlyAppData",
    "TempAppData",
    "make_app_data",
)
