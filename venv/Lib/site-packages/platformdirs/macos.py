"""macOS."""

from __future__ import annotations

import os.path
import sys
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from collections.abc import Iterator

from ._xdg import XDGMixin
from .api import PlatformDirsABC

if TYPE_CHECKING:
    from pathlib import Path


class _MacOSDefaults(PlatformDirsABC):  # ruff:ignore[too-many-public-methods]
    """Default platform directories for macOS without XDG environment variable overrides.

    Follows the guidance from `Apple's File System Programming Guide
    <https://developer.apple.com/library/archive/documentation/FileManagement/Conceptual/FileSystemProgrammingGuide/MacOSXDirectories/MacOSXDirectories.html>`_.
    The XDG env var handling is in :class:`~platformdirs._xdg.XDGMixin`.

    """

    def _base_user_app_support_dir(self) -> str:
        return self._append_app_name_and_version(os.path.expanduser("~/Library/Application Support"))  # ruff:ignore[os-path-expanduser]

    def _base_site_dirs(self) -> list[str]:
        is_homebrew = "/opt/python" in sys.prefix
        homebrew_prefix = sys.prefix.split("/opt/python")[0] if is_homebrew else ""
        path_list = [self._append_app_name_and_version(f"{homebrew_prefix}/share")] if is_homebrew else []
        path_list.append(self._append_app_name_and_version("/Library/Application Support"))
        return path_list

    @property
    def user_data_dir(self) -> str:
        """Data directory tied to the user, e.g. ``~/Library/Application Support/$appname/$version``."""
        return self._base_user_app_support_dir()

    @property
    def _site_data_dirs(self) -> list[str]:
        return self._base_site_dirs()

    @property
    def site_data_path(self) -> Path:
        """Data path shared by users. Only return the first item, even if ``multipath`` is set to ``True``."""
        return self._first_item_as_path_if_multipath(self.site_data_dir)

    @property
    def site_config_path(self) -> Path:
        """Config path shared by users. Only return the first item, even if ``multipath`` is set to ``True``."""
        return self._first_item_as_path_if_multipath(self.site_config_dir)

    @property
    def user_config_dir(self) -> str:
        """Config directory tied to the user, same as `user_data_dir`."""
        return self._base_user_app_support_dir()

    @property
    def _site_config_dirs(self) -> list[str]:
        return self._base_site_dirs()

    @property
    def user_cache_dir(self) -> str:
        """Cache directory tied to the user, e.g. ``~/Library/Caches/$appname/$version``."""
        return self._append_app_name_and_version(os.path.expanduser("~/Library/Caches"))  # ruff:ignore[os-path-expanduser]

    @property
    def site_cache_dir(self) -> str:
        """Cache directory shared by users, e.g. ``/Library/Caches/$appname/$version``. If we're using a Python binary managed by `Homebrew <https://brew.sh>`_, the directory will be under the Homebrew prefix, e.g. ``$homebrew_prefix/var/cache/$appname/$version``. If `multipath <platformdirs.api.PlatformDirsABC.multipath>` is enabled, and we're in Homebrew, the response is a multi-path string separated by ":", e.g. ``$homebrew_prefix/var/cache/$appname/$version:/Library/Caches/$appname/$version``."""
        is_homebrew = "/opt/python" in sys.prefix
        homebrew_prefix = sys.prefix.split("/opt/python")[0] if is_homebrew else ""
        path_list = [self._append_app_name_and_version(f"{homebrew_prefix}/var/cache")] if is_homebrew else []
        path_list.append(self._append_app_name_and_version("/Library/Caches"))
        if self.multipath:
            return os.pathsep.join(path_list)
        return path_list[0]

    @property
    def site_cache_path(self) -> Path:
        """Cache path shared by users. Only return the first item, even if ``multipath`` is set to ``True``."""
        return self._first_item_as_path_if_multipath(self.site_cache_dir)

    @property
    def user_state_dir(self) -> str:
        """State directory tied to the user, same as `user_data_dir`."""
        return self._base_user_app_support_dir()

    @property
    def site_state_dir(self) -> str:
        """State directory shared by users, same as `site_data_dir`."""
        return self._base_site_dirs()[0]

    @property
    def user_log_dir(self) -> str:
        """Log directory tied to the user, e.g. ``~/Library/Logs/$appname/$version``."""
        return self._append_app_name_and_version(os.path.expanduser("~/Library/Logs"))  # ruff:ignore[os-path-expanduser]

    @property
    def site_log_dir(self) -> str:
        """Log directory shared by users, e.g. ``/Library/Logs/$appname/$version``."""
        return self._append_app_name_and_version("/Library/Logs")

    @property
    def user_documents_dir(self) -> str:
        """Documents directory tied to the user, e.g. ``~/Documents``."""
        return os.path.expanduser("~/Documents")  # ruff:ignore[os-path-expanduser]

    @property
    def user_downloads_dir(self) -> str:
        """Downloads directory tied to the user, e.g. ``~/Downloads``."""
        return os.path.expanduser("~/Downloads")  # ruff:ignore[os-path-expanduser]

    @property
    def user_pictures_dir(self) -> str:
        """Pictures directory tied to the user, e.g. ``~/Pictures``."""
        return os.path.expanduser("~/Pictures")  # ruff:ignore[os-path-expanduser]

    @property
    def user_videos_dir(self) -> str:
        """Videos directory tied to the user, e.g. ``~/Movies``."""
        return os.path.expanduser("~/Movies")  # ruff:ignore[os-path-expanduser]

    @property
    def user_music_dir(self) -> str:
        """Music directory tied to the user, e.g. ``~/Music``."""
        return os.path.expanduser("~/Music")  # ruff:ignore[os-path-expanduser]

    @property
    def user_desktop_dir(self) -> str:
        """Desktop directory tied to the user, e.g. ``~/Desktop``."""
        return os.path.expanduser("~/Desktop")  # ruff:ignore[os-path-expanduser]

    @property
    def user_projects_dir(self) -> str:
        """Projects directory tied to the user, e.g. ``~/Projects``."""
        return os.path.expanduser("~/Projects")  # ruff:ignore[os-path-expanduser]

    @property
    def user_publicshare_dir(self) -> str:
        """Public share directory tied to the user, e.g. ``~/Public``."""
        return os.path.expanduser("~/Public")  # ruff:ignore[os-path-expanduser]  # API returns str, not Path

    @property
    def user_templates_dir(self) -> str:
        """Templates directory tied to the user, e.g. ``~/Templates``."""
        return os.path.expanduser("~/Templates")  # ruff:ignore[os-path-expanduser]  # API returns str, not Path

    @property
    def user_fonts_dir(self) -> str:
        """Fonts directory tied to the user, e.g. ``~/Library/Fonts``."""
        return os.path.expanduser("~/Library/Fonts")  # ruff:ignore[os-path-expanduser]  # API returns str, not Path

    @property
    def user_preference_dir(self) -> str:
        """Preference directory tied to the user, e.g. ``~/Library/Preferences/AppName``."""
        return self._append_app_name_and_version(os.path.expanduser("~/Library/Preferences"))  # ruff:ignore[os-path-expanduser]  # API returns str, not Path

    @property
    def user_bin_dir(self) -> str:
        """Bin directory tied to the user, e.g. ``~/.local/bin``."""
        return os.path.expanduser("~/.local/bin")  # ruff:ignore[os-path-expanduser]

    @property
    def site_bin_dir(self) -> str:
        """Bin directory shared by users, e.g. ``/usr/local/bin``."""
        return "/usr/local/bin"

    @property
    def user_applications_dir(self) -> str:
        """Applications directory tied to the user, e.g. ``~/Applications``."""
        return os.path.expanduser("~/Applications")  # ruff:ignore[os-path-expanduser]

    @property
    def _site_applications_dirs(self) -> list[str]:
        return ["/Applications"]

    @property
    def site_applications_dir(self) -> str:
        """Applications directory shared by users, e.g. ``/Applications``."""
        dirs = self._site_applications_dirs
        return os.pathsep.join(dirs) if self.multipath else dirs[0]

    @property
    def user_runtime_dir(self) -> str:
        """Runtime directory tied to the user, e.g. ``~/Library/Caches/TemporaryItems/$appname/$version``."""
        return self._append_app_name_and_version(os.path.expanduser("~/Library/Caches/TemporaryItems"))  # ruff:ignore[os-path-expanduser]

    @property
    def site_runtime_dir(self) -> str:
        """Runtime directory shared by users, same as `user_runtime_dir`."""
        return self.user_runtime_dir

    def iter_config_dirs(self) -> Iterator[str]:
        """:yield: all user and site configuration directories."""
        yield self.user_config_dir
        yield from self._site_config_dirs

    def iter_data_dirs(self) -> Iterator[str]:
        """:yield: all user and site data directories."""
        yield self.user_data_dir
        yield from self._site_data_dirs


class MacOS(XDGMixin, _MacOSDefaults):
    """Platform directories for the macOS operating system.

    Follows the guidance from `Apple documentation
    <https://developer.apple.com/library/archive/documentation/FileManagement/Conceptual/FileSystemProgrammingGuide/MacOSXDirectories/MacOSXDirectories.html>`_.
    Makes use of the `appname <platformdirs.api.PlatformDirsABC.appname>`, `version
    <platformdirs.api.PlatformDirsABC.version>`, `ensure_exists <platformdirs.api.PlatformDirsABC.ensure_exists>`.

    XDG environment variables (e.g. ``$XDG_DATA_HOME``) are supported and take precedence over macOS defaults.

    """


__all__ = [
    "MacOS",
]
