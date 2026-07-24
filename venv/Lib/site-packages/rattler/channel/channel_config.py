from __future__ import annotations
import os
from rattler.rattler import PyChannelConfig


class ChannelConfig:
    def __init__(self, channel_alias: str = "https://conda.anaconda.org/", root_dir: str = os.getcwd()) -> None:
        """
        Create a new channel configuration.

        Examples
        --------
        ```python
        >>> channel_config = ChannelConfig()
        >>> channel_config # doctest: +ELLIPSIS
        ChannelConfig(channel_alias="https://conda.anaconda.org/", ...
        >>> channel_config = ChannelConfig("https://repo.prefix.dev/", "/path/to/root/dir")
        >>> channel_config
        ChannelConfig(channel_alias="https://repo.prefix.dev/", root_dir="/path/to/root/dir")
        >>>
        ```
        """
        if root_dir is None:
            # Use the current working directory as the root directory.
            root_dir = os.getcwd()
        self._channel_configuration = PyChannelConfig(channel_alias, root_dir)

    def __repr__(self) -> str:
        """
        Return a string representation of this channel configuration.

        Examples
        --------
        ```python
        >>> channel_config = ChannelConfig()
        >>> channel_config # doctest: +ELLIPSIS
        ChannelConfig(channel_alias="https://conda.anaconda.org/", root_dir="...
        >>>
        ```
        """
        alias = self._channel_configuration.channel_alias
        root_dir = self._channel_configuration.root_dir
        return f'ChannelConfig(channel_alias="{alias}", root_dir="{root_dir}")'
