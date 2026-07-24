from __future__ import annotations
from typing import Optional, TYPE_CHECKING

if TYPE_CHECKING:
    from rattler.lock.channel import LockChannel

from rattler.rattler import PyChannel
from rattler.channel.channel_config import ChannelConfig


class Channel:
    def __init__(self, name: str, channel_configuration: Optional[ChannelConfig] = None) -> None:
        """
        Create a new channel.

        ```python
        >>> channel = Channel("conda-forge")
        >>> channel
        Channel(name="conda-forge", base_url="https://conda.anaconda.org/conda-forge/")
        >>>
        ```
        """
        if not channel_configuration:
            channel_configuration = ChannelConfig()

        self._channel = PyChannel(name, channel_configuration._channel_configuration)

    @classmethod
    def _from_py_channel(cls, py_channel: PyChannel) -> Channel:
        channel = cls.__new__(cls)
        channel._channel = py_channel
        return channel

    def to_lock_channel(self) -> LockChannel:
        """
        Returns a new [`LockChannel`] from existing channel.

        ```python
        >>> channel = Channel("conda-forge")
        >>> channel.to_lock_channel()
        LockChannel(url="https://conda.anaconda.org/conda-forge/")
        >>>
        ```
        """
        from rattler.lock.channel import LockChannel

        return LockChannel(self.base_url)

    @property
    def name(self) -> Optional[str]:
        """
        Return the name of this channel.

        Examples
        --------
        ```python
        >>> channel = Channel("conda-forge")
        >>> channel.name
        'conda-forge'
        >>>
        ```
        """
        return self._channel.name

    @property
    def base_url(self) -> str:
        """
        Return the base URL of this channel.

        Examples
        --------
        ```python
        >>> channel = Channel("conda-forge")
        >>> channel.base_url
        'https://conda.anaconda.org/conda-forge/'
        >>>
        ```
        """
        return self._channel.base_url

    def __repr__(self) -> str:
        """
        Return a string representation of this channel.

        Examples
        --------
        ```python
        >>> channel = Channel("conda-forge")
        >>> channel
        Channel(name="conda-forge", base_url="https://conda.anaconda.org/conda-forge/")
        >>>
        ```
        """
        return f'Channel(name="{self.name}", base_url="{self.base_url}")'
