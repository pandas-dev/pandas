from __future__ import annotations

from rattler.rattler import PyLockChannel


class LockChannel:
    _channel: PyLockChannel

    def __init__(self, url: str) -> None:
        """
        Create a new channel.
        """
        self._channel = PyLockChannel(url)

    @classmethod
    def _from_py_lock_channel(cls, channel: PyLockChannel) -> LockChannel:
        """
        Construct Rattler LockChannel from FFI PyLockChannel object.
        """
        chan = cls.__new__(cls)
        chan._channel = channel
        return chan

    def __str__(self) -> str:
        """
        Returns a string representation of the LockChannel.

        Examples
        --------
        ```python
        >>> channel = LockChannel("https://conda.anaconda.org/conda-forge/")
        >>> str(channel)
        'https://conda.anaconda.org/conda-forge/'
        >>>
        ```
        """
        return self._channel.as_str()

    def __repr__(self) -> str:
        """
        Returns a representation of the LockChannel.

        Examples
        --------
        ```python
        >>> channel = LockChannel("https://conda.anaconda.org/conda-forge/")
        >>> channel
        LockChannel(url="https://conda.anaconda.org/conda-forge/")
        >>>
        ```
        """
        return f'LockChannel(url="{self._channel.as_str()}")'
