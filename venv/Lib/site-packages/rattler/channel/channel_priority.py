from enum import Enum

from rattler.rattler import PyChannelPriority


class ChannelPriority(Enum):
    """
    Defines how priority of channels functions during solves. If strict, the channel that the package is first
    found in will be used as the only channel for that package. If disabled, then packages can be retrieved from
    any channel as package version takes precedence.
    """

    Strict = PyChannelPriority.Strict
    Disabled = PyChannelPriority.Disabled
