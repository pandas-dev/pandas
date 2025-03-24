"""
Target Descriptors
"""

from abc import ABCMeta, abstractmethod


class TargetDescriptor(metaclass=ABCMeta):

    def __init__(self, target_name):
        self._target_name = target_name

    @property
    @abstractmethod
    def typing_context(self):
        ...

    @property
    @abstractmethod
    def target_context(self):
        ...
