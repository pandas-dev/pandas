import contextlib
#from functools import cached_property
from numba.core.utils import threadsafe_cached_property as cached_property

from numba.core.descriptors import TargetDescriptor
from numba.core import utils, typing, dispatcher, cpu

# -----------------------------------------------------------------------------
# Default CPU target descriptors


class CPUTarget(TargetDescriptor):
    options = cpu.CPUTargetOptions

    @cached_property
    def _toplevel_target_context(self):
        # Lazily-initialized top-level target context, for all threads
        return cpu.CPUContext(self.typing_context, self._target_name)

    @cached_property
    def _toplevel_typing_context(self):
        # Lazily-initialized top-level typing context, for all threads
        return typing.Context()

    @property
    def target_context(self):
        """
        The target context for CPU targets.
        """
        return self._toplevel_target_context

    @property
    def typing_context(self):
        """
        The typing context for CPU targets.
        """
        return self._toplevel_typing_context


# The global CPU target
cpu_target = CPUTarget('cpu')


class CPUDispatcher(dispatcher.Dispatcher):
    targetdescr = cpu_target


class DelayedRegistry(utils.UniqueDict):
    """
    A unique dictionary but with deferred initialisation of the values.

    Attributes
    ----------
    ondemand:

        A dictionary of key -> value, where value is executed
        the first time it is is used.  It is used for part of a deferred
        initialization strategy.
    """
    def __init__(self, *args, **kws):
        self.ondemand = utils.UniqueDict()
        self.key_type = kws.pop('key_type', None)
        self.value_type = kws.pop('value_type', None)
        self._type_check = self.key_type or self.value_type
        super(DelayedRegistry, self).__init__(*args, **kws)

    def __getitem__(self, item):
        if item in self.ondemand:
            self[item] = self.ondemand[item]()
            del self.ondemand[item]
        return super(DelayedRegistry, self).__getitem__(item)

    def __setitem__(self, key, value):
        if self._type_check:
            def check(x, ty_x):
                if isinstance(ty_x, type):
                    assert ty_x in x.__mro__, (x, ty_x)
                else:
                    assert isinstance(x, ty_x), (x, ty_x)
            if self.key_type is not None:
                check(key, self.key_type)
            if self.value_type is not None:
                check(value, self.value_type)
        return super(DelayedRegistry, self).__setitem__(key, value)
