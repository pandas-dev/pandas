from tensorflow.python.trackable.base import Trackable

class _ResourceMetaclass(type): ...

# Internal type that is commonly used as a base class
# it is needed for the public signatures of some APIs.
class CapturableResource(Trackable, metaclass=_ResourceMetaclass): ...

def __getattr__(name: str): ...  # incomplete module
