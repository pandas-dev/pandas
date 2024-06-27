from numba import types
from numba.cuda.stubs import _vector_type_stubs


class SimulatedVectorType:
    attributes = ['x', 'y', 'z', 'w']

    def __init__(self, *args):
        args_flattened = []
        for arg in args:
            if isinstance(arg, SimulatedVectorType):
                args_flattened += arg.as_list()
            else:
                args_flattened.append(arg)
        self._attrs = self.attributes[:len(args_flattened)]
        if not self.num_elements == len(args_flattened):
            raise TypeError(
                f"{self.name} expects {self.num_elements}"
                f" elements, got {len(args_flattened)}"
            )

        for arg, attr in zip(args_flattened, self._attrs):
            setattr(self, attr, arg)

    @property
    def name(self):
        raise NotImplementedError()

    @property
    def num_elements(self):
        raise NotImplementedError()

    def as_list(self):
        return [getattr(self, attr) for attr in self._attrs]


def make_simulated_vector_type(num_elements, name):
    obj = type(name, (SimulatedVectorType,), {
        "num_elements": num_elements,
        "base_type": types.float32,
        "name": name
    })
    obj.user_facing_object = obj
    return obj


def _initialize():
    _simulated_vector_types = {}
    for stub in _vector_type_stubs:
        num_elements = int(stub.__name__[-1])
        _simulated_vector_types[stub.__name__] = (
            make_simulated_vector_type(num_elements, stub.__name__)
        )
        _simulated_vector_types[stub.__name__].aliases = stub.aliases
    return _simulated_vector_types


vector_types = _initialize()
