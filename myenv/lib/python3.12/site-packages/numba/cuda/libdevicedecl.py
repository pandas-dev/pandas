from numba.cuda import libdevice, libdevicefuncs
from numba.core.typing.templates import ConcreteTemplate, Registry

registry = Registry()
register_global = registry.register_global


def libdevice_declare(func, retty, args):
    class Libdevice_function(ConcreteTemplate):
        cases = [libdevicefuncs.create_signature(retty, args)]

    pyfunc = getattr(libdevice, func[5:])
    register_global(pyfunc)(Libdevice_function)


for func, (retty, args) in libdevicefuncs.functions.items():
    libdevice_declare(func, retty, args)
