# compatibility shim for ipykernel < 6.18
import sys
from IPython import get_ipython
import comm


def requires_ipykernel_shim():
    if "ipykernel" in sys.modules:
        import ipykernel

        version = ipykernel.version_info
        return version < (6, 18)
    else:
        return False


def get_comm_manager():
    if requires_ipykernel_shim():
        ip = get_ipython()

        if ip is not None and getattr(ip, "kernel", None) is not None:
            return get_ipython().kernel.comm_manager
    else:
        return comm.get_comm_manager()


def create_comm(*args, **kwargs):
    if requires_ipykernel_shim():
        from ipykernel.comm import Comm

        return Comm(*args, **kwargs)
    else:
        return comm.create_comm(*args, **kwargs)
