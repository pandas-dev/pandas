def initialize_all():
    # Import models to register them with the data model manager
    import numba.cuda.models  # noqa: F401

    from numba.cuda.decorators import jit
    from numba.cuda.dispatcher import CUDADispatcher
    from numba.core.target_extension import (target_registry,
                                             dispatcher_registry,
                                             jit_registry)

    cuda_target = target_registry["cuda"]
    jit_registry[cuda_target] = jit
    dispatcher_registry[cuda_target] = CUDADispatcher
