from grpc_reflection.v1alpha import reflection_pb2_grpc

class BaseReflectionServicer(reflection_pb2_grpc.ServerReflectionServicer):
    def __init__(self, service_names, pool=None) -> None: ...

__all__ = ["BaseReflectionServicer"]
