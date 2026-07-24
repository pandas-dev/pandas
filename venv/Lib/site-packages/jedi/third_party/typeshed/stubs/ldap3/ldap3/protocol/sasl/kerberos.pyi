posix_gssapi_unavailable: bool
windows_gssapi_unavailable: bool
NO_SECURITY_LAYER: int
INTEGRITY_PROTECTION: int
CONFIDENTIALITY_PROTECTION: int

def get_channel_bindings(ssl_socket): ...
def sasl_gssapi(connection, controls): ...
