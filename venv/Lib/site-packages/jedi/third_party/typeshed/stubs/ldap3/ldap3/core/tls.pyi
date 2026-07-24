from _typeshed import Incomplete

use_ssl_context: bool

class Tls:
    ssl_options: Incomplete
    validate: Incomplete
    ca_certs_file: Incomplete
    ca_certs_path: Incomplete
    ca_certs_data: Incomplete
    private_key_password: Incomplete
    version: Incomplete
    private_key_file: Incomplete
    certificate_file: Incomplete
    valid_names: Incomplete
    ciphers: Incomplete
    sni: Incomplete
    def __init__(
        self,
        local_private_key_file=None,
        local_certificate_file=None,
        validate=...,
        version=None,
        ssl_options=None,
        ca_certs_file=None,
        valid_names=None,
        ca_certs_path=None,
        ca_certs_data=None,
        local_private_key_password=None,
        ciphers=None,
        sni=None,
    ) -> None: ...
    def wrap_socket(self, connection, do_handshake: bool = False) -> None: ...
    def start_tls(self, connection): ...

def check_hostname(sock, server_name, additional_names) -> None: ...
