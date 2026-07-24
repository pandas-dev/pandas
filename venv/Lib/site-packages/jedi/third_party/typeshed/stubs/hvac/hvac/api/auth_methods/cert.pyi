from hvac.api.vault_api_base import VaultApiBase

class Cert(VaultApiBase):
    def create_ca_certificate_role(
        self,
        name,
        certificate: str = "",
        certificate_file: str = "",
        allowed_common_names: str = "",
        allowed_dns_sans: str = "",
        allowed_email_sans: str = "",
        allowed_uri_sans: str = "",
        allowed_organizational_units: str = "",
        required_extensions: str = "",
        display_name: str = "",
        token_ttl: int = 0,
        token_max_ttl: int = 0,
        token_policies=[],
        token_bound_cidrs=[],
        token_explicit_max_ttl: int = 0,
        token_no_default_policy: bool = False,
        token_num_uses: int = 0,
        token_period: int = 0,
        token_type: str = "",
        mount_point: str = "cert",
    ): ...
    def read_ca_certificate_role(self, name, mount_point: str = "cert"): ...
    def list_certificate_roles(self, mount_point: str = "cert"): ...
    def delete_certificate_role(self, name, mount_point: str = "cert"): ...
    def configure_tls_certificate(self, mount_point: str = "cert", disable_binding: bool = False): ...
    def login(
        self,
        name: str = "",
        cacert: bool = False,
        cert_pem: str = "",
        key_pem: str = "",
        mount_point: str = "cert",
        use_token: bool = True,
    ): ...

    class CertificateAuthError(Exception): ...
