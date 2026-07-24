from hvac.api.vault_api_base import VaultApiBase

DEFAULT_MOUNT_POINT: str

class Ssh(VaultApiBase):
    def create_or_update_key(self, name: str = "", key: str = "", mount_point: str = "ssh"): ...
    def delete_key(self, name: str = "", mount_point: str = "ssh"): ...
    def create_role(
        self,
        name: str = "",
        key: str = "",
        admin_user: str = "",
        default_user: str = "",
        cidr_list: str = "",
        exclude_cidr_list: str = "",
        port: int = 22,
        key_type: str = "",
        key_bits: int = 1024,
        install_script: str = "",
        allowed_users: str = "",
        allowed_users_template: str = "",
        allowed_domains: str = "",
        key_option_specs: str = "",
        ttl: str = "",
        max_ttl: str = "",
        allowed_critical_options: str = "",
        allowed_extensions: str = "",
        default_critical_options=None,
        default_extensions=None,
        allow_user_certificates: str = "",
        allow_host_certificates: bool = False,
        allow_bare_domains: bool = False,
        allow_subdomains: bool = False,
        allow_user_key_ids: bool = False,
        key_id_format: str = "",
        allowed_user_key_lengths=None,
        algorithm_signer: str = "",
        mount_point="ssh",
    ): ...
    def read_role(self, name: str = "", mount_point: str = "ssh"): ...
    def list_roles(self, mount_point: str = "ssh"): ...
    def delete_role(self, name: str = "", mount_point: str = "ssh"): ...
    def list_zeroaddress_roles(self, mount_point: str = "ssh"): ...
    def configure_zeroaddress_roles(self, roles: str = "", mount_point: str = "ssh"): ...
    def delete_zeroaddress_role(self, mount_point: str = "ssh"): ...
    def generate_ssh_credentials(self, name: str = "", username: str = "", ip: str = "", mount_point: str = "ssh"): ...
    def list_roles_by_ip(self, ip: str = "", mount_point: str = "ssh"): ...
    def verify_ssh_otp(self, otp, mount_point="ssh"): ...
    def submit_ca_information(
        self,
        private_key: str = "",
        public_key: str = "",
        generate_signing_key: bool = True,
        key_type: str = "ssh-rsa",
        key_bits: int = 0,
        mount_point: str = "ssh",
    ): ...
    def delete_ca_information(self, mount_point: str = "ssh"): ...
    def read_public_key(self, mount_point: str = "ssh"): ...
    def sign_ssh_key(
        self,
        name: str = "",
        public_key: str = "",
        ttl: str = "",
        valid_principals: str = "",
        cert_type: str = "user",
        key_id: str = "",
        critical_options=None,
        extensions=None,
        mount_point: str = "ssh",
    ): ...
