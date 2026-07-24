from hvac.api.vault_api_base import VaultApiBase

DEFAULT_MOUNT_POINT: str

class Transit(VaultApiBase):
    def create_key(
        self,
        name,
        convergent_encryption=None,
        derived=None,
        exportable=None,
        allow_plaintext_backup=None,
        key_type=None,
        mount_point="transit",
        auto_rotate_period=None,
    ): ...
    def read_key(self, name, mount_point="transit"): ...
    def list_keys(self, mount_point="transit"): ...
    def delete_key(self, name, mount_point="transit"): ...
    def update_key_configuration(
        self,
        name,
        min_decryption_version=None,
        min_encryption_version=None,
        deletion_allowed=None,
        exportable=None,
        allow_plaintext_backup=None,
        mount_point="transit",
        auto_rotate_period=None,
    ): ...
    def rotate_key(self, name, mount_point="transit"): ...
    def export_key(self, name, key_type, version=None, mount_point="transit"): ...
    def encrypt_data(
        self,
        name,
        plaintext=None,
        context=None,
        key_version=None,
        nonce=None,
        batch_input=None,
        type=None,
        convergent_encryption=None,
        mount_point: str = "transit",
        associated_data: str | None = None,
    ): ...
    def decrypt_data(
        self,
        name,
        ciphertext=None,
        context=None,
        nonce=None,
        batch_input=None,
        mount_point: str = "transit",
        associated_data: str | None = None,
    ): ...
    def rewrap_data(
        self, name, ciphertext, context=None, key_version=None, nonce=None, batch_input=None, mount_point="transit"
    ): ...
    def generate_data_key(self, name, key_type, context=None, nonce=None, bits=None, mount_point="transit"): ...
    def generate_random_bytes(self, n_bytes=None, output_format=None, mount_point="transit"): ...
    def hash_data(self, hash_input, algorithm=None, output_format=None, mount_point="transit"): ...
    def generate_hmac(self, name, hash_input, key_version=None, algorithm=None, mount_point="transit"): ...
    def sign_data(
        self,
        name,
        hash_input=None,
        key_version=None,
        hash_algorithm=None,
        context=None,
        prehashed=None,
        signature_algorithm=None,
        marshaling_algorithm=None,
        salt_length=None,
        mount_point="transit",
        batch_input=None,
    ): ...
    def verify_signed_data(
        self,
        name,
        hash_input,
        signature=None,
        hmac=None,
        hash_algorithm=None,
        context=None,
        prehashed=None,
        signature_algorithm=None,
        salt_length=None,
        marshaling_algorithm=None,
        mount_point="transit",
    ): ...
    def backup_key(self, name, mount_point="transit"): ...
    def restore_key(self, backup, name=None, force=None, mount_point="transit"): ...
    def trim_key(self, name, min_version, mount_point="transit"): ...
