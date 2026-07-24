import logging

from hvac.api.vault_api_base import VaultApiBase

DEFAULT_MOUNT_POINT: str

logger: logging.Logger

class Gcp(VaultApiBase):
    def configure(
        self, credentials=None, google_certs_endpoint="https://www.googleapis.com/oauth2/v3/certs", mount_point="gcp"
    ): ...
    def read_config(self, mount_point="gcp"): ...
    def delete_config(self, mount_point="gcp"): ...
    def create_role(
        self,
        name,
        role_type,
        project_id,
        ttl=None,
        max_ttl=None,
        period=None,
        policies=None,
        bound_service_accounts=None,
        max_jwt_exp=None,
        allow_gce_inference=None,
        bound_zones=None,
        bound_regions=None,
        bound_instance_groups=None,
        bound_labels=None,
        mount_point="gcp",
    ): ...
    def edit_service_accounts_on_iam_role(self, name, add=None, remove=None, mount_point="gcp"): ...
    def edit_labels_on_gce_role(self, name, add=None, remove=None, mount_point="gcp"): ...
    def read_role(self, name, mount_point="gcp"): ...
    def list_roles(self, mount_point="gcp"): ...
    def delete_role(self, role, mount_point="gcp"): ...
    def login(self, role, jwt, use_token: bool = True, mount_point="gcp"): ...
