from hvac.api.vault_api_base import VaultApiBase

class JWT(VaultApiBase):
    DEFAULT_PATH: str
    def resolve_path(self, path): ...
    def configure(
        self,
        oidc_discovery_url=None,
        oidc_discovery_ca_pem=None,
        oidc_client_id=None,
        oidc_client_secret=None,
        oidc_response_mode=None,
        oidc_response_types=None,
        jwks_url=None,
        jwks_ca_pem=None,
        jwt_validation_pubkeys=None,
        bound_issuer=None,
        jwt_supported_algs=None,
        default_role=None,
        provider_config=None,
        path: str | None = None,
        namespace_in_state: bool | None = None,
    ): ...
    def read_config(self, path=None): ...
    def create_role(
        self,
        name,
        user_claim,
        allowed_redirect_uris,
        role_type: str = "jwt",
        bound_audiences=None,
        clock_skew_leeway=None,
        expiration_leeway=None,
        not_before_leeway=None,
        bound_subject=None,
        bound_claims=None,
        groups_claim=None,
        claim_mappings=None,
        oidc_scopes=None,
        bound_claims_type: str = "string",
        verbose_oidc_logging: bool = False,
        token_ttl=None,
        token_max_ttl=None,
        token_policies=None,
        token_bound_cidrs=None,
        token_explicit_max_ttl=None,
        token_no_default_policy=None,
        token_num_uses=None,
        token_period=None,
        token_type=None,
        path=None,
        user_claim_json_pointer=None,
    ): ...
    def read_role(self, name, path=None): ...
    def list_roles(self, path=None): ...
    def delete_role(self, name, path=None): ...
    def oidc_authorization_url_request(self, role, redirect_uri, path=None): ...
    def oidc_callback(self, state, nonce, code, path=None): ...
    def jwt_login(self, role, jwt, use_token: bool = True, path=None): ...
