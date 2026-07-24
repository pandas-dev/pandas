from hvac.api.auth_methods.jwt import JWT

class OIDC(JWT):
    DEFAULT_PATH: str
    def create_role(
        self,
        name,
        user_claim,
        allowed_redirect_uris,
        role_type: str = "oidc",
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
    ) -> None: ...
