from authlib.oidc.core import UserInfo

def is_openid_scope(scope: str | None) -> bool: ...
def validate_request_prompt(grant, redirect_uri, redirect_fragment: bool = False): ...
def validate_nonce(request, exists_nonce, required: bool = False): ...
def generate_id_token(
    token: dict[str, str | int],
    user_info: UserInfo,
    key: str,
    iss: str,
    aud: list[str],
    alg: str = "RS256",
    exp: int = 3600,
    nonce: str | None = None,
    auth_time: int | None = None,
    acr: str | None = None,
    amr: list[str] | None = None,
    code: str | None = None,
    kid: str | None = None,
) -> str: ...
def create_response_mode_response(redirect_uri, params, response_mode): ...
