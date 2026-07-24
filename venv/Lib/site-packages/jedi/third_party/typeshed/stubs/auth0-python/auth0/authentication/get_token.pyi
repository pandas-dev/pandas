from .base import AuthenticationBase

class GetToken(AuthenticationBase):
    def authorization_code(self, code: str, redirect_uri: str | None, grant_type: str = "authorization_code"): ...
    def authorization_code_pkce(
        self, code_verifier: str, code: str, redirect_uri: str | None, grant_type: str = "authorization_code"
    ): ...
    def client_credentials(self, audience: str, grant_type: str = "client_credentials", organization: str | None = None): ...
    def login(
        self,
        username: str,
        password: str,
        scope: str | None = None,
        realm: str | None = None,
        audience: str | None = None,
        grant_type: str = "http://auth0.com/oauth/grant-type/password-realm",
        forwarded_for: str | None = None,
    ): ...
    def refresh_token(self, refresh_token: str, scope: str = "", grant_type: str = "refresh_token"): ...
    def passwordless_login(self, username: str, otp: str, realm: str, scope: str, audience: str): ...
    def backchannel_login(self, auth_req_id: str, grant_type: str = "urn:openid:params:grant-type:ciba"): ...
    def access_token_for_connection(
        self,
        subject_token_type: str,
        subject_token: str,
        requested_token_type: str,
        connection: str | None = None,
        grant_type: str = ...,
    ): ...
