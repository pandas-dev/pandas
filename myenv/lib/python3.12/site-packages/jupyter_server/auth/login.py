"""Tornado handlers for logging into the Jupyter Server."""

# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.
import os
import re
import uuid
from urllib.parse import urlparse

from tornado.escape import url_escape

from ..base.handlers import JupyterHandler
from .decorator import allow_unauthenticated
from .security import passwd_check, set_password


class LoginFormHandler(JupyterHandler):
    """The basic tornado login handler

    accepts login form, passed to IdentityProvider.process_login_form.
    """

    def _render(self, message=None):
        """Render the login form."""
        self.write(
            self.render_template(
                "login.html",
                next=url_escape(self.get_argument("next", default=self.base_url)),
                message=message,
            )
        )

    def _redirect_safe(self, url, default=None):
        """Redirect if url is on our PATH

        Full-domain redirects are allowed if they pass our CORS origin checks.

        Otherwise use default (self.base_url if unspecified).
        """
        if default is None:
            default = self.base_url
        # protect chrome users from mishandling unescaped backslashes.
        # \ is not valid in urls, but some browsers treat it as /
        # instead of %5C, causing `\\` to behave as `//`
        url = url.replace("\\", "%5C")
        # urllib and browsers interpret extra '/' in the scheme separator (`scheme:///host/path`)
        # differently.
        # urllib gives scheme=scheme, netloc='', path='/host/path', while
        # browsers get scheme=scheme, netloc='host', path='/path'
        # so make sure ':///*' collapses to '://' by splitting and stripping any additional leading slash
        # don't allow any kind of `:/` shenanigans by splitting on ':' only
        # and replacing `:/*` with exactly `://`
        if ":" in url:
            scheme, _, rest = url.partition(":")
            url = f"{scheme}://{rest.lstrip('/')}"
        parsed = urlparse(url)
        # full url may be `//host/path` (empty scheme == same scheme as request)
        # or `https://host/path`
        # or even `https:///host/path` (invalid, but accepted and ambiguously interpreted)
        if (parsed.scheme or parsed.netloc) or not (parsed.path + "/").startswith(self.base_url):
            # require that next_url be absolute path within our path
            allow = False
            # OR pass our cross-origin check
            if parsed.scheme or parsed.netloc:
                # if full URL, run our cross-origin check:
                origin = f"{parsed.scheme}://{parsed.netloc}"
                origin = origin.lower()
                if self.allow_origin:
                    allow = self.allow_origin == origin
                elif self.allow_origin_pat:
                    allow = bool(re.match(self.allow_origin_pat, origin))
            if not allow:
                # not allowed, use default
                self.log.warning("Not allowing login redirect to %r" % url)
                url = default
        self.redirect(url)

    @allow_unauthenticated
    def get(self):
        """Get the login form."""
        if self.current_user:
            next_url = self.get_argument("next", default=self.base_url)
            self._redirect_safe(next_url)
        else:
            self._render()

    @allow_unauthenticated
    def post(self):
        """Post a login."""
        user = self.current_user = self.identity_provider.process_login_form(self)
        if user is None:
            self.set_status(401)
            self._render(message={"error": "Invalid credentials"})
            return

        self.log.info(f"User {user.username} logged in.")
        self.identity_provider.set_login_cookie(self, user)
        next_url = self.get_argument("next", default=self.base_url)
        self._redirect_safe(next_url)


class LegacyLoginHandler(LoginFormHandler):
    """Legacy LoginHandler, implementing most custom auth configuration.

    Deprecated in jupyter-server 2.0.
    Login configuration has moved to IdentityProvider.
    """

    @property
    def hashed_password(self):
        return self.password_from_settings(self.settings)

    def passwd_check(self, a, b):
        """Check a passwd."""
        return passwd_check(a, b)

    @allow_unauthenticated
    def post(self):
        """Post a login form."""
        typed_password = self.get_argument("password", default="")
        new_password = self.get_argument("new_password", default="")

        if self.get_login_available(self.settings):
            if self.passwd_check(self.hashed_password, typed_password) and not new_password:
                self.set_login_cookie(self, uuid.uuid4().hex)
            elif self.token and self.token == typed_password:
                self.set_login_cookie(self, uuid.uuid4().hex)
                if new_password and getattr(self.identity_provider, "allow_password_change", False):
                    config_dir = self.settings.get("config_dir", "")
                    config_file = os.path.join(config_dir, "jupyter_server_config.json")
                    if hasattr(self.identity_provider, "hashed_password"):
                        self.identity_provider.hashed_password = self.settings["password"] = (
                            set_password(new_password, config_file=config_file)
                        )
                    self.log.info("Wrote hashed password to %s" % config_file)
            else:
                self.set_status(401)
                self._render(message={"error": "Invalid credentials"})
                return

        next_url = self.get_argument("next", default=self.base_url)
        self._redirect_safe(next_url)

    @classmethod
    def set_login_cookie(cls, handler, user_id=None):
        """Call this on handlers to set the login cookie for success"""
        cookie_options = handler.settings.get("cookie_options", {})
        cookie_options.setdefault("httponly", True)
        # tornado <4.2 has a bug that considers secure==True as soon as
        # 'secure' kwarg is passed to set_secure_cookie
        if handler.settings.get("secure_cookie", handler.request.protocol == "https"):
            cookie_options.setdefault("secure", True)
        cookie_options.setdefault("path", handler.base_url)
        handler.set_secure_cookie(handler.cookie_name, user_id, **cookie_options)
        return user_id

    auth_header_pat = re.compile(r"token\s+(.+)", re.IGNORECASE)

    @classmethod
    def get_token(cls, handler):
        """Get the user token from a request

        Default:

        - in URL parameters: ?token=<token>
        - in header: Authorization: token <token>
        """

        user_token = handler.get_argument("token", "")
        if not user_token:
            # get it from Authorization header
            m = cls.auth_header_pat.match(handler.request.headers.get("Authorization", ""))
            if m:
                user_token = m.group(1)
        return user_token

    @classmethod
    def should_check_origin(cls, handler):
        """DEPRECATED in 2.0, use IdentityProvider API"""
        return not cls.is_token_authenticated(handler)

    @classmethod
    def is_token_authenticated(cls, handler):
        """DEPRECATED in 2.0, use IdentityProvider API"""
        if getattr(handler, "_user_id", None) is None:
            # ensure get_user has been called, so we know if we're token-authenticated
            handler.current_user  # noqa: B018
        return getattr(handler, "_token_authenticated", False)

    @classmethod
    def get_user(cls, handler):
        """DEPRECATED in 2.0, use IdentityProvider API"""
        # Can't call this get_current_user because it will collide when
        # called on LoginHandler itself.
        if getattr(handler, "_user_id", None):
            return handler._user_id
        token_user_id = cls.get_user_token(handler)
        cookie_user_id = cls.get_user_cookie(handler)
        # prefer token to cookie if both given,
        # because token is always explicit
        user_id = token_user_id or cookie_user_id
        if token_user_id:
            # if token-authenticated, persist user_id in cookie
            # if it hasn't already been stored there
            if user_id != cookie_user_id:
                cls.set_login_cookie(handler, user_id)
            # Record that the current request has been authenticated with a token.
            # Used in is_token_authenticated above.
            handler._token_authenticated = True

        if user_id is None:
            # If an invalid cookie was sent, clear it to prevent unnecessary
            # extra warnings. But don't do this on a request with *no* cookie,
            # because that can erroneously log you out (see gh-3365)
            if handler.get_cookie(handler.cookie_name) is not None:
                handler.log.warning("Clearing invalid/expired login cookie %s", handler.cookie_name)
                handler.clear_login_cookie()
            if not handler.login_available:
                # Completely insecure! No authentication at all.
                # No need to warn here, though; validate_security will have already done that.
                user_id = "anonymous"

        # cache value for future retrievals on the same request
        handler._user_id = user_id
        return user_id

    @classmethod
    def get_user_cookie(cls, handler):
        """DEPRECATED in 2.0, use IdentityProvider API"""
        get_secure_cookie_kwargs = handler.settings.get("get_secure_cookie_kwargs", {})
        user_id = handler.get_secure_cookie(handler.cookie_name, **get_secure_cookie_kwargs)
        if user_id:
            user_id = user_id.decode()
        return user_id

    @classmethod
    def get_user_token(cls, handler):
        """DEPRECATED in 2.0, use IdentityProvider API"""
        token = handler.token
        if not token:
            return None
        # check login token from URL argument or Authorization header
        user_token = cls.get_token(handler)
        authenticated = False
        if user_token == token:
            # token-authenticated, set the login cookie
            handler.log.debug(
                "Accepting token-authenticated connection from %s",
                handler.request.remote_ip,
            )
            authenticated = True

        if authenticated:
            # token does not correspond to user-id,
            # which is stored in a cookie.
            # still check the cookie for the user id
            user_id = cls.get_user_cookie(handler)
            if user_id is None:
                # no cookie, generate new random user_id
                user_id = uuid.uuid4().hex
                handler.log.info(
                    f"Generating new user_id for token-authenticated request: {user_id}"
                )
            return user_id
        else:
            return None

    @classmethod
    def validate_security(cls, app, ssl_options=None):
        """DEPRECATED in 2.0, use IdentityProvider API"""
        if not app.ip:
            warning = "WARNING: The Jupyter server is listening on all IP addresses"
            if ssl_options is None:
                app.log.warning(f"{warning} and not using encryption. This is not recommended.")
            if not app.password and not app.token:
                app.log.warning(
                    f"{warning} and not using authentication. "
                    "This is highly insecure and not recommended."
                )
        elif not app.password and not app.token:
            app.log.warning(
                "All authentication is disabled."
                "  Anyone who can connect to this server will be able to run code."
            )

    @classmethod
    def password_from_settings(cls, settings):
        """DEPRECATED in 2.0, use IdentityProvider API"""
        return settings.get("password", "")

    @classmethod
    def get_login_available(cls, settings):
        """DEPRECATED in 2.0, use IdentityProvider API"""

        return bool(cls.password_from_settings(settings) or settings.get("token"))


# deprecated import, so deprecated implementations get the Legacy class instead
LoginHandler = LegacyLoginHandler
