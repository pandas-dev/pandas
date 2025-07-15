"""Identity Provider interface

This defines the _authentication_ layer of Jupyter Server,
to be used in combination with Authorizer for _authorization_.

.. versionadded:: 2.0
"""

from __future__ import annotations

import binascii
import datetime
import json
import os
import re
import sys
import typing as t
import uuid
from dataclasses import asdict, dataclass
from http.cookies import Morsel

from tornado import escape, httputil, web
from traitlets import Bool, Dict, Type, Unicode, default
from traitlets.config import LoggingConfigurable

from jupyter_server.transutils import _i18n

from .security import passwd_check, set_password
from .utils import get_anonymous_username

_non_alphanum = re.compile(r"[^A-Za-z0-9]")


@dataclass
class User:
    """Object representing a User

    This or a subclass should be returned from IdentityProvider.get_user
    """

    username: str  # the only truly required field

    # these fields are filled from username if not specified
    # name is the 'real' name of the user
    name: str = ""
    # display_name is a shorter name for us in UI,
    # if different from name. e.g. a nickname
    display_name: str = ""

    # these fields are left as None if undefined
    initials: str | None = None
    avatar_url: str | None = None
    color: str | None = None

    # TODO: extension fields?
    # ext: Dict[str, Dict[str, Any]] = field(default_factory=dict)

    def __post_init__(self):
        self.fill_defaults()

    def fill_defaults(self):
        """Fill out default fields in the identity model

        - Ensures all values are defined
        - Fills out derivative values for name fields fields
        - Fills out null values for optional fields
        """

        # username is the only truly required field
        if not self.username:
            msg = f"user.username must not be empty: {self}"
            raise ValueError(msg)

        # derive name fields from username -> name -> display name
        if not self.name:
            self.name = self.username
        if not self.display_name:
            self.display_name = self.name


def _backward_compat_user(got_user: t.Any) -> User:
    """Backward-compatibility for LoginHandler.get_user

    Prior to 2.0, LoginHandler.get_user could return anything truthy.

    Typically, this was either a simple string username,
    or a simple dict.

    Make some effort to allow common patterns to keep working.
    """
    if isinstance(got_user, str):
        return User(username=got_user)
    elif isinstance(got_user, dict):
        kwargs = {}
        if "username" not in got_user and "name" in got_user:
            kwargs["username"] = got_user["name"]
        for field in User.__dataclass_fields__:
            if field in got_user:
                kwargs[field] = got_user[field]
        try:
            return User(**kwargs)
        except TypeError:
            msg = f"Unrecognized user: {got_user}"
            raise ValueError(msg) from None
    else:
        msg = f"Unrecognized user: {got_user}"
        raise ValueError(msg)


class IdentityProvider(LoggingConfigurable):
    """
    Interface for providing identity management and authentication.

    Two principle methods:

    - :meth:`~jupyter_server.auth.IdentityProvider.get_user` returns a :class:`~.User` object
      for successful authentication, or None for no-identity-found.
    - :meth:`~jupyter_server.auth.IdentityProvider.identity_model` turns a :class:`~jupyter_server.auth.User` into a JSONable dict.
      The default is to use :py:meth:`dataclasses.asdict`,
      and usually shouldn't need override.

    Additional methods can customize authentication.

    .. versionadded:: 2.0
    """

    cookie_name: str | Unicode[str, str | bytes] = Unicode(
        "",
        config=True,
        help=_i18n("Name of the cookie to set for persisting login. Default: username-${Host}."),
    )

    cookie_options = Dict(
        config=True,
        help=_i18n(
            "Extra keyword arguments to pass to `set_secure_cookie`."
            " See tornado's set_secure_cookie docs for details."
        ),
    )

    secure_cookie: bool | Bool[bool | None, bool | int | None] = Bool(
        None,
        allow_none=True,
        config=True,
        help=_i18n(
            "Specify whether login cookie should have the `secure` property (HTTPS-only)."
            "Only needed when protocol-detection gives the wrong answer due to proxies."
        ),
    )

    get_secure_cookie_kwargs = Dict(
        config=True,
        help=_i18n(
            "Extra keyword arguments to pass to `get_secure_cookie`."
            " See tornado's get_secure_cookie docs for details."
        ),
    )

    token: str | Unicode[str, str | bytes] = Unicode(
        "<generated>",
        help=_i18n(
            """Token used for authenticating first-time connections to the server.

        The token can be read from the file referenced by JUPYTER_TOKEN_FILE or set directly
        with the JUPYTER_TOKEN environment variable.

        When no password is enabled,
        the default is to generate a new, random token.

        Setting to an empty string disables authentication altogether, which is NOT RECOMMENDED.

        Prior to 2.0: configured as ServerApp.token
        """
        ),
    ).tag(config=True)

    login_handler_class = Type(
        default_value="jupyter_server.auth.login.LoginFormHandler",
        klass=web.RequestHandler,
        config=True,
        help=_i18n("The login handler class to use, if any."),
    )

    logout_handler_class = Type(
        default_value="jupyter_server.auth.logout.LogoutHandler",
        klass=web.RequestHandler,
        config=True,
        help=_i18n("The logout handler class to use."),
    )

    token_generated = False

    @default("token")
    def _token_default(self):
        if os.getenv("JUPYTER_TOKEN"):
            self.token_generated = False
            return os.environ["JUPYTER_TOKEN"]
        if os.getenv("JUPYTER_TOKEN_FILE"):
            self.token_generated = False
            with open(os.environ["JUPYTER_TOKEN_FILE"]) as token_file:
                return token_file.read()
        if not self.need_token:
            # no token if password is enabled
            self.token_generated = False
            return ""
        else:
            self.token_generated = True
            return binascii.hexlify(os.urandom(24)).decode("ascii")

    need_token: bool | Bool[bool, t.Union[bool, int]] = Bool(True)

    def get_user(self, handler: web.RequestHandler) -> User | None | t.Awaitable[User | None]:
        """Get the authenticated user for a request

        Must return a :class:`jupyter_server.auth.User`,
        though it may be a subclass.

        Return None if the request is not authenticated.

        _may_ be a coroutine
        """
        return self._get_user(handler)

    # not sure how to have optional-async type signature
    # on base class with `async def` without splitting it into two methods

    async def _get_user(self, handler: web.RequestHandler) -> User | None:
        """Get the user."""
        if getattr(handler, "_jupyter_current_user", None):
            # already authenticated
            return t.cast(User, handler._jupyter_current_user)  # type:ignore[attr-defined]
        _token_user: User | None | t.Awaitable[User | None] = self.get_user_token(handler)
        if isinstance(_token_user, t.Awaitable):
            _token_user = await _token_user
        token_user: User | None = _token_user  # need second variable name to collapse type
        _cookie_user = self.get_user_cookie(handler)
        if isinstance(_cookie_user, t.Awaitable):
            _cookie_user = await _cookie_user
        cookie_user: User | None = _cookie_user
        # prefer token to cookie if both given,
        # because token is always explicit
        user = token_user or cookie_user

        if user is not None and token_user is not None:
            # if token-authenticated, persist user_id in cookie
            # if it hasn't already been stored there
            if user != cookie_user:
                self.set_login_cookie(handler, user)
            # Record that the current request has been authenticated with a token.
            # Used in is_token_authenticated above.
            handler._token_authenticated = True  # type:ignore[attr-defined]

        if user is None:
            # If an invalid cookie was sent, clear it to prevent unnecessary
            # extra warnings. But don't do this on a request with *no* cookie,
            # because that can erroneously log you out (see gh-3365)
            cookie_name = self.get_cookie_name(handler)
            cookie = handler.get_cookie(cookie_name)
            if cookie is not None:
                self.log.warning(f"Clearing invalid/expired login cookie {cookie_name}")
                self.clear_login_cookie(handler)
            if not self.auth_enabled:
                # Completely insecure! No authentication at all.
                # No need to warn here, though; validate_security will have already done that.
                user = self.generate_anonymous_user(handler)
                # persist user on first request
                # so the user data is stable for a given browser session
                self.set_login_cookie(handler, user)

        return user

    def identity_model(self, user: User) -> dict[str, t.Any]:
        """Return a User as an Identity model"""
        # TODO: validate?
        return asdict(user)

    def get_handlers(self) -> list[tuple[str, object]]:
        """Return list of additional handlers for this identity provider

        For example, an OAuth callback handler.
        """
        handlers = []
        if self.login_available:
            handlers.append((r"/login", self.login_handler_class))
        if self.logout_available:
            handlers.append((r"/logout", self.logout_handler_class))
        return handlers

    def user_to_cookie(self, user: User) -> str:
        """Serialize a user to a string for storage in a cookie

        If overriding in a subclass, make sure to define user_from_cookie as well.

        Default is just the user's username.
        """
        # default: username is enough
        cookie = json.dumps(
            {
                "username": user.username,
                "name": user.name,
                "display_name": user.display_name,
                "initials": user.initials,
                "color": user.color,
            }
        )
        return cookie

    def user_from_cookie(self, cookie_value: str) -> User | None:
        """Inverse of user_to_cookie"""
        user = json.loads(cookie_value)
        return User(
            user["username"],
            user["name"],
            user["display_name"],
            user["initials"],
            None,
            user["color"],
        )

    def get_cookie_name(self, handler: web.RequestHandler) -> str:
        """Return the login cookie name

        Uses IdentityProvider.cookie_name, if defined.
        Default is to generate a string taking host into account to avoid
        collisions for multiple servers on one hostname with different ports.
        """
        if self.cookie_name:
            return self.cookie_name
        else:
            return _non_alphanum.sub("-", f"username-{handler.request.host}")

    def set_login_cookie(self, handler: web.RequestHandler, user: User) -> None:
        """Call this on handlers to set the login cookie for success"""
        cookie_options = {}
        cookie_options.update(self.cookie_options)
        cookie_options.setdefault("httponly", True)
        # tornado <4.2 has a bug that considers secure==True as soon as
        # 'secure' kwarg is passed to set_secure_cookie
        secure_cookie = self.secure_cookie
        if secure_cookie is None:
            secure_cookie = handler.request.protocol == "https"
        if secure_cookie:
            cookie_options.setdefault("secure", True)
        cookie_options.setdefault("path", handler.base_url)  # type:ignore[attr-defined]
        cookie_name = self.get_cookie_name(handler)
        handler.set_secure_cookie(cookie_name, self.user_to_cookie(user), **cookie_options)

    def _force_clear_cookie(
        self, handler: web.RequestHandler, name: str, path: str = "/", domain: str | None = None
    ) -> None:
        """Deletes the cookie with the given name.

        Tornado's cookie handling currently (Jan 2018) stores cookies in a dict
        keyed by name, so it can only modify one cookie with a given name per
        response. The browser can store multiple cookies with the same name
        but different domains and/or paths. This method lets us clear multiple
        cookies with the same name.

        Due to limitations of the cookie protocol, you must pass the same
        path and domain to clear a cookie as were used when that cookie
        was set (but there is no way to find out on the server side
        which values were used for a given cookie).
        """
        name = escape.native_str(name)
        expires = datetime.datetime.now(tz=datetime.timezone.utc) - datetime.timedelta(days=365)

        morsel: Morsel[t.Any] = Morsel()
        morsel.set(name, "", '""')
        morsel["expires"] = httputil.format_timestamp(expires)
        morsel["path"] = path
        if domain:
            morsel["domain"] = domain
        handler.add_header("Set-Cookie", morsel.OutputString())

    def clear_login_cookie(self, handler: web.RequestHandler) -> None:
        """Clear the login cookie, effectively logging out the session."""
        cookie_options = {}
        cookie_options.update(self.cookie_options)
        path = cookie_options.setdefault("path", handler.base_url)  # type:ignore[attr-defined]
        cookie_name = self.get_cookie_name(handler)
        handler.clear_cookie(cookie_name, path=path)
        if path and path != "/":
            # also clear cookie on / to ensure old cookies are cleared
            # after the change in path behavior.
            # N.B. This bypasses the normal cookie handling, which can't update
            # two cookies with the same name. See the method above.
            self._force_clear_cookie(handler, cookie_name)

    def get_user_cookie(
        self, handler: web.RequestHandler
    ) -> User | None | t.Awaitable[User | None]:
        """Get user from a cookie

        Calls user_from_cookie to deserialize cookie value
        """
        _user_cookie = handler.get_secure_cookie(
            self.get_cookie_name(handler),
            **self.get_secure_cookie_kwargs,
        )
        if not _user_cookie:
            return None
        user_cookie = _user_cookie.decode()
        # TODO: try/catch in case of change in config?
        try:
            return self.user_from_cookie(user_cookie)
        except Exception as e:
            # log bad cookie itself, only at debug-level
            self.log.debug(f"Error unpacking user from cookie: cookie={user_cookie}", exc_info=True)
            self.log.error(f"Error unpacking user from cookie: {e}")
            return None

    auth_header_pat = re.compile(r"(token|bearer)\s+(.+)", re.IGNORECASE)

    def get_token(self, handler: web.RequestHandler) -> str | None:
        """Get the user token from a request

        Default:

        - in URL parameters: ?token=<token>
        - in header: Authorization: token <token>
        """
        user_token = handler.get_argument("token", "")
        if not user_token:
            # get it from Authorization header
            m = self.auth_header_pat.match(handler.request.headers.get("Authorization", ""))
            if m:
                user_token = m.group(2)
        return user_token

    async def get_user_token(self, handler: web.RequestHandler) -> User | None:
        """Identify the user based on a token in the URL or Authorization header

        Returns:
        - uuid if authenticated
        - None if not
        """
        token = t.cast("str | None", handler.token)  # type:ignore[attr-defined]
        if not token:
            return None
        # check login token from URL argument or Authorization header
        user_token = self.get_token(handler)
        authenticated = False
        if user_token == token:
            # token-authenticated, set the login cookie
            self.log.debug(
                "Accepting token-authenticated request from %s",
                handler.request.remote_ip,
            )
            authenticated = True

        if authenticated:
            # token does not correspond to user-id,
            # which is stored in a cookie.
            # still check the cookie for the user id
            _user = self.get_user_cookie(handler)
            if isinstance(_user, t.Awaitable):
                _user = await _user
            user: User | None = _user
            if user is None:
                user = self.generate_anonymous_user(handler)
            return user
        else:
            return None

    def generate_anonymous_user(self, handler: web.RequestHandler) -> User:
        """Generate a random anonymous user.

        For use when a single shared token is used,
        but does not identify a user.
        """
        user_id = uuid.uuid4().hex
        moon = get_anonymous_username()
        name = display_name = f"Anonymous {moon}"
        initials = f"A{moon[0]}"
        color = None
        handler.log.debug(f"Generating new user for token-authenticated request: {user_id}")  # type:ignore[attr-defined]
        return User(user_id, name, display_name, initials, None, color)

    def should_check_origin(self, handler: web.RequestHandler) -> bool:
        """Should the Handler check for CORS origin validation?

        Origin check should be skipped for token-authenticated requests.

        Returns:
        - True, if Handler must check for valid CORS origin.
        - False, if Handler should skip origin check since requests are token-authenticated.
        """
        return not self.is_token_authenticated(handler)

    def is_token_authenticated(self, handler: web.RequestHandler) -> bool:
        """Returns True if handler has been token authenticated. Otherwise, False.

        Login with a token is used to signal certain things, such as:

        - permit access to REST API
        - xsrf protection
        - skip origin-checks for scripts
        """
        # ensure get_user has been called, so we know if we're token-authenticated
        handler.current_user  # noqa: B018
        return getattr(handler, "_token_authenticated", False)

    def validate_security(
        self,
        app: t.Any,
        ssl_options: dict[str, t.Any] | None = None,
    ) -> None:
        """Check the application's security.

        Show messages, or abort if necessary, based on the security configuration.
        """
        if not app.ip:
            warning = "WARNING: The Jupyter server is listening on all IP addresses"
            if ssl_options is None:
                app.log.warning(f"{warning} and not using encryption. This is not recommended.")
            if not self.auth_enabled:
                app.log.warning(
                    f"{warning} and not using authentication. "
                    "This is highly insecure and not recommended."
                )
        elif not self.auth_enabled:
            app.log.warning(
                "All authentication is disabled."
                "  Anyone who can connect to this server will be able to run code."
            )

    def process_login_form(self, handler: web.RequestHandler) -> User | None:
        """Process login form data

        Return authenticated User if successful, None if not.
        """
        typed_password = handler.get_argument("password", default="")
        user = None
        if not self.auth_enabled:
            self.log.warning("Accepting anonymous login because auth fully disabled!")
            return self.generate_anonymous_user(handler)

        if self.token and self.token == typed_password:
            return t.cast(User, self.user_for_token(typed_password))  # type:ignore[attr-defined]

        return user

    @property
    def auth_enabled(self):
        """Is authentication enabled?

        Should always be True, but may be False in rare, insecure cases
        where requests with no auth are allowed.

        Previously: LoginHandler.get_login_available
        """
        return True

    @property
    def login_available(self):
        """Whether a LoginHandler is needed - and therefore whether the login page should be displayed."""
        return self.auth_enabled

    @property
    def logout_available(self):
        """Whether a LogoutHandler is needed."""
        return True


class PasswordIdentityProvider(IdentityProvider):
    """A password identity provider."""

    hashed_password = Unicode(
        "",
        config=True,
        help=_i18n(
            """
            Hashed password to use for web authentication.

            To generate, type in a python/IPython shell:

                from jupyter_server.auth import passwd; passwd()

            The string should be of the form type:salt:hashed-password.
            """
        ),
    )

    password_required = Bool(
        False,
        config=True,
        help=_i18n(
            """
            Forces users to use a password for the Jupyter server.
            This is useful in a multi user environment, for instance when
            everybody in the LAN can access each other's machine through ssh.

            In such a case, serving on localhost is not secure since
            any user can connect to the Jupyter server via ssh.

            """
        ),
    )

    allow_password_change = Bool(
        True,
        config=True,
        help=_i18n(
            """
            Allow password to be changed at login for the Jupyter server.

            While logging in with a token, the Jupyter server UI will give the opportunity to
            the user to enter a new password at the same time that will replace
            the token login mechanism.

            This can be set to False to prevent changing password from the UI/API.
            """
        ),
    )

    @default("need_token")
    def _need_token_default(self):
        return not bool(self.hashed_password)

    @property
    def login_available(self) -> bool:
        """Whether a LoginHandler is needed - and therefore whether the login page should be displayed."""
        return self.auth_enabled

    @property
    def auth_enabled(self) -> bool:
        """Return whether any auth is enabled"""
        return bool(self.hashed_password or self.token)

    def passwd_check(self, password):
        """Check password against our stored hashed password"""
        return passwd_check(self.hashed_password, password)

    def process_login_form(self, handler: web.RequestHandler) -> User | None:
        """Process login form data

        Return authenticated User if successful, None if not.
        """
        typed_password = handler.get_argument("password", default="")
        new_password = handler.get_argument("new_password", default="")
        user = None
        if not self.auth_enabled:
            self.log.warning("Accepting anonymous login because auth fully disabled!")
            return self.generate_anonymous_user(handler)

        if self.passwd_check(typed_password) and not new_password:
            return self.generate_anonymous_user(handler)
        elif self.token and self.token == typed_password:
            user = self.generate_anonymous_user(handler)
            if new_password and self.allow_password_change:
                config_dir = handler.settings.get("config_dir", "")
                config_file = os.path.join(config_dir, "jupyter_server_config.json")
                self.hashed_password = set_password(new_password, config_file=config_file)
                self.log.info(_i18n("Wrote hashed password to {file}").format(file=config_file))

        return user

    def validate_security(
        self,
        app: t.Any,
        ssl_options: dict[str, t.Any] | None = None,
    ) -> None:
        """Handle security validation."""
        super().validate_security(app, ssl_options)
        if self.password_required and (not self.hashed_password):
            self.log.critical(
                _i18n("Jupyter servers are configured to only be run with a password.")
            )
            self.log.critical(_i18n("Hint: run the following command to set a password"))
            self.log.critical(_i18n("\t$ python -m jupyter_server.auth password"))
            sys.exit(1)


class LegacyIdentityProvider(PasswordIdentityProvider):
    """Legacy IdentityProvider for use with custom LoginHandlers

    Login configuration has moved from LoginHandler to IdentityProvider
    in Jupyter Server 2.0.
    """

    # settings must be passed for
    settings = Dict()

    @default("settings")
    def _default_settings(self):
        return {
            "token": self.token,
            "password": self.hashed_password,
        }

    @default("login_handler_class")
    def _default_login_handler_class(self):
        from .login import LegacyLoginHandler

        return LegacyLoginHandler

    @property
    def auth_enabled(self):
        return self.login_available

    def get_user(self, handler: web.RequestHandler) -> User | None:
        """Get the user."""
        user = self.login_handler_class.get_user(handler)  # type:ignore[attr-defined]
        if user is None:
            return None
        return _backward_compat_user(user)

    @property
    def login_available(self) -> bool:
        return bool(
            self.login_handler_class.get_login_available(  # type:ignore[attr-defined]
                self.settings
            )
        )

    def should_check_origin(self, handler: web.RequestHandler) -> bool:
        """Whether we should check origin."""
        return bool(self.login_handler_class.should_check_origin(handler))  # type:ignore[attr-defined]

    def is_token_authenticated(self, handler: web.RequestHandler) -> bool:
        """Whether we are token authenticated."""
        return bool(self.login_handler_class.is_token_authenticated(handler))  # type:ignore[attr-defined]

    def validate_security(
        self,
        app: t.Any,
        ssl_options: dict[str, t.Any] | None = None,
    ) -> None:
        """Validate security."""
        if self.password_required and (not self.hashed_password):
            self.log.critical(
                _i18n("Jupyter servers are configured to only be run with a password.")
            )
            self.log.critical(_i18n("Hint: run the following command to set a password"))
            self.log.critical(_i18n("\t$ python -m jupyter_server.auth password"))
            sys.exit(1)
        self.login_handler_class.validate_security(  # type:ignore[attr-defined]
            app, ssl_options
        )
