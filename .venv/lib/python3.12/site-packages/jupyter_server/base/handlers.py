"""Base Tornado handlers for the Jupyter server."""

# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.
from __future__ import annotations

import functools
import inspect
import ipaddress
import json
import mimetypes
import os
import re
import types
import warnings
from collections.abc import Awaitable, Coroutine, Sequence
from http.client import responses
from logging import Logger
from typing import TYPE_CHECKING, Any, cast
from urllib.parse import urlparse

import prometheus_client
from jinja2 import TemplateNotFound
from jupyter_core.paths import is_hidden
from tornado import web
from tornado.log import app_log
from traitlets.config import Application

import jupyter_server
from jupyter_server import CallContext
from jupyter_server._sysinfo import get_sys_info
from jupyter_server._tz import utcnow
from jupyter_server.auth.decorator import allow_unauthenticated, authorized
from jupyter_server.auth.identity import User
from jupyter_server.i18n import combine_translations
from jupyter_server.services.security import csp_report_uri
from jupyter_server.utils import (
    ensure_async,
    filefind,
    url_escape,
    url_is_absolute,
    url_path_join,
    urldecode_unix_socket_path,
)

if TYPE_CHECKING:
    from jupyter_client.kernelspec import KernelSpecManager
    from jupyter_events import EventLogger
    from jupyter_server_terminals.terminalmanager import TerminalManager
    from tornado.concurrent import Future

    from jupyter_server.auth.authorizer import Authorizer
    from jupyter_server.auth.identity import IdentityProvider
    from jupyter_server.serverapp import ServerApp
    from jupyter_server.services.config.manager import ConfigManager
    from jupyter_server.services.contents.manager import ContentsManager
    from jupyter_server.services.kernels.kernelmanager import AsyncMappingKernelManager
    from jupyter_server.services.sessions.sessionmanager import SessionManager

# -----------------------------------------------------------------------------
# Top-level handlers
# -----------------------------------------------------------------------------

_sys_info_cache = None


def json_sys_info():
    """Get sys info as json."""
    global _sys_info_cache  # noqa: PLW0603
    if _sys_info_cache is None:
        _sys_info_cache = json.dumps(get_sys_info())
    return _sys_info_cache


def log() -> Logger:
    """Get the application log."""
    if Application.initialized():
        return cast(Logger, Application.instance().log)
    else:
        return app_log


class AuthenticatedHandler(web.RequestHandler):
    """A RequestHandler with an authenticated user."""

    @property
    def base_url(self) -> str:
        return cast(str, self.settings.get("base_url", "/"))

    @property
    def content_security_policy(self) -> str:
        """The default Content-Security-Policy header

        Can be overridden by defining Content-Security-Policy in settings['headers']
        """
        if "Content-Security-Policy" in self.settings.get("headers", {}):
            # user-specified, don't override
            return cast(str, self.settings["headers"]["Content-Security-Policy"])

        return "; ".join(
            [
                "frame-ancestors 'self'",
                # Make sure the report-uri is relative to the base_url
                "report-uri "
                + self.settings.get("csp_report_uri", url_path_join(self.base_url, csp_report_uri)),
            ]
        )

    def set_default_headers(self) -> None:
        """Set the default headers."""
        headers = {}
        headers["X-Content-Type-Options"] = "nosniff"
        headers.update(self.settings.get("headers", {}))

        headers["Content-Security-Policy"] = self.content_security_policy

        # Allow for overriding headers
        for header_name, value in headers.items():
            try:
                self.set_header(header_name, value)
            except Exception as e:
                # tornado raise Exception (not a subclass)
                # if method is unsupported (websocket and Access-Control-Allow-Origin
                # for example, so just ignore)
                self.log.exception(  # type:ignore[attr-defined]
                    "Could not set default headers: %s", e
                )

    @property
    def cookie_name(self) -> str:
        warnings.warn(
            """JupyterHandler.login_handler is deprecated in 2.0,
            use JupyterHandler.identity_provider.
            """,
            DeprecationWarning,
            stacklevel=2,
        )
        return self.identity_provider.get_cookie_name(self)

    def force_clear_cookie(self, name: str, path: str = "/", domain: str | None = None) -> None:
        """Force a cookie clear."""
        warnings.warn(
            """JupyterHandler.login_handler is deprecated in 2.0,
            use JupyterHandler.identity_provider.
            """,
            DeprecationWarning,
            stacklevel=2,
        )
        self.identity_provider._force_clear_cookie(self, name, path=path, domain=domain)

    def clear_login_cookie(self) -> None:
        """Clear a login cookie."""
        warnings.warn(
            """JupyterHandler.login_handler is deprecated in 2.0,
            use JupyterHandler.identity_provider.
            """,
            DeprecationWarning,
            stacklevel=2,
        )
        self.identity_provider.clear_login_cookie(self)

    def get_current_user(self) -> str:
        """Get the current user."""
        clsname = self.__class__.__name__
        msg = (
            f"Calling `{clsname}.get_current_user()` directly is deprecated in jupyter-server 2.0."
            " Use `self.current_user` instead (works in all versions)."
        )
        if hasattr(self, "_jupyter_current_user"):
            # backward-compat: return _jupyter_current_user
            warnings.warn(
                msg,
                DeprecationWarning,
                stacklevel=2,
            )
            return cast(str, self._jupyter_current_user)
        # haven't called get_user in prepare, raise
        raise RuntimeError(msg)

    def skip_check_origin(self) -> bool:
        """Ask my login_handler if I should skip the origin_check

        For example: in the default LoginHandler, if a request is token-authenticated,
        origin checking should be skipped.
        """
        if self.request.method == "OPTIONS":
            # no origin-check on options requests, which are used to check origins!
            return True
        return not self.identity_provider.should_check_origin(self)

    @property
    def token_authenticated(self) -> bool:
        """Have I been authenticated with a token?"""
        return self.identity_provider.is_token_authenticated(self)

    @property
    def logged_in(self) -> bool:
        """Is a user currently logged in?"""
        user = self.current_user
        return bool(user and user != "anonymous")

    @property
    def login_handler(self) -> Any:
        """Return the login handler for this application, if any."""
        warnings.warn(
            """JupyterHandler.login_handler is deprecated in 2.0,
            use JupyterHandler.identity_provider.
            """,
            DeprecationWarning,
            stacklevel=2,
        )
        return self.identity_provider.login_handler_class

    @property
    def token(self) -> str | None:
        """Return the login token for this application, if any."""
        return self.identity_provider.token

    @property
    def login_available(self) -> bool:
        """May a user proceed to log in?

        This returns True if login capability is available, irrespective of
        whether the user is already logged in or not.

        """
        return cast(bool, self.identity_provider.login_available)

    @property
    def authorizer(self) -> Authorizer:
        if "authorizer" not in self.settings:
            warnings.warn(
                "The Tornado web application does not have an 'authorizer' defined "
                "in its settings. In future releases of jupyter_server, this will "
                "be a required key for all subclasses of `JupyterHandler`. For an "
                "example, see the jupyter_server source code for how to "
                "add an authorizer to the tornado settings: "
                "https://github.com/jupyter-server/jupyter_server/blob/"
                "653740cbad7ce0c8a8752ce83e4d3c2c754b13cb/jupyter_server/serverapp.py"
                "#L234-L256",
                stacklevel=2,
            )
            from jupyter_server.auth import AllowAllAuthorizer

            self.settings["authorizer"] = AllowAllAuthorizer(
                config=self.settings.get("config", None),
                identity_provider=self.identity_provider,
            )

        return cast("Authorizer", self.settings.get("authorizer"))

    @property
    def identity_provider(self) -> IdentityProvider:
        if "identity_provider" not in self.settings:
            warnings.warn(
                "The Tornado web application does not have an 'identity_provider' defined "
                "in its settings. In future releases of jupyter_server, this will "
                "be a required key for all subclasses of `JupyterHandler`. For an "
                "example, see the jupyter_server source code for how to "
                "add an identity provider to the tornado settings: "
                "https://github.com/jupyter-server/jupyter_server/blob/v2.0.0/"
                "jupyter_server/serverapp.py#L242",
                stacklevel=2,
            )
            from jupyter_server.auth import IdentityProvider

            # no identity provider set, load default
            self.settings["identity_provider"] = IdentityProvider(
                config=self.settings.get("config", None)
            )
        return cast("IdentityProvider", self.settings["identity_provider"])


class JupyterHandler(AuthenticatedHandler):
    """Jupyter-specific extensions to authenticated handling

    Mostly property shortcuts to Jupyter-specific settings.
    """

    @property
    def config(self) -> dict[str, Any] | None:
        return cast("dict[str, Any] | None", self.settings.get("config", None))

    @property
    def log(self) -> Logger:
        """use the Jupyter log by default, falling back on tornado's logger"""
        return log()

    @property
    def jinja_template_vars(self) -> dict[str, Any]:
        """User-supplied values to supply to jinja templates."""
        return cast("dict[str, Any]", self.settings.get("jinja_template_vars", {}))

    @property
    def serverapp(self) -> ServerApp | None:
        return cast("ServerApp | None", self.settings["serverapp"])

    # ---------------------------------------------------------------
    # URLs
    # ---------------------------------------------------------------

    @property
    def version_hash(self) -> str:
        """The version hash to use for cache hints for static files"""
        return cast(str, self.settings.get("version_hash", ""))

    @property
    def mathjax_url(self) -> str:
        url = cast(str, self.settings.get("mathjax_url", ""))
        if not url or url_is_absolute(url):
            return url
        return url_path_join(self.base_url, url)

    @property
    def mathjax_config(self) -> str:
        return cast(str, self.settings.get("mathjax_config", "TeX-AMS-MML_HTMLorMML-full,Safe"))

    @property
    def default_url(self) -> str:
        return cast(str, self.settings.get("default_url", ""))

    @property
    def ws_url(self) -> str:
        return cast(str, self.settings.get("websocket_url", ""))

    @property
    def contents_js_source(self) -> str:
        self.log.debug(
            "Using contents: %s",
            self.settings.get("contents_js_source", "services/contents"),
        )
        return cast(str, self.settings.get("contents_js_source", "services/contents"))

    # ---------------------------------------------------------------
    # Manager objects
    # ---------------------------------------------------------------

    @property
    def kernel_manager(self) -> AsyncMappingKernelManager:
        return cast("AsyncMappingKernelManager", self.settings["kernel_manager"])

    @property
    def contents_manager(self) -> ContentsManager:
        return cast("ContentsManager", self.settings["contents_manager"])

    @property
    def session_manager(self) -> SessionManager:
        return cast("SessionManager", self.settings["session_manager"])

    @property
    def terminal_manager(self) -> TerminalManager:
        return cast("TerminalManager", self.settings["terminal_manager"])

    @property
    def kernel_spec_manager(self) -> KernelSpecManager:
        return cast("KernelSpecManager", self.settings["kernel_spec_manager"])

    @property
    def config_manager(self) -> ConfigManager:
        return cast("ConfigManager", self.settings["config_manager"])

    @property
    def event_logger(self) -> EventLogger:
        return cast("EventLogger", self.settings["event_logger"])

    # ---------------------------------------------------------------
    # CORS
    # ---------------------------------------------------------------

    @property
    def allow_origin(self) -> str:
        """Normal Access-Control-Allow-Origin"""
        return cast(str, self.settings.get("allow_origin", ""))

    @property
    def allow_origin_pat(self) -> str | None:
        """Regular expression version of allow_origin"""
        return cast("str | None", self.settings.get("allow_origin_pat", None))

    @property
    def allow_credentials(self) -> bool:
        """Whether to set Access-Control-Allow-Credentials"""
        return cast(bool, self.settings.get("allow_credentials", False))

    def set_default_headers(self) -> None:
        """Add CORS headers, if defined"""
        super().set_default_headers()

    def set_cors_headers(self) -> None:
        """Add CORS headers, if defined

        Now that current_user is async (jupyter-server 2.0),
        must be called at the end of prepare(), instead of in set_default_headers.
        """
        if self.allow_origin:
            self.set_header("Access-Control-Allow-Origin", self.allow_origin)
        elif self.allow_origin_pat:
            origin = self.get_origin()
            if origin and re.match(self.allow_origin_pat, origin):
                self.set_header("Access-Control-Allow-Origin", origin)
        elif self.token_authenticated and "Access-Control-Allow-Origin" not in self.settings.get(
            "headers", {}
        ):
            # allow token-authenticated requests cross-origin by default.
            # only apply this exception if allow-origin has not been specified.
            self.set_header("Access-Control-Allow-Origin", self.request.headers.get("Origin", ""))

        if self.allow_credentials:
            self.set_header("Access-Control-Allow-Credentials", "true")

    def set_attachment_header(self, filename: str) -> None:
        """Set Content-Disposition: attachment header

        As a method to ensure handling of filename encoding
        """
        escaped_filename = url_escape(filename)
        self.set_header(
            "Content-Disposition",
            f"attachment; filename*=utf-8''{escaped_filename}",
        )

    def get_origin(self) -> str | None:
        # Handle WebSocket Origin naming convention differences
        # The difference between version 8 and 13 is that in 8 the
        # client sends a "Sec-Websocket-Origin" header and in 13 it's
        # simply "Origin".
        if "Origin" in self.request.headers:
            origin = self.request.headers.get("Origin")
        else:
            origin = self.request.headers.get("Sec-Websocket-Origin", None)
        return origin

    # origin_to_satisfy_tornado is present because tornado requires
    # check_origin to take an origin argument, but we don't use it
    def check_origin(self, origin_to_satisfy_tornado: str = "") -> bool:
        """Check Origin for cross-site API requests, including websockets

        Copied from WebSocket with changes:

        - allow unspecified host/origin (e.g. scripts)
        - allow token-authenticated requests
        """
        if self.allow_origin == "*" or self.skip_check_origin():
            return True

        host = self.request.headers.get("Host")
        origin = self.request.headers.get("Origin")

        # If no header is provided, let the request through.
        # Origin can be None for:
        # - same-origin (IE, Firefox)
        # - Cross-site POST form (IE, Firefox)
        # - Scripts
        # The cross-site POST (XSRF) case is handled by tornado's xsrf_token
        if origin is None or host is None:
            return True

        origin = origin.lower()
        origin_host = urlparse(origin).netloc

        # OK if origin matches host
        if origin_host == host:
            return True

        # Check CORS headers
        if self.allow_origin:
            allow = bool(self.allow_origin == origin)
        elif self.allow_origin_pat:
            allow = bool(re.match(self.allow_origin_pat, origin))
        else:
            # No CORS headers deny the request
            allow = False
        if not allow:
            self.log.warning(
                "Blocking Cross Origin API request for %s.  Origin: %s, Host: %s",
                self.request.path,
                origin,
                host,
            )
        return allow

    def check_referer(self) -> bool:
        """Check Referer for cross-site requests.
        Disables requests to certain endpoints with
        external or missing Referer.
        If set, allow_origin settings are applied to the Referer
        to whitelist specific cross-origin sites.
        Used on GET for api endpoints and /files/
        to block cross-site inclusion (XSSI).
        """
        if self.allow_origin == "*" or self.skip_check_origin():
            return True

        host = self.request.headers.get("Host")
        referer = self.request.headers.get("Referer")

        if not host:
            self.log.warning("Blocking request with no host")
            return False
        if not referer:
            self.log.warning("Blocking request with no referer")
            return False

        referer_url = urlparse(referer)
        referer_host = referer_url.netloc
        if referer_host == host:
            return True

        # apply cross-origin checks to Referer:
        origin = f"{referer_url.scheme}://{referer_url.netloc}"
        if self.allow_origin:
            allow = self.allow_origin == origin
        elif self.allow_origin_pat:
            allow = bool(re.match(self.allow_origin_pat, origin))
        else:
            # No CORS settings, deny the request
            allow = False

        if not allow:
            self.log.warning(
                "Blocking Cross Origin request for %s.  Referer: %s, Host: %s",
                self.request.path,
                origin,
                host,
            )
        return allow

    def check_xsrf_cookie(self) -> None:
        """Bypass xsrf cookie checks when token-authenticated"""
        if not hasattr(self, "_jupyter_current_user"):
            # Called too early, will be checked later
            return None
        if self.token_authenticated or self.settings.get("disable_check_xsrf", False):
            # Token-authenticated requests do not need additional XSRF-check
            # Servers without authentication are vulnerable to XSRF
            return None
        try:
            if not self.check_origin():
                raise web.HTTPError(404)
            return super().check_xsrf_cookie()
        except web.HTTPError as e:
            if self.request.method in {"GET", "HEAD"}:
                # Consider Referer a sufficient cross-origin check for GET requests
                if not self.check_referer():
                    referer = self.request.headers.get("Referer")
                    if referer:
                        msg = f"Blocking Cross Origin request from {referer}."
                    else:
                        msg = "Blocking request from unknown origin"
                    raise web.HTTPError(403, msg) from e
            else:
                raise

    def check_host(self) -> bool:
        """Check the host header if remote access disallowed.

        Returns True if the request should continue, False otherwise.
        """
        if self.settings.get("allow_remote_access", False):
            return True

        # Remove port (e.g. ':8888') from host
        match = re.match(r"^(.*?)(:\d+)?$", self.request.host)
        assert match is not None
        host = match.group(1)

        # Browsers format IPv6 addresses like [::1]; we need to remove the []
        if host.startswith("[") and host.endswith("]"):
            host = host[1:-1]

        # UNIX socket handling
        check_host = urldecode_unix_socket_path(host)
        if check_host.startswith("/") and os.path.exists(check_host):
            allow = True
        else:
            try:
                addr = ipaddress.ip_address(host)
            except ValueError:
                # Not an IP address: check against hostnames
                allow = host in self.settings.get("local_hostnames", ["localhost"])
            else:
                allow = addr.is_loopback

        if not allow:
            self.log.warning(
                (
                    "Blocking request with non-local 'Host' %s (%s). "
                    "If the server should be accessible at that name, "
                    "set ServerApp.allow_remote_access to disable the check."
                ),
                host,
                self.request.host,
            )
        return allow

    async def prepare(self, *, _redirect_to_login=True) -> Awaitable[None] | None:  # type:ignore[override]
        """Prepare a response."""
        # Set the current Jupyter Handler context variable.
        CallContext.set(CallContext.JUPYTER_HANDLER, self)

        if not self.check_host():
            self.current_user = self._jupyter_current_user = None
            raise web.HTTPError(403)

        from jupyter_server.auth import IdentityProvider

        mod_obj = inspect.getmodule(self.get_current_user)
        assert mod_obj is not None
        user: User | None = None

        if type(self.identity_provider) is IdentityProvider and mod_obj.__name__ != __name__:
            # check for overridden get_current_user + default IdentityProvider
            # deprecated way to override auth (e.g. JupyterHub < 3.0)
            # allow deprecated, overridden get_current_user
            warnings.warn(
                "Overriding JupyterHandler.get_current_user is deprecated in jupyter-server 2.0."
                " Use an IdentityProvider class.",
                DeprecationWarning,
                stacklevel=1,
            )
            user = User(self.get_current_user())
        else:
            _user = self.identity_provider.get_user(self)
            if isinstance(_user, Awaitable):
                # IdentityProvider.get_user _may_ be async
                _user = await _user
            user = _user

        # self.current_user for tornado's @web.authenticated
        # self._jupyter_current_user for backward-compat in deprecated get_current_user calls
        # and our own private checks for whether .current_user has been set
        self.current_user = self._jupyter_current_user = user
        # complete initial steps which require auth to resolve first:
        self.set_cors_headers()
        if self.request.method not in {"GET", "HEAD", "OPTIONS"}:
            self.check_xsrf_cookie()

        if not self.settings.get("allow_unauthenticated_access", False):
            if not self.request.method:
                raise HTTPError(403)
            method = getattr(self, self.request.method.lower())
            if not getattr(method, "__allow_unauthenticated", False):
                if _redirect_to_login:
                    # reuse `web.authenticated` logic, which redirects to the login
                    # page on GET and HEAD and otherwise raises 403
                    return web.authenticated(lambda _: super().prepare())(self)
                else:
                    # raise 403 if user is not known without redirecting to login page
                    user = self.current_user
                    if user is None:
                        self.log.warning(
                            f"Couldn't authenticate {self.__class__.__name__} connection"
                        )
                        raise web.HTTPError(403)

        return super().prepare()

    # ---------------------------------------------------------------
    # template rendering
    # ---------------------------------------------------------------

    def get_template(self, name):
        """Return the jinja template object for a given name"""
        return self.settings["jinja2_env"].get_template(name)

    def render_template(self, name, **ns):
        """Render a template by name."""
        ns.update(self.template_namespace)
        template = self.get_template(name)
        return template.render(**ns)

    @property
    def template_namespace(self) -> dict[str, Any]:
        return dict(
            base_url=self.base_url,
            default_url=self.default_url,
            ws_url=self.ws_url,
            logged_in=self.logged_in,
            allow_password_change=getattr(self.identity_provider, "allow_password_change", False),
            auth_enabled=self.identity_provider.auth_enabled,
            login_available=self.identity_provider.login_available,
            token_available=bool(self.token),
            static_url=self.static_url,
            sys_info=json_sys_info(),
            contents_js_source=self.contents_js_source,
            version_hash=self.version_hash,
            xsrf_form_html=self.xsrf_form_html,
            token=self.token,
            xsrf_token=self.xsrf_token.decode("utf8"),
            nbjs_translations=json.dumps(
                combine_translations(self.request.headers.get("Accept-Language", ""))
            ),
            **self.jinja_template_vars,
        )

    def get_json_body(self) -> dict[str, Any] | None:
        """Return the body of the request as JSON data."""
        if not self.request.body:
            return None
        # Do we need to call body.decode('utf-8') here?
        body = self.request.body.strip().decode("utf-8")
        try:
            model = json.loads(body)
        except Exception as e:
            self.log.debug("Bad JSON: %r", body)
            self.log.error("Couldn't parse JSON", exc_info=True)
            raise web.HTTPError(400, "Invalid JSON in body of request") from e
        return cast("dict[str, Any]", model)

    def write_error(self, status_code: int, **kwargs: Any) -> None:
        """render custom error pages"""
        exc_info = kwargs.get("exc_info")
        message = ""
        status_message = responses.get(status_code, "Unknown HTTP Error")

        if exc_info:
            exception = exc_info[1]
            # get the custom message, if defined
            try:
                message = exception.log_message % exception.args
            except Exception:
                pass

            # construct the custom reason, if defined
            reason = getattr(exception, "reason", "")
            if reason:
                status_message = reason
        else:
            exception = "(unknown)"

        # build template namespace
        ns = {
            "status_code": status_code,
            "status_message": status_message,
            "message": message,
            "exception": exception,
        }

        self.set_header("Content-Type", "text/html")
        # render the template
        try:
            html = self.render_template("%s.html" % status_code, **ns)
        except TemplateNotFound:
            html = self.render_template("error.html", **ns)

        self.write(html)


class APIHandler(JupyterHandler):
    """Base class for API handlers"""

    async def prepare(self) -> None:  # type:ignore[override]
        """Prepare an API response."""
        await super().prepare()
        if not self.check_origin():
            raise web.HTTPError(404)

    def write_error(self, status_code: int, **kwargs: Any) -> None:
        """APIHandler errors are JSON, not human pages"""
        self.set_header("Content-Type", "application/json")
        message = responses.get(status_code, "Unknown HTTP Error")
        reply: dict[str, Any] = {
            "message": message,
        }
        exc_info = kwargs.get("exc_info")
        if exc_info:
            e = exc_info[1]
            if isinstance(e, HTTPError):
                reply["message"] = e.log_message or message
                reply["reason"] = e.reason
            else:
                reply["message"] = "Unhandled error"
                reply["reason"] = None
                # backward-compatibility: traceback field is present,
                # but always empty
                reply["traceback"] = ""
        self.log.warning("wrote error: %r", reply["message"], exc_info=True)
        self.finish(json.dumps(reply))

    def get_login_url(self) -> str:
        """Get the login url."""
        # if get_login_url is invoked in an API handler,
        # that means @web.authenticated is trying to trigger a redirect.
        # instead of redirecting, raise 403 instead.
        if not self.current_user:
            raise web.HTTPError(403)
        return super().get_login_url()

    @property
    def content_security_policy(self) -> str:
        csp = "; ".join(  # noqa: FLY002
            [
                super().content_security_policy,
                "default-src 'none'",
            ]
        )
        return csp

    # set _track_activity = False on API handlers that shouldn't track activity
    _track_activity = True

    def update_api_activity(self) -> None:
        """Update last_activity of API requests"""
        # record activity of authenticated requests
        if (
            self._track_activity
            and getattr(self, "_jupyter_current_user", None)
            and self.get_argument("no_track_activity", None) is None
        ):
            self.settings["api_last_activity"] = utcnow()

    def finish(self, *args: Any, **kwargs: Any) -> Future[Any]:
        """Finish an API response."""
        self.update_api_activity()
        # Allow caller to indicate content-type...
        set_content_type = kwargs.pop("set_content_type", "application/json")
        self.set_header("Content-Type", set_content_type)
        return super().finish(*args, **kwargs)

    @allow_unauthenticated
    def options(self, *args: Any, **kwargs: Any) -> None:
        """Get the options."""
        if "Access-Control-Allow-Headers" in self.settings.get("headers", {}):
            self.set_header(
                "Access-Control-Allow-Headers",
                self.settings["headers"]["Access-Control-Allow-Headers"],
            )
        else:
            self.set_header(
                "Access-Control-Allow-Headers",
                "accept, content-type, authorization, x-xsrftoken",
            )
        self.set_header("Access-Control-Allow-Methods", "GET, PUT, POST, PATCH, DELETE, OPTIONS")

        # if authorization header is requested,
        # that means the request is token-authenticated.
        # avoid browser-side rejection of the preflight request.
        # only allow this exception if allow_origin has not been specified
        # and Jupyter server authentication is enabled.
        # If the token is not valid, the 'real' request will still be rejected.
        requested_headers = self.request.headers.get("Access-Control-Request-Headers", "").split(
            ","
        )
        if (
            requested_headers
            and any(h.strip().lower() == "authorization" for h in requested_headers)
            and (
                # FIXME: it would be even better to check specifically for token-auth,
                # but there is currently no API for this.
                self.login_available
            )
            and (
                self.allow_origin
                or self.allow_origin_pat
                or "Access-Control-Allow-Origin" in self.settings.get("headers", {})
            )
        ):
            self.set_header("Access-Control-Allow-Origin", self.request.headers.get("Origin", ""))


class Template404(JupyterHandler):
    """Render our 404 template"""

    async def prepare(self) -> None:  # type:ignore[override]
        """Prepare a 404 response."""
        await super().prepare()
        raise web.HTTPError(404)


class AuthenticatedFileHandler(JupyterHandler, web.StaticFileHandler):
    """static files should only be accessible when logged in"""

    auth_resource = "contents"

    @property
    def content_security_policy(self) -> str:
        # In case we're serving HTML/SVG, confine any Javascript to a unique
        # origin so it can't interact with the Jupyter server.
        return super().content_security_policy + "; sandbox allow-scripts"

    @web.authenticated
    @authorized
    def head(self, path: str) -> Awaitable[None]:  # type:ignore[override]
        """Get the head response for a path."""
        self.check_xsrf_cookie()
        return super().head(path)

    @web.authenticated
    @authorized
    def get(  # type:ignore[override]
        self, path: str, **kwargs: Any
    ) -> Awaitable[None]:
        """Get a file by path."""
        self.check_xsrf_cookie()
        if os.path.splitext(path)[1] == ".ipynb" or self.get_argument("download", None):
            name = path.rsplit("/", 1)[-1]
            self.set_attachment_header(name)

        return web.StaticFileHandler.get(self, path, **kwargs)

    def get_content_type(self) -> str:
        """Get the content type."""
        assert self.absolute_path is not None
        path = self.absolute_path.strip("/")
        if "/" in path:
            _, name = path.rsplit("/", 1)
        else:
            name = path
        if name.endswith(".ipynb"):
            return "application/x-ipynb+json"
        else:
            cur_mime = mimetypes.guess_type(name)[0]
            if cur_mime == "text/plain":
                return "text/plain; charset=UTF-8"
            else:
                return super().get_content_type()

    def set_headers(self) -> None:
        """Set the headers."""
        super().set_headers()
        # disable browser caching, rely on 304 replies for savings
        if "v" not in self.request.arguments:
            self.add_header("Cache-Control", "no-cache")

    def compute_etag(self) -> str | None:
        """Compute the etag."""
        return None

    def validate_absolute_path(self, root: str, absolute_path: str) -> str:
        """Validate and return the absolute path.

        Requires tornado 3.1

        Adding to tornado's own handling, forbids the serving of hidden files.
        """
        abs_path = super().validate_absolute_path(root, absolute_path)
        abs_root = os.path.abspath(root)
        assert abs_path is not None
        if not self.contents_manager.allow_hidden and is_hidden(abs_path, abs_root):
            self.log.info(
                "Refusing to serve hidden file, via 404 Error, use flag 'ContentsManager.allow_hidden' to enable"
            )
            raise web.HTTPError(404)
        return abs_path


def json_errors(method: Any) -> Any:  # pragma: no cover
    """Decorate methods with this to return GitHub style JSON errors.

    This should be used on any JSON API on any handler method that can raise HTTPErrors.

    This will grab the latest HTTPError exception using sys.exc_info
    and then:

    1. Set the HTTP status code based on the HTTPError
    2. Create and return a JSON body with a message field describing
       the error in a human readable form.
    """
    warnings.warn(
        "@json_errors is deprecated in notebook 5.2.0. Subclass APIHandler instead.",
        DeprecationWarning,
        stacklevel=2,
    )

    @functools.wraps(method)
    def wrapper(self, *args, **kwargs):
        self.write_error = types.MethodType(APIHandler.write_error, self)
        return method(self, *args, **kwargs)

    return wrapper


# -----------------------------------------------------------------------------
# File handler
# -----------------------------------------------------------------------------

# to minimize subclass changes:
HTTPError = web.HTTPError


class FileFindHandler(JupyterHandler, web.StaticFileHandler):
    """subclass of StaticFileHandler for serving files from a search path

    The setting "static_immutable_cache" can be set up to serve some static
    file as immutable (e.g. file name containing a hash). The setting is a
    list of base URL, every static file URL starting with one of those will
    be immutable.
    """

    # cache search results, don't search for files more than once
    _static_paths: dict[str, str] = {}
    root: tuple[str]  # type:ignore[assignment]

    def set_headers(self) -> None:
        """Set the headers."""
        super().set_headers()

        immutable_paths = self.settings.get("static_immutable_cache", [])

        # allow immutable cache for files
        if any(self.request.path.startswith(path) for path in immutable_paths):
            self.set_header("Cache-Control", "public, max-age=31536000, immutable")

        # disable browser caching, rely on 304 replies for savings
        elif "v" not in self.request.arguments or any(
            self.request.path.startswith(path) for path in self.no_cache_paths
        ):
            self.set_header("Cache-Control", "no-cache")

    def initialize(
        self,
        path: str | list[str],
        default_filename: str | None = None,
        no_cache_paths: list[str] | None = None,
    ) -> None:
        """Initialize the file find handler."""
        self.no_cache_paths = no_cache_paths or []

        if isinstance(path, str):
            path = [path]

        self.root = tuple(os.path.abspath(os.path.expanduser(p)) + os.sep for p in path)  # type:ignore[assignment]
        self.default_filename = default_filename

    def compute_etag(self) -> str | None:
        """Compute the etag."""
        return None

    # access is allowed as this class is used to serve static assets on login page
    # TODO: create an allow-list of files used on login page and remove this decorator
    @allow_unauthenticated
    def get(self, path: str, include_body: bool = True) -> Coroutine[Any, Any, None]:
        return super().get(path, include_body)

    # access is allowed as this class is used to serve static assets on login page
    # TODO: create an allow-list of files used on login page and remove this decorator
    @allow_unauthenticated
    def head(self, path: str) -> Awaitable[None]:
        return super().head(path)

    @classmethod
    def get_absolute_path(cls, roots: Sequence[str], path: str) -> str:
        """locate a file to serve on our static file search path"""
        with cls._lock:
            if path in cls._static_paths:
                return cls._static_paths[path]
            try:
                abspath = os.path.abspath(filefind(path, roots))
            except OSError:
                # IOError means not found
                return ""

            cls._static_paths[path] = abspath

            log().debug(f"Path {path} served from {abspath}")
            return abspath

    def validate_absolute_path(self, root: str, absolute_path: str) -> str | None:
        """check if the file should be served (raises 404, 403, etc.)"""
        if not absolute_path:
            raise web.HTTPError(404)

        for root in self.root:
            if (absolute_path + os.sep).startswith(root):
                break

        return super().validate_absolute_path(root, absolute_path)


class APIVersionHandler(APIHandler):
    """An API handler for the server version."""

    _track_activity = False

    @allow_unauthenticated
    def get(self) -> None:
        """Get the server version info."""
        # not authenticated, so give as few info as possible
        self.finish(json.dumps({"version": jupyter_server.__version__}))


class TrailingSlashHandler(web.RequestHandler):
    """Simple redirect handler that strips trailing slashes

    This should be the first, highest priority handler.
    """

    @allow_unauthenticated
    def get(self) -> None:
        """Handle trailing slashes in a get."""
        assert self.request.uri is not None
        path, *rest = self.request.uri.partition("?")
        # trim trailing *and* leading /
        # to avoid misinterpreting repeated '//'
        path = "/" + path.strip("/")
        new_uri = "".join([path, *rest])
        self.redirect(new_uri)

    post = put = get


class MainHandler(JupyterHandler):
    """Simple handler for base_url."""

    @allow_unauthenticated
    def get(self) -> None:
        """Get the main template."""
        html = self.render_template("main.html")
        self.write(html)

    post = put = get


class FilesRedirectHandler(JupyterHandler):
    """Handler for redirecting relative URLs to the /files/ handler"""

    @staticmethod
    async def redirect_to_files(self: Any, path: str) -> None:
        """make redirect logic a reusable static method

        so it can be called from other handlers.
        """
        cm = self.contents_manager
        if await ensure_async(cm.dir_exists(path)):
            # it's a *directory*, redirect to /tree
            url = url_path_join(self.base_url, "tree", url_escape(path))
        else:
            orig_path = path
            # otherwise, redirect to /files
            parts = path.split("/")

            if not await ensure_async(cm.file_exists(path=path)) and "files" in parts:
                # redirect without files/ iff it would 404
                # this preserves pre-2.0-style 'files/' links
                self.log.warning("Deprecated files/ URL: %s", orig_path)
                parts.remove("files")
                path = "/".join(parts)

            if not await ensure_async(cm.file_exists(path=path)):
                raise web.HTTPError(404)

            url = url_path_join(self.base_url, "files", url_escape(path))
        self.log.debug("Redirecting %s to %s", self.request.path, url)
        self.redirect(url)

    @allow_unauthenticated
    async def get(self, path: str = "") -> None:
        return await self.redirect_to_files(self, path)


class RedirectWithParams(web.RequestHandler):
    """Same as web.RedirectHandler, but preserves URL parameters"""

    def initialize(self, url: str, permanent: bool = True) -> None:
        """Initialize a redirect handler."""
        self._url = url
        self._permanent = permanent

    @allow_unauthenticated
    def get(self) -> None:
        """Get a redirect."""
        sep = "&" if "?" in self._url else "?"
        url = sep.join([self._url, self.request.query])
        self.redirect(url, permanent=self._permanent)


class PrometheusMetricsHandler(JupyterHandler):
    """
    Return prometheus metrics for this server
    """

    @allow_unauthenticated
    def get(self) -> None:
        """Get prometheus metrics."""
        if self.settings["authenticate_prometheus"] and not self.logged_in:
            raise web.HTTPError(403)

        self.set_header("Content-Type", prometheus_client.CONTENT_TYPE_LATEST)
        self.write(prometheus_client.generate_latest(prometheus_client.REGISTRY))


class PublicStaticFileHandler(web.StaticFileHandler):
    """Same as web.StaticFileHandler, but decorated to acknowledge that auth is not required."""

    @allow_unauthenticated
    def head(self, path: str) -> Awaitable[None]:
        return super().head(path)

    @allow_unauthenticated
    def get(self, path: str, include_body: bool = True) -> Coroutine[Any, Any, None]:
        return super().get(path, include_body)


# -----------------------------------------------------------------------------
# URL pattern fragments for reuse
# -----------------------------------------------------------------------------

# path matches any number of `/foo[/bar...]` or just `/` or ''
path_regex = r"(?P<path>(?:(?:/[^/]+)+|/?))"

# -----------------------------------------------------------------------------
# URL to handler mappings
# -----------------------------------------------------------------------------


default_handlers = [
    (r".*/", TrailingSlashHandler),
    (r"api", APIVersionHandler),
    (r"/(robots\.txt|favicon\.ico)", PublicStaticFileHandler),
    (r"/metrics", PrometheusMetricsHandler),
]
