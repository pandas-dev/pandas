"""Tornado handlers for the sessions web service.

Preliminary documentation at https://github.com/ipython/ipython/wiki/IPEP-16%3A-Notebook-multi-directory-dashboard-and-URL-mapping#sessions-api
"""

# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.
import asyncio
import json

try:
    from jupyter_client.jsonutil import json_default
except ImportError:
    from jupyter_client.jsonutil import date_default as json_default

from jupyter_client.kernelspec import NoSuchKernel
from jupyter_core.utils import ensure_async
from tornado import web

from jupyter_server.auth.decorator import authorized
from jupyter_server.utils import url_path_join

from ...base.handlers import APIHandler

AUTH_RESOURCE = "sessions"


class SessionsAPIHandler(APIHandler):
    """A Sessions API handler."""

    auth_resource = AUTH_RESOURCE


class SessionRootHandler(SessionsAPIHandler):
    """A Session Root API handler."""

    @web.authenticated
    @authorized
    async def get(self):
        """Get a list of running sessions."""
        sm = self.session_manager
        sessions = await ensure_async(sm.list_sessions())
        self.finish(json.dumps(sessions, default=json_default))

    @web.authenticated
    @authorized
    async def post(self):
        """Create a new session."""
        # (unless a session already exists for the named session)
        sm = self.session_manager

        model = self.get_json_body()
        if model is None:
            raise web.HTTPError(400, "No JSON data provided")

        if "notebook" in model:
            self.log.warning("Sessions API changed, see updated swagger docs")
            model["type"] = "notebook"
            if "name" in model["notebook"]:
                model["path"] = model["notebook"]["name"]
            elif "path" in model["notebook"]:
                model["path"] = model["notebook"]["path"]

        try:
            # There is a high chance here that `path` is not a path but
            # a unique session id
            path = model["path"]
        except KeyError as e:
            raise web.HTTPError(400, "Missing field in JSON data: path") from e

        try:
            mtype = model["type"]
        except KeyError as e:
            raise web.HTTPError(400, "Missing field in JSON data: type") from e

        name = model.get("name", None)
        kernel = model.get("kernel", {})
        kernel_name = kernel.get("name", None)
        kernel_id = kernel.get("id", None)

        if not kernel_id and not kernel_name:
            self.log.debug("No kernel specified, using default kernel")
            kernel_name = None

        exists = await ensure_async(sm.session_exists(path=path))
        if exists:
            s_model = await sm.get_session(path=path)
        else:
            try:
                s_model = await sm.create_session(
                    path=path,
                    kernel_name=kernel_name,
                    kernel_id=kernel_id,
                    name=name,
                    type=mtype,
                )
            except NoSuchKernel:
                msg = (
                    "The '%s' kernel is not available. Please pick another "
                    "suitable kernel instead, or install that kernel." % kernel_name
                )
                status_msg = "%s not found" % kernel_name
                self.log.warning("Kernel not found: %s" % kernel_name)
                self.set_status(501)
                self.finish(json.dumps({"message": msg, "short_message": status_msg}))
                return
            except Exception as e:
                raise web.HTTPError(500, str(e)) from e

        location = url_path_join(self.base_url, "api", "sessions", s_model["id"])
        self.set_header("Location", location)
        self.set_status(201)
        self.finish(json.dumps(s_model, default=json_default))


class SessionHandler(SessionsAPIHandler):
    """A handler for a single session."""

    @web.authenticated
    @authorized
    async def get(self, session_id):
        """Get the JSON model for a single session."""
        sm = self.session_manager
        model = await sm.get_session(session_id=session_id)
        self.finish(json.dumps(model, default=json_default))

    @web.authenticated
    @authorized
    async def patch(self, session_id):
        """Patch updates sessions:

        - path updates session to track renamed paths
        - kernel.name starts a new kernel with a given kernelspec
        """
        sm = self.session_manager
        km = self.kernel_manager
        model = self.get_json_body()
        if model is None:
            raise web.HTTPError(400, "No JSON data provided")

        # get the previous session model
        before = await sm.get_session(session_id=session_id)

        changes = {}
        if "notebook" in model and "path" in model["notebook"]:
            self.log.warning("Sessions API changed, see updated swagger docs")
            model["path"] = model["notebook"]["path"]
            model["type"] = "notebook"
        if "path" in model:
            changes["path"] = model["path"]
        if "name" in model:
            changes["name"] = model["name"]
        if "type" in model:
            changes["type"] = model["type"]
        if "kernel" in model:
            # Kernel id takes precedence over name.
            if model["kernel"].get("id") is not None:
                kernel_id = model["kernel"]["id"]
                if kernel_id not in km:
                    raise web.HTTPError(400, "No such kernel: %s" % kernel_id)
                changes["kernel_id"] = kernel_id
            elif model["kernel"].get("name") is not None:
                kernel_name = model["kernel"]["name"]
                kernel_id = await sm.start_kernel_for_session(
                    session_id,
                    kernel_name=kernel_name,
                    name=before["name"],
                    path=before["path"],
                    type=before["type"],
                )
                changes["kernel_id"] = kernel_id

        await sm.update_session(session_id, **changes)
        s_model = await sm.get_session(session_id=session_id)

        if s_model["kernel"]["id"] != before["kernel"]["id"]:
            # kernel_id changed because we got a new kernel
            # shutdown the old one
            fut = asyncio.ensure_future(ensure_async(km.shutdown_kernel(before["kernel"]["id"])))
            # If we are not using pending kernels, wait for the kernel to shut down
            if not getattr(km, "use_pending_kernels", None):
                await fut
        self.finish(json.dumps(s_model, default=json_default))

    @web.authenticated
    @authorized
    async def delete(self, session_id):
        """Delete the session with given session_id."""
        sm = self.session_manager
        try:
            await sm.delete_session(session_id)
        except KeyError as e:
            # the kernel was deleted but the session wasn't!
            raise web.HTTPError(410, "Kernel deleted before session") from e
        self.set_status(204)
        self.finish()


# -----------------------------------------------------------------------------
# URL to handler mappings
# -----------------------------------------------------------------------------

_session_id_regex = r"(?P<session_id>\w+-\w+-\w+-\w+-\w+)"

default_handlers = [
    (r"/api/sessions/%s" % _session_id_regex, SessionHandler),
    (r"/api/sessions", SessionRootHandler),
]
