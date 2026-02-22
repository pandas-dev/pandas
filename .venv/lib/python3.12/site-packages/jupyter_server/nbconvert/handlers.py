"""Tornado handlers for nbconvert."""

# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.
import io
import os
import sys
import zipfile

from anyio.to_thread import run_sync
from jupyter_core.utils import ensure_async
from nbformat import from_dict
from tornado import web
from tornado.log import app_log

from jupyter_server.auth.decorator import authorized

from ..base.handlers import FilesRedirectHandler, JupyterHandler, path_regex

AUTH_RESOURCE = "nbconvert"

# datetime.strftime date format for jupyter
# inlined from ipython_genutils
if sys.platform == "win32":
    date_format = "%B %d, %Y"
else:
    date_format = "%B %-d, %Y"


def find_resource_files(output_files_dir):
    """Find the resource files in a directory."""
    files = []
    for dirpath, _, filenames in os.walk(output_files_dir):
        files.extend([os.path.join(dirpath, f) for f in filenames])
    return files


def respond_zip(handler, name, output, resources):
    """Zip up the output and resource files and respond with the zip file.

    Returns True if it has served a zip file, False if there are no resource
    files, in which case we serve the plain output file.
    """
    # Check if we have resource files we need to zip
    output_files = resources.get("outputs", None)
    if not output_files:
        return False

    # Headers
    zip_filename = os.path.splitext(name)[0] + ".zip"
    handler.set_attachment_header(zip_filename)
    handler.set_header("Content-Type", "application/zip")
    handler.set_header("Cache-Control", "no-store, no-cache, must-revalidate, max-age=0")

    # Prepare the zip file
    buffer = io.BytesIO()
    zipf = zipfile.ZipFile(buffer, mode="w", compression=zipfile.ZIP_DEFLATED)
    output_filename = os.path.splitext(name)[0] + resources["output_extension"]
    zipf.writestr(output_filename, output.encode("utf-8"))
    for filename, data in output_files.items():
        zipf.writestr(os.path.basename(filename), data)
    zipf.close()

    handler.finish(buffer.getvalue())
    return True


def get_exporter(format, **kwargs):
    """get an exporter, raising appropriate errors"""
    # if this fails, will raise 500
    try:
        from nbconvert.exporters.base import get_exporter
    except ImportError as e:
        raise web.HTTPError(500, "Could not import nbconvert: %s" % e) from e

    try:
        exporter = get_exporter(format)
    except KeyError as e:
        # should this be 400?
        raise web.HTTPError(404, "No exporter for format: %s" % format) from e

    try:
        return exporter(**kwargs)
    except Exception as e:
        app_log.exception("Could not construct Exporter: %s", exporter)
        raise web.HTTPError(500, "Could not construct Exporter: %s" % e) from e


class NbconvertFileHandler(JupyterHandler):
    """An nbconvert file handler."""

    auth_resource = AUTH_RESOURCE
    SUPPORTED_METHODS = ("GET",)

    @web.authenticated
    @authorized
    async def get(self, format, path):
        """Get a notebook file in a desired format."""
        self.check_xsrf_cookie()
        exporter = get_exporter(format, config=self.config, log=self.log)

        path = path.strip("/")
        # If the notebook relates to a real file (default contents manager),
        # give its path to nbconvert.
        if hasattr(self.contents_manager, "_get_os_path"):
            os_path = self.contents_manager._get_os_path(path)
            ext_resources_dir, basename = os.path.split(os_path)
        else:
            ext_resources_dir = None

        model = await ensure_async(self.contents_manager.get(path=path))
        name = model["name"]
        if model["type"] != "notebook":
            # not a notebook, redirect to files
            return FilesRedirectHandler.redirect_to_files(self, path)

        nb = model["content"]

        self.set_header("Last-Modified", model["last_modified"])

        # create resources dictionary
        mod_date = model["last_modified"].strftime(date_format)
        nb_title = os.path.splitext(name)[0]

        resource_dict = {
            "metadata": {"name": nb_title, "modified_date": mod_date},
            "config_dir": self.application.settings["config_dir"],
        }

        if ext_resources_dir:
            resource_dict["metadata"]["path"] = ext_resources_dir

        # Exporting can take a while, delegate to a thread so we don't block the event loop
        try:
            output, resources = await run_sync(
                lambda: exporter.from_notebook_node(nb, resources=resource_dict)
            )
        except Exception as e:
            self.log.exception("nbconvert failed: %r", e)
            raise web.HTTPError(500, "nbconvert failed: %s" % e) from e

        if respond_zip(self, name, output, resources):
            return None

        # Force download if requested
        if self.get_argument("download", "false").lower() == "true":
            filename = os.path.splitext(name)[0] + resources["output_extension"]
            self.set_attachment_header(filename)

        # MIME type
        if exporter.output_mimetype:
            self.set_header("Content-Type", "%s; charset=utf-8" % exporter.output_mimetype)

        self.set_header("Cache-Control", "no-store, no-cache, must-revalidate, max-age=0")
        self.finish(output)


class NbconvertPostHandler(JupyterHandler):
    """An nbconvert post handler."""

    SUPPORTED_METHODS = ("POST",)
    auth_resource = AUTH_RESOURCE

    @web.authenticated
    @authorized
    async def post(self, format):
        """Convert a notebook file to a desired format."""
        exporter = get_exporter(format, config=self.config)

        model = self.get_json_body()
        assert model is not None
        name = model.get("name", "notebook.ipynb")
        nbnode = from_dict(model["content"])

        try:
            output, resources = await run_sync(
                lambda: exporter.from_notebook_node(
                    nbnode,
                    resources={
                        "metadata": {"name": name[: name.rfind(".")]},
                        "config_dir": self.application.settings["config_dir"],
                    },
                )
            )
        except Exception as e:
            raise web.HTTPError(500, "nbconvert failed: %s" % e) from e

        if respond_zip(self, name, output, resources):
            return

        # MIME type
        if exporter.output_mimetype:
            self.set_header("Content-Type", "%s; charset=utf-8" % exporter.output_mimetype)

        self.finish(output)


# -----------------------------------------------------------------------------
# URL to handler mappings
# -----------------------------------------------------------------------------

_format_regex = r"(?P<format>\w+)"


default_handlers = [
    (r"/nbconvert/%s" % _format_regex, NbconvertPostHandler),
    (rf"/nbconvert/{_format_regex}{path_regex}", NbconvertFileHandler),
]
