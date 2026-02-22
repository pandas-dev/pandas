"""Kernelspecs API Handlers."""

import mimetypes

from jupyter_core.utils import ensure_async
from tornado import web

from jupyter_server.auth.decorator import authorized

from ..base.handlers import JupyterHandler
from ..services.kernelspecs.handlers import kernel_name_regex

AUTH_RESOURCE = "kernelspecs"


class KernelSpecResourceHandler(web.StaticFileHandler, JupyterHandler):
    """A Kernelspec resource handler."""

    SUPPORTED_METHODS = ("GET", "HEAD")
    auth_resource = AUTH_RESOURCE

    def initialize(self):
        """Initialize a kernelspec resource handler."""
        web.StaticFileHandler.initialize(self, path="")

    @web.authenticated
    @authorized
    async def get(self, kernel_name, path, include_body=True):
        """Get a kernelspec resource."""
        ksm = self.kernel_spec_manager
        if path.lower().endswith(".png"):
            self.set_header("Cache-Control", f"max-age={60*60*24*30}")
        ksm = self.kernel_spec_manager
        if hasattr(ksm, "get_kernel_spec_resource"):
            # If the kernel spec manager defines a method to get kernelspec resources,
            # then use that instead of trying to read from disk.
            kernel_spec_res = await ksm.get_kernel_spec_resource(kernel_name, path)
            if kernel_spec_res is not None:
                # We have to explicitly specify the `absolute_path` attribute so that
                # the underlying StaticFileHandler methods can calculate an etag.
                self.absolute_path = path
                mimetype: str = mimetypes.guess_type(path)[0] or "text/plain"
                self.set_header("Content-Type", mimetype)
                self.finish(kernel_spec_res)
                return None
            else:
                self.log.warning(
                    f"Kernelspec resource '{path}' for '{kernel_name}' not found.  Kernel spec manager may"
                    " not support resource serving. Falling back to reading from disk"
                )
        try:
            kspec = await ensure_async(ksm.get_kernel_spec(kernel_name))
            self.root = kspec.resource_dir
        except KeyError as e:
            raise web.HTTPError(404, "Kernel spec %s not found" % kernel_name) from e
        self.log.debug("Serving kernel resource from: %s", self.root)
        return await web.StaticFileHandler.get(self, path, include_body=include_body)

    @web.authenticated
    @authorized
    async def head(self, kernel_name, path):
        """Get the head info for a kernel resource."""
        return await ensure_async(self.get(kernel_name, path, include_body=False))


default_handlers = [
    (r"/kernelspecs/%s/(?P<path>.*)" % kernel_name_regex, KernelSpecResourceHandler),
]
