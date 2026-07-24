"""Tornado handlers for kernel specifications.

Preliminary documentation at https://github.com/ipython/ipython/wiki/IPEP-25%3A-Registry-of-installed-kernels#rest-api
"""

# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.
from __future__ import annotations

import asyncio
import glob
import json
import os
from typing import Any

pjoin = os.path.join

from jupyter_core.utils import ensure_async
from tornado import web

from jupyter_server.auth.decorator import authorized

from ...base.handlers import APIHandler
from ...utils import url_path_join, url_unescape

AUTH_RESOURCE = "kernelspecs"


def kernelspec_model(handler, name, spec_dict, resource_dir):
    """Load a KernelSpec by name and return the REST API model"""
    d = {"name": name, "spec": spec_dict, "resources": {}}

    # Add resource files if they exist
    for resource in ["kernel.js", "kernel.css"]:
        if os.path.exists(pjoin(resource_dir, resource)):
            d["resources"][resource] = url_path_join(
                handler.base_url, "kernelspecs", name, resource
            )
    for logo_file in glob.glob(pjoin(resource_dir, "logo-*")):
        fname = os.path.basename(logo_file)
        no_ext, _ = os.path.splitext(fname)
        d["resources"][no_ext] = url_path_join(handler.base_url, "kernelspecs", name, fname)
    return d


def is_kernelspec_model(spec_dict):
    """Returns True if spec_dict is already in proper form.  This will occur when using a gateway."""
    return (
        isinstance(spec_dict, dict)
        and "name" in spec_dict
        and "spec" in spec_dict
        and "resources" in spec_dict
    )


class KernelSpecsAPIHandler(APIHandler):
    """A kernel spec API handler."""

    auth_resource = AUTH_RESOURCE


class MainKernelSpecHandler(KernelSpecsAPIHandler):
    """The root kernel spec handler."""

    @web.authenticated
    @authorized
    async def get(self):
        """Get the list of kernel specs."""
        ksm = self.kernel_spec_manager
        km = self.kernel_manager
        model: dict[str, Any] = {}
        model["default"] = km.default_kernel_name
        model["kernelspecs"] = specs = {}
        kspecs = await ensure_async(ksm.get_all_specs())
        tasks = {}
        for kernel_name, kernel_info in kspecs.items():
            if is_kernelspec_model(kernel_info):
                specs[kernel_name] = kernel_info
            else:
                kernel_spec = kernel_info.get("spec")
                kernel_resource_dir = kernel_info.get("resource_dir")
                if kernel_spec is None:
                    self.log.error("Kernel spec is missing for %s", kernel_name)
                    continue

                if kernel_resource_dir is None:
                    self.log.error("Kernel resource_dir is missing for %s", kernel_name)
                    continue
                tasks[kernel_name] = asyncio.create_task(
                    asyncio.to_thread(
                        kernelspec_model,
                        self,
                        kernel_name,
                        kernel_spec,
                        kernel_resource_dir,
                    )
                )
        for kernel_name, task in tasks.items():
            try:
                specs[kernel_name] = await task
            except Exception:
                self.log.error("Failed to load kernel spec: '%s'", kernel_name, exc_info=True)
        self.set_header("Content-Type", "application/json")
        self.finish(json.dumps(model))


class KernelSpecHandler(KernelSpecsAPIHandler):
    """A handler for an individual kernel spec."""

    @web.authenticated
    @authorized
    async def get(self, kernel_name):
        """Get a kernel spec model."""
        ksm = self.kernel_spec_manager
        kernel_name = url_unescape(kernel_name)
        try:
            spec = await ensure_async(ksm.get_kernel_spec(kernel_name))
        except KeyError as e:
            raise web.HTTPError(404, "Kernel spec %s not found" % kernel_name) from e
        if is_kernelspec_model(spec):
            model = spec
        else:
            model = await asyncio.to_thread(
                kernelspec_model,
                self,
                kernel_name,
                spec.to_dict(),
                spec.resource_dir,
            )
        self.set_header("Content-Type", "application/json")
        self.finish(json.dumps(model))


# URL to handler mappings

kernel_name_regex = r"(?P<kernel_name>[\w\.\-%]+)"

default_handlers = [
    (r"/api/kernelspecs", MainKernelSpecHandler),
    (r"/api/kernelspecs/%s" % kernel_name_regex, KernelSpecHandler),
]
