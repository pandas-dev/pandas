"""JupyterLab Server config."""

# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.
from __future__ import annotations

import json
from itertools import chain
from pathlib import Path
from typing import Any

DEFAULT_TEMPLATE_PATH = str(Path(__file__).parent / "templates")


def get_package_url(data: dict[str, Any]) -> str:
    """Get the url from the extension data."""
    # homepage, repository  are optional
    if "homepage" in data:
        return str(data["homepage"])
    if "repository" in data and isinstance(data["repository"], dict):
        return str(data["repository"].get("url", ""))
    return ""


def get_federated_extensions(labextensions_path: list[str]) -> dict[str, Any]:
    """Get the metadata about federated extensions."""
    federated_extensions = {}
    for ext_dir in labextensions_path:
        ext_dir_path = Path(ext_dir)
        # extensions are either top-level directories, or two-deep in @org directories
        for pkg_file in chain(
            ext_dir_path.glob("[!@]*/package.json"),
            ext_dir_path.glob("@*/*/package.json"),
        ):
            with pkg_file.open(encoding="utf-8") as fid:
                pkgdata = json.load(fid)
            if pkgdata["name"] not in federated_extensions:
                data: dict[str, Any] = {
                    "name": pkgdata["name"],
                    "version": pkgdata["version"],
                    "description": pkgdata.get("description", ""),
                    "url": get_package_url(pkgdata),
                    "ext_dir": ext_dir,
                    "ext_path": str(pkg_file.parent),
                    "is_local": False,
                    "dependencies": pkgdata.get("dependencies", {}),
                    "jupyterlab": pkgdata.get("jupyterlab", {}),
                }

                # Add repository info if available
                if "repository" in pkgdata and "url" in pkgdata.get("repository", {}):
                    data["repository"] = {"url": pkgdata.get("repository").get("url")}

                install_path = pkg_file.parent / "install.json"
                if install_path.exists():
                    with install_path.open(encoding="utf-8") as fid:
                        data["install"] = json.load(fid)
                federated_extensions[data["name"]] = data
    return federated_extensions
