"""Pytest Fixtures exported by Jupyter Server."""

# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.
import json
from pathlib import Path

import pytest

from jupyter_server.services.contents.filemanager import AsyncFileContentsManager
from jupyter_server.services.contents.largefilemanager import AsyncLargeFileManager

pytest_plugins = ["pytest_jupyter.jupyter_server"]

some_resource = "The very model of a modern major general"
sample_kernel_json = {
    "argv": ["cat", "{connection_file}"],
    "display_name": "Test kernel",
}


@pytest.fixture  # type:ignore[misc]
def jp_kernelspecs(jp_data_dir: Path) -> None:
    """Configures some sample kernelspecs in the Jupyter data directory."""
    spec_names = ["sample", "sample2", "bad"]
    for name in spec_names:
        sample_kernel_dir = jp_data_dir.joinpath("kernels", name)
        sample_kernel_dir.mkdir(parents=True)
        # Create kernel json file
        sample_kernel_file = sample_kernel_dir.joinpath("kernel.json")
        kernel_json = sample_kernel_json.copy()
        if name == "bad":
            kernel_json["argv"] = ["non_existent_path"]
        sample_kernel_file.write_text(json.dumps(kernel_json))
        # Create resources text
        sample_kernel_resources = sample_kernel_dir.joinpath("resource.txt")
        sample_kernel_resources.write_text(some_resource)


@pytest.fixture(params=[True, False])
def jp_contents_manager(request, tmp_path):
    """Returns an AsyncFileContentsManager instance based on the use_atomic_writing parameter value."""
    return AsyncFileContentsManager(root_dir=str(tmp_path), use_atomic_writing=request.param)


@pytest.fixture
def jp_large_contents_manager(tmp_path):
    """Returns an AsyncLargeFileManager instance."""
    return AsyncLargeFileManager(root_dir=str(tmp_path))
