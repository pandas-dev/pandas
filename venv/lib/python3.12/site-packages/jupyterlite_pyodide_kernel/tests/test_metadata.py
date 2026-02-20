from pathlib import Path

import jupyterlite_pyodide_kernel
from jupyterlite_pyodide_kernel.constants import PYODIDE_VERSION

import pytest
import re

RE_CDN_ROOT = r"https://cdn.jsdelivr.net/pyodide/v(.*)/full"
RE_CDN_JS_URL = rf"{RE_CDN_ROOT}/pyodide.js"

VERSION_PATH_PATTERNS = {
    "packages/pyodide-kernel/package.json": r"\"pyodide\": \"(.*)\"",
    "examples/jupyter-lite.json": rf"{RE_CDN_ROOT}/pyodide-lock.json",
    "packages/pyodide-kernel-extension/src/index.ts": RE_CDN_JS_URL,
    "packages/pyodide-kernel-extension/schema/kernel.v0.schema.json": RE_CDN_JS_URL,
}


def test_labextension_meta():
    paths = jupyterlite_pyodide_kernel._jupyter_labextension_paths()
    assert len(paths) == 1


@pytest.mark.parametrize("repo_path", sorted(VERSION_PATH_PATTERNS))
def test_pyodide_versions(repo_path):
    path = Path(repo_path)
    if not path.exists():
        pytest.skip("not in full repo checkout")

    text = path.read_text(encoding="utf-8")
    pattern = VERSION_PATH_PATTERNS[repo_path]
    version = re.search(pattern, text)
    assert version, f"no version of pyodide found in {repo_path}: {pattern}"
    assert version[1] == PYODIDE_VERSION, f"unexpected pyodide version in {repo_path}"
