import os
import shutil
import sys
from pathlib import Path
from subprocess import check_output
from typing import List, Text, Union

from ..schema import SPEC_VERSION
from ..types import (
    KeyedLanguageServerSpecs,
    LanguageServerManagerAPI,
    LanguageServerSpec,
    SpecBase,
    Token,
)

# helper scripts for known tricky language servers
HELPERS = Path(__file__).parent / "helpers"

# when building docs, let all specs go through
BUILDING_DOCS = os.environ.get("JUPYTER_LSP_BUILDING_DOCS") is not None


class ShellSpec(SpecBase):  # pragma: no cover
    """Helper for a language server spec for executables on $PATH in the
    notebook server environment.
    """

    cmd = ""

    # [optional] arguments passed to `cmd` which upon execution should print
    # out a non-empty string if the the required language server package
    # is installed, or nothing if it is missing and user action is required.
    is_installed_args: List[Token] = []

    def is_installed(self, mgr: LanguageServerManagerAPI) -> bool:
        cmd = self.solve()

        if not cmd:
            return False

        if not self.is_installed_args:
            return bool(cmd)
        else:
            check_result = check_output([cmd, *self.is_installed_args]).decode(
                encoding="utf-8"
            )
            return check_result != ""

    def solve(self) -> Union[str, None]:
        for ext in ["", ".cmd", ".bat", ".exe"]:
            cmd = shutil.which(self.cmd + ext)
            if cmd:
                break
        return cmd

    def __call__(self, mgr: LanguageServerManagerAPI) -> KeyedLanguageServerSpecs:
        cmd = self.solve()

        spec = dict(self.spec)

        if not cmd:
            troubleshooting = [f"{self.cmd} not found."]
            if "troubleshoot" in spec:
                troubleshooting.append(spec["troubleshoot"])
            spec["troubleshoot"] = "\n\n".join(troubleshooting)

        if not cmd and BUILDING_DOCS:  # pragma: no cover
            cmd = self.cmd

        return {
            self.key: {
                "argv": [cmd, *self.args] if cmd else [self.cmd, *self.args],
                "languages": self.languages,
                "version": SPEC_VERSION,
                **spec,
            }
        }


class PythonModuleSpec(SpecBase):
    """Helper for a python-based language server spec in the notebook server
    environment
    """

    python_module = ""

    def is_installed(self, mgr: LanguageServerManagerAPI) -> bool:
        spec = self.solve()

        if not spec:
            return False

        if not spec.origin:  # pragma: no cover
            return False

        return True

    def solve(self):
        return __import__("importlib").util.find_spec(self.python_module)

    def __call__(self, mgr: LanguageServerManagerAPI) -> KeyedLanguageServerSpecs:
        is_installed = self.is_installed(mgr)

        return {
            self.key: {
                "argv": (
                    [sys.executable, "-m", self.python_module, *self.args]
                    if is_installed
                    else []
                ),
                "languages": self.languages,
                "version": SPEC_VERSION,
                **self.spec,
            }
        }


class NodeModuleSpec(SpecBase):
    """Helper for a nodejs-based language server spec in one of several
    node_modules
    """

    node_module = ""
    script: List[Text] = []

    def is_installed(self, mgr: LanguageServerManagerAPI) -> bool:
        node_module = self.solve(mgr)
        return bool(node_module)

    def solve(self, mgr: LanguageServerManagerAPI):
        return mgr.find_node_module(self.node_module, *self.script)

    def __call__(self, mgr: LanguageServerManagerAPI) -> KeyedLanguageServerSpecs:
        node_module = self.solve(mgr)

        spec = dict(self.spec)

        troubleshooting = ["Node.js is required to install this server."]
        if "troubleshoot" in spec:  # pragma: no cover
            troubleshooting.append(spec["troubleshoot"])
        spec["troubleshoot"] = "\n\n".join(troubleshooting)

        is_installed = self.is_installed(mgr)

        return {
            self.key: {
                "argv": ([mgr.nodejs, node_module, *self.args] if is_installed else []),
                "languages": self.languages,
                "version": SPEC_VERSION,
                **spec,
            }
        }


# these are not desirable to publish to the frontend
# and will be replaced with the simplest schema-compliant values
SKIP_JSON_SPEC = {"argv": [""], "debug_argv": [""], "env": {}}


def censored_spec(spec: LanguageServerSpec) -> LanguageServerSpec:
    return {k: SKIP_JSON_SPEC.get(k, v) for k, v in spec.items()}
