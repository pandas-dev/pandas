from ..types import LanguageServerManagerAPI
from .config import load_config_schema
from .utils import NodeModuleSpec


class BashLanguageServer(NodeModuleSpec):
    node_module = key = "bash-language-server"
    script = ["out", "cli.js"]
    fallback_script = ["bin", "main.js"]
    args = ["start"]
    languages = ["bash", "sh"]
    spec = dict(
        display_name=key,
        mime_types=["text/x-sh", "application/x-sh"],
        urls=dict(
            home="https://github.com/bash-lsp/{}".format(key),
            issues="https://github.com/bash-lsp/{}/issues".format(key),
        ),
        install=dict(
            npm="npm install --save-dev {}".format(key),
            yarn="yarn add --dev {}".format(key),
            jlpm="jlpm add --dev {}".format(key),
        ),
        config_schema=load_config_schema(key),
    )

    def solve(self, mgr: LanguageServerManagerAPI):
        new_path = mgr.find_node_module(self.node_module, *self.script)
        if new_path:
            return new_path
        return mgr.find_node_module(
            self.node_module, *self.fallback_script
        )  # pragma: no cover
