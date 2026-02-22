from .config import load_config_schema
from .utils import NodeModuleSpec


class TypescriptLanguageServer(NodeModuleSpec):
    node_module = key = "typescript-language-server"
    script = ["lib", "cli.mjs"]
    args = ["--stdio"]
    languages = [
        "javascript",
        "jsx",
        "typescript",
        "typescript-jsx",
        "typescriptreact",
        "javascriptreact",
    ]
    spec = dict(
        display_name=key,
        mime_types=[
            "application/typescript",
            "text/typescript-jsx",
            "text/javascript",
            "text/ecmascript",
            "application/javascript",
            "application/x-javascript",
            "application/ecmascript",
            "text/jsx",
        ],
        urls=dict(
            home="https://github.com/typescript-language-server/{}".format(key),
            issues="https://github.com/typescript-language-server/{}/issues".format(
                key
            ),
        ),
        install=dict(
            npm="npm install --save-dev {}".format(key),
            yarn="yarn add --dev {}".format(key),
            jlpm="jlpm add --dev {}".format(key),
        ),
        config_schema=load_config_schema(key),
    )
