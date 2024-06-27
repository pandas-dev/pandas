from .utils import NodeModuleSpec


class JavascriptTypescriptLanguageServer(NodeModuleSpec):
    node_module = key = "javascript-typescript-langserver"
    script = ["lib", "language-server-stdio.js"]
    languages = [
        "javascript",
        "jsx",
        "typescript",
        "typescript-jsx",
        "typescriptreact",
        "javascriptreact",
    ]
    spec = dict(
        display_name=key + " (deprecated)",
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
            home="https://github.com/sourcegraph/{}".format(key),
            issues="https://github.com/sourcegraph/{}/issues".format(key),
        ),
        install=dict(
            npm="npm install --save-dev {}".format(key),
            yarn="yarn add --dev {}".format(key),
            jlpm="jlpm add --dev {}".format(key),
        ),
    )
