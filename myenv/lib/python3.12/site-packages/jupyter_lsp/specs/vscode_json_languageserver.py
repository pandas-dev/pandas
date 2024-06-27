from .utils import NodeModuleSpec


class VSCodeJSONLanguageServer(NodeModuleSpec):
    node_module = key = "vscode-json-languageserver-bin"
    script = ["jsonServerMain.js"]
    args = ["--stdio"]
    languages = ["json"]
    spec = dict(
        display_name=key,
        mime_types=["application/json", "application/x-json", "application/ld+json"],
        urls=dict(
            home="https://github.com/vscode-langservers/{}".format(key),
            issues="https://github.com/vscode-langservers/{}/issues".format(key),
        ),
        install=dict(
            npm="npm install --save-dev {}".format(key),
            yarn="yarn add --dev {}".format(key),
            jlpm="jlpm add --dev {}".format(key),
        ),
    )
