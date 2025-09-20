from .utils import NodeModuleSpec


class VSCodeCSSLanguageServer(NodeModuleSpec):
    node_module = key = "vscode-css-languageserver-bin"
    script = ["cssServerMain.js"]
    args = ["--stdio"]
    languages = ["css", "less", "scss"]
    spec = dict(
        display_name=key,
        mime_types=["text/x-scss", "text/css", "text/x-less"],
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
