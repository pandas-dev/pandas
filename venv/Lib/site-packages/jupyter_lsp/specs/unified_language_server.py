from .utils import NodeModuleSpec


class UnifiedLanguageServer(NodeModuleSpec):
    node_module = key = "unified-language-server"
    script = ["src", "server.js"]
    args = ["--parser=remark-parse", "--stdio"]
    languages = ["markdown", "ipythongfm", "gfm"]
    spec = dict(
        display_name=key,
        mime_types=["text/x-gfm", "text/x-ipythongfm", "text/x-markdown"],
        urls=dict(
            home="https://github.com/unifiedjs/{}".format(key),
            issues="https://github.com/unifiedjs/{}/issues".format(key),
        ),
        install=dict(
            npm="npm install --save-dev {}".format(key),
            yarn="yarn add --dev {}".format(key),
            jlpm="jlpm add --dev {}".format(key),
        ),
    )
