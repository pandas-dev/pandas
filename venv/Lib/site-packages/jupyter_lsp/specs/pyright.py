from .config import load_config_schema
from .utils import NodeModuleSpec


class PyrightLanguageServer(NodeModuleSpec):
    node_module = key = "pyright"
    script = ["langserver.index.js"]
    args = ["--stdio"]
    languages = ["python"]
    spec = dict(
        display_name=key,
        mime_types=["text/python", "text/x-ipython"],
        urls=dict(
            home="https://github.com/microsoft/pyright",
            issues="https://github.com/microsoft/pyright/issues",
        ),
        install=dict(
            npm="npm install --save-dev {}".format(key),
            yarn="yarn add --dev {}".format(key),
            jlpm="jlpm add --dev {}".format(key),
        ),
        config_schema=load_config_schema(key),
        requires_documents_on_disk=False,
    )
