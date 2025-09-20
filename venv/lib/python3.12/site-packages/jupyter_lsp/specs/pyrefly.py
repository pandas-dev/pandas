from .config import load_config_schema
from .utils import ShellSpec


class PyreflyLanguageServer(ShellSpec):
    key = cmd = "pyrefly"
    args = ["lsp"]
    languages = ["python"]
    spec = dict(
        display_name="Pyrefly",
        mime_types=["text/python", "text/x-ipython"],
        urls=dict(
            home="https://github.com/facebook/pyrefly",
            issues="https://github.com/facebook/pyrefly/issues",
        ),
        install=dict(
            pip="pip install pyrefly",
            uv="uv add pyrefly",
            conda="conda install -c conda-forge pyrefly",
        ),
        config_schema=load_config_schema(key),
        requires_documents_on_disk=False,
    )
