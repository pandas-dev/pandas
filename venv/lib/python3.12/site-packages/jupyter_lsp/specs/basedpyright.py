from .config import load_config_schema
from .utils import ShellSpec


class BasedPyrightLanguageServer(ShellSpec):
    key = "basedpyright"
    cmd = "basedpyright-langserver"
    args = ["--stdio"]
    languages = ["python"]
    spec = dict(
        display_name=key,
        mime_types=["text/python", "text/x-ipython"],
        urls=dict(
            home="https://github.com/DetachHead/basedpyright",
            issues="https://github.com/DetachHead/basedpyright/issues",
        ),
        install=dict(
            pip="pip install basedpyright",
            conda="conda install -c conda-forge basedpyright",
        ),
        config_schema=load_config_schema(key),
        requires_documents_on_disk=False,
    )
