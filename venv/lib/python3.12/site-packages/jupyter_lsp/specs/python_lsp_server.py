from .config import load_config_schema
from .utils import PythonModuleSpec


class PythonLSPServer(PythonModuleSpec):
    python_module = key = "pylsp"
    languages = ["python"]
    spec = dict(
        display_name="python-lsp-server (pylsp)",
        mime_types=["text/python", "text/x-ipython"],
        urls=dict(
            home="https://github.com/python-lsp/python-lsp-server",
            issues="https://github.com/python-lsp/python-lsp-server/issues",
        ),
        install=dict(
            pip="pip install 'python-lsp-server[all]'",
            conda="conda install -c conda-forge python-lsp-server",
        ),
        extend=[
            dict(
                display_name="pyls-mypy",
                install=dict(
                    pip="pip install pyls-mypy", conda="conda install pyls-mypy"
                ),
            ),
            dict(
                display_name="pyls-black",
                install=dict(
                    pip="pip install pyls-black", conda="conda install pyls-black"
                ),
            ),
            dict(
                display_name="pyls-isort",
                install=dict(
                    pip="pip install pyls-isort",
                    conda="conda install pyls-isort",
                ),
            ),
            dict(
                display_name="pyls-memestra",
                install=dict(
                    pip="pip install pyls-memestra",
                    conda="conda install pyls-memestra",
                ),
            ),
        ],
        config_schema=load_config_schema(key),
        env=dict(PYTHONUNBUFFERED="1"),
    )
