from .config import load_config_schema
from .utils import PythonModuleSpec


class PalantirPythonLanguageServer(PythonModuleSpec):
    python_module = key = "pyls"
    languages = ["python"]
    spec = dict(
        display_name="pyls",
        mime_types=["text/python", "text/x-ipython"],
        urls=dict(
            home="https://github.com/palantir/python-language-server",
            issues="https://github.com/palantir/python-language-server/issues",
        ),
        install=dict(
            pip="pip install 'python-language-server[all]'",
            conda="conda install -c conda-forge python-language-server",
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
            dict(display_name="pyls-isort", install=dict(pip="pip install pyls-isort")),
        ],
        config_schema=load_config_schema(key),
        env=dict(PYTHONUNBUFFERED="1"),
    )
