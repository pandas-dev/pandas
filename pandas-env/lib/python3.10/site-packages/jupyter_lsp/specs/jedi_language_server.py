from .utils import ShellSpec


class JediLanguageServer(ShellSpec):
    key = cmd = "jedi-language-server"
    languages = ["python"]
    spec = dict(
        display_name="jedi-language-server",
        mime_types=["text/python", "text/x-ipython"],
        urls=dict(
            home="https://github.com/pappasam/jedi-language-server",
            issues="https://github.com/pappasam/jedi-language-server/issues",
        ),
        install=dict(
            pip="pip install -U jedi-language-server",
            conda="conda install -c conda-forge jedi-language-server",
        ),
        env=dict(PYTHONUNBUFFERED="1"),
    )
