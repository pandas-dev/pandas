from .config import load_config_schema
from .utils import ShellSpec

TROUBLESHOOT = """\
Please ensure that texlab executable is in the PATH; \
this should happen automatically when installing texlab from Conda, \
but may require manual configuration of PATH environment variable \
if you compiled texlab from source.

You can ensure check if texlab is in the PATH, by running:

  which texlab

which should return the path to the executable (if found).
"""


class Texlab(ShellSpec):
    cmd = key = "texlab"
    languages = ["tex", "latex"]
    spec = dict(
        display_name="texlab",
        mime_types=["text/x-latex", "text/x-tex"],
        urls=dict(
            home="https://texlab.netlify.app",
            issues="https://github.com/latex-lsp/texlab/issues",
        ),
        install=dict(conda="conda install -c conda-forge texlab chktex"),
        config_schema=load_config_schema(key),
        env=dict(RUST_BACKTRACE="1"),
        troubleshoot=TROUBLESHOOT,
    )
