from .config import load_config_schema
from .utils import ShellSpec

TROUBLESHOOT = """\
Please ensure that RScript executable is in the PATH; \
this should happen automatically when using Linux, Mac OS or Conda, \
but will require manual configuration when using the default R installer on Windows.

For more details please consult documentation:
https://cran.r-project.org/bin/windows/base/rw-FAQ.html#Rcmd-is-not-found-in-my-PATH_0021

If Rscript is already in the PATH, you can check whether \
the language server package is properly installed with:

  Rscript -e "cat(system.file(package='languageserver'))"

which should return the path to the installed package.
"""


class RLanguageServer(ShellSpec):
    package = "languageserver"
    key = "r-languageserver"
    cmd = "Rscript"

    @property
    def args(self):
        return ["--slave", "-e", f"{self.package}::run()"]

    @property
    def is_installed_args(self):
        return ["-e", f"cat(system.file(package='{self.package}'))"]

    languages = ["r"]
    spec = dict(
        display_name=key,
        mime_types=["text/x-rsrc"],
        urls=dict(
            home="https://github.com/REditorSupport/languageserver",
            issues="https://github.com/REditorSupport/languageserver/issues",
        ),
        install=dict(
            cran=f'install.packages("{package}")',
            conda="conda install -c conda-forge r-languageserver",
        ),
        config_schema=load_config_schema(key),
        troubleshoot=TROUBLESHOOT,
    )
