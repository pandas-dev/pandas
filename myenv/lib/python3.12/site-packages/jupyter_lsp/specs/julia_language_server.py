from .config import load_config_schema
from .utils import ShellSpec


class JuliaLanguageServer(ShellSpec):
    key = "julia-language-server"
    languages = ["julia"]
    cmd = "julia"
    args = [
        "--project=.",
        "-e",
        "using LanguageServer, LanguageServer.SymbolServer; runserver()",
        ".",
    ]
    is_installed_args = [
        "-e",
        'print(if (Base.find_package("LanguageServer") === nothing) "" else "yes" end)',
    ]
    spec = dict(
        display_name="LanguageServer.jl",
        mime_types=["text/julia", "text/x-julia", "application/julia"],
        urls=dict(
            home="https://github.com/julia-vscode/LanguageServer.jl",
            issues="https://github.com/julia-vscode/LanguageServer.jl/issues",
        ),
        install=dict(julia='using Pkg; Pkg.add("LanguageServer")'),
        config_schema=load_config_schema(key),
    )
