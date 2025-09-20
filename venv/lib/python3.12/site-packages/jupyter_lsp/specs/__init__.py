""" default specs
"""

# flake8: noqa: F401

from .basedpyright import BasedPyrightLanguageServer
from .bash_language_server import BashLanguageServer
from .dockerfile_language_server_nodejs import DockerfileLanguageServerNodeJS
from .javascript_typescript_langserver import JavascriptTypescriptLanguageServer
from .jedi_language_server import JediLanguageServer
from .julia_language_server import JuliaLanguageServer
from .pyls import PalantirPythonLanguageServer
from .pyrefly import PyreflyLanguageServer
from .pyright import PyrightLanguageServer
from .python_lsp_server import PythonLSPServer
from .r_languageserver import RLanguageServer
from .sql_language_server import SQLLanguageServer
from .texlab import Texlab
from .typescript_language_server import TypescriptLanguageServer
from .unified_language_server import UnifiedLanguageServer
from .vscode_css_languageserver import VSCodeCSSLanguageServer
from .vscode_html_languageserver import VSCodeHTMLLanguageServer
from .vscode_json_languageserver import VSCodeJSONLanguageServer
from .yaml_language_server import YAMLLanguageServer

basedpyright = BasedPyrightLanguageServer()
bash = BashLanguageServer()
css = VSCodeCSSLanguageServer()
dockerfile = DockerfileLanguageServerNodeJS()
html = VSCodeHTMLLanguageServer()
jedi = JediLanguageServer()
json = VSCodeJSONLanguageServer()
julia = JuliaLanguageServer()
md = UnifiedLanguageServer()
py_palantir = PalantirPythonLanguageServer()
py_lsp_server = PythonLSPServer()
pyrefly = PyreflyLanguageServer()
pyright = PyrightLanguageServer()
r = RLanguageServer()
tex = Texlab()
ts_old = JavascriptTypescriptLanguageServer()
ts = TypescriptLanguageServer()
sql = SQLLanguageServer()
yaml = YAMLLanguageServer()
