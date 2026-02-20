from .config import load_config_schema
from .utils import NodeModuleSpec


class SQLLanguageServer(NodeModuleSpec):
    """Supports mysql, postgres and sqlite3"""

    node_module = key = "sql-language-server"
    script = ["dist", "bin", "cli.js"]
    languages = [
        "sql",
    ]
    args = ["up", "--method", "stdio"]
    spec = dict(
        display_name=key,
        mime_types=[
            "application/sql",
            "text/sql",
            "text/x-sql",
            "text/x-mysql",
            "text/x-mariadb",
            "text/x-pgsql",
        ],
        urls=dict(
            home="https://github.com/joe-re/{}".format(key),
            issues="https://github.com/joe-re/{}/issues".format(key),
        ),
        install=dict(
            npm="npm install --save-dev {}".format(key),
            yarn="yarn add --dev {}".format(key),
            jlpm="jlpm add --dev {}".format(key),
        ),
        config_schema=load_config_schema(key),
    )
