import logging
import types
from _typeshed import Incomplete

import click
import flask

# This should be jinja2.Environment, but it does not have stubs and forbidden for requires in METADATA.toml
jinja_env: Incomplete
script_file: str | None
click_root_cmd: str | None
OUTPUT_FOLDER: str
_flask_app: flask.Flask | None
logger: logging.Logger | None

def create_click_web_app(module: types.ModuleType, command: click.Command, root: str = "/") -> flask.Flask: ...
