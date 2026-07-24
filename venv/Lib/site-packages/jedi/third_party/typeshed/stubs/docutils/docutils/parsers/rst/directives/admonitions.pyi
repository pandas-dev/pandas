from typing import Final

from docutils import nodes
from docutils.parsers.rst import Directive

__docformat__: Final = "reStructuredText"

class BaseAdmonition(Directive):
    node_class: type[nodes.Admonition]  # Subclasses must set this to the appropriate admonition node class.

class Admonition(BaseAdmonition):
    node_class: type[nodes.admonition]

class Attention(BaseAdmonition):
    node_class: type[nodes.attention]

class Caution(BaseAdmonition):
    node_class: type[nodes.caution]

class Danger(BaseAdmonition):
    node_class: type[nodes.danger]

class Error(BaseAdmonition):
    node_class: type[nodes.error]

class Hint(BaseAdmonition):
    node_class: type[nodes.hint]

class Important(BaseAdmonition):
    node_class: type[nodes.important]

class Note(BaseAdmonition):
    node_class: type[nodes.note]

class Tip(BaseAdmonition):
    node_class: type[nodes.tip]

class Warning(BaseAdmonition):
    node_class: type[nodes.warning]
