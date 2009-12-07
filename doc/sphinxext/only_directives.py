#
# A pair of directives for inserting content that will only appear in
# either html or latex.
#

from docutils.nodes import Body, Element
from docutils.writers.html4css1 import HTMLTranslator
try:
    from sphinx.latexwriter import LaTeXTranslator
except ImportError:
    from sphinx.writers.latex import LaTeXTranslator

    import warnings
    warnings.warn("The numpydoc.only_directives module is deprecated;"
                  "please use the only:: directive available in Sphinx >= 0.6",
                  DeprecationWarning, stacklevel=2)

from docutils.parsers.rst import directives

class html_only(Body, Element):
    pass

class latex_only(Body, Element):
    pass

def run(content, node_class, state, content_offset):
    text = '\n'.join(content)
    node = node_class(text)
    state.nested_parse(content, content_offset, node)
    return [node]

try:
    from docutils.parsers.rst import Directive
except ImportError:
    from docutils.parsers.rst.directives import _directives

    def html_only_directive(name, arguments, options, content, lineno,
                            content_offset, block_text, state, state_machine):
        return run(content, html_only, state, content_offset)

    def latex_only_directive(name, arguments, options, content, lineno,
                             content_offset, block_text, state, state_machine):
        return run(content, latex_only, state, content_offset)

    for func in (html_only_directive, latex_only_directive):
        func.content = 1
        func.options = {}
        func.arguments = None

    _directives['htmlonly'] = html_only_directive
    _directives['latexonly'] = latex_only_directive
else:
    class OnlyDirective(Directive):
        has_content = True
        required_arguments = 0
        optional_arguments = 0
        final_argument_whitespace = True
        option_spec = {}

        def run(self):
            self.assert_has_content()
            return run(self.content, self.node_class,
                       self.state, self.content_offset)

    class HtmlOnlyDirective(OnlyDirective):
        node_class = html_only

    class LatexOnlyDirective(OnlyDirective):
        node_class = latex_only

    directives.register_directive('htmlonly', HtmlOnlyDirective)
    directives.register_directive('latexonly', LatexOnlyDirective)

def setup(app):
    app.add_node(html_only)
    app.add_node(latex_only)

    # Add visit/depart methods to HTML-Translator:
    def visit_perform(self, node):
        pass
    def depart_perform(self, node):
        pass
    def visit_ignore(self, node):
        node.children = []
    def depart_ignore(self, node):
        node.children = []

    HTMLTranslator.visit_html_only = visit_perform
    HTMLTranslator.depart_html_only = depart_perform
    HTMLTranslator.visit_latex_only = visit_ignore
    HTMLTranslator.depart_latex_only = depart_ignore

    LaTeXTranslator.visit_html_only = visit_ignore
    LaTeXTranslator.depart_html_only = depart_ignore
    LaTeXTranslator.visit_latex_only = visit_perform
    LaTeXTranslator.depart_latex_only = depart_perform
