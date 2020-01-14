"""Sphinx extension for collecting the titles in the rst files and validating
if they follow the capitalization convention.  This sends a warning if
it is not followed.

Usage::

   .. contents::

This will be replaced with nothing (or hello world haha)
"""

from docutils import nodes
from docutils.parsers.rst import Directive
from sphinx.builders import Builder

# Use spelling builder as an example!!!!!


class HeadingCapitalizationValidator(Builder):

    def init(self):
        z = 0
        if z == 0:
            self.warning("Warning: capitalization not followed")
        paragraph_node = nodes.paragraph(text='Hello World!')
        return [paragraph_node]

    def write_doc(self, docname, doctree):
        for node in doctree.traverse(nodes.Text):
            if (node.tagname == '#title'):
                print("Bro")

    def finish(self):
        z = 1
        return


def setup(app):
    app.info("Initializing capitalization validator")
    app.

    return {
        'version': '0.1',
        'parallel_read_safe': True,
        'parallel_write_safe': True,
    }
