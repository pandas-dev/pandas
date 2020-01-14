"""Sphinx extension for collecting the titles in the rst files and validating
if they follow the capitalization convention.  This sends a warning if
it is not followed.

Usage::

   .. contents::

This will be replaced with nothing (or hello world haha)
"""

from docutils import nodes
from docutils.parsers.rst import Directive


class HelloWorld(Directive):

    def run(self):
        paragraph_node = nodes.paragraph(text='Hello World!')
        return [paragraph_node]


def setup(app):
    app.add_directive("helloworld", HelloWorld)

    return {
        'version': '0.1',
        'parallel_read_safe': True,
        'parallel_write_safe': True,
    }
