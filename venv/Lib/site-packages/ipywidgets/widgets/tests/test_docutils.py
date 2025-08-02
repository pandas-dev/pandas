# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.

from unittest import TestCase

from ipywidgets.widgets.docutils import doc_subst


class TestDocSubst(TestCase):

    def test_substitution(self):
        snippets = {'key': '62'}

        @doc_subst(snippets)
        def f():
            """ Docstring with value {key} """

        assert "Docstring with value 62" in f.__doc__

    def test_unused_keys(self):
        snippets = {'key': '62', 'other-key': 'unused'}

        @doc_subst(snippets)
        def f():
            """ Docstring with value {key} """

        assert "Docstring with value 62" in f.__doc__
