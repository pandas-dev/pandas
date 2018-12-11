import pytest

import docutils.nodes
import docutils.utils

import validate_documentation
DocumentChecker = validate_documentation.DocumentChecker


@pytest.fixture()
def doctree():
    return docutils.utils.new_document('test.rst')


@pytest.mark.parametrize('raw_doc,errs', [
    (['   Indented Line\n'], {}),
    (['\n'], {}),
    (['\tIndented Line\n'], {'WS01': [((0, 1), {})]}),
    (['Line   \n'], {'WS02': [((0, 1), {})]}),
    (['\t\n'], {'WS03': [((0, 1), {})]}),
    ([' \n'], {'WS03': [((0, 1), {})]}),
])
def test_line_checkers(doctree, raw_doc, errs):
    checker = DocumentChecker('', raw_doc, doctree)
    checker.validate()

    assert checker.errs == errs


def test_doctree_with_list_in_block_quote(doctree):
    bullet_list = docutils.nodes.bullet_list()
    block_quote = docutils.nodes.block_quote('BlockQuote', bullet_list)
    doctree.append(block_quote)

    checker = DocumentChecker('', [''], doctree)
    checker.validate()

    assert checker.errs == {'DT01': [(None, {})]}


def test_doctree_with_block_quote_in_list(doctree):
    bullet_list = docutils.nodes.bullet_list()
    block_quote = docutils.nodes.block_quote('BlockQuote', bullet_list)
    top_list = docutils.nodes.bullet_list('', block_quote)
    doctree.append(top_list)

    checker = DocumentChecker('', [''], doctree)
    result = checker.validate()

    assert result == {'DT01': [(None, {})]}
