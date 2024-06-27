"""Tests for nbformat corpus"""

# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.
from __future__ import annotations

from nbformat.corpus import words


def test_generate_corpus_id(recwarn):
    """Test generating a corpus id."""
    assert len(words.generate_corpus_id()) > 7
    # 1 in 4294967296 (2^32) times this will fail
    assert words.generate_corpus_id() != words.generate_corpus_id()
    assert len(recwarn) == 0
