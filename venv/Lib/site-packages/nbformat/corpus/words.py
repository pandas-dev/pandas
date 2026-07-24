"""Generate a corpus id."""

from __future__ import annotations

import uuid


def generate_corpus_id():
    """Generate a corpus id."""
    return uuid.uuid4().hex[:8]
