"""Private utilities for consistent node ID handling across SCM backends."""

from __future__ import annotations

# Standard node ID length used across all SCM backends
_NODE_ID_LENGTH = 10


def _slice_node_id(node_id: str) -> str:
    """
    Slice a node ID to a consistent length.

    This ensures that all SCM backends (git, mercurial, archival)
    return the same length node IDs for consistency.

    Args:
        node_id: The full node ID/hash from the SCM

    Returns:
        The node ID sliced to the standard length
    """
    return node_id[:_NODE_ID_LENGTH]


def _format_node_for_output(node_id: str | None) -> str | None:
    """
    Format a node ID for output, applying consistent slicing.

    Args:
        node_id: The full node ID/hash from the SCM or None

    Returns:
        The node ID sliced to standard length for output, or None if input was None
    """
    if node_id is None:
        return None

    # Handle mercurial nodes with 'h' prefix
    if node_id.startswith("h"):
        # For mercurial nodes, slice the part after 'h' and reconstruct
        hg_hash = node_id[1:]  # Remove 'h' prefix
        sliced_hash = _slice_node_id(hg_hash)
        return "h" + sliced_hash

    # For git nodes (with or without 'g' prefix) and others
    return _slice_node_id(node_id)
