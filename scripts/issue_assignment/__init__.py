"""Automation for issue assignment and the pull request lifecycle.

Decision logic lives in :mod:`core` (pure functions, unit tested), the thin
GitHub API layer in :mod:`client`, and contributor-facing text in
:mod:`messages`. The per-workflow entry points (:mod:`gate`,
:mod:`label_awaiting_review`, :mod:`unassign_inactive`) wire them together.
"""

from __future__ import annotations
