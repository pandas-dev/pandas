"""Pure decision logic for issue-assignment automation.

Every function here takes plain Python data and performs no I/O, so the rules
can be exercised directly in ``scripts/tests/test_issue_assignment.py``. The
``TypedDict``\\ s below document the record shapes produced by :mod:`client` and
consumed here.
"""

from __future__ import annotations

from datetime import (
    datetime,
    timedelta,
)
from typing import TypedDict

STALE_ASSIGNEE_DAYS = 14

EXEMPT_ASSOCIATIONS = {"OWNER", "MEMBER", "COLLABORATOR"}
BLOCKING_LABELS = ("Needs Triage", "Needs Discussion")

GATE_LABEL = "Needs Issue Assignment"
AWAITING_REVIEW_LABEL = "Awaiting Review"
STALE_LABEL = "Stale"


class LinkedIssue(TypedDict):
    number: int
    assignees: list[str]


class Review(TypedDict):
    state: str | None
    submitted_at: datetime | None


class Comment(TypedDict):
    author: str | None
    created_at: datetime | None


class IssueActivity(TypedDict):
    assignees: list[str]
    comments: list[Comment]
    open_pr_authors: list[str]


class OpenPRState(TypedDict):
    number: int
    is_draft: bool
    reviews: list[Review]
    last_commit_at: datetime | None
    labels: list[str]


class GateDecision(TypedDict, total=False):
    action: str  # "skip" | "pass" | "flag"; always present
    reason: str  # "exempt" | "no_linked_issue" (skip only)
    variant: str  # "unassigned" | "assigned_other" (flag only)
    issue: int  # flag only
    assignee: str  # assigned_other only


def is_exempt(author_association: str | None, author_is_bot: bool) -> bool:
    """Maintainers, collaborators, and bots are never gated."""
    return bool(author_is_bot) or (author_association or "") in EXEMPT_ASSOCIATIONS


def latest_changes_requested_at(reviews: list[Review]) -> datetime | None:
    """Newest ``CHANGES_REQUESTED`` review time, or ``None``."""
    times = [
        submitted
        for r in reviews
        if r.get("state") == "CHANGES_REQUESTED"
        and (submitted := r["submitted_at"]) is not None
    ]
    return max(times) if times else None


def latest_assignee_comment_at(
    comments: list[Comment], assignees: list[str]
) -> datetime | None:
    """Newest comment time authored by one of ``assignees``, or ``None``."""
    assignee_set = set(assignees)
    times = [
        created
        for c in comments
        if c.get("author") in assignee_set and (created := c["created_at"]) is not None
    ]
    return max(times) if times else None


def awaiting_contributor(
    changes_requested_at: datetime | None, last_commit_at: datetime | None
) -> bool:
    """True when changes were requested and the author hasn't pushed since."""
    if changes_requested_at is None:
        return False
    if last_commit_at is None:
        return True
    return changes_requested_at > last_commit_at


def should_label_awaiting_review(
    is_open: bool,
    is_draft: bool,
    changes_requested_at: datetime | None,
    last_commit_at: datetime | None,
) -> bool:
    """An open, non-draft PR with the ball in the maintainers' court."""
    if not is_open or is_draft:
        return False
    return not awaiting_contributor(changes_requested_at, last_commit_at)


def gate_decision(
    author: str,
    author_association: str | None,
    author_is_bot: bool,
    linked_issues: list[LinkedIssue],
) -> GateDecision:
    """Decide what the assignment gate should do for a pull request.

    Returns a dict with an ``action`` of ``skip``, ``pass``, or ``flag``; a
    ``flag`` also carries a ``variant`` (``unassigned`` / ``assigned_other``),
    the ``issue`` number to reference, and (for ``assigned_other``) the
    ``assignee`` holding it.
    """
    if is_exempt(author_association, author_is_bot):
        return {"action": "skip", "reason": "exempt"}
    if not linked_issues:
        return {"action": "skip", "reason": "no_linked_issue"}
    if any(author in issue["assignees"] for issue in linked_issues):
        return {"action": "pass"}
    unassigned = [issue for issue in linked_issues if not issue["assignees"]]
    if unassigned:
        return {
            "action": "flag",
            "variant": "unassigned",
            "issue": unassigned[0]["number"],
        }
    issue = linked_issues[0]
    return {
        "action": "flag",
        "variant": "assigned_other",
        "issue": issue["number"],
        "assignee": issue["assignees"][0],
    }


def issue_is_active(
    now: datetime,
    has_open_linked_pr_by_assignee: bool,
    last_assignee_comment_at: datetime | None,
    stale_days: int,
) -> bool:
    """Whether a claimed issue should keep its assignment.

    An open linked PR by an assignee keeps the claim regardless of age; failing
    that, the issue is active only if an assignee commented within
    ``stale_days``.
    """
    if has_open_linked_pr_by_assignee:
        return True
    if last_assignee_comment_at is None:
        return False
    return now - last_assignee_comment_at < timedelta(days=stale_days)
