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


class ReviewRequest(TypedDict):
    actor: str | None
    requested_at: datetime | None


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
    author: str | None
    reviews: list[Review]
    review_requests: list[ReviewRequest]
    labels: list[str]


class GateDecision(TypedDict, total=False):
    # always present:
    outcome: str  # "not_in_scope" | "valid_assignment" | "invalid_assignment"
    reason: str  # "exempt" | "no_linked_issue" (not_in_scope only)
    variant: str  # "unassigned" | "assigned_other" (invalid_assignment only)
    issue: int  # invalid_assignment only
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


def latest_rereview_request_at(
    review_requests: list[ReviewRequest], contributors: set[str]
) -> datetime | None:
    """Newest review (re-)request initiated by one of ``contributors``, or ``None``.

    Scoped to the contributor's own action: a maintainer wrangling reviewers
    doesn't count as the contributor signalling "ready for another look."
    """
    times = [
        requested
        for r in review_requests
        if r.get("actor") in contributors
        and (requested := r["requested_at"]) is not None
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
    changes_requested_at: datetime | None, rereview_requested_at: datetime | None
) -> bool:
    """True when changes were requested and no re-review has been asked for since.

    A mere push doesn't move the ball back — WIP commits aren't a claim of
    readiness. The contributor has to explicitly re-request review (newer than
    the changes-requested review) to put the PR back in the maintainers' court.
    """
    if changes_requested_at is None:
        return False
    if rereview_requested_at is None:
        return True
    return changes_requested_at > rereview_requested_at


def should_label_awaiting_review(
    is_open: bool,
    is_draft: bool,
    changes_requested_at: datetime | None,
    rereview_requested_at: datetime | None,
) -> bool:
    """An open, non-draft PR with the ball in the maintainers' court."""
    if not is_open or is_draft:
        return False
    return not awaiting_contributor(changes_requested_at, rereview_requested_at)


def gate_decision(
    author: str,
    author_association: str | None,
    author_is_bot: bool,
    linked_issues: list[LinkedIssue],
) -> GateDecision:
    """Decide what the assignment gate should do for a pull request.

    Returns a dict whose ``outcome`` is one of ``not_in_scope`` (gate doesn't
    apply), ``valid_assignment`` (author is an assignee — clear the flag), or
    ``invalid_assignment`` (author isn't validly assigned — flag it). An
    ``invalid_assignment`` also carries a ``variant`` (``unassigned`` /
    ``assigned_other``), the ``issue`` number to reference, and (for
    ``assigned_other``) the ``assignee`` holding it.
    """
    if is_exempt(author_association, author_is_bot):
        return {"outcome": "not_in_scope", "reason": "exempt"}
    if not linked_issues:
        return {"outcome": "not_in_scope", "reason": "no_linked_issue"}
    if any(author in issue["assignees"] for issue in linked_issues):
        return {"outcome": "valid_assignment"}
    unassigned = [issue for issue in linked_issues if not issue["assignees"]]
    if unassigned:
        return {
            "outcome": "invalid_assignment",
            "variant": "unassigned",
            "issue": unassigned[0]["number"],
        }
    issue = linked_issues[0]
    return {
        "outcome": "invalid_assignment",
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
