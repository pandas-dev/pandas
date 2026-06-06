"""Pure decision logic for issue-assignment automation.

Every function here takes plain Python data and performs no I/O, so the rules
can be exercised directly in ``scripts/tests/test_issue_assignment.py``.
"""

from __future__ import annotations

from datetime import timedelta

STALE_ASSIGNEE_DAYS = 14
TAKEOVER_DAYS = 14

EXEMPT_ASSOCIATIONS = {"OWNER", "MEMBER", "COLLABORATOR"}
BLOCKING_LABELS = ("Needs Triage", "Needs Discussion")

GATE_LABEL = "Needs Issue Assignment"
AWAITING_REVIEW_LABEL = "Awaiting Review"
STALE_LABEL = "Stale"


def is_exempt(author_association, author_is_bot):
    """Maintainers, collaborators, and bots are never gated."""
    return bool(author_is_bot) or (author_association or "") in EXEMPT_ASSOCIATIONS


def latest_changes_requested_at(reviews):
    """Newest ``CHANGES_REQUESTED`` review time, or ``None``.

    ``reviews`` is a list of ``{"state": str, "submitted_at": datetime}``.
    """
    times = [
        r["submitted_at"]
        for r in reviews
        if r.get("state") == "CHANGES_REQUESTED" and r.get("submitted_at")
    ]
    return max(times) if times else None


def latest_assignee_comment_at(comments, assignees):
    """Newest comment time authored by one of ``assignees``, or ``None``.

    ``comments`` is a list of ``{"author": str, "created_at": datetime}``.
    """
    assignees = set(assignees)
    times = [
        c["created_at"]
        for c in comments
        if c.get("author") in assignees and c.get("created_at")
    ]
    return max(times) if times else None


def awaiting_contributor(changes_requested_at, last_commit_at):
    """True when changes were requested and the author hasn't pushed since."""
    if changes_requested_at is None:
        return False
    if last_commit_at is None:
        return True
    return changes_requested_at > last_commit_at


def should_label_awaiting_review(
    is_open, is_draft, changes_requested_at, last_commit_at
):
    """An open, non-draft PR with the ball in the maintainers' court."""
    if not is_open or is_draft:
        return False
    return not awaiting_contributor(changes_requested_at, last_commit_at)


def gate_decision(author, author_association, author_is_bot, linked_issues):
    """Decide what the assignment gate should do for a pull request.

    ``linked_issues`` is a list of ``{"number": int, "assignees": [login]}``.
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
    now, has_open_linked_pr_by_assignee, last_assignee_comment_at, stale_days
):
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
