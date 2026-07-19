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
    timezone,
)
from typing import (
    Required,
    TypedDict,
)

STALE_ASSIGNEE_DAYS = 7
PR_STALE_DAYS = 14
PR_CLOSE_DAYS = 7

# Update to the merge date when this lands: assignments and inactivity that
# predate the automation are never judged by rules nobody was operating under —
# every pre-existing assignment gets a full STALE_ASSIGNEE_DAYS from this date.
ROLLOUT_CUTOFF = datetime(2026, 7, 27, tzinfo=timezone.utc)

EXEMPT_ASSOCIATIONS = {"OWNER", "MEMBER", "COLLABORATOR"}
REVIEW_BLOCKING_ASSOCIATIONS = {"OWNER", "MEMBER"}
# Review states that express a reviewer's standing verdict on the PR;
# COMMENTED and PENDING reviews leave their previous verdict in place.
STANCE_STATES = {"APPROVED", "CHANGES_REQUESTED", "DISMISSED"}
BLOCKING_LABELS = ("Needs Triage", "Needs Discussion")

GATE_LABEL = "Needs Issue Assignment"
AWAITING_REVIEW_LABEL = "Awaiting Review"
STALE_LABEL = "Stale"


class LinkedIssue(TypedDict):
    number: int
    assignees: list[str]


class Review(TypedDict):
    author: str | None
    state: str | None
    submitted_at: datetime | None
    author_association: str | None


class ReviewRequest(TypedDict):
    actor: str | None
    requested_at: datetime | None


class Comment(TypedDict):
    author: str | None
    created_at: datetime | None


class IssueActivity(TypedDict):
    assignees: list[str]
    assigned_at: datetime | None
    comments: list[Comment]
    comments_truncated: bool
    open_pr_authors: list[str]
    crossrefs_truncated: bool


class OpenPRState(TypedDict):
    number: int
    is_draft: bool
    author: str | None
    author_association: str | None
    reviews: list[Review]
    review_requests: list[ReviewRequest]
    last_commit_at: datetime | None
    comments: list[Comment]
    comments_truncated: bool
    reopened_events: list[Comment]
    labels: list[str]
    stale_marked_at: datetime | None


class GateDecision(TypedDict, total=False):
    outcome: Required[str]  # "not_in_scope" | "valid_assignment" | "invalid_assignment"
    reason: str  # "exempt" | "no_linked_issue" (not_in_scope only)
    variant: str  # "unassigned" | "assigned_other" (invalid_assignment only)
    issue: int  # invalid_assignment only
    assignee: str  # assigned_other only


def is_exempt(author_association: str | None, author_is_bot: bool) -> bool:
    """Maintainers, collaborators, and bots are never gated."""
    return bool(author_is_bot) or (author_association or "") in EXEMPT_ASSOCIATIONS


def outstanding_changes_requested_at(reviews: list[Review]) -> datetime | None:
    """Newest maintainer changes-request that is still outstanding, or ``None``.

    Only an ``OWNER``/``MEMBER`` review moves the ball into the contributor's
    court; a drive-by "request changes" from an outside account doesn't count.
    A changes-request is outstanding only while it is that reviewer's *current
    stance*: their own later ``APPROVED`` (or a dismissal, which rewrites the
    review's state to ``DISMISSED``) supersedes it, while ``COMMENTED`` reviews
    leave it standing — matching GitHub's own merge-box semantics. Another
    reviewer's approval does not clear it.
    """
    stances: dict[str | None, tuple[datetime, str | None]] = {}
    for r in reviews:
        submitted = r["submitted_at"]
        if (
            r.get("state") in STANCE_STATES
            and (r.get("author_association") or "") in REVIEW_BLOCKING_ASSOCIATIONS
            and submitted is not None
        ):
            current = stances.get(r.get("author"))
            if current is None or submitted > current[0]:
                stances[r.get("author")] = (submitted, r.get("state"))
    times = [t for t, state in stances.values() if state == "CHANGES_REQUESTED"]
    return max(times) if times else None


def latest_review_submitted_at(
    reviews: list[Review], contributors: set[str]
) -> datetime | None:
    """Newest review submitted by one of ``contributors``, or ``None``.

    Any state counts: replying inline to a review thread creates an implicit
    ``COMMENTED`` review, so this is what registers an author who responds to
    feedback without pushing or leaving a top-level comment.
    """
    times = [
        submitted
        for r in reviews
        if r.get("author") in contributors
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


def latest_event_at(events: list[Comment]) -> datetime | None:
    """Newest event time regardless of actor, or ``None``."""
    times = [created for e in events if (created := e["created_at"]) is not None]
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


def _within(now: datetime, ts: datetime | None, days: int) -> bool:
    """True when ``ts`` is set and less than ``days`` before ``now``."""
    return ts is not None and now - ts < timedelta(days=days)


def issue_is_active(
    now: datetime,
    has_open_linked_pr_by_assignee: bool,
    last_activity_at: datetime | None,
    stale_days: int,
) -> bool:
    """Whether a claimed issue should keep its assignment.

    An open linked PR by an assignee keeps the claim regardless of age; failing
    that, the issue is active only if ``last_activity_at`` — the newest of the
    latest assignee comment, the assignment itself, and ``ROLLOUT_CUTOFF`` —
    falls within ``stale_days``. Anchoring on the assignment means a fresh
    assignee who hasn't commented yet still gets the full window, and the
    cutoff keeps pre-automation claims from being expired retroactively.
    """
    if has_open_linked_pr_by_assignee:
        return True
    return _within(now, last_activity_at, stale_days)


def pr_subject_to_stale(
    is_exempt_author: bool,
    is_draft: bool,
    changes_requested_at: datetime | None,
    rereview_requested_at: datetime | None,
) -> bool:
    """Whether a PR is in the contributor's court and may go stale.

    Exempt authors (owners/members/collaborators) and drafts are never subject;
    otherwise the PR is subject only while ``awaiting_contributor`` (a maintainer
    requested changes and the author hasn't re-requested review since).
    """
    if is_exempt_author or is_draft:
        return False
    return awaiting_contributor(changes_requested_at, rereview_requested_at)


def pr_stale_action(
    now: datetime,
    subject_to_stale: bool,
    label_present: bool,
    stale_marked_at: datetime | None,
    last_author_action_at: datetime | None,
    stale_days: int,
    close_days: int,
) -> str:
    """Decide the stale action for an open PR.

    Returns ``"none"``, ``"mark_stale"``, ``"clear_stale"``, or ``"close"``.
    Whether the PR is stale at all is driven purely by the **author's own**
    activity (``last_author_action_at`` — commit, PR/review/linked-issue
    comment, reopen); the countdown from warning to close is driven by
    ``stale_marked_at``, the time the *engine* applied the ``Stale`` label. A
    close therefore always follows the engine's own warning by a full
    ``close_days``, no matter how long the author was already inactive when the
    PR was first marked. A ``Stale`` label the engine itself never applied
    (e.g. added by a human) triggers a fresh ``mark_stale`` rather than a
    close, so the warning comment is never skipped.
    """
    if not subject_to_stale or _within(now, last_author_action_at, stale_days):
        return "clear_stale" if label_present else "none"
    if not label_present or stale_marked_at is None:
        return "mark_stale"
    if _within(now, stale_marked_at, close_days):
        return "none"
    return "close"


def gate_action(
    decision: GateDecision, label_present: bool, close_enabled: bool
) -> str:
    """What the gate should actually do, given the decision and current labels.

    Returns ``"none"``, ``"clear_label"``, ``"flag"``, or ``"flag_and_close"``.
    In warn-only mode an already-flagged PR is left alone — reopening without
    fixing the assignment shouldn't repost the same comment. In close mode a
    repeat still comments and closes: silently re-closing a reopened PR would
    be far worse than repeating the explanation.
    """
    if decision["outcome"] == "not_in_scope":
        return "none"
    if decision["outcome"] == "valid_assignment":
        return "clear_label" if label_present else "none"
    if close_enabled:
        return "flag_and_close"
    return "none" if label_present else "flag"
