"""Contributor-facing comment text for the Python automation jobs.

Single source of truth for the messages posted by the gate, unassign, and PR
stale jobs. The ``/take`` and ``/untake`` replies live inline in
``.github/workflows/comment-commands.yml``.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from scripts.issue_assignment.core import (
    PR_CLOSE_DAYS,
    PR_STALE_DAYS,
    STALE_ASSIGNEE_DAYS,
)

if TYPE_CHECKING:
    from scripts.issue_assignment.core import GateDecision

# Point to the dev docs so that changes are reflected immediately.
DOCS_URL = (
    "https://pandas.pydata.org/docs/dev/development/contributing.html"
    "#issue-assignment-and-the-pull-request-lifecycle"
)


def _gate_review_note(issue: int) -> str:
    return (
        f"This pull request is unlikely to be reviewed until #{issue} is assigned "
        f"to you."
    )


def gate_flagged(author: str, decision: GateDecision) -> str:
    """The gate comment for an ``invalid_assignment`` decision."""
    issue = decision["issue"]
    if decision["variant"] == "unassigned":
        return gate_unassigned(author, issue)
    return gate_assigned_other(author, issue, decision["assignee"])


def gate_unassigned(author: str, issue: int) -> str:
    return (
        f"Thanks for the pull request, @{author}! It's linked to #{issue}, but "
        f"that issue isn't assigned to you yet. To make sure two people don't "
        f"unknowingly work on the same thing, we ask contributors to claim an "
        f"issue first. You can comment `/take` on #{issue} to claim it, but be "
        f"aware that certain issues cannot be taken if they are labeled e.g. "
        f"`Needs Triage` — see the [contributing guide]({DOCS_URL}) for more "
        f"details.\n\n{_gate_review_note(issue)}"
    )


def gate_assigned_other(author: str, issue: int, assignee: str) -> str:
    return (
        f"Thanks for the pull request, @{author}! It's linked to #{issue}, which "
        f"is currently assigned to @{assignee}, who's already working on it. We "
        f"keep one contributor per issue to avoid duplicated effort. If "
        f"@{assignee} has had no activity for **{STALE_ASSIGNEE_DAYS} days**, the "
        f"issue is released automatically and you'll then be able to claim it "
        f"with `/take` on #{issue} — for now, please coordinate with them on the "
        f"issue. See the [contributing guide]({DOCS_URL}) for details.\n\n"
        f"{_gate_review_note(issue)}"
    )


def _mention_list(assignees: list[str]) -> str:
    return ", ".join(f"@{a}" for a in assignees)


def issue_unassigned_inactive(assignees: list[str]) -> str:
    mentions = _mention_list(assignees)
    return (
        f"This issue has been automatically unassigned from {mentions} because "
        f"there's been no linked pull request or activity for "
        f"**{STALE_ASSIGNEE_DAYS} days**. No worries at all if other things came "
        f"up! {mentions} — you're welcome to comment `/take` to pick it back up "
        f"anytime, and it's now open for anyone else to claim too. See the "
        f"[contributing guide]({DOCS_URL}) for how this works."
    )


def pr_marked_stale(gate_issue: int | None = None) -> str:
    """Stale warning; ``gate_issue`` is the linked issue when the PR is gated.

    A gated PR (``Needs Issue Assignment``) is stale because it can't be
    reviewed, so the fix is claiming the issue — re-requesting review does
    nothing for it.
    """
    if gate_issue is None:
        next_step = (
            "If you've already addressed the feedback, **re-request a review** "
            "(the ↻ next to the reviewer) to move it back into the review "
            "queue. "
        )
    else:
        next_step = (
            f"This pull request is labeled `Needs Issue Assignment` and won't "
            f"be reviewed until #{gate_issue} is assigned to you: comment "
            f"`/take` on #{gate_issue} to claim it, and the label clears via "
            f"the daily job. "
        )
    return (
        f"This pull request has had no activity from its author for "
        f"**{PR_STALE_DAYS} days**, so I've marked it **stale**. If you're still "
        f"on it, just push a commit, reply to a review comment, or leave a "
        f"comment — here or on the linked issue — and the label clears itself. "
        f"{next_step}"
        f"Otherwise it'll be closed in "
        f"**{PR_CLOSE_DAYS} days** to keep the queue manageable. Your branch "
        f"still remains, and you can ask a maintainer to reopen this PR to "
        f"continue. See the [contributing guide]({DOCS_URL}) for how the pull "
        f"request lifecycle works.\n\n"
        f"_Note: labels update via a once-a-day job, so `Stale` may take up to a "
        f"day to clear._"
    )


def pr_closed_stale() -> str:
    return (
        f"Closing this pull request after **{PR_CLOSE_DAYS} days** stale with no "
        f"activity from its author. Thank you for the work you put into it! "
        f"**Nothing is lost** — your branch still remains, and you can ask a "
        f"maintainer to reopen this PR whenever you're ready to continue. If the "
        f"linked issue has since been claimed by someone else, just leave a "
        f"comment there to coordinate. See the [contributing guide]({DOCS_URL})."
    )


def issue_freed_stale_pr() -> str:
    return (
        "The pull request linked to this issue was closed after going stale, so "
        "I've unassigned it — this issue is available again. Anyone interested "
        f"can comment `/take` to claim it. See the [contributing guide]"
        f"({DOCS_URL}) for details."
    )
