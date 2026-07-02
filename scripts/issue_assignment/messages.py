"""Contributor-facing comment text for the Python automation jobs.

Single source of truth for the messages posted by the gate, unassign, and PR
stale jobs. The ``/take`` and ``/untake`` replies live inline in
``.github/workflows/comment-commands.yml``.
"""

from __future__ import annotations

from scripts.issue_assignment.core import (
    PR_CLOSE_DAYS,
    PR_STALE_DAYS,
    STALE_ASSIGNEE_DAYS,
)

DOCS_URL = (
    "https://pandas.pydata.org/docs/development/contributing.html"
    "#issue-assignment-and-the-pull-request-lifecycle"
)


def gate_unassigned(author: str, issue: int) -> str:
    return (
        f"Thanks for the pull request, @{author}! It's linked to #{issue}, but "
        f"that issue isn't assigned to you yet. To make sure two people don't "
        f"unknowingly work on the same thing, we ask contributors to claim an "
        f"issue first. Just comment `/take` on #{issue} to claim it, and you're "
        f"good to go. See the [contributing guide]({DOCS_URL}) for the full flow."
    )


def gate_assigned_other(author: str, issue: int, assignee: str) -> str:
    return (
        f"Thanks for the pull request, @{author}! It's linked to #{issue}, which "
        f"is currently assigned to @{assignee}, who's already working on it. We "
        f"keep one contributor per issue to avoid duplicated effort. If "
        f"@{assignee} has had no activity for **{STALE_ASSIGNEE_DAYS} days**, the "
        f"issue is released automatically and you'll then be able to claim it "
        f"with `/take` on #{issue} — for now, please coordinate with them on the "
        f"issue. See the [contributing guide]({DOCS_URL}) for details."
    )


def gate_close_addendum(issue: int) -> str:
    return (
        "I've closed this PR for now to keep the queue tidy — **none of your "
        "work is lost.** To pick it back up: **1)** comment `/take` on "
        f"#{issue} to claim it, then **2)** reopen this PR with the button "
        "below. Thanks!"
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


def pr_marked_stale() -> str:
    return (
        f"This pull request has had no activity from its author for "
        f"**{PR_STALE_DAYS} days**, so I've marked it **stale**. If you're still "
        f"on it, just push a commit or leave a comment — here or on the linked "
        f"issue — and the label clears itself. If you've already addressed the "
        f"feedback, **re-request a review** (the ↻ next to the reviewer) to move "
        f"it back into the review queue. Otherwise it'll be closed in "
        f"**{PR_CLOSE_DAYS} days** to keep the queue manageable — you can always "
        f"reopen it later to continue. See the [contributing guide]({DOCS_URL}) "
        f"for how the pull request lifecycle works.\n\n"
        f"_Note: labels update via a once-a-day job, so `Stale` may take up to a "
        f"day to clear._"
    )


def pr_closed_stale() -> str:
    return (
        f"Closing this pull request after **{PR_CLOSE_DAYS} days** stale with no "
        f"activity from its author. Thank you for the work you put into it! "
        f"**Nothing is lost** — you can reopen this PR anytime to continue. If "
        f"the linked issue has since been claimed by someone else, just leave a "
        f"comment there to coordinate. See the [contributing guide]({DOCS_URL})."
    )


def issue_freed_stale_pr() -> str:
    return (
        "The pull request linked to this issue was closed after going stale, so "
        "I've unassigned it — this issue is available again. Anyone interested "
        f"can comment `/take` to claim it. See the [contributing guide]"
        f"({DOCS_URL}) for details."
    )
