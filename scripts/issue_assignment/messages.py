"""Contributor-facing comment text for the Python automation jobs.

Single source of truth for the messages posted by the gate and unassign jobs.
The ``/take`` and ``/untake`` replies live inline in
``.github/workflows/comment-commands.yml``; the stale/close PR messages live in
``.github/workflows/stale-pr.yml``.
"""

from __future__ import annotations

from scripts.issue_assignment.core import (
    STALE_ASSIGNEE_DAYS,
    TAKEOVER_DAYS,
)

DOCS_URL = (
    "https://pandas.pydata.org/docs/development/contributing.html"
    "#issue-assignment-and-the-pull-request-lifecycle"
)


def gate_unassigned(author, issue):
    return (
        f"Thanks for the pull request, @{author}! It's linked to #{issue}, but "
        f"that issue isn't assigned to you yet. To make sure two people don't "
        f"unknowingly work on the same thing, we ask contributors to claim an "
        f"issue first. Just comment `/take` on #{issue} to claim it, and you're "
        f"good to go. See the [contributing guide]({DOCS_URL}) for the full flow."
    )


def gate_assigned_other(author, issue, assignee):
    return (
        f"Thanks for the pull request, @{author}! It's linked to #{issue}, which "
        f"is currently assigned to @{assignee}, who's already working on it. We "
        f"keep one contributor per issue to avoid duplicated effort. If "
        f"@{assignee} has been inactive for **{TAKEOVER_DAYS} days or more**, "
        f"you're welcome to comment `/take` on #{issue} to take it over — "
        f"otherwise, please coordinate with them on the issue. See the "
        f"[contributing guide]({DOCS_URL}) for details."
    )


def gate_close_addendum(issue):
    return (
        "I've closed this PR for now to keep the queue tidy — **none of your "
        "work is lost.** To pick it back up: **1)** comment `/take` on "
        f"#{issue} to claim it, then **2)** reopen this PR with the button "
        "below. Thanks!"
    )


def _mention_list(assignees):
    return ", ".join(f"@{a}" for a in assignees)


def issue_unassigned_inactive(assignees):
    mentions = _mention_list(assignees)
    return (
        f"This issue has been automatically unassigned from {mentions} because "
        f"there's been no linked pull request or activity for "
        f"**{STALE_ASSIGNEE_DAYS} days**. No worries at all if other things came "
        f"up! {mentions} — you're welcome to comment `/take` to pick it back up "
        f"anytime, and it's now open for anyone else to claim too. See the "
        f"[contributing guide]({DOCS_URL}) for how this works."
    )


def issue_freed_stale_pr():
    return (
        "The pull request linked to this issue was closed after going stale, so "
        "I've unassigned it — this issue is available again. Anyone interested "
        f"can comment `/take` to claim it. See the [contributing guide]"
        f"({DOCS_URL}) for details."
    )
