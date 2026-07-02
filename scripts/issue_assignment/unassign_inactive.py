"""Daily issue/PR maintenance, plus the PR stale engine.

Two triggers:

* ``schedule`` — (1) sweep open assigned issues and unassign any that are
  inactive (no open linked PR by an assignee, no assignee comment within
  ``STALE_ASSIGNEE_DAYS``); (2) run the PR stale engine (``run_pr_stale_sweep``)
  that replaces ``actions/stale``: mark/clear ``Stale`` and auto-close purely on
  the **author's own** activity, with owners/members/collaborators exempt.
* ``pull_request_target`` ``closed`` — when a *human* closes a ``Stale``-labeled
  PR, free its linked issues. (Auto-closes from the engine free issues inline,
  since a ``GITHUB_TOKEN`` close doesn't re-trigger this event.)
"""

from __future__ import annotations

from datetime import (
    datetime,
    timedelta,
    timezone,
)
import json
import os
from typing import Any

from scripts.issue_assignment import (
    core,
    messages,
)
from scripts.issue_assignment.client import GitHubClient


def free_linked_issues(client: GitHubClient, number: int, author: str | None) -> None:
    """Unassign ``author`` from every linked issue they hold + comment on each."""
    if author is None:
        return
    for issue in client.linked_issues_for_pr(number):
        if author in issue["assignees"]:
            client.remove_assignees(issue["number"], [author])
            client.comment(issue["number"], messages.issue_freed_stale_pr())


def run_sweep(client: GitHubClient) -> None:
    now = datetime.now(timezone.utc)
    for number in client.list_open_assigned_issue_numbers():
        activity = client.issue_activity(number)
        assignees = activity["assignees"]
        if not assignees:
            continue
        assignee_set = set(assignees)
        has_open_pr = any(a in assignee_set for a in activity["open_pr_authors"])
        last_comment_at = core.latest_assignee_comment_at(
            activity["comments"], assignees
        )
        if not core.issue_is_active(
            now, has_open_pr, last_comment_at, core.STALE_ASSIGNEE_DAYS
        ):
            client.remove_assignees(number, assignees)
            client.comment(number, messages.issue_unassigned_inactive(assignees))


def _stale_clock_anchor(
    client: GitHubClient,
    now: datetime,
    pr: core.OpenPRState,
    changes_requested_at: datetime | None,
) -> datetime | None:
    """When the PR's stale clock started: the later of ``changes_requested_at``
    (the ball entering the contributor's court) and the author's last action
    (commit, PR comment, or linked-issue comment).

    Flooring at ``changes_requested_at`` is what keeps a long review wait from
    eating into the stale window — the author gets a full ``PR_STALE_DAYS`` from
    the review, never from whenever they last happened to act before it. Linked
    issues are only read when nothing already shows the author active inside
    ``PR_STALE_DAYS`` — that read can't change the outcome, so it's skipped.
    """
    author = pr["author"]
    candidates = [
        changes_requested_at,
        pr["last_commit_at"],
        core.latest_assignee_comment_at(pr["comments"], [author] if author else []),
    ]
    anchor = max((t for t in candidates if t is not None), default=None)
    if anchor is not None and now - anchor < timedelta(days=core.PR_STALE_DAYS):
        return anchor
    for issue in client.linked_issues_for_pr(pr["number"]):
        comment_at = core.latest_assignee_comment_at(
            client.issue_activity(issue["number"])["comments"],
            [author] if author else [],
        )
        if comment_at is not None and (anchor is None or comment_at > anchor):
            anchor = comment_at
    return anchor


def run_pr_stale_sweep(client: GitHubClient) -> None:
    """Mark/clear ``Stale`` and auto-close open PRs based on author inactivity."""
    now = datetime.now(timezone.utc)
    for pr in client.iter_open_pull_requests_review_state():
        contributors = {pr["author"]} - {None}
        changes_requested_at = core.latest_changes_requested_at(pr["reviews"])
        subject = core.pr_subject_to_stale(
            core.is_exempt(pr["author_association"], False),
            pr["is_draft"],
            changes_requested_at,
            core.latest_rereview_request_at(pr["review_requests"], contributors),
        )
        anchor = (
            _stale_clock_anchor(client, now, pr, changes_requested_at)
            if subject
            else None
        )
        action = core.pr_stale_action(
            now,
            subject,
            core.STALE_LABEL in pr["labels"],
            anchor,
            core.PR_STALE_DAYS,
            core.PR_CLOSE_DAYS,
        )
        number = pr["number"]
        if action == "mark_stale":
            client.add_labels(number, [core.STALE_LABEL])
            client.comment(number, messages.pr_marked_stale())
        elif action == "clear_stale":
            client.remove_label(number, core.STALE_LABEL)
        elif action == "close":
            client.comment(number, messages.pr_closed_stale())
            free_linked_issues(client, number, pr["author"])
            client.close_pull_request(number)


def run_pr_closed(client: GitHubClient, event: dict[str, Any]) -> None:
    pr = event["pull_request"]
    if pr.get("merged"):
        return
    if core.is_exempt(pr.get("author_association"), pr["user"].get("type") == "Bot"):
        return
    labels = [label["name"] for label in pr.get("labels", [])]
    if core.STALE_LABEL not in labels:
        return
    free_linked_issues(client, pr["number"], pr["user"]["login"])


def main() -> None:
    repo = os.environ["GITHUB_REPOSITORY"]
    client = GitHubClient(repo)
    if os.environ.get("GITHUB_EVENT_NAME") == "schedule":
        run_sweep(client)
        run_pr_stale_sweep(client)
        return
    with open(os.environ["GITHUB_EVENT_PATH"]) as fh:
        event = json.load(fh)
    run_pr_closed(client, event)


if __name__ == "__main__":
    main()
