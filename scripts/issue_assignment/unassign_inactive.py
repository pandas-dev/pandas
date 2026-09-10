"""Daily issue/PR maintenance, plus the PR stale engine.

Two commands:

* ``sweep`` (the ``schedule`` trigger) — (1) sweep open assigned issues and
  unassign any non-maintainer assignees that are inactive (no open linked PR by
  an assignee, no assignee comment or fresh assignment within
  ``STALE_ASSIGNEE_DAYS``); (2) run the PR stale engine
  (``run_pr_stale_sweep``) that replaces ``actions/stale``. The sweep first
  reconciles the ``Needs Issue Assignment`` label on every non-exempt PR —
  applied (with the gate comment) whenever the author doesn't hold the
  linked issue, whether they never did or lost it since, and dropped once
  they do — then a PR is subject to stale while it carries that label or
  while a maintainer's changes-request is outstanding; whether it is stale is
  decided purely by the **author's own** activity, the warning-to-close
  countdown by the engine's own ``Stale`` label event. PRs from
  owners/members/collaborators are never touched.
* ``pr-closed`` (the ``pull_request_target`` ``closed`` trigger) — when a
  *human* closes a ``Stale``-labeled PR, free its linked issues. (Auto-closes
  from the engine free issues inline, since a ``GITHUB_TOKEN`` close doesn't
  re-trigger this event.)
"""

from __future__ import annotations

import argparse
from datetime import (
    UTC,
    datetime,
    timedelta,
)
import json
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
    now = datetime.now(UTC)  # noqa: TID251
    window = timedelta(days=core.STALE_ASSIGNEE_DAYS)
    for number in client.list_open_assigned_issue_numbers():
        activity = client.issue_activity(number)
        assignees = activity["assignees"]
        if not assignees:
            continue
        assignee_set = set(assignees)
        has_open_pr = any(a in assignee_set for a in activity["open_pr_authors"])
        candidates = [
            activity["assigned_at"],
            core.latest_assignee_comment_at(activity["comments"], assignees),
            core.ROLLOUT_CUTOFF,
        ]
        last_activity_at = max(t for t in candidates if t is not None)
        if core.issue_is_active(
            now, has_open_pr, last_activity_at, core.STALE_ASSIGNEE_DAYS
        ):
            continue
        # Unassigning is destructive, so an absence of protecting signals in
        # the bounded windows has to be re-checked exactly before acting.
        if activity["comments_truncated"] and (
            client.latest_comment_at_since(number, assignee_set, now - window)
            is not None
        ):
            continue
        if activity["crossrefs_truncated"] and any(
            a in assignee_set for a in client.open_pr_authors_exact(number)
        ):
            continue
        to_unassign = [a for a in assignees if not client.is_exempt_collaborator(a)]
        if not to_unassign:
            continue
        client.remove_assignees(number, to_unassign)
        client.comment(number, messages.issue_unassigned_inactive(to_unassign))


def _stale_clock_anchor(
    client: GitHubClient,
    now: datetime,
    pr: core.OpenPRState,
    changes_requested_at: datetime | None,
    gate_marked_at: datetime | None,
) -> datetime | None:
    """When the PR's stale clock started.

    The latest of: the changes-requested review, the gate label application,
    the author's own activity (commit, PR comment, review reply, or
    linked-issue comment), and any reopen. The clock never starts before the
    changes-requested review or the gate label, so time spent waiting on
    review never counts against the author and a gated PR's clock runs from
    the bot's comment rather than from an older commit. The reads escalate in
    cost, and each is skipped once a cheaper one already shows the author
    active within ``PR_STALE_DAYS``.
    """
    author = pr["author"]
    authors = [author] if author else []
    window = timedelta(days=core.PR_STALE_DAYS)
    candidates = [
        changes_requested_at,
        gate_marked_at,
        pr["last_commit_at"],
        core.latest_assignee_comment_at(pr["comments"], authors),
        core.latest_review_submitted_at(pr["reviews"], set(authors)),
        core.latest_event_at(pr["reopened_events"]),
    ]
    anchor = max((t for t in candidates if t is not None), default=None)
    if anchor is not None and now - anchor < window:
        return anchor
    if pr["comments_truncated"]:
        comment_at = client.latest_comment_at_since(
            pr["number"], set(authors), now - window
        )
        if comment_at is not None and (anchor is None or comment_at > anchor):
            anchor = comment_at
            if now - anchor < window:
                return anchor
    for issue in pr["linked_issues"]:
        activity = client.issue_activity(issue["number"])
        comment_at = core.latest_assignee_comment_at(activity["comments"], authors)
        if comment_at is None and activity["comments_truncated"]:
            comment_at = client.latest_comment_at_since(
                issue["number"], set(authors), now - window
            )
        if comment_at is not None and (anchor is None or comment_at > anchor):
            anchor = comment_at
    return anchor


def run_pr_stale_sweep(client: GitHubClient) -> None:
    """Mark/clear ``Stale`` and auto-close open PRs based on author inactivity."""
    now = datetime.now(UTC)  # noqa: TID251
    for pr in client.iter_open_pull_requests_review_state():
        number = pr["number"]
        author = pr["author"]
        contributors = {author} - {None}
        exempt = core.is_exempt(pr["author_association"], pr["author_is_bot"])
        gated = core.GATE_LABEL in pr["labels"]
        gate_marked_at = pr["gate_marked_at"]
        gate_issue: int | None = None
        if not exempt:
            # The gate workflow only runs when a PR is (re)opened, so the
            # label is reconciled here daily: dropped once the author holds
            # the assignment (or the PR no longer links an issue), and applied
            # when the author doesn't — including losing it after the gate
            # last passed (``/untake``, an inactivity release, or a closing
            # keyword added after opening). Exempt authors are never touched.
            decision = core.gate_decision(
                author or "",
                pr["author_association"],
                pr["author_is_bot"],
                pr["linked_issues"],
            )
            gate_action = core.gate_action(decision, gated)
            if gate_action == "clear_label":
                client.remove_label(number, core.GATE_LABEL)
                gated = False
            elif gate_action == "flag":
                client.add_labels(number, [core.GATE_LABEL])
                client.comment(number, messages.gate_flagged(author or "", decision))
                gated = True
                # The stale clock runs from this comment, not older activity.
                gate_marked_at = now
            if gated:
                gate_issue = decision["issue"]
        changes_requested_at = core.outstanding_changes_requested_at(pr["reviews"])
        subject = core.pr_subject_to_stale(
            exempt,
            pr["is_draft"],
            gated,
            changes_requested_at,
            core.latest_rereview_request_at(pr["review_requests"], contributors),
        )
        anchor = (
            _stale_clock_anchor(
                client,
                now,
                pr,
                changes_requested_at,
                gate_marked_at if gated else None,
            )
            if subject
            else None
        )
        label_present = core.STALE_LABEL in pr["labels"]
        action = core.pr_stale_action(
            now,
            subject,
            label_present,
            pr["stale_marked_at"],
            anchor,
            core.PR_STALE_DAYS,
            core.PR_CLOSE_DAYS,
        )
        if action == "mark_stale":
            if label_present:
                # A label someone else applied carries no engine label-event;
                # cycle it so the close countdown has a fresh, own anchor.
                client.remove_label(number, core.STALE_LABEL)
            client.add_labels(number, [core.STALE_LABEL])
            client.comment(number, messages.pr_marked_stale(gate_issue))
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


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=["sweep", "pr-closed"])
    parser.add_argument("--repo", required=True, help="owner/name to operate on")
    parser.add_argument(
        "--event-path", help="path to the event payload (pr-closed only)"
    )
    args = parser.parse_args(argv)
    if args.command == "pr-closed" and args.event_path is None:
        parser.error("pr-closed requires --event-path")

    client = GitHubClient(args.repo)
    if args.command == "sweep":
        run_sweep(client)
        run_pr_stale_sweep(client)
        return
    with open(args.event_path) as fh:
        event = json.load(fh)
    run_pr_closed(client, event)


if __name__ == "__main__":
    main()
