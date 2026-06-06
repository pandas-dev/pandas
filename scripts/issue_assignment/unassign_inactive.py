"""Free issues whose claim has lapsed.

Two triggers:

* ``schedule`` — sweep open assigned issues and unassign any that are inactive
  (no open linked PR by an assignee, no assignee comment within
  ``STALE_ASSIGNEE_DAYS``).
* ``pull_request_target`` ``closed`` — when a ``Stale``-labeled PR is closed,
  free its linked issues immediately (the ``Stale`` label is a reliable proxy
  for "contributor's court, abandoned" because review-blocked PRs are exempt
  from going stale).
"""

from __future__ import annotations

from datetime import (
    datetime,
    timezone,
)
import json
import os

from scripts.issue_assignment import (
    core,
    messages,
)
from scripts.issue_assignment.client import GitHubClient


def run_sweep(client):
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


def run_pr_closed(client, event):
    pr = event["pull_request"]
    if pr.get("merged"):
        return
    labels = [label["name"] for label in pr.get("labels", [])]
    if core.STALE_LABEL not in labels:
        return
    author = pr["user"]["login"]
    for issue in client.linked_issues_for_pr(pr["number"]):
        if author in issue["assignees"]:
            client.remove_assignees(issue["number"], [author])
            client.comment(issue["number"], messages.issue_freed_stale_pr())


def main():
    repo = os.environ["GITHUB_REPOSITORY"]
    client = GitHubClient(repo)
    if os.environ.get("GITHUB_EVENT_NAME") == "schedule":
        run_sweep(client)
        return
    with open(os.environ["GITHUB_EVENT_PATH"]) as fh:
        event = json.load(fh)
    run_pr_closed(client, event)


if __name__ == "__main__":
    main()
