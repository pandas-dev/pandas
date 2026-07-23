"""Assignment gate: warn (and optionally close) PRs from non-assignees.

Triggered by ``pull_request_target`` on ``opened`` / ``reopened``.
"""

from __future__ import annotations

import json
import os
import time

from scripts.issue_assignment import (
    core,
    messages,
)
from scripts.issue_assignment.client import GitHubClient


def main() -> None:
    with open(os.environ["GITHUB_EVENT_PATH"]) as fh:
        event = json.load(fh)
    repo = os.environ["GITHUB_REPOSITORY"]
    close_enabled = os.environ.get("CLOSE_ENABLED", "false").lower() == "true"

    pr = event["pull_request"]
    number = pr["number"]
    author = pr["user"]["login"]
    author_is_bot = pr["user"].get("type") == "Bot"
    label_present = core.GATE_LABEL in [label["name"] for label in pr.get("labels", [])]

    client = GitHubClient(repo)
    # GitHub resolves closing keywords into linked issues asynchronously after
    # a PR is created, so an immediate read can come back empty; poll briefly
    # before concluding the PR really has no linked issue.
    linked_issues = client.linked_issues_for_pr(number)
    for _ in range(11):
        if linked_issues:
            break
        time.sleep(5)
        linked_issues = client.linked_issues_for_pr(number)

    decision = core.gate_decision(
        author, pr.get("author_association"), author_is_bot, linked_issues
    )
    action = core.gate_action(decision, label_present, close_enabled)

    if action == "none":
        return
    if action == "clear_label":
        client.remove_label(number, core.GATE_LABEL)
        return

    issue = decision["issue"]
    if decision["variant"] == "unassigned":
        body = messages.gate_unassigned(author, issue)
    else:
        body = messages.gate_assigned_other(author, issue, decision["assignee"])
    if action == "flag_and_close":
        body = f"{body}\n\n{messages.gate_close_addendum(issue)}"

    client.add_labels(number, [core.GATE_LABEL])
    client.comment(number, body)
    if action == "flag_and_close":
        client.close_pull_request(number)


if __name__ == "__main__":
    main()
