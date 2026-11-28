"""Assignment gate: flag PRs whose author doesn't hold the linked issue.

A PR against an unclaimed issue is labeled and left open; one against an
issue someone else already holds is labeled, commented on, and closed.

Triggered by ``pull_request_target`` on ``opened`` / ``reopened``.
"""

from __future__ import annotations

import argparse
import json
import time
from typing import Any

from scripts.issue_assignment import (
    core,
    messages,
)
from scripts.issue_assignment.client import GitHubClient


def run_gate(client: GitHubClient, pr: dict[str, Any]) -> None:
    number = pr["number"]
    author = pr["user"]["login"]
    author_is_bot = pr["user"].get("type") == "Bot"
    label_present = core.GATE_LABEL in [label["name"] for label in pr.get("labels", [])]

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
    action = core.gate_action(decision, label_present, close_assigned_other=True)

    if action == "none":
        return
    if action == "clear_label":
        client.remove_label(number, core.GATE_LABEL)
        return
    # GITHUB_TOKEN reports a private org member as CONTRIBUTOR, so
    # ``gate_decision`` cannot exempt them; check their role before flagging.
    if client.is_exempt_collaborator(author):
        return

    closing = action == "flag_and_close"
    client.add_labels(number, [core.GATE_LABEL])
    client.comment(number, messages.gate_flagged(author, decision, closing))
    if closing:
        client.close_pull_request(number)


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", required=True, help="owner/name to operate on")
    parser.add_argument(
        "--event-path", required=True, help="path to the pull_request event payload"
    )
    args = parser.parse_args(argv)

    with open(args.event_path) as fh:
        event = json.load(fh)

    run_gate(GitHubClient(args.repo), event["pull_request"])


if __name__ == "__main__":
    main()
