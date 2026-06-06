"""Reconcile the ``Awaiting Review`` label across all open pull requests.

Run once a day from the scheduled maintenance job (folded in to share a single
runner boot rather than firing on every push). The label exempts review-blocked
PRs from going stale (see ``stale-pr.yml``); lagging reality by up to a day is
harmless against the 14-day stale window. Current labels are read alongside the
review state so PRs already in the right state cost no write request.
"""

from __future__ import annotations

import os

from scripts.issue_assignment import core
from scripts.issue_assignment.client import GitHubClient


def reconcile_all(client):
    label = core.AWAITING_REVIEW_LABEL
    for pr in client.iter_open_pull_requests_review_state():
        changes_requested_at = core.latest_changes_requested_at(pr["reviews"])
        want = core.should_label_awaiting_review(
            True, pr["is_draft"], changes_requested_at, pr["last_commit_at"]
        )
        has = label in pr["labels"]
        if want and not has:
            client.add_labels(pr["number"], [label])
        elif has and not want:
            client.remove_label(pr["number"], label)


def main():
    reconcile_all(GitHubClient(os.environ["GITHUB_REPOSITORY"]))


if __name__ == "__main__":
    main()
