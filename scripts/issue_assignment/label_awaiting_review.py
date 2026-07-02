"""Reconcile the ``Awaiting Review`` label across all open pull requests.

Run once a day from the scheduled maintenance job (folded in to share a single
runner boot rather than firing on every push). The label is informational: it
flags PRs sitting in the maintainers' court. (The matching state — not the label
— is what exempts a PR from the stale engine in ``unassign_inactive``.) Exempt
authors and gated PRs are skipped. Lagging reality by up to a day is harmless
against the 14-day stale window; current labels are read alongside the review
state so PRs already in the right state cost no write request.
"""

from __future__ import annotations

import os
from typing import (
    TYPE_CHECKING,
    Protocol,
)

from scripts.issue_assignment import core
from scripts.issue_assignment.client import GitHubClient

if TYPE_CHECKING:
    from collections.abc import Iterator

    from scripts.issue_assignment.core import OpenPRState


class SupportsLabelReconcile(Protocol):
    """The slice of :class:`GitHubClient` that ``reconcile_all`` needs."""

    def iter_open_pull_requests_review_state(self) -> Iterator[OpenPRState]: ...

    def add_labels(self, number: int, labels: list[str]) -> None: ...

    def remove_label(self, number: int, label: str) -> None: ...


def reconcile_all(client: SupportsLabelReconcile) -> None:
    label = core.AWAITING_REVIEW_LABEL
    for pr in client.iter_open_pull_requests_review_state():
        if core.is_exempt(pr["author_association"], False):
            continue
        changes_requested_at = core.latest_changes_requested_at(pr["reviews"])
        contributors = {pr["author"]} - {None}
        rereview_requested_at = core.latest_rereview_request_at(
            pr["review_requests"], contributors
        )
        want = core.GATE_LABEL not in pr[
            "labels"
        ] and core.should_label_awaiting_review(
            True, pr["is_draft"], changes_requested_at, rereview_requested_at
        )
        has = label in pr["labels"]
        if want and not has:
            client.add_labels(pr["number"], [label])
        elif has and not want:
            client.remove_label(pr["number"], label)


def main() -> None:
    reconcile_all(GitHubClient(os.environ["GITHUB_REPOSITORY"]))


if __name__ == "__main__":
    main()
