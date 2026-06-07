from __future__ import annotations

from datetime import (
    datetime,
    timedelta,
    timezone,
)
from typing import TYPE_CHECKING

import pytest

from scripts.issue_assignment import (
    core,
    label_awaiting_review,
)

if TYPE_CHECKING:
    from collections.abc import Iterator

NOW = datetime(2026, 6, 1, tzinfo=timezone.utc)


def dt(days_ago: int) -> datetime:
    return NOW - timedelta(days=days_ago)


class TestIsExempt:
    @pytest.mark.parametrize("association", ["OWNER", "MEMBER", "COLLABORATOR"])
    def test_exempt_associations(self, association: str) -> None:
        assert core.is_exempt(association, False)

    @pytest.mark.parametrize(
        "association", ["CONTRIBUTOR", "NONE", "FIRST_TIMER", None]
    )
    def test_non_exempt_associations(self, association: str | None) -> None:
        assert not core.is_exempt(association, False)

    def test_bot_is_exempt_regardless(self) -> None:
        assert core.is_exempt("NONE", True)


class TestAwaitingContributor:
    def test_no_changes_requested(self) -> None:
        assert core.awaiting_contributor(None, dt(1)) is False

    def test_changes_requested_no_commits(self) -> None:
        assert core.awaiting_contributor(dt(2), None) is True

    def test_pushed_after_changes_requested(self) -> None:
        assert core.awaiting_contributor(dt(5), dt(2)) is False

    def test_changes_requested_after_last_push(self) -> None:
        assert core.awaiting_contributor(dt(2), dt(5)) is True


class TestShouldLabelAwaitingReview:
    def test_open_no_review(self) -> None:
        assert core.should_label_awaiting_review(True, False, None, dt(1)) is True

    def test_open_awaiting_contributor(self) -> None:
        assert core.should_label_awaiting_review(True, False, dt(1), dt(3)) is False

    def test_draft_never_labeled(self) -> None:
        assert core.should_label_awaiting_review(True, True, None, dt(1)) is False

    def test_closed_never_labeled(self) -> None:
        assert core.should_label_awaiting_review(False, False, None, dt(1)) is False


class TestLatestSelectors:
    def test_latest_changes_requested_at(self) -> None:
        reviews: list[core.Review] = [
            {"state": "COMMENTED", "submitted_at": dt(1)},
            {"state": "CHANGES_REQUESTED", "submitted_at": dt(5)},
            {"state": "CHANGES_REQUESTED", "submitted_at": dt(3)},
            {"state": "APPROVED", "submitted_at": dt(0)},
        ]
        assert core.latest_changes_requested_at(reviews) == dt(3)

    def test_latest_changes_requested_none(self) -> None:
        reviews: list[core.Review] = [{"state": "APPROVED", "submitted_at": dt(1)}]
        assert core.latest_changes_requested_at(reviews) is None

    def test_latest_assignee_comment_at(self) -> None:
        comments: list[core.Comment] = [
            {"author": "alice", "created_at": dt(10)},
            {"author": "bob", "created_at": dt(1)},
            {"author": "alice", "created_at": dt(4)},
        ]
        assert core.latest_assignee_comment_at(comments, ["alice"]) == dt(4)

    def test_latest_assignee_comment_ignores_non_assignees(self) -> None:
        comments: list[core.Comment] = [{"author": "bob", "created_at": dt(1)}]
        assert core.latest_assignee_comment_at(comments, ["alice"]) is None


class TestGateDecision:
    def test_exempt_author_skipped(self) -> None:
        out = core.gate_decision("maint", "MEMBER", False, [])
        assert out["action"] == "skip" and out["reason"] == "exempt"

    def test_no_linked_issue_skipped(self) -> None:
        out = core.gate_decision("newbie", "NONE", False, [])
        assert out["action"] == "skip" and out["reason"] == "no_linked_issue"

    def test_author_is_assignee_passes(self) -> None:
        linked: list[core.LinkedIssue] = [{"number": 5, "assignees": ["newbie"]}]
        assert core.gate_decision("newbie", "NONE", False, linked)["action"] == "pass"

    def test_author_assignee_on_one_of_several_passes(self) -> None:
        linked: list[core.LinkedIssue] = [
            {"number": 5, "assignees": ["someone"]},
            {"number": 6, "assignees": ["newbie"]},
        ]
        assert core.gate_decision("newbie", "NONE", False, linked)["action"] == "pass"

    def test_unassigned_issue_flagged(self) -> None:
        linked: list[core.LinkedIssue] = [{"number": 7, "assignees": []}]
        out = core.gate_decision("newbie", "NONE", False, linked)
        assert out == {"action": "flag", "variant": "unassigned", "issue": 7}

    def test_assigned_to_other_flagged(self) -> None:
        linked: list[core.LinkedIssue] = [{"number": 8, "assignees": ["someone"]}]
        out = core.gate_decision("newbie", "NONE", False, linked)
        assert out == {
            "action": "flag",
            "variant": "assigned_other",
            "issue": 8,
            "assignee": "someone",
        }

    def test_unassigned_preferred_over_assigned_other(self) -> None:
        linked: list[core.LinkedIssue] = [
            {"number": 8, "assignees": ["someone"]},
            {"number": 9, "assignees": []},
        ]
        out = core.gate_decision("newbie", "NONE", False, linked)
        assert out["variant"] == "unassigned" and out["issue"] == 9


class TestIssueIsActive:
    def test_open_pr_keeps_claim_even_when_old(self) -> None:
        assert core.issue_is_active(NOW, True, dt(99), 14) is True

    def test_recent_assignee_comment_active(self) -> None:
        assert core.issue_is_active(NOW, False, dt(13), 14) is True

    def test_stale_comment_inactive(self) -> None:
        assert core.issue_is_active(NOW, False, dt(15), 14) is False

    def test_no_pr_no_comment_inactive(self) -> None:
        assert core.issue_is_active(NOW, False, None, 14) is False

    def test_boundary_exactly_stale_days_inactive(self) -> None:
        assert core.issue_is_active(NOW, False, dt(14), 14) is False


class FakeClient:
    def __init__(self, prs: list[core.OpenPRState]) -> None:
        self._prs = prs
        self.added: list[tuple[int, tuple[str, ...]]] = []
        self.removed: list[tuple[int, str]] = []

    def iter_open_pull_requests_review_state(self) -> Iterator[core.OpenPRState]:
        return iter(self._prs)

    def add_labels(self, number: int, labels: list[str]) -> None:
        self.added.append((number, tuple(labels)))

    def remove_label(self, number: int, label: str) -> None:
        self.removed.append((number, label))


def pr(
    number: int,
    *,
    draft: bool = False,
    reviews: list[core.Review] | None = None,
    last_commit_at: datetime | None = None,
    labels: list[str] | None = None,
) -> core.OpenPRState:
    return {
        "number": number,
        "is_draft": draft,
        "reviews": reviews if reviews is not None else [],
        "last_commit_at": last_commit_at,
        "labels": labels if labels is not None else [],
    }


AWAITING = core.AWAITING_REVIEW_LABEL


class TestReconcileAll:
    def test_adds_when_wanted_and_absent(self) -> None:
        client = FakeClient([pr(1, last_commit_at=dt(1))])
        label_awaiting_review.reconcile_all(client)
        assert client.added == [(1, (AWAITING,))]
        assert client.removed == []

    def test_removes_when_unwanted_and_present(self) -> None:
        client = FakeClient([pr(2, draft=True, labels=[AWAITING])])
        label_awaiting_review.reconcile_all(client)
        assert client.removed == [(2, AWAITING)]
        assert client.added == []

    def test_noop_when_already_correct(self) -> None:
        client = FakeClient(
            [
                pr(3, last_commit_at=dt(1), labels=[AWAITING]),
                pr(4, draft=True),
            ]
        )
        label_awaiting_review.reconcile_all(client)
        assert client.added == []
        assert client.removed == []

    def test_awaiting_contributor_pr_gets_label_removed(self) -> None:
        client = FakeClient(
            [
                pr(
                    5,
                    reviews=[{"state": "CHANGES_REQUESTED", "submitted_at": dt(1)}],
                    last_commit_at=dt(3),
                    labels=[AWAITING],
                )
            ]
        )
        label_awaiting_review.reconcile_all(client)
        assert client.removed == [(5, AWAITING)]
        assert client.added == []
