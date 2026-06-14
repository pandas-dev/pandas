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
    messages,
    unassign_inactive,
)

if TYPE_CHECKING:
    from collections.abc import Iterator

NOW = datetime(2026, 6, 1, tzinfo=timezone.utc)


def dt(days_ago: int) -> datetime:
    return NOW - timedelta(days=days_ago)


def review(state: str, days_ago: int, association: str = "MEMBER") -> core.Review:
    return {
        "state": state,
        "submitted_at": dt(days_ago),
        "author_association": association,
    }


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

    def test_changes_requested_no_rereview(self) -> None:
        assert core.awaiting_contributor(dt(2), None) is True

    def test_rereview_requested_after_changes(self) -> None:
        # re-review (dt(2), newer) postdates changes-requested (dt(5)) -> ball moves
        assert core.awaiting_contributor(dt(5), dt(2)) is False

    def test_rereview_predates_latest_changes(self) -> None:
        # changes requested again (dt(2)) after the last re-review (dt(5))
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
            review("COMMENTED", 1),
            review("CHANGES_REQUESTED", 5),
            review("CHANGES_REQUESTED", 3),
            review("APPROVED", 0),
        ]
        assert core.latest_changes_requested_at(reviews) == dt(3)

    def test_latest_changes_requested_none(self) -> None:
        reviews: list[core.Review] = [review("APPROVED", 1)]
        assert core.latest_changes_requested_at(reviews) is None

    def test_changes_requested_by_outsider_ignored(self) -> None:
        # A "request changes" from a non-member doesn't move the ball.
        reviews: list[core.Review] = [
            review("CHANGES_REQUESTED", 3, association="NONE"),
            review("CHANGES_REQUESTED", 2, association="CONTRIBUTOR"),
        ]
        assert core.latest_changes_requested_at(reviews) is None

    def test_changes_requested_owner_counts(self) -> None:
        reviews: list[core.Review] = [
            review("CHANGES_REQUESTED", 3, association="NONE"),
            review("CHANGES_REQUESTED", 2, association="OWNER"),
        ]
        assert core.latest_changes_requested_at(reviews) == dt(2)

    def test_latest_rereview_request_at(self) -> None:
        requests: list[core.ReviewRequest] = [
            {"actor": "alice", "requested_at": dt(8)},
            {"actor": "maintainer", "requested_at": dt(1)},
            {"actor": "alice", "requested_at": dt(4)},
        ]
        assert core.latest_rereview_request_at(requests, {"alice"}) == dt(4)

    def test_latest_rereview_request_none_for_other_actor(self) -> None:
        requests: list[core.ReviewRequest] = [{"actor": "bob", "requested_at": dt(1)}]
        assert core.latest_rereview_request_at(requests, {"alice"}) is None

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
    def test_exempt_author_not_in_scope(self) -> None:
        out = core.gate_decision("maint", "MEMBER", False, [])
        assert out["outcome"] == "not_in_scope" and out["reason"] == "exempt"

    def test_no_linked_issue_not_in_scope(self) -> None:
        out = core.gate_decision("newbie", "NONE", False, [])
        assert out["outcome"] == "not_in_scope" and out["reason"] == "no_linked_issue"

    def test_author_is_assignee_valid(self) -> None:
        linked: list[core.LinkedIssue] = [{"number": 5, "assignees": ["newbie"]}]
        out = core.gate_decision("newbie", "NONE", False, linked)
        assert out["outcome"] == "valid_assignment"

    def test_author_assignee_on_one_of_several_valid(self) -> None:
        linked: list[core.LinkedIssue] = [
            {"number": 5, "assignees": ["someone"]},
            {"number": 6, "assignees": ["newbie"]},
        ]
        out = core.gate_decision("newbie", "NONE", False, linked)
        assert out["outcome"] == "valid_assignment"

    def test_unassigned_issue_invalid(self) -> None:
        linked: list[core.LinkedIssue] = [{"number": 7, "assignees": []}]
        out = core.gate_decision("newbie", "NONE", False, linked)
        assert out == {
            "outcome": "invalid_assignment",
            "variant": "unassigned",
            "issue": 7,
        }

    def test_assigned_to_other_invalid(self) -> None:
        linked: list[core.LinkedIssue] = [{"number": 8, "assignees": ["someone"]}]
        out = core.gate_decision("newbie", "NONE", False, linked)
        assert out == {
            "outcome": "invalid_assignment",
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


class TestPrSubjectToStale:
    def test_exempt_author_never_subject(self) -> None:
        assert core.pr_subject_to_stale(True, False, dt(5), None) is False

    def test_draft_never_subject(self) -> None:
        assert core.pr_subject_to_stale(False, True, dt(5), None) is False

    def test_awaiting_contributor_is_subject(self) -> None:
        # changes requested, no re-review since
        assert core.pr_subject_to_stale(False, False, dt(5), None) is True

    def test_awaiting_review_not_subject(self) -> None:
        # author re-requested review after changes -> back in maintainers' court
        assert core.pr_subject_to_stale(False, False, dt(5), dt(2)) is False

    def test_no_changes_requested_not_subject(self) -> None:
        assert core.pr_subject_to_stale(False, False, None, None) is False


class TestPrStaleAction:
    def test_not_subject_no_label_noop(self) -> None:
        assert core.pr_stale_action(NOW, False, False, dt(99), 14, 7) == "none"

    def test_not_subject_clears_existing_label(self) -> None:
        assert core.pr_stale_action(NOW, False, True, dt(99), 14, 7) == "clear_stale"

    def test_active_author_clears_label(self) -> None:
        assert core.pr_stale_action(NOW, True, True, dt(3), 14, 7) == "clear_stale"

    def test_active_author_no_label_noop(self) -> None:
        assert core.pr_stale_action(NOW, True, False, dt(3), 14, 7) == "none"

    def test_inactive_marks_stale(self) -> None:
        assert core.pr_stale_action(NOW, True, False, dt(15), 14, 7) == "mark_stale"

    def test_stale_within_close_window_waits(self) -> None:
        # stale, inactive 18d: past 14 but inside the 14+7 close window
        assert core.pr_stale_action(NOW, True, True, dt(18), 14, 7) == "none"

    def test_stale_past_close_window_closes(self) -> None:
        assert core.pr_stale_action(NOW, True, True, dt(22), 14, 7) == "close"

    def test_no_action_ever_treated_inactive(self) -> None:
        assert core.pr_stale_action(NOW, True, False, None, 14, 7) == "mark_stale"


class FakeClient:
    def __init__(
        self,
        prs: list[core.OpenPRState],
        linked: dict[int, list[core.LinkedIssue]] | None = None,
        activity: dict[int, core.IssueActivity] | None = None,
    ) -> None:
        self._prs = prs
        self._linked = linked or {}
        self._activity = activity or {}
        self.added: list[tuple[int, tuple[str, ...]]] = []
        self.removed: list[tuple[int, str]] = []
        self.comments: list[tuple[int, str]] = []
        self.unassigned: list[tuple[int, tuple[str, ...]]] = []
        self.closed: list[int] = []

    def iter_open_pull_requests_review_state(self) -> Iterator[core.OpenPRState]:
        return iter(self._prs)

    def linked_issues_for_pr(self, number: int) -> list[core.LinkedIssue]:
        return self._linked.get(number, [])

    def issue_activity(self, number: int) -> core.IssueActivity:
        return self._activity.get(
            number, {"assignees": [], "comments": [], "open_pr_authors": []}
        )

    def add_labels(self, number: int, labels: list[str]) -> None:
        self.added.append((number, tuple(labels)))

    def remove_label(self, number: int, label: str) -> None:
        self.removed.append((number, label))

    def comment(self, number: int, body: str) -> None:
        self.comments.append((number, body))

    def remove_assignees(self, number: int, assignees: list[str]) -> None:
        self.unassigned.append((number, tuple(assignees)))

    def close_pull_request(self, number: int) -> None:
        self.closed.append(number)


def pr(
    number: int,
    *,
    draft: bool = False,
    author: str | None = "alice",
    author_association: str | None = "NONE",
    reviews: list[core.Review] | None = None,
    review_requests: list[core.ReviewRequest] | None = None,
    last_commit_at: datetime | None = None,
    comments: list[core.Comment] | None = None,
    labels: list[str] | None = None,
) -> core.OpenPRState:
    return {
        "number": number,
        "is_draft": draft,
        "author": author,
        "author_association": author_association,
        "reviews": reviews if reviews is not None else [],
        "review_requests": review_requests if review_requests is not None else [],
        "last_commit_at": last_commit_at,
        "comments": comments if comments is not None else [],
        "labels": labels if labels is not None else [],
    }


def rereq(actor: str, days_ago: int) -> core.ReviewRequest:
    return {"actor": actor, "requested_at": dt(days_ago)}


AWAITING = core.AWAITING_REVIEW_LABEL
GATE = core.GATE_LABEL


class TestReconcileAll:
    def test_adds_when_wanted_and_absent(self) -> None:
        client = FakeClient([pr(1)])
        label_awaiting_review.reconcile_all(client)
        assert client.added == [(1, (AWAITING,))]
        assert client.removed == []

    def test_removes_when_unwanted_and_present(self) -> None:
        client = FakeClient([pr(2, draft=True, labels=[AWAITING])])
        label_awaiting_review.reconcile_all(client)
        assert client.removed == [(2, AWAITING)]
        assert client.added == []

    def test_gated_pr_not_labeled(self) -> None:
        client = FakeClient([pr(20, labels=[GATE])])
        label_awaiting_review.reconcile_all(client)
        assert client.added == []
        assert client.removed == []

    def test_gated_pr_loses_awaiting_review(self) -> None:
        client = FakeClient([pr(21, labels=[GATE, AWAITING])])
        label_awaiting_review.reconcile_all(client)
        assert client.removed == [(21, AWAITING)]
        assert client.added == []

    def test_exempt_author_skipped(self) -> None:
        client = FakeClient([pr(22, author_association="MEMBER")])
        label_awaiting_review.reconcile_all(client)
        assert client.added == []
        assert client.removed == []

    def test_noop_when_already_correct(self) -> None:
        client = FakeClient(
            [
                pr(3, labels=[AWAITING]),
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
                    reviews=[review("CHANGES_REQUESTED", 3)],
                    labels=[AWAITING],
                )
            ]
        )
        label_awaiting_review.reconcile_all(client)
        assert client.removed == [(5, AWAITING)]
        assert client.added == []

    def test_push_alone_does_not_label(self) -> None:
        # A commit after changes-requested is no longer enough; without a
        # re-review request the PR stays in the contributor's court.
        client = FakeClient([pr(6, reviews=[review("CHANGES_REQUESTED", 3)])])
        label_awaiting_review.reconcile_all(client)
        assert client.added == []
        assert client.removed == []

    def test_rereview_request_relabels(self) -> None:
        client = FakeClient(
            [
                pr(
                    7,
                    reviews=[review("CHANGES_REQUESTED", 3)],
                    review_requests=[rereq("alice", 1)],
                )
            ]
        )
        label_awaiting_review.reconcile_all(client)
        assert client.added == [(7, (AWAITING,))]

    def test_rereview_request_by_non_contributor_ignored(self) -> None:
        # A maintainer wrangling reviewers isn't the contributor signalling ready.
        client = FakeClient(
            [
                pr(
                    8,
                    author="alice",
                    reviews=[review("CHANGES_REQUESTED", 3)],
                    review_requests=[rereq("maintainer", 1)],
                )
            ]
        )
        label_awaiting_review.reconcile_all(client)
        assert client.added == []
        assert client.removed == []


STALE = core.STALE_LABEL


def activity(comments: list[core.Comment]) -> core.IssueActivity:
    return {"assignees": [], "comments": comments, "open_pr_authors": []}


def ago(days: int) -> datetime:
    # run_pr_stale_sweep reads the wall clock, so fixtures are relative to real now.
    return datetime.now(timezone.utc) - timedelta(days=days)


CHANGES = [review("CHANGES_REQUESTED", 5)]  # a member changes-requested review


class TestRunPrStaleSweep:
    def test_marks_stale_when_author_inactive(self) -> None:
        client = FakeClient([pr(1, reviews=CHANGES, last_commit_at=ago(20))])
        unassign_inactive.run_pr_stale_sweep(client)
        assert client.added == [(1, (STALE,))]
        assert client.comments == [(1, messages.pr_marked_stale())]
        assert client.closed == []

    def test_active_author_clears_stale(self) -> None:
        client = FakeClient(
            [pr(2, reviews=CHANGES, last_commit_at=ago(2), labels=[STALE])]
        )
        unassign_inactive.run_pr_stale_sweep(client)
        assert client.removed == [(2, STALE)]
        assert client.added == []

    def test_linked_issue_comment_prevents_stale(self) -> None:
        # PR untouched 20d, but the author commented on a linked issue 2d ago.
        client = FakeClient(
            [pr(3, reviews=CHANGES, last_commit_at=ago(20))],
            linked={3: [{"number": 30, "assignees": ["alice"]}]},
            activity={30: activity([{"author": "alice", "created_at": ago(2)}])},
        )
        unassign_inactive.run_pr_stale_sweep(client)
        assert client.added == []
        assert client.removed == []

    def test_linked_issue_comment_clears_existing_stale(self) -> None:
        client = FakeClient(
            [pr(4, reviews=CHANGES, last_commit_at=ago(20), labels=[STALE])],
            linked={4: [{"number": 40, "assignees": ["alice"]}]},
            activity={40: activity([{"author": "alice", "created_at": ago(2)}])},
        )
        unassign_inactive.run_pr_stale_sweep(client)
        assert client.removed == [(4, STALE)]

    def test_exempt_author_never_stale(self) -> None:
        client = FakeClient(
            [
                pr(
                    5,
                    author_association="MEMBER",
                    reviews=CHANGES,
                    last_commit_at=ago(99),
                )
            ]
        )
        unassign_inactive.run_pr_stale_sweep(client)
        assert client.added == []
        assert client.closed == []

    def test_awaiting_review_not_swept(self) -> None:
        # No changes requested -> maintainers' court -> never stale, even if old.
        client = FakeClient([pr(6, last_commit_at=ago(99))])
        unassign_inactive.run_pr_stale_sweep(client)
        assert client.added == []
        assert client.closed == []

    def test_closes_and_unassigns_after_window(self) -> None:
        client = FakeClient(
            [pr(7, reviews=CHANGES, last_commit_at=ago(22), labels=[STALE])],
            linked={7: [{"number": 70, "assignees": ["alice", "bob"]}]},
        )
        unassign_inactive.run_pr_stale_sweep(client)
        assert client.closed == [7]
        assert client.unassigned == [(70, ("alice",))]
        assert client.comments == [
            (7, messages.pr_closed_stale()),
            (70, messages.issue_freed_stale_pr()),
        ]

    def test_stale_within_close_window_waits(self) -> None:
        client = FakeClient(
            [pr(8, reviews=CHANGES, last_commit_at=ago(18), labels=[STALE])]
        )
        unassign_inactive.run_pr_stale_sweep(client)
        assert client.closed == []
        assert client.removed == []
        assert client.added == []
