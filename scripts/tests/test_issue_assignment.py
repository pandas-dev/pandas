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


def review(
    state: str,
    days_ago: int,
    association: str = "MEMBER",
    author: str = "maint",
) -> core.Review:
    return {
        "author": author,
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


class TestOutstandingChangesRequested:
    def test_outstanding_changes_requested_at(self) -> None:
        reviews: list[core.Review] = [
            review("CHANGES_REQUESTED", 5),
            review("CHANGES_REQUESTED", 3),
        ]
        assert core.outstanding_changes_requested_at(reviews) == dt(3)

    def test_none_without_changes_requested(self) -> None:
        reviews: list[core.Review] = [review("APPROVED", 1)]
        assert core.outstanding_changes_requested_at(reviews) is None

    def test_changes_requested_by_outsider_ignored(self) -> None:
        # A "request changes" from a non-member doesn't move the ball.
        reviews: list[core.Review] = [
            review("CHANGES_REQUESTED", 3, association="NONE", author="alice"),
            review("CHANGES_REQUESTED", 2, association="CONTRIBUTOR", author="bob"),
        ]
        assert core.outstanding_changes_requested_at(reviews) is None

    def test_changes_requested_owner_counts(self) -> None:
        reviews: list[core.Review] = [
            review("CHANGES_REQUESTED", 3, association="NONE", author="alice"),
            review("CHANGES_REQUESTED", 2, association="OWNER"),
        ]
        assert core.outstanding_changes_requested_at(reviews) == dt(2)

    def test_own_approval_supersedes_changes_requested(self) -> None:
        # The reviewer approving clears their outstanding request — the PR
        # must not sit in the contributor's court after being approved.
        reviews: list[core.Review] = [
            review("CHANGES_REQUESTED", 5),
            review("APPROVED", 2),
        ]
        assert core.outstanding_changes_requested_at(reviews) is None

    def test_comment_review_leaves_request_standing(self) -> None:
        reviews: list[core.Review] = [
            review("CHANGES_REQUESTED", 5),
            review("COMMENTED", 2),
        ]
        assert core.outstanding_changes_requested_at(reviews) == dt(5)

    def test_dismissed_review_does_not_count(self) -> None:
        # Dismissal rewrites the review's state in place.
        reviews: list[core.Review] = [review("DISMISSED", 5)]
        assert core.outstanding_changes_requested_at(reviews) is None

    def test_other_reviewer_approval_does_not_clear(self) -> None:
        # Reviewer beta approving leaves alpha's request outstanding.
        reviews: list[core.Review] = [
            review("CHANGES_REQUESTED", 5, author="alpha"),
            review("APPROVED", 2, author="beta"),
        ]
        assert core.outstanding_changes_requested_at(reviews) == dt(5)

    def test_renewed_changes_request_counts_again(self) -> None:
        reviews: list[core.Review] = [
            review("CHANGES_REQUESTED", 9),
            review("APPROVED", 6),
            review("CHANGES_REQUESTED", 2),
        ]
        assert core.outstanding_changes_requested_at(reviews) == dt(2)


class TestLatestSelectors:
    def test_latest_review_submitted_at(self) -> None:
        reviews: list[core.Review] = [
            review("COMMENTED", 9, association="NONE", author="alice"),
            review("COMMENTED", 4, association="NONE", author="alice"),
            review("APPROVED", 1),
        ]
        assert core.latest_review_submitted_at(reviews, {"alice"}) == dt(4)

    def test_latest_review_submitted_none_for_others(self) -> None:
        reviews: list[core.Review] = [review("APPROVED", 1)]
        assert core.latest_review_submitted_at(reviews, {"alice"}) is None

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


INVALID_DECISION: core.GateDecision = {
    "outcome": "invalid_assignment",
    "variant": "unassigned",
    "issue": 7,
}


class TestGateAction:
    def test_not_in_scope_none(self) -> None:
        decision: core.GateDecision = {"outcome": "not_in_scope", "reason": "exempt"}
        assert core.gate_action(decision, False, False) == "none"

    def test_valid_clears_present_label(self) -> None:
        decision: core.GateDecision = {"outcome": "valid_assignment"}
        assert core.gate_action(decision, True, False) == "clear_label"

    def test_valid_without_label_noop(self) -> None:
        decision: core.GateDecision = {"outcome": "valid_assignment"}
        assert core.gate_action(decision, False, False) == "none"

    def test_invalid_flags(self) -> None:
        assert core.gate_action(INVALID_DECISION, False, False) == "flag"

    def test_invalid_already_flagged_not_recommented(self) -> None:
        # Reopening without fixing the assignment mustn't repost the comment.
        assert core.gate_action(INVALID_DECISION, True, False) == "none"

    @pytest.mark.parametrize("label_present", [False, True])
    def test_close_mode_always_comments_and_closes(self, label_present: bool) -> None:
        out = core.gate_action(INVALID_DECISION, label_present, True)
        assert out == "flag_and_close"


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
        assert core.pr_stale_action(NOW, False, False, None, dt(99), 14, 7) == "none"

    def test_not_subject_clears_existing_label(self) -> None:
        out = core.pr_stale_action(NOW, False, True, dt(20), dt(99), 14, 7)
        assert out == "clear_stale"

    def test_active_author_clears_label(self) -> None:
        out = core.pr_stale_action(NOW, True, True, dt(20), dt(3), 14, 7)
        assert out == "clear_stale"

    def test_active_author_no_label_noop(self) -> None:
        assert core.pr_stale_action(NOW, True, False, None, dt(3), 14, 7) == "none"

    def test_inactive_marks_stale(self) -> None:
        out = core.pr_stale_action(NOW, True, False, None, dt(15), 14, 7)
        assert out == "mark_stale"

    def test_recent_mark_waits(self) -> None:
        # warned 4d ago: inside the 7-day close window
        assert core.pr_stale_action(NOW, True, True, dt(4), dt(18), 14, 7) == "none"

    def test_mark_older_than_close_window_closes(self) -> None:
        out = core.pr_stale_action(NOW, True, True, dt(8), dt(22), 14, 7)
        assert out == "close"

    def test_close_counts_from_mark_not_from_inactivity(self) -> None:
        # Author inactive 40d but the engine only warned 2d ago: the close
        # countdown runs from the warning, so it waits.
        assert core.pr_stale_action(NOW, True, True, dt(2), dt(40), 14, 7) == "none"

    def test_foreign_label_remarked_not_closed(self) -> None:
        # Stale label present but the engine never applied it (no label event):
        # warn first, never close cold.
        out = core.pr_stale_action(NOW, True, True, None, dt(40), 14, 7)
        assert out == "mark_stale"

    def test_no_action_ever_treated_inactive(self) -> None:
        out = core.pr_stale_action(NOW, True, False, None, None, 14, 7)
        assert out == "mark_stale"


class FakeClient:
    def __init__(
        self,
        prs: list[core.OpenPRState],
        linked: dict[int, list[core.LinkedIssue]] | None = None,
        activity: dict[int, core.IssueActivity] | None = None,
        exempt: set[str] | None = None,
        exact_comments: dict[int, datetime] | None = None,
        exact_pr_authors: dict[int, list[str]] | None = None,
    ) -> None:
        self._prs = prs
        self._linked = linked or {}
        self._activity = activity or {}
        self._exempt = exempt or set()
        self._exact_comments = exact_comments or {}
        self._exact_pr_authors = exact_pr_authors or {}
        self.added: list[tuple[int, tuple[str, ...]]] = []
        self.removed: list[tuple[int, str]] = []
        self.comments: list[tuple[int, str]] = []
        self.unassigned: list[tuple[int, tuple[str, ...]]] = []
        self.closed: list[int] = []
        self.exact_comment_reads: list[int] = []

    def iter_open_pull_requests_review_state(self) -> Iterator[core.OpenPRState]:
        return iter(self._prs)

    def linked_issues_for_pr(self, number: int) -> list[core.LinkedIssue]:
        return self._linked.get(number, [])

    def list_open_assigned_issue_numbers(self) -> list[int]:
        return sorted(self._activity)

    def issue_activity(self, number: int) -> core.IssueActivity:
        return self._activity.get(
            number,
            {
                "assignees": [],
                "assigned_at": None,
                "comments": [],
                "comments_truncated": False,
                "open_pr_authors": [],
                "crossrefs_truncated": False,
            },
        )

    def latest_comment_at_since(
        self, number: int, logins: set[str], since: datetime
    ) -> datetime | None:
        self.exact_comment_reads.append(number)
        found = self._exact_comments.get(number)
        return found if found is not None and found >= since else None

    def open_pr_authors_exact(self, number: int) -> list[str]:
        return self._exact_pr_authors.get(number, [])

    def is_exempt_collaborator(self, username: str) -> bool:
        return username in self._exempt

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
    comments_truncated: bool = False,
    reopened_events: list[core.Comment] | None = None,
    labels: list[str] | None = None,
    stale_marked_at: datetime | None = None,
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
        "comments_truncated": comments_truncated,
        "reopened_events": reopened_events if reopened_events is not None else [],
        "labels": labels if labels is not None else [],
        "stale_marked_at": stale_marked_at,
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


def activity(
    comments: list[core.Comment],
    assignees: list[str] | None = None,
    assigned_at: datetime | None = None,
    comments_truncated: bool = False,
    open_pr_authors: list[str] | None = None,
    crossrefs_truncated: bool = False,
) -> core.IssueActivity:
    return {
        "assignees": assignees if assignees is not None else [],
        "assigned_at": assigned_at,
        "comments": comments,
        "comments_truncated": comments_truncated,
        "open_pr_authors": open_pr_authors if open_pr_authors is not None else [],
        "crossrefs_truncated": crossrefs_truncated,
    }


def ago(days: int) -> datetime:
    # The sweeps read the wall clock, so fixtures are relative to real now.
    return datetime.now(timezone.utc) - timedelta(days=days)


def changes(days_ago: int) -> list[core.Review]:
    # A member changes-requested review, dated relative to real now (the sweep's
    # clock anchors on this time, so it must line up with the ``ago``-based fixtures).
    return [
        {
            "author": "maint",
            "state": "CHANGES_REQUESTED",
            "submitted_at": ago(days_ago),
            "author_association": "MEMBER",
        }
    ]


class TestRunPrStaleSweep:
    def test_marks_stale_when_author_inactive(self) -> None:
        client = FakeClient([pr(1, reviews=changes(30), last_commit_at=ago(20))])
        unassign_inactive.run_pr_stale_sweep(client)
        assert client.added == [(1, (STALE,))]
        assert client.comments == [(1, messages.pr_marked_stale())]
        assert client.closed == []

    def test_active_author_clears_stale(self) -> None:
        client = FakeClient(
            [pr(2, reviews=changes(30), last_commit_at=ago(2), labels=[STALE])]
        )
        unassign_inactive.run_pr_stale_sweep(client)
        assert client.removed == [(2, STALE)]
        assert client.added == []

    def test_recent_changes_request_floors_clock(self) -> None:
        # Author pushed 30d ago, then review took a while and changes were
        # requested only 2d ago. The clock runs from the review, not the push,
        # so the PR is NOT stale — review latency never penalizes the author.
        client = FakeClient([pr(3, reviews=changes(2), last_commit_at=ago(30))])
        unassign_inactive.run_pr_stale_sweep(client)
        assert client.added == []
        assert client.closed == []

    def test_linked_issue_comment_prevents_stale(self) -> None:
        # PR untouched 20d, but the author commented on a linked issue 2d ago.
        client = FakeClient(
            [pr(4, reviews=changes(30), last_commit_at=ago(20))],
            linked={4: [{"number": 40, "assignees": ["alice"]}]},
            activity={40: activity([{"author": "alice", "created_at": ago(2)}])},
        )
        unassign_inactive.run_pr_stale_sweep(client)
        assert client.added == []
        assert client.removed == []

    def test_linked_issue_comment_clears_existing_stale(self) -> None:
        client = FakeClient(
            [pr(5, reviews=changes(30), last_commit_at=ago(20), labels=[STALE])],
            linked={5: [{"number": 50, "assignees": ["alice"]}]},
            activity={50: activity([{"author": "alice", "created_at": ago(2)}])},
        )
        unassign_inactive.run_pr_stale_sweep(client)
        assert client.removed == [(5, STALE)]

    def test_exempt_author_never_stale(self) -> None:
        client = FakeClient(
            [
                pr(
                    6,
                    author_association="MEMBER",
                    reviews=changes(30),
                    last_commit_at=ago(99),
                )
            ]
        )
        unassign_inactive.run_pr_stale_sweep(client)
        assert client.added == []
        assert client.closed == []

    def test_awaiting_review_not_swept(self) -> None:
        # No changes requested -> maintainers' court -> never stale, even if old.
        client = FakeClient([pr(7, last_commit_at=ago(99))])
        unassign_inactive.run_pr_stale_sweep(client)
        assert client.added == []
        assert client.closed == []

    def test_closes_and_unassigns_after_window(self) -> None:
        client = FakeClient(
            [
                pr(
                    8,
                    reviews=changes(30),
                    last_commit_at=ago(22),
                    labels=[STALE],
                    stale_marked_at=ago(8),
                )
            ],
            linked={8: [{"number": 80, "assignees": ["alice", "bob"]}]},
        )
        unassign_inactive.run_pr_stale_sweep(client)
        assert client.closed == [8]
        assert client.unassigned == [(80, ("alice",))]
        assert client.comments == [
            (8, messages.pr_closed_stale()),
            (80, messages.issue_freed_stale_pr()),
        ]

    def test_stale_within_close_window_waits(self) -> None:
        client = FakeClient(
            [
                pr(
                    9,
                    reviews=changes(30),
                    last_commit_at=ago(18),
                    labels=[STALE],
                    stale_marked_at=ago(4),
                )
            ]
        )
        unassign_inactive.run_pr_stale_sweep(client)
        assert client.closed == []
        assert client.removed == []
        assert client.added == []

    def test_close_counts_from_mark_not_from_inactivity(self) -> None:
        # Author inactive far past 14+7, but the engine warned only 2d ago.
        client = FakeClient(
            [
                pr(
                    10,
                    reviews=changes(60),
                    last_commit_at=ago(50),
                    labels=[STALE],
                    stale_marked_at=ago(2),
                )
            ]
        )
        unassign_inactive.run_pr_stale_sweep(client)
        assert client.closed == []
        assert client.added == []
        assert client.removed == []

    def test_foreign_stale_label_recycled_with_warning(self) -> None:
        # Label present but the engine never applied it (human / old bot):
        # cycle the label for a fresh event and warn — never close cold.
        client = FakeClient(
            [pr(11, reviews=changes(60), last_commit_at=ago(50), labels=[STALE])]
        )
        unassign_inactive.run_pr_stale_sweep(client)
        assert client.closed == []
        assert client.removed == [(11, STALE)]
        assert client.added == [(11, (STALE,))]
        assert client.comments == [(11, messages.pr_marked_stale())]

    @pytest.mark.parametrize("actor", ["alice", "maint"])
    def test_reopen_counts_as_activity(self, actor: str) -> None:
        # Reopening — whether by the author or a maintainer — is an explicit
        # act of keeping the PR alive; it must not be insta-reclosed.
        client = FakeClient(
            [
                pr(
                    12,
                    reviews=changes(30),
                    last_commit_at=ago(20),
                    reopened_events=[{"author": actor, "created_at": ago(2)}],
                )
            ]
        )
        unassign_inactive.run_pr_stale_sweep(client)
        assert client.added == []
        assert client.closed == []

    def test_author_review_reply_counts_as_activity(self) -> None:
        # An inline reply creates an implicit COMMENTED review by the author.
        reply: core.Review = {
            "author": "alice",
            "state": "COMMENTED",
            "submitted_at": ago(2),
            "author_association": "NONE",
        }
        client = FakeClient(
            [pr(13, reviews=[*changes(30), reply], last_commit_at=ago(20))]
        )
        unassign_inactive.run_pr_stale_sweep(client)
        assert client.added == []
        assert client.closed == []

    def test_truncated_comments_checked_exactly_before_marking(self) -> None:
        client = FakeClient(
            [
                pr(
                    14,
                    reviews=changes(30),
                    last_commit_at=ago(20),
                    comments_truncated=True,
                )
            ],
            exact_comments={14: ago(2)},
        )
        unassign_inactive.run_pr_stale_sweep(client)
        assert client.exact_comment_reads == [14]
        assert client.added == []

    def test_untruncated_comments_skip_exact_read(self) -> None:
        client = FakeClient([pr(15, reviews=changes(30), last_commit_at=ago(20))])
        unassign_inactive.run_pr_stale_sweep(client)
        assert client.exact_comment_reads == []
        assert client.added == [(15, (STALE,))]


class TestRunSweep:
    @pytest.fixture(autouse=True)
    def _rollout_long_past(self, monkeypatch: pytest.MonkeyPatch) -> None:
        monkeypatch.setattr(core, "ROLLOUT_CUTOFF", ago(100))

    def test_unassigns_inactive_assignee(self) -> None:
        client = FakeClient(
            [], activity={1: activity([], assignees=["alice"], assigned_at=ago(30))}
        )
        unassign_inactive.run_sweep(client)
        assert client.unassigned == [(1, ("alice",))]
        assert client.comments == [(1, messages.issue_unassigned_inactive(["alice"]))]

    def test_fresh_assignment_without_comment_kept(self) -> None:
        # Assigned 2d ago, never commented: the clock runs from the assignment.
        client = FakeClient(
            [], activity={2: activity([], assignees=["alice"], assigned_at=ago(2))}
        )
        unassign_inactive.run_sweep(client)
        assert client.unassigned == []

    def test_recent_assignee_comment_kept(self) -> None:
        client = FakeClient(
            [],
            activity={
                3: activity(
                    [{"author": "alice", "created_at": ago(3)}],
                    assignees=["alice"],
                    assigned_at=ago(30),
                )
            },
        )
        unassign_inactive.run_sweep(client)
        assert client.unassigned == []

    def test_open_pr_keeps_claim(self) -> None:
        client = FakeClient(
            [],
            activity={
                4: activity(
                    [],
                    assignees=["alice"],
                    assigned_at=ago(90),
                    open_pr_authors=["alice"],
                )
            },
        )
        unassign_inactive.run_sweep(client)
        assert client.unassigned == []

    def test_pre_rollout_assignment_gets_grace(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        # Assignment long predates the cutoff: the clock runs from the cutoff,
        # so nothing is expired retroactively at rollout.
        monkeypatch.setattr(core, "ROLLOUT_CUTOFF", ago(5))
        client = FakeClient(
            [], activity={5: activity([], assignees=["alice"], assigned_at=ago(400))}
        )
        unassign_inactive.run_sweep(client)
        assert client.unassigned == []

    def test_exempt_assignee_never_unassigned(self) -> None:
        client = FakeClient(
            [],
            activity={6: activity([], assignees=["maint"], assigned_at=ago(90))},
            exempt={"maint"},
        )
        unassign_inactive.run_sweep(client)
        assert client.unassigned == []
        assert client.comments == []

    def test_only_non_exempt_co_assignee_unassigned(self) -> None:
        client = FakeClient(
            [],
            activity={
                7: activity([], assignees=["maint", "alice"], assigned_at=ago(90))
            },
            exempt={"maint"},
        )
        unassign_inactive.run_sweep(client)
        assert client.unassigned == [(7, ("alice",))]
        assert client.comments == [(7, messages.issue_unassigned_inactive(["alice"]))]

    def test_truncated_comments_checked_exactly(self) -> None:
        client = FakeClient(
            [],
            activity={
                8: activity(
                    [],
                    assignees=["alice"],
                    assigned_at=ago(30),
                    comments_truncated=True,
                )
            },
            exact_comments={8: ago(2)},
        )
        unassign_inactive.run_sweep(client)
        assert client.exact_comment_reads == [8]
        assert client.unassigned == []

    def test_truncated_crossrefs_checked_exactly(self) -> None:
        client = FakeClient(
            [],
            activity={
                9: activity(
                    [],
                    assignees=["alice"],
                    assigned_at=ago(30),
                    crossrefs_truncated=True,
                )
            },
            exact_pr_authors={9: ["alice"]},
        )
        unassign_inactive.run_sweep(client)
        assert client.unassigned == []
