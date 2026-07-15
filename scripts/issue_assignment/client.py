"""Thin GitHub API layer built on the pre-authenticated ``gh`` CLI.

This is the only module that performs I/O, so it is deliberately small; all
judgement lives in :mod:`core`. ``gh`` reads ``GH_TOKEN`` / ``GITHUB_TOKEN``
from the environment, so no explicit auth wiring is needed.
"""

from __future__ import annotations

from datetime import (
    datetime,
    timezone,
)
import json
import subprocess
from typing import (
    TYPE_CHECKING,
    Any,
)
from urllib.parse import quote

from scripts.issue_assignment.core import STALE_LABEL

if TYPE_CHECKING:
    from collections.abc import Iterator

    from scripts.issue_assignment.core import (
        Comment,
        IssueActivity,
        LinkedIssue,
        OpenPRState,
        Review,
        ReviewRequest,
    )

_LINKED_ISSUES_QUERY = """
query($owner: String!, $name: String!, $number: Int!) {
  repository(owner: $owner, name: $name) {
    pullRequest(number: $number) {
      closingIssuesReferences(first: 20) {
        nodes {
          number
          repository { nameWithOwner }
          assignees(first: 20) { nodes { login } }
        }
      }
    }
  }
}
"""

_ISSUE_ACTIVITY_QUERY = """
query($owner: String!, $name: String!, $number: Int!) {
  repository(owner: $owner, name: $name) {
    issue(number: $number) {
      assignees(first: 20) { nodes { login } }
      comments(last: 50) { totalCount nodes { author { login } createdAt } }
      crossrefs: timelineItems(last: 50, itemTypes: [CROSS_REFERENCED_EVENT]) {
        totalCount
        nodes {
          ... on CrossReferencedEvent {
            source { ... on PullRequest { state author { login } } }
          }
        }
      }
      assignedEvents: timelineItems(last: 50, itemTypes: [ASSIGNED_EVENT]) {
        nodes { ... on AssignedEvent { createdAt } }
      }
    }
  }
}
"""

_OPEN_PRS_QUERY = """
query($owner: String!, $name: String!, $cursor: String) {
  repository(owner: $owner, name: $name) {
    pullRequests(states: OPEN, first: 50, after: $cursor) {
      pageInfo { hasNextPage endCursor }
      nodes {
        number
        isDraft
        author { login }
        authorAssociation
        reviews(last: 50) {
          nodes { author { login } state submittedAt authorAssociation }
        }
        reviewRequests: timelineItems(
          last: 50, itemTypes: [REVIEW_REQUESTED_EVENT]
        ) {
          nodes { ... on ReviewRequestedEvent { createdAt actor { login } } }
        }
        labelEvents: timelineItems(last: 50, itemTypes: [LABELED_EVENT]) {
          nodes {
            ... on LabeledEvent { createdAt actor { login } label { name } }
          }
        }
        reopenedEvents: timelineItems(last: 10, itemTypes: [REOPENED_EVENT]) {
          nodes { ... on ReopenedEvent { createdAt actor { login } } }
        }
        commits(last: 1) { nodes { commit { committedDate } } }
        comments(last: 50) { totalCount nodes { author { login } createdAt } }
        labels(first: 50) { nodes { name } }
      }
    }
  }
}
"""

# Login the LABELED_EVENT actor reports for GITHUB_TOKEN-driven workflow runs;
# how the engine recognizes its own ``Stale`` marks (finding the warning time
# that anchors the close countdown).
_BOT_LOGIN = "github-actions"

# Repository roles that correspond to the OWNER/MEMBER/COLLABORATOR author
# associations; ``read`` is what everyone has on a public repository.
_EXEMPT_ROLES = {"admin", "maintain", "write", "triage"}


def parse_dt(value: str | None) -> datetime | None:
    if not value:
        return None
    return datetime.fromisoformat(value.replace("Z", "+00:00"))


class GitHubClient:
    def __init__(self, repo: str) -> None:
        self.repo = repo
        self.owner, self.name = repo.split("/")

    def _gh(self, args: list[str]) -> str:
        return subprocess.run(
            ["gh", *args], check=True, capture_output=True, text=True
        ).stdout

    def _graphql(self, query: str, **variables: str | int | None) -> dict[str, Any]:
        args = [
            "api",
            "graphql",
            "-f",
            f"query={query}",
            "-f",
            f"owner={self.owner}",
            "-f",
            f"name={self.name}",
        ]
        for key, value in variables.items():
            if value is None:
                continue
            if isinstance(value, bool):
                args += ["-F", f"{key}={'true' if value else 'false'}"]
            elif isinstance(value, int):
                args += ["-F", f"{key}={value}"]
            else:
                args += ["-f", f"{key}={value}"]
        return json.loads(self._gh(args))["data"]["repository"]

    # --- reads ---------------------------------------------------------------

    def linked_issues_for_pr(self, number: int) -> list[LinkedIssue]:
        nodes = self._graphql(_LINKED_ISSUES_QUERY, number=number)["pullRequest"][
            "closingIssuesReferences"
        ]["nodes"]
        return [
            {
                "number": node["number"],
                "assignees": [a["login"] for a in node["assignees"]["nodes"]],
            }
            for node in nodes
            # Closing keywords can reference issues in other repositories;
            # their numbers mean nothing here.
            if node["repository"]["nameWithOwner"].lower() == self.repo.lower()
        ]

    def iter_open_pull_requests_review_state(self) -> Iterator[OpenPRState]:
        """Yield lifecycle state for every open PR, paginating in batches of 50.

        Each item carries enough for the daily ``Awaiting Review`` reconcile and
        the PR stale engine to decide without per-PR follow-up requests: author +
        association, reviews (with reviewer and association), re-review requests,
        the engine's own ``Stale`` label-event time, reopens, the latest commit
        time, the author's PR comments, and current labels. The comment window
        is bounded, so ``comments_truncated`` flags when an absence of author
        comments must be re-checked exactly before acting on it; linked-issue
        comments are fetched separately, only for PRs the PR-side activity can't
        already rule active.
        """
        cursor: str | None = None
        while True:
            connection = self._graphql(_OPEN_PRS_QUERY, cursor=cursor)["pullRequests"]
            for node in connection["nodes"]:
                reviews: list[Review] = [
                    {
                        "author": (r.get("author") or {}).get("login"),
                        "state": r["state"],
                        "submitted_at": parse_dt(r["submittedAt"]),
                        "author_association": r["authorAssociation"],
                    }
                    for r in node["reviews"]["nodes"]
                ]
                review_requests: list[ReviewRequest] = [
                    {
                        "actor": (r.get("actor") or {}).get("login"),
                        "requested_at": parse_dt(r.get("createdAt")),
                    }
                    for r in node["reviewRequests"]["nodes"]
                ]
                stale_marks = [
                    marked
                    for event in node["labelEvents"]["nodes"]
                    if event["label"]["name"] == STALE_LABEL
                    and (event.get("actor") or {}).get("login") == _BOT_LOGIN
                    and (marked := parse_dt(event["createdAt"])) is not None
                ]
                reopened_events: list[Comment] = [
                    {
                        "author": (event.get("actor") or {}).get("login"),
                        "created_at": parse_dt(event["createdAt"]),
                    }
                    for event in node["reopenedEvents"]["nodes"]
                ]
                comments: list[Comment] = [
                    {
                        "author": (c.get("author") or {}).get("login"),
                        "created_at": parse_dt(c["createdAt"]),
                    }
                    for c in node["comments"]["nodes"]
                ]
                commit_nodes = node["commits"]["nodes"]
                last_commit_at = (
                    parse_dt(commit_nodes[0]["commit"]["committedDate"])
                    if commit_nodes
                    else None
                )
                yield {
                    "number": node["number"],
                    "is_draft": node["isDraft"],
                    "author": (node.get("author") or {}).get("login"),
                    "author_association": node.get("authorAssociation"),
                    "reviews": reviews,
                    "review_requests": review_requests,
                    "last_commit_at": last_commit_at,
                    "comments": comments,
                    "comments_truncated": node["comments"]["totalCount"]
                    > len(comments),
                    "reopened_events": reopened_events,
                    "labels": [label["name"] for label in node["labels"]["nodes"]],
                    "stale_marked_at": max(stale_marks) if stale_marks else None,
                }
            page = connection["pageInfo"]
            if not page["hasNextPage"]:
                break
            cursor = page["endCursor"]

    def list_open_assigned_issue_numbers(self) -> list[int]:
        # The REST list endpoint rather than a search: search results are
        # hard-capped at 1000, which pagination would eventually trip over.
        out = self._gh(
            [
                "api",
                "--paginate",
                "--slurp",
                f"repos/{self.repo}/issues?state=open&assignee=*&per_page=100",
            ]
        )
        return [
            item["number"]
            for page in json.loads(out)
            for item in page
            # This endpoint returns pull requests too; the sweep is issues-only.
            if "pull_request" not in item
        ]

    def issue_activity(self, number: int) -> IssueActivity:
        issue = self._graphql(_ISSUE_ACTIVITY_QUERY, number=number)["issue"]
        assignees = [a["login"] for a in issue["assignees"]["nodes"]]
        assigned_times = [
            assigned
            for event in issue["assignedEvents"]["nodes"]
            if (assigned := parse_dt(event["createdAt"])) is not None
        ]
        comments: list[Comment] = [
            {
                "author": (c["author"] or {}).get("login"),
                "created_at": parse_dt(c["createdAt"]),
            }
            for c in issue["comments"]["nodes"]
        ]
        crossrefs = issue["crossrefs"]
        open_pr_authors: list[str] = []
        for node in crossrefs["nodes"]:
            source = node.get("source") or {}
            if source.get("state") == "OPEN":
                author = (source.get("author") or {}).get("login")
                if author:
                    open_pr_authors.append(author)
        return {
            "assignees": assignees,
            "assigned_at": max(assigned_times) if assigned_times else None,
            "comments": comments,
            "comments_truncated": issue["comments"]["totalCount"] > len(comments),
            "open_pr_authors": open_pr_authors,
            "crossrefs_truncated": crossrefs["totalCount"] > len(crossrefs["nodes"]),
        }

    def latest_comment_at_since(
        self, number: int, logins: set[str], since: datetime
    ) -> datetime | None:
        """Exact newest comment by one of ``logins``, via paginated REST.

        The fallback behind the bounded comment windows above: called before a
        destructive decision when ``comments_truncated`` says the batched read
        can't prove the absence of activity. ``since`` bounds the fetch to the
        window that matters (a comment *created* after ``since`` is necessarily
        also updated after it, so the filter can't miss one).
        """
        since_arg = since.astimezone(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
        out = self._gh(
            [
                "api",
                "--paginate",
                "--slurp",
                f"repos/{self.repo}/issues/{number}/comments"
                f"?since={since_arg}&per_page=100",
            ]
        )
        times = [
            created
            for page in json.loads(out)
            for c in page
            if (c.get("user") or {}).get("login") in logins
            and (created := parse_dt(c.get("created_at"))) is not None
        ]
        return max(times) if times else None

    def open_pr_authors_exact(self, number: int) -> list[str]:
        """Exact authors of open PRs cross-referencing an issue, via REST.

        The fallback behind the bounded cross-reference window: an assignee's
        open PR scrolling out of the batched read must not get them unassigned.
        """
        out = self._gh(
            [
                "api",
                "--paginate",
                "--slurp",
                f"repos/{self.repo}/issues/{number}/timeline?per_page=100",
            ]
        )
        authors = []
        for page in json.loads(out):
            for event in page:
                if event.get("event") != "cross-referenced":
                    continue
                source = (event.get("source") or {}).get("issue") or {}
                if "pull_request" in source and source.get("state") == "open":
                    login = (source.get("user") or {}).get("login")
                    if login:
                        authors.append(login)
        return authors

    def is_exempt_collaborator(self, username: str) -> bool:
        """Whether ``username`` has triage-or-higher access to the repository.

        The sweep-side equivalent of the OWNER/MEMBER/COLLABORATOR author
        associations, for users where no association is available (issue
        assignees). ``role_name`` rather than the legacy ``permission`` field,
        which collapses to ``read`` for everyone on a public repository.
        """
        try:
            out = self._gh(
                ["api", f"repos/{self.repo}/collaborators/{quote(username)}/permission"]
            )
        except subprocess.CalledProcessError as err:
            if "404" in (err.stderr or ""):
                return False
            raise
        return json.loads(out).get("role_name") in _EXEMPT_ROLES

    # --- writes --------------------------------------------------------------

    def add_labels(self, number: int, labels: list[str]) -> None:
        args = ["api", "--method", "POST", f"repos/{self.repo}/issues/{number}/labels"]
        for label in labels:
            args += ["-f", f"labels[]={label}"]
        self._gh(args)

    def remove_label(self, number: int, label: str) -> None:
        path = f"repos/{self.repo}/issues/{number}/labels/{quote(label)}"
        try:
            self._gh(["api", "--method", "DELETE", path])
        except subprocess.CalledProcessError as err:
            if "404" not in (err.stderr or ""):
                raise

    def comment(self, number: int, body: str) -> None:
        self._gh(
            [
                "api",
                "--method",
                "POST",
                f"repos/{self.repo}/issues/{number}/comments",
                "-f",
                f"body={body}",
            ]
        )

    def remove_assignees(self, number: int, assignees: list[str]) -> None:
        args = [
            "api",
            "--method",
            "DELETE",
            f"repos/{self.repo}/issues/{number}/assignees",
        ]
        for assignee in assignees:
            args += ["-f", f"assignees[]={assignee}"]
        self._gh(args)

    def close_pull_request(self, number: int) -> None:
        self._gh(
            [
                "api",
                "--method",
                "PATCH",
                f"repos/{self.repo}/pulls/{number}",
                "-f",
                "state=closed",
            ]
        )
