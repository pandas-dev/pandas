"""Thin GitHub API layer built on the pre-authenticated ``gh`` CLI.

This is the only module that performs I/O, so it is deliberately small; all
judgement lives in :mod:`core`. ``gh`` reads ``GH_TOKEN`` / ``GITHUB_TOKEN``
from the environment.
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

# Used to determine whether it was the bot that marked a PR as stale, rather than
# manually done.
_BOT_LOGIN = "github-actions"

# Repository roles that correspond to the OWNER/MEMBER/COLLABORATOR author
# associations.
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
            # Ignore linked issues to other repositories.
            if node["repository"]["nameWithOwner"].lower() == self.repo.lower()
        ]

    def iter_open_pull_requests_review_state(self) -> Iterator[OpenPRState]:
        """Yield an ``OpenPRState`` for each open PR, paginating in batches.

        Each yielded state carries enough for the daily jobs to decide without
        per-PR follow-up requests. The comment window is bounded, so
        ``comments_truncated`` means the batch cannot prove an absence of
        author comments — re-check with the more thorough ``latest_comment_at_since``
        before acting on one.
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
                # Use only the newest Stale label application, and only if
                # the bot applied it.
                stale_events = [
                    (marked, (event.get("actor") or {}).get("login"))
                    for event in node["labelEvents"]["nodes"]
                    if event["label"]["name"] == STALE_LABEL
                    and (marked := parse_dt(event["createdAt"])) is not None
                ]
                newest_stale = max(stale_events, default=None)
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
                    "stale_marked_at": (
                        newest_stale[0]
                        if newest_stale is not None and newest_stale[1] == _BOT_LOGIN
                        else None
                    ),
                }
            page = connection["pageInfo"]
            if not page["hasNextPage"]:
                break
            cursor = page["endCursor"]

    def list_open_assigned_issue_numbers(self) -> list[int]:
        # Purposely avoiding search: search results are hard-capped at 1000.
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
            # This endpoint returns PRs too.
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
        """Exact newest comment by one of ``logins``.

        Used when a batched read sets ``comments_truncated``. The REST
        ``since`` parameter filters on update time, so it may return comments
        created before ``since`` but can never miss one created after it.
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
        """Exact authors of open PRs cross-referencing an issue.

        Used when a batched read sets ``crossrefs_truncated``.
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

        Used where no GitHub provided author association is available (issue
        assignees). Reads ``role_name``; the legacy ``permission`` field
        reports ``read`` for everyone on a public repository.
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
