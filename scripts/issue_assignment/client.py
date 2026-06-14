"""Thin GitHub API layer built on the pre-authenticated ``gh`` CLI.

This is the only module that performs I/O, so it is deliberately small; all
judgement lives in :mod:`core`. ``gh`` reads ``GH_TOKEN`` / ``GITHUB_TOKEN``
from the environment, so no explicit auth wiring is needed.
"""

from __future__ import annotations

from datetime import datetime
import json
import subprocess
from typing import (
    TYPE_CHECKING,
    Any,
)
from urllib.parse import quote

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
        nodes { number assignees(first: 20) { nodes { login } } }
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
      comments(last: 30) { nodes { author { login } createdAt } }
      timelineItems(last: 50, itemTypes: [CROSS_REFERENCED_EVENT]) {
        nodes {
          ... on CrossReferencedEvent {
            source { ... on PullRequest { state author { login } } }
          }
        }
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
        reviews(last: 20) { nodes { state submittedAt } }
        timelineItems(last: 30, itemTypes: [REVIEW_REQUESTED_EVENT]) {
          nodes { ... on ReviewRequestedEvent { createdAt actor { login } } }
        }
        labels(first: 50) { nodes { name } }
      }
    }
  }
}
"""


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
        ]

    def iter_open_pull_requests_review_state(self) -> Iterator[OpenPRState]:
        """Yield review state for every open PR, paginating in batches of 50.

        Each item has ``number``, ``is_draft``, ``author``, ``reviews``,
        ``review_requests``, and current ``labels`` — enough for the daily
        ``Awaiting Review`` reconciliation to decide without any per-PR
        follow-up requests.
        """
        cursor: str | None = None
        while True:
            connection = self._graphql(_OPEN_PRS_QUERY, cursor=cursor)["pullRequests"]
            for node in connection["nodes"]:
                reviews: list[Review] = [
                    {"state": r["state"], "submitted_at": parse_dt(r["submittedAt"])}
                    for r in node["reviews"]["nodes"]
                ]
                review_requests: list[ReviewRequest] = [
                    {
                        "actor": (r.get("actor") or {}).get("login"),
                        "requested_at": parse_dt(r.get("createdAt")),
                    }
                    for r in node["timelineItems"]["nodes"]
                ]
                yield {
                    "number": node["number"],
                    "is_draft": node["isDraft"],
                    "author": (node.get("author") or {}).get("login"),
                    "reviews": reviews,
                    "review_requests": review_requests,
                    "labels": [label["name"] for label in node["labels"]["nodes"]],
                }
            page = connection["pageInfo"]
            if not page["hasNextPage"]:
                break
            cursor = page["endCursor"]

    def list_open_assigned_issue_numbers(self) -> list[int]:
        out = self._gh(
            [
                "issue",
                "list",
                "--repo",
                self.repo,
                "--search",
                "is:open assignee:*",
                "--json",
                "number",
                "--limit",
                "1500",
            ]
        )
        return [item["number"] for item in json.loads(out)]

    def issue_activity(self, number: int) -> IssueActivity:
        issue = self._graphql(_ISSUE_ACTIVITY_QUERY, number=number)["issue"]
        assignees = [a["login"] for a in issue["assignees"]["nodes"]]
        comments: list[Comment] = [
            {
                "author": (c["author"] or {}).get("login"),
                "created_at": parse_dt(c["createdAt"]),
            }
            for c in issue["comments"]["nodes"]
        ]
        open_pr_authors: list[str] = []
        for node in issue["timelineItems"]["nodes"]:
            source = node.get("source") or {}
            if source.get("state") == "OPEN":
                author = (source.get("author") or {}).get("login")
                if author:
                    open_pr_authors.append(author)
        return {
            "assignees": assignees,
            "comments": comments,
            "open_pr_authors": open_pr_authors,
        }

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
