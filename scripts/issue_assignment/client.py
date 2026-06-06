"""Thin GitHub API layer built on the pre-authenticated ``gh`` CLI.

This is the only module that performs I/O, so it is deliberately small; all
judgement lives in :mod:`core`. ``gh`` reads ``GH_TOKEN`` / ``GITHUB_TOKEN``
from the environment, so no explicit auth wiring is needed.
"""

from __future__ import annotations

from datetime import datetime
import json
import subprocess
from urllib.parse import quote

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
        reviews(last: 20) { nodes { state submittedAt } }
        commits(last: 1) { nodes { commit { committedDate } } }
        labels(first: 50) { nodes { name } }
      }
    }
  }
}
"""


def parse_dt(value):
    if not value:
        return None
    return datetime.fromisoformat(value.replace("Z", "+00:00"))


class GitHubClient:
    def __init__(self, repo):
        self.repo = repo
        self.owner, self.name = repo.split("/")

    def _gh(self, args):
        return subprocess.run(
            ["gh", *args], check=True, capture_output=True, text=True
        ).stdout

    def _graphql(self, query, **variables):
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

    def linked_issues_for_pr(self, number):
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

    def iter_open_pull_requests_review_state(self):
        """Yield review state for every open PR, paginating in batches of 50.

        Each item has ``number``, ``is_draft``, ``reviews``, ``last_commit_at``,
        and current ``labels`` — enough for the daily ``Awaiting Review``
        reconciliation to decide without any per-PR follow-up requests.
        """
        cursor = None
        while True:
            connection = self._graphql(_OPEN_PRS_QUERY, cursor=cursor)["pullRequests"]
            for node in connection["nodes"]:
                reviews = [
                    {"state": r["state"], "submitted_at": parse_dt(r["submittedAt"])}
                    for r in node["reviews"]["nodes"]
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
                    "reviews": reviews,
                    "last_commit_at": last_commit_at,
                    "labels": [label["name"] for label in node["labels"]["nodes"]],
                }
            page = connection["pageInfo"]
            if not page["hasNextPage"]:
                break
            cursor = page["endCursor"]

    def list_open_assigned_issue_numbers(self):
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
                "500",
            ]
        )
        return [item["number"] for item in json.loads(out)]

    def issue_activity(self, number):
        issue = self._graphql(_ISSUE_ACTIVITY_QUERY, number=number)["issue"]
        assignees = [a["login"] for a in issue["assignees"]["nodes"]]
        comments = [
            {
                "author": (c["author"] or {}).get("login"),
                "created_at": parse_dt(c["createdAt"]),
            }
            for c in issue["comments"]["nodes"]
        ]
        open_pr_authors = []
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

    def add_labels(self, number, labels):
        args = ["api", "--method", "POST", f"repos/{self.repo}/issues/{number}/labels"]
        for label in labels:
            args += ["-f", f"labels[]={label}"]
        self._gh(args)

    def remove_label(self, number, label):
        path = f"repos/{self.repo}/issues/{number}/labels/{quote(label)}"
        try:
            self._gh(["api", "--method", "DELETE", path])
        except subprocess.CalledProcessError as err:
            if "404" not in (err.stderr or ""):
                raise

    def comment(self, number, body):
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

    def remove_assignees(self, number, assignees):
        args = [
            "api",
            "--method",
            "DELETE",
            f"repos/{self.repo}/issues/{number}/assignees",
        ]
        for assignee in assignees:
            args += ["-f", f"assignees[]={assignee}"]
        self._gh(args)

    def close_pull_request(self, number):
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
