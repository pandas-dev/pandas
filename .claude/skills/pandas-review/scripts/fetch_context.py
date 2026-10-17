"""
Gather what pandas-review reads from git and GitHub into a directory of files.

Usage, from the root of a pandas clone::

    python .claude/skills/pandas-review/scripts/fetch_context.py [TARGET] [--out DIR]

TARGET is a PR number or URL, a branch name (``owner:branch`` for a fork's PR), or
omitted to review the current branch. Prints a summary naming the files written.
Read-only: changes nothing in the clone, its git state, or on GitHub.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import re
import subprocess
import sys
import tempfile
import time
from urllib.parse import quote

REPO = "pandas-dev/pandas"
MAX_HEAD_FILES = 100
MAX_CI_LOGS = 6
# lines kept before each "##[error]" marker in a failed job's log
LOG_WINDOW = 100
# color codes, raw or as the literal "^[[" that gh writes
ANSI = re.compile(r"(?:\x1b|\^\[)\[[0-9;]*[A-Za-z]")
# `gh run view --log-failed` prefixes each line with "job<TAB>step<TAB>timestamp "
LOG_PREFIX = re.compile(r"^[^\t]*\t[^\t]*\t\d{4}-\d\d-\d\dT[\d:.]+Z ?")
NOT_FOUND = ("(HTTP 404)", "Could not resolve", "no pull requests found")
FAILED = {"FAILURE", "TIMED_OUT", "STARTUP_FAILURE", "ERROR"}


class FetchError(Exception):
    pass


def run(*args: str, retries: int = 3) -> str:
    for attempt in range(retries):
        proc = subprocess.run(args, capture_output=True, check=False)
        if proc.returncode == 0:
            return proc.stdout.decode("utf-8", "replace")
        err = proc.stderr.decode("utf-8", "replace").strip()
        # retrying won't help a missing PR, file, or branch
        if any(msg in err for msg in NOT_FOUND):
            break
        if attempt < retries - 1:
            time.sleep(2 * (attempt + 1))
    raise FetchError(f"`{' '.join(args)}` failed: {err}")


def git(*args: str) -> str:
    return run("git", *args, retries=1).strip()


def gh_json(*args: str):
    return json.loads(run("gh", *args))


def gh_paged(endpoint: str) -> list[dict]:
    pages = gh_json(
        "api", "-X", "GET", "--paginate", "--slurp", f"repos/{REPO}/{endpoint}"
    )
    return [item for page in pages for item in page]


def upstream_remote() -> str | None:
    for line in git("remote", "-v").splitlines():
        name, url = line.split()[:2]
        if re.search(r"[/:]pandas-dev/pandas(\.git)?/?$", url):
            return name
    return None


def resolve_pr(target: str, branch: str) -> tuple[int | None, bool]:
    """Return the PR number (None if there is none) and whether this is local mode."""
    if not target or target == branch:
        if not branch:
            # detached HEAD
            return None, True
        hits = gh_json(
            "pr",
            "list",
            "--repo",
            REPO,
            "--head",
            branch,
            "--author",
            "@me",
            "--json",
            "number",
        )
        # --head matches the branch name in any fork
        return (hits[0]["number"] if len(hits) == 1 else None), True
    number = target.lstrip("#")
    if number.isdigit() or target.startswith("http") or ":" in target:
        meta = gh_json(
            "pr",
            "view",
            number if number.isdigit() else target,
            "--repo",
            REPO,
            "--json",
            "number",
        )
        return meta["number"], False
    hits = gh_json("pr", "list", "--repo", REPO, "--head", target, "--json", "number")
    if len(hits) != 1:
        raise FetchError(
            f"{len(hits)} open PRs have head branch {target!r}; "
            "pass the PR number or owner:branch instead"
        )
    return hits[0]["number"], False


def who(item: dict, pr_author: str) -> str:
    login = (item.get("user") or {}).get("login", "ghost")
    tags = []
    if login == pr_author:
        tags.append("PR author")
    assoc = item.get("author_association", "NONE")
    if assoc != "NONE":
        tags.append(assoc.lower())
    return f"@{login}" + (f" ({', '.join(tags)})" if tags else "")


def render_thread(pr: int, pr_author: str) -> tuple[list[str], dict[str, int]]:
    reviews = gh_paged(f"pulls/{pr}/reviews")
    inline = gh_paged(f"pulls/{pr}/comments")
    convo = gh_paged(f"issues/{pr}/comments")
    lines = ["## Reviews", ""]
    for rev in reviews:
        # an empty COMMENTED review only wraps inline comments, listed below
        if rev["state"] == "COMMENTED" and not (rev.get("body") or "").strip():
            continue
        lines += [
            f"### {who(rev, pr_author)} {rev['state']} {rev.get('submitted_at')}",
            "",
            rev.get("body") or "",
            "",
        ]
    lines += ["## Inline comments", ""]
    threads: dict[int, list[dict]] = {}
    for com in inline:
        threads.setdefault(com.get("in_reply_to_id") or com["id"], []).append(com)
    for thread in threads.values():
        root = thread[0]
        line = root.get("line")
        where = (
            f"{root['path']}:{line}"
            if line
            else (f"{root['path']}:{root.get('original_line')} (outdated)")
        )
        lines += [f"### {where}", ""]
        for com in sorted(thread, key=lambda c: c["created_at"]):
            lines += [
                f"**{who(com, pr_author)}** {com['created_at']}:",
                "",
                com["body"],
                "",
            ]
    lines += ["## Conversation", ""]
    for com in convo:
        lines += [f"### {who(com, pr_author)} {com['created_at']}", "", com["body"], ""]
    counts = {
        "reviews": len(reviews),
        "inline": len(inline),
        "conversation": len(convo),
    }
    return lines, counts


def render_issue(number: int, pr_author: str) -> list[str]:
    issue = gh_json("api", "-X", "GET", f"repos/{REPO}/issues/{number}")
    comments = gh_paged(f"issues/{number}/comments")
    reactions = gh_paged(f"issues/{number}/reactions")
    lines = [
        f"## Linked issue #{number}: {issue['title']} ({issue['state']})",
        "",
        f"Opened by {who(issue, pr_author)}",
        "",
        issue.get("body") or "",
        "",
    ]
    reacts = ", ".join(f"@{r['user']['login']} {r['content']}" for r in reactions)
    lines += [f"Reactions on the issue: {reacts or 'none'}", ""]
    for com in comments:
        lines += [f"### {who(com, pr_author)} {com['created_at']}", "", com["body"], ""]
    return lines


def log_excerpt(log: str) -> str:
    rows = [ANSI.sub("", LOG_PREFIX.sub("", row)) for row in log.splitlines()]
    keep: set[int] = set()
    for i, row in enumerate(rows):
        if "##[error]" in row:
            keep.update(range(max(0, i - LOG_WINDOW), i + 1))
    if not keep:
        keep = set(range(max(0, len(rows) - LOG_WINDOW), len(rows)))
    out, prev = [], -1
    for i in sorted(keep):
        if i != prev + 1:
            out.append("[...]")
        out.append(rows[i])
        prev = i
    return "\n".join(out)


def render_ci(checks: list[dict], out: Path) -> tuple[list[str], list[str]]:
    failed = [
        chk
        for chk in checks
        if chk.get("conclusion") in FAILED or chk.get("state") in FAILED
    ]
    names = [chk.get("name") or chk.get("context") or "?" for chk in failed]
    lines = ["## Failing CI checks", ""]
    if not failed:
        return [*lines, "None.", ""], names
    for n, chk in enumerate(failed):
        url = chk.get("detailsUrl") or chk.get("targetUrl") or ""
        lines += [f"### {names[n]}", "", url, ""]
        job = re.search(r"/job/(\d+)", url)
        if job is None or n >= MAX_CI_LOGS:
            lines += ["(log not fetched)", ""]
            continue
        try:
            log = run(
                "gh", "run", "view", "--repo", REPO, "--job", job[1], "--log-failed"
            )
        except FetchError as err:
            lines += [f"(log unavailable: {err})", ""]
            continue
        path = out / "ci" / f"{job[1]}.log"
        path.parent.mkdir(exist_ok=True)
        path.write_text(ANSI.sub("", log), encoding="utf-8")
        lines += [f"Full log: {path}", "", "```", log_excerpt(log), "```", ""]
    return lines, names


def fetch_head_files(meta: dict, out: Path) -> list[str]:
    notes = []
    for n, file in enumerate(meta["files"]):
        if n >= MAX_HEAD_FILES:
            notes.append(f"only the first {MAX_HEAD_FILES} changed files were saved")
            break
        endpoint = (
            f"repos/{REPO}/contents/{quote(file['path'])}?ref={meta['headRefOid']}"
        )
        try:
            text = run(
                "gh",
                "api",
                "-X",
                "GET",
                endpoint,
                "-H",
                "Accept: application/vnd.github.raw",
            )
        except FetchError as err:
            if "(HTTP 404)" not in str(err):
                raise
            notes.append(f"{file['path']} is not at the PR head (deleted)")
            continue
        dest = out / "head" / file["path"]
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_text(text, encoding="utf-8")
    return notes


def collect(target: str, out: Path) -> str:
    remote = upstream_remote()
    branch = git("branch", "--show-current")
    pr, local = resolve_pr(target, branch)
    summary, context = [], []

    meta = None
    if pr is not None:
        meta = gh_json(
            "pr",
            "view",
            str(pr),
            "--repo",
            REPO,
            "--json",
            "number,url,title,author,state,isDraft,headRefOid,files,body,"
            "closingIssuesReferences,statusCheckRollup",
        )
        context += [
            f"# PR #{pr}: {meta['title']}",
            "",
            f"{meta['url']} by @{meta['author']['login']}, {meta['state']}"
            f"{', draft' if meta['isDraft'] else ''}, head {meta['headRefOid']}",
            "",
            "## Description",
            "",
            meta["body"] or "",
            "",
        ]
        summary.append(f"PR: #{pr} {meta['title']} ({meta['url']})")

    head = git("rev-parse", "HEAD")
    if local:
        if remote is None:
            raise FetchError(
                "no git remote points at pandas-dev/pandas; add one with "
                "`git remote add upstream https://github.com/pandas-dev/pandas`"
            )
        base = git("merge-base", "HEAD", f"{remote}/main")
        diff = run("git", "diff", base, retries=1)
        untracked = [
            row[3:]
            for row in run(
                "git", "status", "--porcelain", "-uall", retries=1
            ).splitlines()
            if row.startswith("?? ")
        ]
        summary.append(
            f"Mode: local, branch {branch or '(detached HEAD)'}"
            + ("" if pr is not None else ", no PR yet")
        )
        summary.append("Code under review: the working tree")
        if untracked:
            summary.append(
                "Untracked files (not in diff.patch; Read them): "
                + ", ".join(untracked)
            )
        if meta is not None and meta["headRefOid"] != head:
            summary.append(
                f"The PR head on GitHub is {meta['headRefOid']}, not local HEAD; "
                "thread and CI refer to the pushed head"
            )
    else:
        assert meta is not None
        diff = run("gh", "pr", "diff", str(pr), "--repo", REPO)
        summary.append("Mode: PR")
        if meta["headRefOid"] == head:
            summary.append(
                "Code under review: the local checkout (HEAD is the PR head)"
            )
            local = True
        else:
            notes = fetch_head_files(meta, out)
            summary.append(
                f"Code under review: the PR head on GitHub; changed files saved "
                f"under {out / 'head'}"
            )
            summary += notes
    (out / "diff.patch").write_text(diff, encoding="utf-8")
    summary.append(f"Diff: {out / 'diff.patch'} ({len(diff.splitlines())} lines)")

    if local and remote is not None:
        behind = int(git("rev-list", "--count", f"HEAD..{remote}/main"))
        if behind:
            summary.append(
                f"Behind {remote}/main by {behind} commits (as of the last fetch)"
            )

    if meta is not None:
        author = meta["author"]["login"]
        thread, counts = render_thread(pr, author)
        context += thread
        summary.append("Thread: " + ", ".join(f"{v} {k}" for k, v in counts.items()))
        issues = [
            ref["number"]
            for ref in meta["closingIssuesReferences"]
            if ref["repository"]["owner"]["login"] == "pandas-dev"
            and ref["repository"]["name"] == "pandas"
        ]
        for number in issues:
            context += render_issue(number, author)
        summary.append(
            "Linked issues: " + (", ".join(f"#{n}" for n in issues) or "none")
        )
        ci, failing = render_ci(meta["statusCheckRollup"] or [], out)
        context += ci
        summary.append("Failing CI checks: " + (", ".join(failing) or "none"))
        (out / "context.md").write_text("\n".join(context), encoding="utf-8")
        summary.append(f"Description, thread, linked issues, CI: {out / 'context.md'}")
    return "\n".join(summary)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    parser.add_argument("target", nargs="?", default="")
    parser.add_argument(
        "--out", help="directory to write into (default: a new temp dir)"
    )
    args = parser.parse_args(argv)
    out = (
        Path(args.out) if args.out else Path(tempfile.mkdtemp(prefix="pandas-review-"))
    )
    out.mkdir(parents=True, exist_ok=True)
    try:
        print(collect(args.target.strip(), out))
    except FetchError as err:
        print(f"ERROR: {err}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
