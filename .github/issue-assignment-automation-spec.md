# Issue-Assignment & PR Lifecycle Automation — Spec

## Goal

Reduce the poor contributor experience caused by premature and competing PRs:
gate PRs on issue assignment, give contributors a self-service `/take` flow, and
clean up abandoned work — **without ever penalizing a contributor for the team's
review latency**.

---

## High-level summary of features

1. **Self-service issue assignment** — contributors claim an issue by commenting
   `/take`, and release it with `/untake`. No maintainer action required.
2. **PR assignment gate** — a PR linked to an issue whose author isn't an
   assignee is warned (and, once enabled, closed) with clear recovery steps.
   PRs with no linked issue are left alone. Maintainers and collaborators are
   exempt.
3. **`Awaiting Review` auto-labeling** — while the ball is in the maintainers'
   court, the PR is labeled `Awaiting Review` and made exempt from going stale,
   so a contributor is never penalized for slow review.
4. **Automatic unassignment of inactive claims** — an issue claimed but left
   inactive for 14 days (no PR, no assignee comment) is automatically freed, with
   a comment inviting anyone to reclaim it.
5. **Tightened stale-PR lifecycle** — a PR with the ball in the contributor's
   court is labeled `Stale` after 14 days of inactivity and auto-closed 7 days
   later; the linked issue is freed when that stale PR closes.

The guiding principle throughout: **a contributor's claim — and their PR — should
only ever lapse due to the contributor's own inactivity, never the team's.**

---

## The contributor flow

1. **Find an available issue.** Available = not labeled `Needs Triage` or
   `Needs Discussion`, and not already assigned to someone else.
2. **Claim it:** comment `/take`. The bot assigns you and reacts 👍.
   - If the issue is still `Needs Triage` / `Needs Discussion`, the bot declines
     and explains it isn't ready yet.
   - If someone else already holds it, the bot declines and points you to the
     takeover policy (you may take over after 14 days of their inactivity).
3. **Open a PR** that links the issue (e.g. `closes #1234`). Because you're an
   assignee, the assignment gate passes silently.
   - If you open a PR **without** being assigned, the bot labels it
     `Needs Issue Assignment` and comments (and, once close-mode is enabled,
     closes it) with the recovery steps: *comment `/take` on the issue to get
     assigned, then reopen this PR.* You can reopen your own PR yourself.
4. **While waiting on review,** your PR carries `Awaiting Review` and is exempt
   from going stale. You keep the issue no matter how long review takes.
5. **If changes are requested and you go quiet for 14 days,** the PR is labeled
   `Stale`, then auto-closed 7 days later. The linked issue is freed when the PR
   closes — so it and the PR release together.
6. **Change your mind early?** Comment `/untake` to release the issue immediately
   so someone else can pick it up.
7. **Claimed an issue but never opened a PR and went quiet for 14 days?** The
   issue is automatically unassigned, with a comment inviting anyone to `/take`
   it again.

---

## Labels

| Label | Status | Owner | Purpose |
|---|---|---|---|
| `Needs Issue Assignment` | **create** | bot | PR linked to an issue whose author isn't an assignee |
| `Awaiting Review` | **create** | bot | PR where the ball is with maintainers; exempts it from stale |
| `Needs Triage` | exists | human | blocks `/take` |
| `Needs Discussion` | exists | human | blocks `/take` |
| `Stale`, `Blocked` | exist | — | unchanged |

---

## Component 1 — `/take` and `/untake` (github-script, *not* Python)

In `comment-commands.yml`; latency-sensitive and trivial logic. Issue comments
only (`!issue.pull_request`).

**`/take`:**
- Issue has `Needs Triage` or `Needs Discussion` → refuse with a comment
  ("not available until triaged"); no assignment.
- Issue unassigned → assign commenter, react 👍.
- Issue already assigned to someone else → refuse, point to the 14-day
  takeover policy (`TAKEOVER_DAYS`).

**`/untake`:**
- Remove the commenter from assignees, react 👍. (Maintainers unassign others via
  the GitHub UI — not via command.)

Multi-contributor collaboration: no special automation. A maintainer manually
adds co-assignees; the gate checks "author is *among* assignees."

---

## Component 2 — PR assignment gate (Python, `scripts/`)

- **Trigger:** `pull_request_target`, types `[opened, reopened]`.
- **Exemptions (skip entirely):** `author_association ∈ {OWNER, MEMBER,
  COLLABORATOR}` (covers all org/team members incl. pandas-core/pandas-triage,
  plus hand-added collaborators) **and** bots (`user.type == 'Bot'`).
- **Linked issues:** resolved via GraphQL `closingIssuesReferences` (authoritative
  Development links).
- **Decisions:**
  - No linked issue → **do nothing** (rule removed entirely).
  - Author is an assignee of ≥1 linked issue → pass (clear `Needs Issue
    Assignment` if present).
  - Author is *not* an assignee of any linked issue → **warn**: add `Needs Issue
    Assignment` + comment. **Close only if `CLOSE_ENABLED == true`** (default
    `false`).
- This single rule subsumes the Needs-Triage/Discussion case (those issues can't
  be taken → no assignee → flagged).
- **Close comment states recovery order explicitly:** "comment `/take` on issue
  #N to get assigned, *then* reopen this PR." Gate re-fires on `reopened`, so the
  author's self-service reopen works.

---

## Component 3 — `Awaiting Review` label (Python, `scripts/`)

- **`awaiting_contributor`** ≝ latest `CHANGES_REQUESTED` review is newer than the
  latest commit (author hasn't pushed since changes were requested). No
  changes-requested, or pushed-after → *not* awaiting_contributor.
- **Apply** `Awaiting Review` when the PR is open, non-draft, and **not**
  `awaiting_contributor`. **Remove** otherwise (awaiting_contributor, draft,
  closed).
- **Trigger:** a **daily reconciliation** over all open PRs, folded into the
  scheduled maintenance job (Component 4) so it shares a single runner boot
  rather than firing on every push. The batched GraphQL read also returns
  current labels, so PRs already in the right state cost no write. Up to ~24h of
  label lag is harmless against the 14-day stale window.
- Wired into `stale-pr.yml` `exempt-pr-labels` → an awaiting-review PR **never
  goes Stale**.

---

## Component 4 — Issue unassign / inactivity (Python, `scripts/`)

- **Tunable:** `STALE_ASSIGNEE_DAYS = 14` (single top-level constant).
- **Issue stays claimed** if an assignee has **any open linked PR** (Stale or not)
  **or** an assignee commented within 14 days (assignee-scoped — third-party
  chatter doesn't count).
- **Issue is freed when:*\
  1. A `Stale`-labeled PR by the assignee is **closed** → free immediately (on the
     close event). `Stale`-at-close is a reliable proxy for "contributor's court,
     abandoned," because awaiting-review PRs are exempt from Stale.
  2. Assignee has **no open PR** and no comment within 14 days → scheduled job
     unassigns + comments "freed up — comment `/take` to reclaim."
  3. `/untake`.
- PR closed **without** `Stale` → no auto-free (deliberate close); the 14-day timer
  handles it only if they then vanish.
- **Triggers:** daily `schedule` + `pull_request_target` `[closed]`.

---

## Component 5 — `stale-pr.yml` changes

- `days-before-pr-stale: 14` (down from 30) → applies `Stale`.
- `days-before-close: 7` → auto-close 7 days after `Stale` (was `-1`/never).
- **`remove-stale-when-updated: true`** (was `false`) — mandatory, so a push
  clears `Stale` and an active contributor can't be auto-closed mid-fix.
- `exempt-pr-labels`: replace the **dead `Needs Review`** with **`Awaiting
  Review`**; keep `Blocked`, `Needs Discussion`. (`Needs Issue Assignment` is
  *not* exempt — gated PRs may go stale.)

---

## Cross-cutting: architecture & conventions

- **github-script** for `/take`/`/untake`; **Python in `scripts/`** for gate,
  Awaiting Review, unassign — split into an API-client layer + **pure decision
  functions**.
- **Unit tests** in `scripts/tests/` covering the pure logic
  (awaiting_contributor, activity rules, gate decision, exemptions) with mocked
  API data.
- **Thin YAML:** checkout (base repo only) → `setup-python` (pip-cached) →
  `python scripts/…`, event data passed via env.
- **GitHub API:** minimal deps — `gh` CLI (pre-authed) or a thin `requests`
  wrapper for REST + the one GraphQL call. No PyGithub.
- All actions **SHA-pinned**; `permissions: {}` top-level + granular per-job;
  `if: github.repository_owner == 'pandas-dev'`; `ubuntu-24.04`.
  `pull_request_target` is safe here — we only read metadata and call the API,
  never check out or execute PR code.

---

## Files

| File | Action |
|---|---|
| `.github/workflows/comment-commands.yml` | edit — add `/take`, `/untake` |
| `.github/workflows/pr-issue-gate.yml` | new |
| `.github/workflows/unassign-inactive.yml` | new — daily sweep + `Awaiting Review` reconcile + PR-close unassign |
| `.github/workflows/stale-pr.yml` | edit — thresholds, `remove-stale-when-updated`, exempt-label swap |
| `scripts/issue_assignment/…` | new — Python logic |
| `scripts/tests/test_issue_assignment*.py` | new — unit tests |
| `doc/source/development/contributing.rst` | edit — document new flow |

---

## Docs (`contributing.rst`, ~lines 51–64)

"Leave a comment with your intention" must change: a comment no longer grants a
claim — **only `/take` does**. Document `/take`/`/untake`, the assignment gate,
Needs-Triage/Discussion blocking, and the 14-day auto-unassign.

---

## Rollout sequencing (warm-up, not phases)

1. Ship everything with gate `CLOSE_ENABLED = false` (warn-only).
2. Swap `exempt-pr-labels` to `Awaiting Review` **before** enabling stale
   auto-close (else nothing protects awaiting-review PRs).
3. Let `Awaiting Review` labeling run and prove correct, **then** flip stale
   auto-close on.
4. Then flip gate `CLOSE_ENABLED = true`.
5. Migration: warn-only (+ optional "PRs opened after cutoff") spares in-flight
   PRs and old-style intent-commenters from retroactive closes.

---

## Tunable parameters (single constants)

`STALE_ASSIGNEE_DAYS = 14` · `days-before-pr-stale = 14` · `days-before-close = 7`
· `TAKEOVER_DAYS = 14` · gate `CLOSE_ENABLED = false`.
