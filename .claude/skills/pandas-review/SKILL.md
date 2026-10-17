---
name: pandas-review
description: Self-review your own pandas PR (or your local branch before opening one) the way a pandas maintainer would. Prints findings and a verdict in the terminal; never posts to GitHub or edits files.
disable-model-invocation: true
argument-hint: "[PR number, URL, or branch]"
allowed-tools: Bash(python .claude/skills/pandas-review/scripts/fetch_context.py), Bash(python .claude/skills/pandas-review/scripts/fetch_context.py *), Bash(gh pr view *), Bash(gh api -X GET *), Bash(git diff *), Bash(git rev-parse *), Bash(git show *), Bash(git log *), Bash(git status *), Bash(pytest *), Bash(python -m pytest *), Read, Grep, Glob
---

Review a pandas pull request the way a pandas maintainer would, and report to its author.

Target: $ARGUMENTS

The target is a PR number or URL, a branch name, or nothing (including when the line above
reads literally `$ARGUMENTS`) to review the current branch with its uncommitted changes.

**Review as a maintainer who did not write this code**, even if you did. Verify what the PR
description, commit messages, and code comments claim instead of taking them on trust.
Report only: don't edit files, and don't fix anything during the review. The author decides
what to change afterwards.

## 1. Gather the change and its thread

From the repo root, run

```
python .claude/skills/pandas-review/scripts/fetch_context.py <target>
```

with no argument for the current branch, and `owner:branch` for a branch in a fork. It
prints the PR it found (if any), where the code under review is, and the files it wrote.
Read all of them in full, never truncated: `diff.patch`, any untracked files it lists, and
`context.md` (the PR description, reviews, inline comments, conversation, linked issues with
their reactions, and failing CI jobs with log excerpts).

- **If the script fails, stop and report its error.** A failed fetch is unknown, not empty;
  never review as if the thread were silent.
- When the code under review is the PR head on GitHub, its changed files are saved under
  `head/`. Read any other file at the head with
  `gh api -X GET "repos/pandas-dev/pandas/contents/<path>?ref=<head sha>" -H "Accept: application/vnd.github.raw"`.
- If it reports the branch is behind `main`, say so. CI tests the PR merged into `main`, so
  a recent change on `main` can break it there even when it passes locally.
- Review verdicts matter most: a maintainer's approval, requested changes, or "I don't
  think this is worth doing" usually lives there, not inline.
- **The thread outranks this file.** Where a maintainer on the thread has asked for something
  that conflicts with the guidance below, follow the maintainer.

## 2. Investigate

Read `AGENTS.md` at the repo root. Consult `doc/source/development/contributing_codebase.rst`
when a question about tests, typing, or backward compatibility comes up. Use Grep/Glob/Read
on the surrounding code to check the change against how the codebase already does things.

**Run code when the local build is the code under review.** Run the PR's new and changed
tests, the existing test files covering the code it touches, and every repro you report,
and show real output. If pandas fails to import or build, say so and review from reading
alone; don't try to repair the environment.

## 3. What to look for

In priority order:

1. **Maintainer position.** If a maintainer has said the PR shouldn't proceed, the verdict
   is `CLOSE`; if they want a different approach, `WRONG APPROACH`. Skip to section 4 and
   report only that. If they've said it's blocked on a decision, the verdict is `BLOCKED`,
   but review the rest to give the verdict it would get once unblocked.
   A core maintainer's 👍 reaction on the linked issue counts as buy-in. A `/take` with no
   activity for several days doesn't block anyone.
2. **Unaddressed review comments.** For each reviewer request in the thread, check whether
   the current diff addresses it or the author replied with a reason. List any that are
   neither. Don't draft replies: pandas asks that authors speak for themselves on the
   thread (see `doc/source/development/contributing.rst`). If the author has already
   explained on the thread why one of your findings doesn't apply, say so and don't count it
   toward the verdict, unless the explanation is factually wrong (then say what is wrong with
   it).
3. **Correctness.** Bugs in the diff, and above all **regressions**: inputs that worked on
   `main` and break with this PR.
   - A fix that covers some code paths but leaves sibling paths unfixed is **not** a defect.
     At most mention it as an optional follow-up. Only regressions and wrong results on the
     paths the PR claims to fix drive `NEEDS WORK`.
   - A failing CI job the diff causes is a correctness finding; quote the log line.
4. **Necessity, simplicity, brevity.** Flag only concrete excess you can point at:
   - code, parameters, or tests the fix doesn't need (a guard no input reaches, with the
     reason; a refactor riding along; a test duplicating an existing one);
   - a new helper duplicating existing code (name it at `file:line`), a one-caller
     abstraction that inlines with no loss, or a materially smaller equally-correct approach
     (sketch it);
   - the same rationale pasted at several sites, or comments narrating code the PR removed.
   Don't nit comment length otherwise. Equal-size alternative designs and naming taste are
   not findings.
5. **Conventions that commonly come up in review:**
   - A whatsnew entry describes the user-visible change, not the implementation, and names
     only public API. A pure refactor gets no entry.
   - Tests: module-level functions or an existing class are both fine; so is a new class
     with a reasonable scope.

Do **not** report:
- things that are fine ("tests look good", "no whatsnew needed"); silence means no concern;
- PR process or hygiene: template checkboxes, labels, merge readiness;
- CI failures the diff didn't cause (infrastructure, flaky tests, a break already on `main`);
- drift between the PR description and the diff; review the diff;
- predictions of what a pre-commit hook or CI will catch (the author should run
  `pre-commit` separately).

## 4. Output

Open with one or two sentences on what the PR does, then the findings, most important first.

**Anchor each localized finding** at `path:line` on its own line, using the line number at
the head of the code under review (not a hunk offset). Use a range when it spans lines. For
diffuse findings (wrong approach, missing coverage) use prose; don't invent an anchor.
Every finding ends with a concrete fix.

```
pandas/core/indexes/multi.py:1487
`codes` is reused after the early `continue`, so the vectorized path sees the previous
level's array. Fix: move the assignment above the branch.
```

**One repro per distinct failure mode**: runnable Python, imports included, actual next to
expected. Different inputs that fail the same way share one repro. Say in one line whether
you executed it; never present predicted output as observed.

```python
import pandas as pd

idx = pd.MultiIndex.from_product([["a", "b"], [1, 2]])
print(idx.get_locs([["a"], [1, 2]]))
# actual:   [0]
# expected: [0, 1]
```

Keep it concise: findings, not a checklist of what you inspected.

**End with the verdict**, the last thing in the output, in a fenced box:

```
┏━━━━━━━━━━━━━━━━━━━━━━━┓
┃  VERDICT: NEEDS WORK  ┃
┗━━━━━━━━━━━━━━━━━━━━━━━┛
```

| Verdict | Means |
|---|---|
| `LGTM` | Ready for maintainer review as-is. |
| `NITS ONLY` | Nothing blocking; the findings are optional. |
| `NEEDS WORK` | Real defects; the approach is sound. |
| `WRONG APPROACH` | Can't be patched into correctness as written. |
| `BLOCKED` | Waiting on a maintainer decision or another PR, not on the author. Add the verdict it would get once unblocked: `BLOCKED → NITS ONLY`. |
| `CLOSE` | Unlikely to be merged in any form; say why. |

## For the author

- `NITS ONLY` clears the bar for review, but `LGTM` gets reviewed sooner. Fix the nits you
  agree with and rerun once; then stop, since further reruns surface new optional nits, not
  blockers.
- If you disagree with a finding, say why on the thread in your own words rather than
  changing code to satisfy this review.
- A green verdict doesn't replace your own review of the change, which pandas requires, and
  your PR description still needs the tool, model, and effort you used (see the automated
  contributions policy in `contributing.rst`).

## Hard rules

- Everything stays in the terminal. Never post, comment, review, react, or push on GitHub.
- Don't modify files or git state.
