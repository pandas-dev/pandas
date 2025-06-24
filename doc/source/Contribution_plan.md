# CONTRIBUTION_PLAN

## 1. Project
- pandas (📘 documentation)

## 2. Issue
- #60647: “DOC: Simplify pandas theme footer”

## 3. Why
- Beginner-friendly, docs-only, no code.
- Helps modernize the docs and maintain consistency.

## 4. Steps
1. Fork & clone the repo.
2. Create branch `docs/simplify-footer`.
3. Edit `doc/source/conf.py`, simplify/remove custom footer.
4. Commit: `git commit -m "DOC: simplify theme footer (fix #60647)"`.
5. Push to branch, open PR.

## 5. PR Expectations
- Title: `DOC: simplify theme footer (closes #60647)`
- Body: Briefly explain change and reference issue.
- CI in pandas checks: doc build only, so no tests should fail.
- Respond to maintainers' feedback.
