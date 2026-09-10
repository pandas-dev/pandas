#!/bin/bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage: scripts/push_tag_for_release.sh <version> <branch> [options]

Create the release commit and tag for a new pandas version and push them to the
remote, which triggers the wheels/PyPI, docs and GitHub release automation.

Arguments:
  <version>              Release version without the leading "v" (e.g. 1.5.2, 1.4.0rc0)
  <branch>               Branch to release from (e.g. 3.0.x, main)

Options:
  --remote NAME          Remote to push to (default: upstream)
  --rc-branch X.Y.x      New maintenance branch to create (release candidate only)
  --next-version X.Y.Z   Start the next dev cycle on <branch> (release candidate only)
  --yes                  Skip the interactive confirmation prompt
  -h, --help             Show this help and exit
USAGE
}

die() {
  echo "Error: $*" >&2
  exit 1
}

VERSION=""
BRANCH=""
REMOTE="upstream"
RC_BRANCH=""
NEXT_VERSION=""
ASSUME_YES="false"

POSITIONAL=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    --remote)
      [[ $# -ge 2 ]] || die "--remote requires a value"
      REMOTE="$2"
      shift 2
      ;;
    --rc-branch)
      [[ $# -ge 2 ]] || die "--rc-branch requires a value"
      RC_BRANCH="$2"
      shift 2
      ;;
    --next-version)
      [[ $# -ge 2 ]] || die "--next-version requires a value"
      NEXT_VERSION="$2"
      shift 2
      ;;
    --yes)
      ASSUME_YES="true"
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    -*)
      die "Unknown option: $1"
      ;;
    *)
      POSITIONAL+=("$1")
      shift
      ;;
  esac
done

[[ ${#POSITIONAL[@]} -eq 2 ]] || { usage >&2; die "Expected exactly 2 positional arguments, got ${#POSITIONAL[@]}"; }
VERSION="${POSITIONAL[0]}"
BRANCH="${POSITIONAL[1]}"

version_pat='^[0-9]+\.[0-9]+\.[0-9]+(rc[0-9]+)?$'
next_version_pat='^[0-9]+\.[0-9]+\.[0-9]+$'

[[ "$VERSION" =~ $version_pat ]] || die "Invalid version '$VERSION' (expected e.g. 1.5.2 or 1.4.0rc0)"

IS_RC="false"
if [[ "$VERSION" == *rc* ]]; then
  IS_RC="true"
fi

if [[ "$IS_RC" == "true" ]]; then
  [[ -n "$RC_BRANCH" && -n "$NEXT_VERSION" ]] || die "Release candidate requires both --rc-branch and --next-version"
  [[ "$NEXT_VERSION" =~ $next_version_pat ]] || die "Invalid --next-version '$NEXT_VERSION' (expected e.g. 1.5.0)"
else
  [[ -z "$RC_BRANCH" && -z "$NEXT_VERSION" ]] || die "--rc-branch/--next-version are only valid for a release candidate"
fi

git rev-parse --is-inside-work-tree >/dev/null 2>&1 || die "Not inside a git work tree"
git remote get-url "$REMOTE" >/dev/null 2>&1 || die "Remote '$REMOTE' does not exist"

git ls-remote --exit-code --heads "$REMOTE" "$BRANCH" >/dev/null 2>&1 || die "Branch '$BRANCH' does not exist on remote '$REMOTE'"

TAG="v$VERSION"
git rev-parse -q --verify "refs/tags/$TAG" >/dev/null 2>&1 && die "Tag '$TAG' already exists locally"
if git ls-remote --exit-code --tags "$REMOTE" "$TAG" >/dev/null 2>&1; then
  die "Tag '$TAG' already exists on remote '$REMOTE'"
fi

echo "About to release pandas $VERSION from branch '$BRANCH' (remote '$REMOTE')."
echo "  - 'git clean -xdf' will DELETE all untracked files in the working tree."
echo "  - An empty 'RLS: $VERSION' commit and tag '$TAG' will be pushed to '$REMOTE/$BRANCH'."
if [[ "$IS_RC" == "true" ]]; then
  echo "  - Maintenance branch '$RC_BRANCH' will be created and pushed."
  echo "  - The next dev cycle '$NEXT_VERSION' (tag 'v$NEXT_VERSION.dev0') will be started on '$BRANCH'."
fi

if [[ "$ASSUME_YES" != "true" ]]; then
  read -r -p "Proceed? [y/N] " reply
  case "$reply" in
    y|Y|yes|Yes) ;;
    *) die "Aborted by user" ;;
  esac
fi

git checkout "$BRANCH"
git pull --ff-only "$REMOTE" "$BRANCH"
git clean -xdf
git commit --allow-empty --author="pandas Development Team <pandas-dev@python.org>" -m "RLS: $VERSION"
git tag -a "$TAG" -m "Version $VERSION"
git push "$REMOTE" "$BRANCH" --follow-tags

if [[ "$IS_RC" == "true" ]]; then
  git checkout -b "$RC_BRANCH"
  git push "$REMOTE" "$RC_BRANCH"
  git checkout "$BRANCH"
  git commit --allow-empty -m "Start $NEXT_VERSION"
  git tag -a "v$NEXT_VERSION.dev0" -m "DEV: Start $NEXT_VERSION"
  git push "$REMOTE" "$BRANCH" --follow-tags
fi

echo "Done. Pushed $TAG to $REMOTE/$BRANCH."
