"""Git utilities."""

# Used also from setup.py, so don't pull in anything additional here (like mypy or typing):
from __future__ import annotations

import os
import subprocess


def is_git_repo(dir: str) -> bool:
    """Is the given directory version-controlled with git?"""
    return os.path.exists(os.path.join(dir, ".git"))


def have_git() -> bool:
    """Can we run the git executable?"""
    try:
        subprocess.check_output(["git", "--help"])
        return True
    except subprocess.CalledProcessError:
        return False
    except OSError:
        return False


def git_revision(dir: str) -> bytes:
    """Get the SHA-1 of the HEAD of a git repository."""
    return subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=dir).strip()


def git_revision_no_subprocess(dir: str) -> bytes | None:
    """Get the SHA-1 of HEAD by reading git files directly, without a subprocess.

    Returns None if the revision cannot be determined this way (e.g. unusual
    git state), in which case the caller should fall back to git_revision().
    """
    try:
        head_path = os.path.join(dir, ".git", "HEAD")
        with open(head_path, "rb") as f:
            head = f.read().strip()
        if head.startswith(b"ref: "):
            ref = head[5:]  # e.g. b"refs/heads/main"
            ref_path = os.path.join(dir, ".git", os.fsdecode(ref))
            if os.path.exists(ref_path):
                with open(ref_path, "rb") as f:
                    return f.read().strip()
            # The ref may be in packed-refs instead of a loose file.
            packed_refs_path = os.path.join(dir, ".git", "packed-refs")
            if os.path.exists(packed_refs_path):
                with open(packed_refs_path, "rb") as f:
                    for line in f:
                        if line.startswith(b"#"):
                            continue
                        parts = line.strip().split()
                        if len(parts) >= 2 and parts[1] == ref:
                            return parts[0]
            return None
        # Detached HEAD: content is the SHA itself.
        if len(head) == 40 and all(chr(c) in "0123456789abcdef" for c in head):
            return head
        return None
    except OSError:
        return None


def is_dirty(dir: str) -> bool:
    """Check whether a git repository has uncommitted changes."""
    output = subprocess.check_output(["git", "status", "-uno", "--porcelain"], cwd=dir)
    return output.strip() != b""
