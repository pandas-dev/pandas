"""Integration tests for misc/diff-cache.py and misc/apply-cache-diff.py."""

from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
import tempfile
import time
import unittest

from mypy.defaults import SQLITE_NUM_SHARDS
from mypy.test.config import PREFIX

_MISC_DIR = os.path.join(PREFIX, "misc")
_DIFF_CACHE_PATH = os.path.join(_MISC_DIR, "diff-cache.py")
_APPLY_CACHE_DIFF_PATH = os.path.join(_MISC_DIR, "apply-cache-diff.py")


class DiffCacheIntegrationTests(unittest.TestCase):
    """Run mypy twice with different sources, diff the caches, and apply the diff."""

    def test_diff_cache_produces_valid_json(self) -> None:
        # Use a single source directory with two cache directories so that
        # source paths in the cache metadata are identical between runs.
        # b.py is modified and c.py is added in the second run.
        src_dir = tempfile.mkdtemp()
        output_file = os.path.join(tempfile.mkdtemp(), "diff.json")
        env = os.environ.copy()
        env["PYTHONPATH"] = PREFIX
        try:
            cache1 = os.path.join(src_dir, "cache1")
            cache2 = os.path.join(src_dir, "cache2")

            # Write sources and run mypy for cache1 (using sqlite cache)
            with open(os.path.join(src_dir, "a.py"), "w") as f:
                f.write("x: int = 1\n")
            with open(os.path.join(src_dir, "b.py"), "w") as f:
                f.write("import a\ndef foo() -> int:\n    return 1\n")
            result = subprocess.run(
                [
                    sys.executable,
                    "-m",
                    "mypy",
                    "--sqlite-cache",
                    "--cache-fine-grained",
                    "--cache-dir",
                    cache1,
                    "b.py",
                ],
                cwd=src_dir,
                capture_output=True,
                text=True,
                env=env,
            )
            assert result.returncode == 0, f"mypy run 1 failed: {result.stderr}"

            # Sleep so that mtimes will be different between runs
            time.sleep(1)

            # Touch a.py to change its mtime without modifying content
            os.utime(os.path.join(src_dir, "a.py"))

            # Modify b.py to access a.x (adding a fine-grained dependency),
            # and add a new c.py, then run mypy for cache2
            with open(os.path.join(src_dir, "b.py"), "w") as f:
                f.write("import a\ndef foo() -> str:\n    return str(a.x)\n")
            with open(os.path.join(src_dir, "c.py"), "w") as f:
                f.write("import a\ny: str = 'world'\n")
            result = subprocess.run(
                [
                    sys.executable,
                    "-m",
                    "mypy",
                    "--sqlite-cache",
                    "--cache-fine-grained",
                    "--cache-dir",
                    cache2,
                    "b.py",
                    "c.py",
                ],
                cwd=src_dir,
                capture_output=True,
                text=True,
                env=env,
            )
            assert result.returncode == 0, f"mypy run 2 failed: {result.stderr}"

            # Find the Python version subdirectory (e.g. "3.14")
            subdirs = [
                e
                for e in os.listdir(cache1)
                if os.path.isdir(os.path.join(cache1, e)) and e[0].isdigit()
            ]
            assert len(subdirs) == 1, f"Expected one version subdir, got {subdirs}"
            ver = subdirs[0]

            # Run diff-cache.py with --sqlite
            result = subprocess.run(
                [
                    sys.executable,
                    _DIFF_CACHE_PATH,
                    "--sqlite",
                    os.path.join(cache1, ver),
                    os.path.join(cache2, ver),
                    output_file,
                ],
                capture_output=True,
                text=True,
                env=env,
            )
            assert result.returncode == 0, f"diff-cache.py failed: {result.stderr}"

            # Verify the output is valid JSON
            with open(output_file) as f:
                data = json.load(f)
            assert isinstance(data, dict)
            assert len(data) > 0, "Expected non-empty diff"

            # Only modified or new files should appear in the diff.
            # b.py changed and c.py is new, so both should be present.
            # a.py did not change, so no a.* keys should appear.
            keys = set(data.keys())
            b_keys = {k for k in keys if "/b." in k or k.startswith("b.")}
            c_keys = {k for k in keys if "/c." in k or k.startswith("c.")}
            a_keys = {k for k in keys if "/a." in k or k.startswith("a.")}
            assert len(a_keys) == 0, f"Unexpected a.* entries in diff: {a_keys}"
            assert len(b_keys) == 3, f"Expected 3 b.* entries in diff, got: {b_keys}"
            assert len(c_keys) == 4, f"Expected 4 c.* entries in diff, got: {c_keys}"

            # The new access to a.x in b.py should create a fine-grained
            # dependency recorded in @root.deps.json.
            assert "@root.deps.json" in keys
            root_deps = json.loads(data["@root.deps.json"])
            assert set(root_deps.keys()) == {
                "<a.x>",
                "<a>",
            }, f"Unexpected root deps keys: {sorted(root_deps.keys())}"
            assert sorted(root_deps["<a.x>"]) == ["b.foo"]
            assert sorted(root_deps["<a>"]) == ["b.foo", "c"]

            # Apply the diff to a copy of cache1 and verify the result.
            cache2_ver = os.path.join(cache2, ver)
            patched = os.path.join(src_dir, "patched")
            patched_ver = os.path.join(patched, ver)
            shutil.copytree(cache1, patched)

            # Snapshot cache entries before applying the diff
            from mypy.metastore import SqliteMetadataStore

            def read_all(cache_dir: str) -> dict[str, bytes]:
                store = SqliteMetadataStore(cache_dir, num_shards=SQLITE_NUM_SHARDS)
                result = {name: store.read(name) for name in store.list_all()}
                store.close()
                return result

            before = read_all(patched_ver)

            # Apply the diff
            result = subprocess.run(
                [sys.executable, _APPLY_CACHE_DIFF_PATH, "--sqlite", patched_ver, output_file],
                capture_output=True,
                text=True,
                env=env,
            )
            assert result.returncode == 0, f"apply-cache-diff.py failed: {result.stderr}"

            after = read_all(patched_ver)

            # a.py entries should be unchanged
            for name in before:
                if name.startswith("a.") or "/a." in name:
                    assert name in after, f"{name} missing after apply"
                    assert before[name] == after[name], f"{name} changed after apply"

            # b.py and c.py entries should match cache2 after applying the diff.
            # Skip .meta.ff files since they contain mtimes that legitimately differ.
            target = read_all(cache2_ver)
            for prefix in ("b.", "c."):
                for name in target:
                    if not (name.startswith(prefix) or f"/{prefix}" in name):
                        continue
                    assert name in after, f"{name} missing after apply"
                    if name.endswith(".meta.ff"):
                        # mtimes legitimately differ, but content should not be identical
                        # to the pre-apply version (it was updated by the diff)
                        assert after[name] != before.get(name), f"{name} unchanged after apply"
                    else:
                        assert after[name] == target[name], f"{name} differs from target"

            # Verify fine-grained deps were applied correctly
            from mypy.util import json_loads

            applied_root_deps = json_loads(after["@root.deps.json"])
            assert set(applied_root_deps.keys()) == {
                "<a.x>",
                "<a>",
            }, f"Unexpected applied root deps keys: {sorted(applied_root_deps.keys())}"
            assert sorted(applied_root_deps["<a.x>"]) == ["b.foo"]
            assert sorted(applied_root_deps["<a>"]) == ["b.foo", "c"]
        finally:
            shutil.rmtree(src_dir, ignore_errors=True)
            shutil.rmtree(os.path.dirname(output_file), ignore_errors=True)


if __name__ == "__main__":
    unittest.main()
