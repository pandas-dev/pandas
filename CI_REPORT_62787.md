# CI Status Update for GH#62787

- Problem: Many unit test jobs failed across platforms after initial fix. Likely due to fragile CoW cleanup popping by stale indices and encountering dead weakrefs.
- Action: Strengthened CoW cleanup in core/internals/blocks.py to rebuild referenced_blocks by object identity, skipping dead weakrefs and avoiding index-based pops. This is more robust across multiprocessing/freethreading and platform variations.
- Tests: Existing regression tests retained; pre-commit clean. Could not run full pytest locally due to binary build requirements, relying on CI for full matrix.
- Next: Wait for CI re-run. If failures persist, share the specific traceback from a failing job to iterate quickly.
