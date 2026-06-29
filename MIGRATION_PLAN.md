# pandas 2.x private optimization port to pandas 3.0.1

## Baselines and protected references

- Target branch: `migration/pandas-3.0.1`.
- Target baseline: official `v3.0.1` at `e04b26f375`.
- Prior port baseline: official `v3.0.3` at `9297b0f74b`.
- Prior port source: `migration/pandas2-private-optimizations`.
- Protected backups: `migration/pandas-3.0.3-backup`,
  `migration-pandas-3.0.3-backup`, and `v3.0.3-optimized`.
- Existing uncommitted `AGENTS.md` changes are unrelated and must not be
  included in migration commits.

The authoritative private-change intent remains
`../pandas_git_export_20260626_zenB_huang`. The reconstructed `../pandas`
repository and the pandas 3.0.3 port are supporting evidence.

## Commit inventory

There are 38 non-merge commits reachable from the prior port after
`v3.0.3`. The merge `117d4e984f` is not replayed; its only independent
change is `3fa2758641` (ASV configuration), which is handled as its own
batch. The final PR merge `08f4d0eedb` is also excluded.

| Order | Prior pandas 3.0.3 commit | Theme | pandas 3.0.1 disposition |
| ---: | --- | --- | --- |
| 1 | `d978bd13ee` | nancorr | adapted |
| 2 | `28df2030b5` | take helpers | pending |
| 3 | `f6d861f8a6` | SwissTable | pending |
| 4 | `8e0304f403` | maybe_convert_objects | pending |
| 5 | `f39ba34d1d` | groupby loops | pending |
| 6 | `3fa2758641` | ASV configuration | pending |
| 7 | `ee85531203` | duplicated nogil repair | pending |
| 8 | `2bfbc167cb` | migration inventory | folded into these documents |
| 9-38 | `2f489b472c` through `4ceff7983d` | remaining audited batches | pending |

The detailed 89-row private commit inventory remains available at
`migration/pandas2-private-optimizations:MIGRATION_PLAN.md`. It will be
ported into this document as each corresponding 3.0.1 batch is resolved,
so statuses describe the 3.0.1 implementation rather than the 3.0.3
implementation.

## Batch order

1. Existing early ports: nancorr, take, SwissTable, object conversion,
   groupby, ASV config, and duplicated repair.
2. Index and join primitives.
3. Low-level hash, sorting, scalar algorithms, and object helpers.
4. Reshape, rolling, skiplist, and benchmark coverage.
5. Datetime and offset fast paths.
6. Python-layer apply, value_counts, take, fillna, and astype paths.
7. Final index/groupby completion and migration audit.

Each code batch is replayed separately. A clean cherry-pick is accepted
only after comparing the touched implementation with both pandas 3.0.1
and pandas 3.0.3. Conflicts are resolved from exported intent, not by
choosing an entire side.

## Batch 1: nancorr

- Original private commit: `eba6d76f8c`.
- Prior pandas 3.0.3 port: `d978bd13ee`.
- pandas 3.0.1 result: adapted to the pre-`no_nans` Welford
  implementation in `pandas/_libs/algos.pyx`.
- Compatibility change: the fully-valid pair fast path requires `N > 0`
  before reading `mat[0]`, preserving empty-frame behavior that the
  exported patch did not guard explicitly.
- Risk: the two-pass path trades Welford's incremental arithmetic for
  fewer divisions and mask branches. Numerical behavior needs the
  existing corr suite and performance workloads on a built extension.
