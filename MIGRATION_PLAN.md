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
| 2 | `28df2030b5` | pad and take helpers | split and adapted |
| 3 | `f6d861f8a6` | SwissTable | directly applied and audited |
| 4 | `8e0304f403` | maybe_convert_objects | directly applied and audited |
| 5 | `f39ba34d1d` | groupby loops | reimplemented from exported commits |
| 6 | `3fa2758641` | ASV configuration | no direct migration |
| 7 | `ee85531203` | pandas 3.0.3 nogil repair | not applicable |
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

## Batch 2a: pad inplace loops

- Original private commits: `706991b28f`, `38f97b58aa`, and
  `6fd0359851`.
- Prior pandas 3.0.3 port: bundled into `28df2030b5` despite that
  commit's take-only subject.
- pandas 3.0.1 result: split into a dedicated batch and adapted without
  importing unrelated pandas 3.0.2/3.0.3 Cython changes.
- Compatibility change: retain the 3.0.1 GIL structure while adding
  no-limit paths, leading-missing-prefix handling, and the four-way 1-D
  loop unroll.
- Risk: mask mutation for leading missing values and limited fills must
  remain identical to the original implementation.

## Batch 2b: take helpers

- Original private commits: `2e964bd967` and `171fe464a6`.
- Prior pandas 3.0.3 port: `28df2030b5`.
- pandas 3.0.1 result: replayed from the authoritative exported API
  shape rather than the pandas 3.0.3 `allow_fill` adaptation.
- Compatibility change: pandas 3.0.1 take helpers always interpret `-1`
  as a fill position, so no later `allow_fill` parameter is introduced.
- Object helpers replace the generated object/object variants and use
  explicit INCREF/XDECREF around raw `PyObject **` assignments.
- The unused exported `Py_DECREF` cimport is omitted.
- Risk: object reference ownership and arbitrary row/column strides
  require a Cython build plus object take/reindex tests.

## Batch 3: SwissTable

- Original private commits: the 24-commit SwissTable series from
  `26f530f608` through the final integration commits
  `12c293f745`, `8c1c10d14f`, `7d2d23951e`, `4423bd227e`,
  `a985978612`, `56996122c8`, and `108a06f167`.
- Prior pandas 3.0.3 port: `f6d861f8a6`.
- pandas 3.0.1 result: the prior port applies cleanly because the
  extension is additive and its four Python integration files retain
  compatible APIs.
- Compatibility checks cover the generated integer maps/factorizers,
  explicit float/complex implementations, `.pxd` declarations, C++17
  Meson target, `compute.use_swisstable`, factorize/unique/safe_sort,
  and merge factorizer selection.
- The later `ee85531203` change is not part of this batch: despite its
  subject, its diff only restructures `pad_inplace` returns inside a
  pandas 3.0.3 `nogil` block. The pandas 3.0.1 adaptation deliberately
  retained a GIL-held pad implementation, so that fix is inapplicable.
- Risk: this is the largest native-code batch and requires full Cython,
  C++ compiler, ABI, dtype, NA/mask, and merge correctness validation.

## Batch 4: maybe_convert_objects

- Original private commits: `844af97539` and `aae84a43d5`.
- Prior pandas 3.0.3 port: `8e0304f403`.
- pandas 3.0.1 result: directly applied; pandas 3.0.1 already uses the
  same single-array return contract as the prior port.
- The exact-builtin helper handles homogeneous bool, float, or int
  families and falls back for mixed families, `pd.NA`, unsafe mode,
  empty/all-missing input, and integers outside uint64/int64 bounds.
- The integer general path computes one complex conversion and reuses
  its real component for the float scratch array.
- Risk: nullable uint selection, mixed signed/large-unsigned fallback,
  and object subclasses require runtime inference tests.

## Batch 5: GroupBy hot loops

- Original private commits: `598e709a9b`, `e38ed47c0e`,
  `d019bc4d8e`, `63b8520f8c`, `82d3f0d7c5`, `eaa0ec8472`,
  `c60600ce36`, and `829e42a602`.
- Prior pandas 3.0.3 port: `f39ba34d1d`.
- pandas 3.0.1 result: reimplemented function by function from the
  exported patches because a whole-commit cherry-pick crossed later
  upstream GroupBy behavior fixes.
- The port covers `group_cummin_max`, `group_any_all`,
  `group_cumprod`, `group_shift_indexer`, `group_last`,
  `group_mean`, `group_min_max`, and `group_prod`.
- pandas 3.0.1 compatibility:
  - arbitrary-stride `values` and masks remain accepted;
  - `skipna=False` poisoning is retained for mean/min/max/prod;
  - Kahan compensation is retained for mean;
  - `lab < 0` cumulative handling stays at pandas 3.0.1 behavior
    instead of importing the later 3.0.2 nullable-dtype bug fix;
  - empty/all-NA shift inputs avoid zero-length pointer extraction.
- Risk: duplicated loop bodies increase review surface; every
  mask/skipna/datetime/min_count branch needs native runtime coverage.

## Batch 6: environment-only follow-ups

- Prior ASV commit `3fa2758641`: no direct migration.
  - It is an environment-specific change rather than a private pandas
    optimization.
  - It hard-codes `HEAD` and Python 3.14, removes the documented config
    form, and contains a duplicate `versioneer[toml]` matrix key.
  - pandas 3.0.1's existing ASV configuration is retained.
- Prior repair `ee85531203`: not applicable.
  - Its subject mentions duplicated, but its actual diff only moves
    returns inside pandas 3.0.3 `pad_inplace`/`pad_2d_inplace` `nogil`
    blocks.
  - The pandas 3.0.1 pad adaptation intentionally remains GIL-held, so
    the Cython return restriction does not exist.

## Batch 7: index search primitives and groupsort

- Original private commits: `0744c0555f`, duplicate export
  `df0c75b9fd`, `53b373262d`, and the groupsort portion of
  `403df2a143`.
- Prior pandas 3.0.3 port: `2f489b472c`.
- pandas 3.0.1 result: directly adapted because the affected engine and
  generated-template APIs are unchanged from the target baseline.
- The port vectorizes `BaseMultiIndexCodesEngine` code preparation,
  adds typed left/right binary-search hooks for non-complex unmasked
  engines, and unrolls the two `groupsort_indexer` row loops.
- pandas 3.0.1 compatibility:
  - `ObjectEngine` keeps its tuple-aware `_bin_search` implementation;
  - complex and masked engines retain NumPy search behavior;
  - generated numeric methods rely on the existing `_check_type`
    contract before converting the lookup key;
  - empty input remains valid because all unrolled loops use a
    four-element-aligned limit.
- The join take-helper portion of `403df2a143` remains for the join
  batch; it is not mixed into this index-engine commit.
- Risk: generated Cython methods and MultiIndex code-width handling
  require a native build and focused index tests.

## Batch 8: RangeIndex concat planning

- Original private commit: `1b285b9697`.
- Prior pandas 3.0.3 port: `79ba386d4a`.
- pandas 3.0.1 result: adapted after a narrow conflict in
  `RangeIndex._concat`; the conflict was caused by nearby upstream type
  changes, not by different concat semantics.
- `pandas._libs.lib.concat_range_indexes` now performs the all-RangeIndex
  planning loop and returns a compact tag plus construction arguments.
- The Python layer retains pandas 3.0.1 Index construction, integer
  fallback, name propagation, repeated-range tiling, and empty-range
  behavior.
- Risk: the helper accepts Python objects and accesses `_range`
  directly; native tests must cover mixed Index inputs, one-element
  ranges, repeated ranges, empty ranges, and discontinuous ranges.

## Batch 9: khash insertion macros

- Original private commits: `2777612697` and `42dc387d2c`.
- Prior pandas 3.0.3 port: `e23914d92f`.
- pandas 3.0.1 result: directly applied to the unchanged vendored khash
  insertion macro.
- Portable likely/unlikely macros annotate the common empty-slot path,
  and local key/flag pointers avoid repeated structure member loads.
- Pointer caches are initialized after any resize, so they cannot refer
  to pre-resize storage.
- Risk: this vendored header is expanded across many generated hash
  tables and needs compiler coverage on GCC/Clang and MSVC.

## Batch 10: stable safe_sort controls

- Original private commits: `51a2b98159`, `31a63abb20`, and
  `1c4c300a36`.
- Prior pandas 3.0.3 port: `dded84d01c`.
- pandas 3.0.1 result: adapted with only the local `SortKind` type
  import; the broader 3.0.3 typing-import reorganization was excluded.
- `safe_sort` threads an explicit NumPy sort kind through normal,
  mixed, and reverse-code sorting. `factorize(sort=True)` selects
  stable ordering to preserve equal-value order.
- The existing pandas 3.0.1 SwissTable-aware factorization path remains
  unchanged.
- Risk: stable remapping must preserve NA sentinels, mixed integer/string
  ordering, tuple fallback, and ExtensionArray behavior.

## Batch 11: scalar Cython algorithm hot paths

- Original private commits: `4c4a5096dc` and `0e2677777a`.
- Prior pandas 3.0.3 port: `28f873ba8e`.
- pandas 3.0.1 result: directly applied to the matching Cython
  implementations.
- `kth_smallest_c` advances typed pointers instead of repeatedly
  indexing from the array base. `checknull` converts exact float objects
  once and calls C `isnan`, while complex values retain self-comparison.
- The removed pandas 2.x `inf_as_na` signature is not restored; infinity
  remains non-missing in pandas 3.0.1.
- Risk: quickselect pointer bounds and NumPy float scalar conversion
  require native compilation and runtime coverage.

## Batch 12: object array construction

- Original private commits: `fefbc7b5a4` and `d72753640f`.
- Prior pandas 3.0.3 port: `0672469c2a`.
- pandas 3.0.1 result: directly applied to compatible cast and lib
  interfaces.
- Flat list/tuple input uses NumPy's direct object-array construction;
  nested list-like and ndarray input uses a Cython top-level copy helper
  so the result remains one-dimensional.
- The helper uses `PySequence_Fast`, explicitly transfers references
  into NumPy's initialized object slots, and releases the temporary
  sequence.
- Risk: reference ownership, sequence subclasses, nested arrays, and
  exception paths require a native extension build.

## Batch 13: lib object pointer helpers

- Original private commits: `a60133b265`, `c064f8cd15`,
  `defb42e92e`, and `454f5e27f1`.
- Prior pandas 3.0.3 port: `5efbbfabd2`.
- pandas 3.0.1 result: directly applied to compatible `_libs.lib`
  functions and consolidated Python C-API declarations.
- Contiguous object-array equivalence uses flat pointers; `fast_zip`
  writes tuple items through result pointers; `eq_NA_compat` traverses
  the one-dimensional IndexEngine values buffer using its stride.
- Non-contiguous object equality retains the existing MultiIter path.
- Risk: borrowed/new/stolen reference rules and non-contiguous object
  arrays need native memory-safety and behavior tests.

## Batch 14: reshape writers

- Original private commits: `15dba60291` and `8232deb8a0`.
- Prior pandas 3.0.3 port: `b545a39fc2`.
- pandas 3.0.1 result: dense get_dummies applies directly; numeric
  unstack required a target-specific reimplementation.
- Wide dense dummies use a fused Cython writer over the Fortran-order
  buffer, with column-wise comparison for small cardinality and the
  original NumPy fallback for unsupported output/code dtypes.
- Unlike pandas 3.0.3, pandas 3.0.1 still passes `new_mask` to
  `_libs.reshape.unstack`, so exported commit `8232deb8a0` is applicable.
  Numeric specializations write the mask in a separate pass; object
  specialization retains the original single pass and reference writes.
- This is a documented divergence from the prior 3.0.3 audit, which
  correctly skipped the optimization only because 3.0.3 removed
  `new_mask`.
- Risk: fused-type specialization, Fortran buffer indexing, NA codes,
  and object-reference writes require a native build.

## Batch 15: rolling sum loop structure

- Original private commit: `66863407cc`.
- Prior pandas 3.0.3 port: `e7b5a3234a`.
- pandas 3.0.1 result: directly applied to the matching
  `roll_sum` implementation.
- Monotonic windows initialize once, cache prior bounds, and update by
  removing/adding only changed ranges. Non-monotonic windows use a
  separate full-recomputation loop without redundant cleanup stores.
- Kahan add/remove compensation and consecutive-identical-value
  handling are preserved.
- Risk: empty/disjoint/overlapping windows, non-monotonic bounds,
  infinities, NaNs, and minimum-period behavior need native tests.

## Batch 16: index behavior no-port audit

- Original private commits: `eaaffa04b4`, its corrective revert
  `0b958ce8fd`, and `72efd58130`.
- Prior pandas 3.0.3 audit commits: `f001d1bdb8` and `4c339e66c5`.
- `is_monotonic`: no code migration. The later private commit reverted
  the block optimization after a Xiecheng MergeJoin correctness failure;
  pandas 3.0.1 already retains the restored scalar upstream loop.
- `ObjectEngine._get_loc_duplicates`: no duplicate migration. pandas
  3.0.1 already has a monotonic duplicate override using tuple-safe
  `_bin_search` and `_bin_search_right`.
- Risk: the original MergeJoin failure should be reproduced before any
  future monotonic replacement is considered.
