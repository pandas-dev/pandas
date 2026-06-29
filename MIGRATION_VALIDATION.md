# pandas 3.0.1 private optimization migration validation

Full ASV is not available in this environment. No ASV result is claimed.

## Repository protection checks

Executed:

- `git status --short --branch`
- `git branch -a`
- `git tag --list`
- `git log --all --decorate --oneline --graph -100`
- `git show --no-patch --decorate HEAD`
- ancestry and first-parent commit-range checks for the protected 3.0.3
  branches and tags

Findings:

- `HEAD` started at official `v3.0.1` on
  `migration/pandas-3.0.1`.
- The prior port and its backup refs remain untouched.
- `AGENTS.md` was already modified and is excluded from every migration
  commit.
- The active Python is 3.14.5 with NumPy 2.5.0.
- `python -m cython --version` fails because Cython is not installed.

## Batch 1: nancorr

Sources inspected:

- Exported private commit `eba6d76f8cb5e66cc30a97d85569066302da8ec5`
  in `all_commits_with_file_diffs.patch`.
- Reconstructed pandas 2.x matching history.
- Prior pandas 3.0.3 port `d978bd13ee`.
- pandas 3.0.1 and pandas 3.0.3 `nancorr` implementations.
- Existing DataFrame correlation tests for missing values, constant
  columns, small numbers, and `min_periods`.

Adaptation checks:

- Column validity is computed once per column.
- The two-pass path is used only when both columns are fully finite and
  `N > 0`.
- Missing-data pairs retain pandas 3.0.1's Welford loop and pairwise
  finite mask.
- Constant fully-valid columns return NaN without the centered pass.
- Empty input cannot read `mat[0]`.
- Covariance divisor and correlation clipping remain unchanged.

Executed before commit:

- `git diff --check`
- conflict-marker scan
- `python -m compileall -q pandas`
- declaration/use inspection for the added Cython memory view and scalar
  variables
- `python -m pytest pandas/tests/frame/methods/test_cov_corr.py -q`
  stopped while loading `pandas/conftest.py` because `dateutil` is not
  installed; no tests were collected.

Not executed:

- Cython compilation: Cython is not installed.
- Runtime pandas tests: blocked by the missing `dateutil` dependency and,
  after that is resolved, may require a pandas extension build compatible
  with this checkout.
- ASV: unavailable in this environment.

Follow-up validation:

- `python -m pytest pandas/tests/frame/methods/test_cov_corr.py`
- `python -m pytest pandas/tests/test_nanops.py -k "nancorr or nancov"`
- Build the Cython extensions before the runtime tests.
- Run the DataFrame correlation ASV cases with fully valid columns,
  sparse missing columns, constant columns, and varying `min_periods`.

## Batch 2a: pad inplace loops

Sources inspected:

- Exported commits `706991b28f`, `38f97b58aa`, and `6fd0359851`.
- Prior pandas 3.0.3 port `28df2030b5`.
- pandas 3.0.1 `pad_inplace` and `pad_2d_inplace`.

Adaptation checks:

- Leading missing prefixes remain missing because there is no prior
  value to propagate.
- The no-limit 1-D path clears masks only after the first valid value.
- The four-way unroll updates the carried value in source order.
- Limited fills reset their count after every valid value.
- The 2-D no-limit path applies the same rules independently per row.
- The pandas 3.0.3-only surrounding `nogil` and scalar-fill helpers are
  not imported.

Executed:

- `git diff --check`
- `python -m compileall -q pandas`
- `python -m pytest pandas/tests/libs/test_libalgos.py -q` stopped
  during conftest import because `dateutil` is missing.

Not executed:

- Cython compilation: Cython is not installed.
- Runtime pad/backfill tests: pytest is blocked by missing `dateutil`.
- ASV: not run in this environment.

Follow-up validation:

- `python -m pytest pandas/tests/libs/test_libalgos.py -k "pad or backfill"`
- `python -m pytest pandas/tests/frame/methods/test_fillna.py -k "pad or backfill"`
- Run pad/backfill ASV cases for unlimited and limited 1-D/2-D inputs,
  including leading, trailing, and all-missing masks.

## Batch 2b: take helpers

Sources inspected:

- Exported commits `2e964bd967` and `171fe464a6`.
- Reconstructed commits `f0fd052` and `529be9b`.
- Prior pandas 3.0.3 port `28df2030b5`.
- pandas 3.0.1 take template signatures and generated dispatch list.

Adaptation checks:

- The pandas 3.0.3-only `allow_fill` parameter is not introduced.
- Numeric contiguous-slice copying requires matching input/output
  element types and unit inner strides.
- Non-contiguous numeric paths retain `-1` fill handling in both the
  unrolled body and tail.
- Object helpers honor source and output strides in pointer units.
- Every raw object store increments the incoming reference before
  decrementing the replaced output reference.
- Empty 2-D axis-1 inputs return before taking row pointers.
- The generated object/object variants are disabled exactly once.

Executed:

- Applied the reconstructed commits' exact diffs with Git's patch
  machinery after confirming they match the export.
- `git diff --check`
- declaration/use scan for `PyObject`, `Py_INCREF`, and `Py_XDECREF`
- `python -m compileall -q pandas`
- `python -m pytest pandas/tests/test_take.py -q` stopped during
  conftest import because `dateutil` is missing.

Not executed:

- Cython compilation: Cython is not installed.
- Runtime take tests: pytest is blocked by missing `dateutil`.
- ASV: not run in this environment.

Follow-up validation:

- `python -m pytest pandas/tests/test_take.py`
- `python -m pytest pandas/tests/frame/methods/test_reindex.py`
- Run ASV for numeric `take_2d_axis1` contiguous slices, irregular
  indexers, fill indexers, and object 1-D/2-D take variants.

## Batch 3: SwissTable

Sources inspected:

- All 24 SwissTable export commits identified by the 89-row inventory.
- Prior pandas 3.0.3 port `f6d861f8a6`.
- Follow-up commit `ee85531203`, which does not touch SwissTable.
- pandas 3.0.1 Meson, algorithms, config, and merge integration points.

Static consistency checks:

- `swisstable.pyx` includes the generated
  `swisstable_class_helper.pxi`.
- The `.pxd` declares all integer, float, and complex Map types used by
  Python integration.
- The helper generates integer Map, membership, value-count,
  duplicated, and Factorizer implementations.
- Explicit float/complex implementations expose the same operations.
- Empty contiguous duplicated inputs are guarded by `n > 0` before
  taking `&values[0]` or `&result[0]`.
- Meson generates the helper, builds the extension as C++17, includes
  pandas/NumPy headers, and adds ARM CRC flags only on AArch64.
- `algorithms.py` and `merge.py` fall back to khash for unsupported
  dtypes and object keys.
- Masked and Arrow numeric merge keys select their NumPy dtype before
  choosing a Swiss Factorizer.

Executed:

- `git diff --check`
- `python -m compileall -q pandas asv_bench/benchmarks/swisstable.py`
- `python -m pytest pandas/tests/libs/test_swisstable.py -q` stopped
  during conftest import because `dateutil` is missing.
- Native tool probes found no `g++`, Meson, or Ninja executables.

Not executed:

- Cython compilation: Cython is not installed.
- Runtime SwissTable tests: pytest is blocked by missing `dateutil`.
- Full Meson/native build and ABI validation: Cython, a C++ compiler,
  Meson, and Ninja are unavailable.
- ASV: not run in this environment.

Follow-up validation:

- `python -m pytest pandas/tests/libs/test_swisstable.py`
- `python -m pytest pandas/tests/test_algorithms.py -k "unique or factorize or safe_sort"`
- `python -m pytest pandas/tests/reshape/merge`
- Run `asv_bench/benchmarks/swisstable.py` with the option enabled and
  disabled, including strided arrays, masks, NaN/NA, complex values,
  and all merge factorizer dtypes.

## Batch 4: maybe_convert_objects

Sources inspected:

- Exported commits `844af97539` and `aae84a43d5`.
- Prior pandas 3.0.3 port `8e0304f403`.
- pandas 3.0.1 `maybe_convert_objects` signature and return contract.

Static consistency checks:

- Fast-path entry conditions preserve `safe` and non-numeric conversion
  behavior.
- Missing bool values require nullable conversion; otherwise the
  general path remains responsible.
- Missing integer values produce nullable integer arrays only when
  requested and float64 otherwise.
- Mixed negative and greater-than-int64 positive integers fall back.
- Values outside int64/uint64 bounds fall back.
- The general path retains all datetime, timedelta, period, string, and
  object detection.

Executed:

- `git diff --check`
- `python -m compileall -q pandas`
- `python -m pytest pandas/tests/dtypes/test_inference.py -q -k
  maybe_convert_objects` stopped during conftest import because
  `dateutil` is missing.

Not executed:

- Cython compilation: Cython is not installed.
- Runtime inference tests: pytest is blocked by missing `dateutil`.
- ASV: not run in this environment.

Follow-up validation:

- `python -m pytest pandas/tests/dtypes/test_inference.py -k maybe_convert_objects`
- Run object-conversion ASV for homogeneous numeric blocks, nullable
  results, mixed numeric families, and fallback-heavy object arrays.

## Batch 5: GroupBy hot loops

Sources inspected:

- Eight authoritative export commits and reconstructed commits
  `41e01ce`, `c742f34`, `e942322`, `7adf4d5`, `292604f`,
  `1e905a4`, `7f0119c`, and `9740cc5`.
- Prior pandas 3.0.3 port `f39ba34d1d`.
- pandas 3.0.1 implementations and the pandas 3.0.2/3.0.3 GroupBy
  changes touching the same functions.

Adaptation checks:

- Loop-invariant mask, skipna, direction, and min/max branches move
  outside inner column loops.
- Pointer row caches are used only for internally allocated contiguous
  outputs/accumulators; input values and masks keep arbitrary strides.
- Mean retains Kahan compensation and resets non-finite compensation.
- Mean, min/max, and prod preserve `skipna=False` poisoned results.
- Cumulative min/max and prod retain pandas 3.0.1 negative-label
  behavior rather than importing later upstream fixes.
- Shift handles positive/negative periods without `sign * i`, and
  returns safely for zero rows or zero groups.
- Source-only non-ASCII and decorative comments were removed.
- No conflict markers remain.

Executed:

- `git diff --check`
- conflict-marker and non-ASCII scan on `groupby.pyx`
- `python -m compileall -q pandas`
- `python -m pytest pandas/tests/groupby/test_reductions.py -q`
  stopped during conftest import because `dateutil` is missing.

Not executed:

- Cython compilation: Cython is not installed.
- Runtime GroupBy tests: pytest is blocked by missing `dateutil`.
- ASV: not run in this environment.

Follow-up validation:

- `python -m pytest pandas/tests/groupby/test_groupby.py`
- `python -m pytest pandas/tests/groupby/test_reductions.py`
- `python -m pytest pandas/tests/groupby -k "cummin or cummax or cumprod or shift or first or last or mean or min or max or prod"`
- Run the corresponding GroupBy ASV matrix for mask/no-mask,
  skipna true/false, datetime, nullable, negative labels, and empty
  groups.

## Batch 6: environment-only follow-ups

Inspected:

- Prior ASV config commit `3fa2758641`.
- Current pandas 3.0.1 `asv_bench/asv.conf.json`.
- Prior repair `ee85531203` and its full `algos.pyx` diff.

Decision:

- No ASV config file change. The prior commit is machine-specific and
  contains a duplicate matrix key; retaining the 3.0.1 config avoids
  changing benchmark environment policy as part of a code port.
- No `ee85531203` code change. Its `nogil` return fix targets a later
  implementation structure absent from this branch.

Validation:

- Compared both prior diffs with the current 3.0.1 files.
- ASV was not run.

## Batch 7: index search primitives and groupsort

Sources inspected:

- Exported commits `0744c0555f`, `df0c75b9fd`, `53b373262d`, and the
  groupsort portion of `403df2a143`.
- Prior pandas 3.0.3 port `2f489b472c`.
- pandas 3.0.1 index engine, generated engine template, and
  `groupsort_indexer` implementations.

Static consistency checks:

- `_searchsorted_left` and `_searchsorted_right` share the same
  caller-side `_check_type` contract.
- Generated binary search is excluded for complex and masked engines.
- `ObjectEngine` retains tuple-key-safe search behavior.
- `level_has_nans` remains a one-dimensional per-level boolean
  container and has one indexed consumer.
- MultiIndex codes are shifted as signed `int64` before being viewed as
  `uint64`, preserving the missing-code representation.
- Four-way groupsort loops handle empty and non-multiple-of-four input
  through the scalar tail.
- The template imports match every generated Cython scalar type.

Executed:

- `git diff --cached --check`
- Static call-site and declaration inspection with `rg`.
- `python -m compileall -q pandas`
- `python -m pytest pandas/tests/indexes/test_engines.py -q` stopped
  while importing `pandas/conftest.py` because `dateutil` is missing.

Not executed:

- Cython compilation: Cython is not installed.
- Runtime index assertions: pytest did not reach collection.
- ASV: not run in this environment.

Follow-up validation:

- `python -m pytest pandas/tests/indexes/multi`
- `python -m pytest pandas/tests/indexes/test_engines.py`
- `python -m pytest pandas/tests/reshape/merge`
- Run ASV for MultiIndex engine construction, duplicate monotonic
  lookup, and `groupsort_indexer` consumers.

## Batch 8: RangeIndex concat planning

Sources inspected:

- Exported commit `1b285b9697`.
- Prior pandas 3.0.3 port `79ba386d4a`.
- pandas 3.0.1 `RangeIndex._concat` fallback and construction paths.

Adaptation checks:

- Non-RangeIndex input uses the unchanged parent fallback.
- A single RangeIndex is returned with the requested name behavior.
- Repeated identical non-empty ranges use `np.tile`.
- Discontinuous ranges concatenate materialized values.
- A single non-empty range preserves its original step.
- All-empty input constructs the standard empty RangeIndex.
- The `.pyi` declaration matches the four-item planning tuple.

Executed:

- `git diff --check`
- `python -m py_compile pandas/core/indexes/range.py`
- `python -m compileall -q pandas`
- `python -m pytest pandas/tests/indexes/ranges -q` stopped while
  importing `pandas/conftest.py` because `dateutil` is missing.

Not executed:

- Cython compilation: Cython is not installed.
- Runtime assertions: pytest did not reach collection.
- ASV: not run in this environment.

Follow-up validation:

- `python -m pytest pandas/tests/indexes/ranges`
- Run RangeIndex concat ASV for consecutive, discontinuous, repeated,
  singleton, empty, and mixed-index inputs.

## Batch 9: khash insertion macros

Sources inspected:

- Exported commits `2777612697` and `42dc387d2c`.
- Prior pandas 3.0.3 port `e23914d92f`.
- pandas 3.0.1 vendored `kh_put` macro expansion.

Static consistency checks:

- GNU/Clang branch prediction builtins have a portable fallback.
- Key and flag pointers are cached only after possible table resize.
- Probe, deleted-slot, occupancy, and return-code behavior is unchanged.
- Every subsequent insertion access uses the cached pointers.

Executed:

- `git diff --cached --check`
- Static macro expansion inspection.

Not executed:

- C/C++ compilation: no native compiler toolchain is available.
- Runtime hash-table tests: edited extensions cannot be rebuilt.
- ASV: not run in this environment.

Follow-up validation:

- Build generated hash tables with GCC/Clang and MSVC.
- Run hashtable, factorize, unique, merge, and SwissTable-off tests.
- Benchmark collision-heavy and mostly-new-key insertion workloads.

## Batch 10: stable safe_sort controls

Sources inspected:

- Exported commits `51a2b98159`, `31a63abb20`, and `1c4c300a36`.
- Prior pandas 3.0.3 port `dded84d01c`.
- pandas 3.0.1 `factorize`, `safe_sort`, and `_sort_mixed`.

Adaptation checks:

- Only `SortKind` is added to the 3.0.1 `TYPE_CHECKING` imports.
- Stable sort reaches values, mixed numeric/string subsets, and reverse
  code remapping.
- Default callers retain quicksort behavior.
- `factorize(sort=True)` explicitly requests stable ordering.
- NA sentinel verification and tuple fallback are unchanged.

Executed:

- `git diff --cached --check`
- `python -m py_compile pandas/core/algorithms.py`
- `python -m pytest pandas/tests/test_algos.py -q -k
  'safe_sort or factorize'` stopped while importing
  `pandas/conftest.py` because `dateutil` is missing.

Not executed:

- Runtime assertions: pytest did not reach collection.
- ASV: not run in this environment.

Follow-up validation:

- `python -m pytest pandas/tests/test_algos.py -k
  "safe_sort or factorize"`
- Benchmark sorted factorization for numeric, object, mixed, duplicate,
  and missing-value inputs.

## Batch 11: scalar Cython algorithm hot paths

Sources inspected:

- Exported commits `4c4a5096dc` and `0e2677777a`.
- Prior pandas 3.0.3 port `28f873ba8e`.
- pandas 3.0.1 quickselect and missing-value signatures.

Adaptation checks:

- Pointer quickselect preserves pivot, partition, and target movement.
- Existing callers provide non-empty arrays and valid `k`.
- Float `checknull` calls C `isnan` after one `float64_t` conversion.
- Complex, decimal, datetime64, timedelta64, `None`, `NaT`, and `pd.NA`
  paths are unchanged.
- NumPy float scalar regression coverage includes NaN, finite, and
  infinity values.

Executed:

- `git diff --cached --check`
- `python -m py_compile pandas/tests/dtypes/test_missing.py`
- `python -m compileall -q pandas`
- `python -m pytest pandas/tests/dtypes/test_missing.py -q -k
  checknull` stopped during conftest import because `dateutil` is
  missing.

Not executed:

- Cython compilation: Cython is not installed.
- Runtime assertions: pytest did not reach collection.
- ASV: not run in this environment.

Follow-up validation:

- `python -m pytest pandas/tests/dtypes/test_missing.py -k checknull`
- Run selection/groupby quantile tests for all fused numeric types.
- Benchmark scalar missing checks and kth-selection workloads.

## Batch 12: object array construction

Sources inspected:

- Exported commits `fefbc7b5a4` and `d72753640f`.
- Prior pandas 3.0.3 port `0672469c2a`.
- pandas 3.0.1 cast helper call sites and lib declarations.

Adaptation checks:

- Flat list/tuple input cannot be expanded into extra dimensions.
- Nested list-like input remains a one-dimensional object array.
- Two-dimensional ndarray rows remain top-level objects.
- Empty input returns an empty object array.
- `PySequence_Fast` ownership and item INCREF/object-slot DECREF are
  paired, and the `.pyi` declaration matches the callable.

Executed:

- `git diff --cached --check`
- `python -m py_compile pandas/core/dtypes/cast.py
  pandas/tests/dtypes/cast/test_construct_object_arr.py`
- `python -m compileall -q pandas`
- `python -m pytest
  pandas/tests/dtypes/cast/test_construct_object_arr.py -q` stopped
  during conftest import because `dateutil` is missing.

Not executed:

- Cython compilation: Cython is not installed.
- Runtime assertions: pytest did not reach collection.
- ASV: not run in this environment.

Follow-up validation:

- Run the focused construction test file after rebuilding extensions.
- Exercise lists, tuples, generators with length, nested arrays,
  sequence subclasses, empty input, and object lifetimes.
- Benchmark Series/DataFrame object construction from flat Python
  sequences.

## Batch 13: lib object pointer helpers

Sources inspected:

- Exported commits `a60133b265`, `c064f8cd15`, `defb42e92e`, and
  `454f5e27f1`.
- Prior pandas 3.0.3 port `5efbbfabd2`.
- pandas 3.0.1 callers in IndexEngine, MultiIndex, sparse construction,
  and object-array equivalence.

Static consistency checks:

- Contiguous equality flattens equal-shaped arrays and recursively
  preserves nested ndarray handling.
- Non-contiguous equality retains the original broadcast iterator.
- `eq_NA_compat` is called with the one-dimensional engine values
  buffer and honors arbitrary first-axis stride.
- `fast_zip` preserves `PyTuple_SET_ITEM` stolen-reference ownership by
  incrementing each Python value before insertion.
- Existing `.pxd` declaration for `eq_NA_compat` remains compatible.

Executed:

- `git diff --cached --check`
- Static call-site and Python C-API ownership inspection.
- `python -m compileall -q pandas`
- Focused pytest invocation stopped during conftest import because
  `dateutil` is missing.

Not executed:

- Cython/native compilation: Cython and a compiler are unavailable.
- Runtime assertions and reference-leak checks.
- ASV: not run in this environment.

Follow-up validation:

- Run missing-value, index engine, object equivalence, MultiIndex, and
  sparse construction tests after rebuilding extensions.
- Add debug-build/refcount coverage for nested object arrays and tuple
  construction.
- Benchmark contiguous/non-contiguous object equality, NA lookup, and
  MultiIndex tuple materialization.

## Batch 14: reshape writers

Sources inspected:

- Exported commits `15dba60291` and `8232deb8a0`.
- Prior pandas 3.0.3 port `b545a39fc2`.
- pandas 2.x numeric/object unstack branches and pandas 3.0.1/3.0.3
  unstack signatures.

Adaptation checks:

- Dense dummy writes use `code * n_rows + row` against a Fortran-order
  flattened buffer.
- Bool output is viewed as uint8; uint8, int64, and float64 have fused
  writers; unsupported dtypes use the original NumPy fallback.
- Negative codes are skipped, preserving non-dummy-NA behavior.
- Small-cardinality column comparisons preserve all dtype behavior.
- Numeric unstack pre-materializes every mask position, then writes only
  valid values.
- Object unstack retains the original single-pass assignment and mask
  write, avoiding an extra reference traversal.
- `.pyi` helper declaration matches the Python-visible Cython function.

Executed:

- `git diff --cached --check`
- `python -m py_compile pandas/core/reshape/encoding.py
  pandas/tests/reshape/test_get_dummies.py`
- `python -m compileall -q pandas`
- Focused get_dummies/unstack pytest invocation stopped during conftest
  import because `dateutil` is missing.

Not executed:

- Cython compilation: Cython is not installed.
- Runtime assertions: pytest did not reach collection.
- ASV: not run in this environment.

Follow-up validation:

- Run get_dummies and unstack test suites after rebuilding extensions.
- Cover bool/uint8/int64/float64/unsupported dtypes, NA codes,
  drop-first, dummy-NA, object values, and empty shapes.
- Benchmark wide dense get_dummies and numeric unstack mask-heavy and
  mask-light workloads.

## Batch 15: rolling sum loop structure

Sources inspected:

- Exported commit `66863407cc`.
- Prior pandas 3.0.3 port `e7b5a3234a`.
- pandas 3.0.1 rolling sum state and bounds handling.

Static consistency checks:

- The monotonic branch is entered only for non-empty bounds arrays.
- Prior start/end bounds are updated after every output.
- Disjoint windows reset all rolling state before recomputation.
- Overlapping windows remove expired values and add new values.
- Non-monotonic windows recompute every output from an empty state.
- Kahan add/remove compensation and `calc_sum` inputs are unchanged.

Executed:

- `git diff --cached --check`
- Static control-flow inspection.
- `python -m compileall -q pandas`
- `python -m pytest pandas/tests/window/test_rolling.py -q -k sum`
  stopped during conftest import because `dateutil` is missing.

Not executed:

- Cython compilation: Cython is not installed.
- Runtime assertions: pytest did not reach collection.
- ASV: not run in this environment.

Follow-up validation:

- Run rolling/window sum tests after rebuilding extensions.
- Cover expanding, fixed, variable, empty, disjoint, overlapping, and
  non-monotonic bounds with NaN/inf and different `min_periods`.
- Run rolling sum ASV for monotonic and variable-window workloads.

## Batch 16: index behavior no-port audit

Sources inspected:

- Exported `is_monotonic` commits `eaaffa04b4` and `0b958ce8fd`.
- Exported ObjectEngine commit `72efd58130`.
- Prior audit commits `f001d1bdb8` and `4c339e66c5`.
- pandas 3.0.1 scalar monotonic loop and ObjectEngine duplicate lookup.

Decision:

- Do not restore the reverted block-based monotonic implementation; the
  export history records a correctness regression.
- Do not duplicate ObjectEngine code; the target already contains the
  tuple-safe equivalent optimization.

Executed:

- Static comparison of export history and target implementations.
- `git diff --check`

Not executed:

- Runtime regression tests: pytest remains blocked by missing
  `dateutil`.
- ASV: not run in this environment.

Follow-up validation:

- `python -m pytest pandas/tests/indexes/test_monotonic.py`
- Run object-index `get_loc` tests for monotonic duplicates and tuple
  keys.
- Reproduce the Xiecheng MergeJoin failure before designing another
  monotonic fast path.
- Benchmark ObjectEngine duplicate lookup; no monotonic benchmark claim
  is made for the reverted optimization.

## Batch 17: object join indexers

Sources inspected:

- Exported object-indexer commit `52f9c792ec`.
- Exported many-to-many commit `07ccc6e28e`.
- Prior pandas 3.0.3 port `6e766c4650`.
- pandas 3.0.1 grouped joins and all four monotonic indexer loops.

Adaptation checks:

- Raw pointers are entered only for C-contiguous object ndarrays.
- Generic fused loops remain the fallback for every other input.
- Count and fill passes use the same duplicate-run advancement.
- Empty left/right paths allocate correctly sized output/indexers.
- Python comparisons use exception-aware `PyObject_RichCompareBool`.
- Added tests cover unique, left, inner, and outer object indexers.
- Existing `sort=False` grouped join logic already covers
  `07ccc6e28e`.

Executed:

- `git diff --cached --check`
- Static loop/count/fill and call-guard inspection.
- `python -m py_compile pandas/tests/libs/test_join.py`
- `python -m compileall -q pandas`
- `python -m pytest pandas/tests/libs/test_join.py -q -k object`
  stopped during conftest import because `dateutil` is missing.

Not executed:

- Cython compilation: Cython is not installed.
- Runtime assertions: pytest did not reach collection.
- ASV: not run in this environment.

Follow-up validation:

- Run low-level join and merge test suites after rebuilding extensions.
- Cover empty, duplicate-heavy, non-contiguous, mixed-comparison-error,
  tuple, and missing object values.
- Benchmark monotonic object indexers and many-to-many `sort=False`
  grouped joins.

## Batch 18: skiplist allocation and traversal

Sources inspected:

- Exported commit `91019df988`.
- Prior pandas 3.0.3 port `8ece4c5407`.
- pandas 3.0.1 and 3.0.3 skiplist header differences.

Adaptation checks:

- Flexible storage offsets follow `node_t`, pointer-array, then
  int-array alignment requirements.
- Zero-level NIL nodes use null next/width pointers.
- Destruction recursively follows levels and frees each node block once
  according to the existing reference count.
- Allocation-failure cleanup still accepts partially initialized
  structures.
- Width/rank updates and duplicate ordering are preserved.
- Prefetch calls use reachable NIL/node pointers and compile to no-ops
  outside GCC/Clang.
- No conflict markers remain.

Executed:

- `git diff --cached --check`
- Allocation/free, reference-count, traversal, and conflict-marker
  inspection.

Not executed:

- C syntax/native build: no C compiler is available.
- Rolling quantile/median and rank runtime tests.
- ASV: not run in this environment.

Follow-up validation:

- Build with GCC/Clang and MSVC, including debug allocators/sanitizers.
- Run rolling quantile/median and ranking tests for empty, singleton,
  duplicate-heavy, NaN, insert/remove, and allocation-failure cases.
- Benchmark skiplist-heavy rolling quantile/median and rank workloads.

## Batch 19: Xiecheng ASV workloads

Sources inspected:

- Exported commits `a3739ab1e4`, `42b376312f`, `baca9d6de2`, and
  `fcdd01e0a4`.
- Prior pandas 3.0.3 port `473b0e68ff`.
- Final benchmark method names and setup mutations.

Executed:

- `git diff --cached --check`
- `python -m py_compile asv_bench/benchmarks/xiecheng.py`
- Static check that benchmark method names are unique and categorical
  conversion does not overwrite source object columns.

Not executed:

- ASV was not run in this environment.
- Product runtime validation is not applicable to this benchmark-only
  batch.

Follow-up validation:

- Build pandas 3.0.1 in the ASV environment.
- Run `asv run -b xiecheng` or the repository-supported equivalent.
- Start with the smallest parameter case, then run the full rows/users/
  categories matrix while monitoring memory.

## Batch 20: add_overflowsafe fast paths

Sources inspected:

- Exported commits `9e640a5d2c`, `3f77b3b691`, and `82081460b3`.
- Prior pandas 3.0.3 port `e48041d9cc`.
- pandas 3.0.1 overflow decorator, shape contract, and MultiIter path.

Adaptation checks:

- Pointer paths require C-contiguous storage.
- Array/array path also requires identical ndim, size, and shape.
- Scalar right reads one int64 value.
- NaT on either side produces NaT before addition.
- Cython overflow checking wraps all three addition loops.
- Non-contiguous/broadcast inputs retain MultiIter.
- Tests cover scalar, same-shape NaT, non-contiguous, and overflow.

Executed:

- `git diff --cached --check`
- Conflict-marker scan.
- `python -m py_compile pandas/tests/tslibs/test_np_datetime.py`
- `python -m compileall -q pandas`
- Focused pytest invocation stopped during conftest import because
  `dateutil` is missing.

Not executed:

- Cython compilation: Cython is not installed.
- Runtime assertions: pytest did not reach collection.
- ASV: not run in this environment.

Follow-up validation:

- Run `TestAddOverflowSafe` after rebuilding extensions.
- Cover 0-D/1-D/N-D, C/F/non-contiguous, broadcast, empty, NaT, positive
  and negative overflow.
- Benchmark datetime/timedelta scalar and equal-shape arithmetic.

## Batch 21: datetime object and name construction

Sources inspected:

- Exported commits `e951e5f721` and `a22b61a327`.
- Prior pandas 3.0.3 port `766dc263de`.
- pandas 3.0.1/3.0.3 `fields.pyx` difference and upstream commits.

Adaptation checks:

- `import_datetime()` is called during module initialization.
- C constructors receive the same year/month/day/time/tz/fold fields as
  the previous Python constructors.
- Existing timezone conversion branches remain unchanged.
- Default day/month constants are already correctly cased.
- Locale-provided names are capitalized once before the output loop.
- NaT output remains `np.nan`.

Executed:

- `git diff --cached --check`
- Static C-constructor argument and locale-loop inspection.
- `python -m compileall -q pandas`
- Focused datetime pytest invocation stopped during conftest import
  because `dateutil` is missing.

Not executed:

- Cython compilation: Cython is not installed.
- Runtime datetime/locale assertions.
- ASV: not run in this environment.

Follow-up validation:

- Run datetime accessor, `to_pydatetime`, day-name, and month-name tests.
- Cover timezone-naive/aware, fold 0/1, date/time boxing, all supported
  resolutions, NaT, default locale, and non-English locales.
- Benchmark DatetimeAccessor day/month names and object conversion.

## Batch 22: month and business-day shifts

Sources inspected:

- Exported commits `9686250639` and `42c76d772f`.
- Reconstructed commits `cd5c272` and `1ab19b1`.
- Prior pandas 3.0.3 port `4d00dfadf4`.
- pandas 3.0.1/3.0.3 offsets architecture delta.

Adaptation checks:

- Direct month shifting is limited to contiguous one-dimensional input.
- Days up to 28 avoid month-length lookup; later days still clamp.
- NaT and all datetime resolutions retain conversion behavior.
- MultiIter remains for non-contiguous/N-D/named-day inputs.
- BusinessDay uses floor epoch days, corrects negative modulo, and adds
  only a whole-day delta to preserve time-of-day.
- Positive/negative periods and weekend adjustment reproduce
  `_adjust_ndays`.
- Added tests cover early/late days, NaT, non-contiguous/N-D input,
  negative epoch, positive/negative BusinessDay, and intraday values.

Executed:

- `git diff --cached --check`
- Conflict-marker and reconstructed-export comparison.
- `python -m py_compile` on both modified test files.
- `python -m compileall -q pandas`
- Focused pytest invocation stopped during conftest import because
  `dateutil` is missing.

Not executed:

- Cython compilation: Cython is not installed.
- Runtime assertions: pytest did not reach collection.
- ASV: not run in this environment.

Follow-up validation:

- Run datetime arithmetic and BusinessDay offset tests after rebuild.
- Cover all resolutions, pre/post epoch, weekend, leap day, month end,
  NaT, normalization, strides, and N-D arrays.
- Benchmark `shift_months` contiguous/fallback and BusinessDay array
  application.

## Batch 23: SemiMonth offset fast paths

Sources inspected:

- Exported commit `d1f4ed4720`.
- Prior pandas 3.0.3 port `de87482916`.
- pandas 3.0.1 SemiMonth shared array path and calendar helper.

Adaptation checks:

- Gregorian leap-year calculation matches the existing helper.
- SemiMonthBegin never needs source/destination month-end lookup.
- SemiMonthEnd checks source month end for positive roll convention and
  clamps only destination day 31.
- NaT and datetime resolution conversion remain unchanged.
- Array/scalar regression cases cover begin/end, positive/negative,
  January/February month ends, and intraday values.
- The object-construction `cast.py` rollback from the export is omitted
  and documented.

Executed:

- `git diff --cached --check`
- Static branch and calendar arithmetic inspection.
- `python -m py_compile pandas/tests/tseries/offsets/test_month.py`
- `python -m compileall -q pandas`
- Focused SemiMonth pytest invocation stopped during conftest import
  because `dateutil` is missing.

Not executed:

- Cython compilation: Cython is not installed.
- Runtime assertions: pytest did not reach collection.
- ASV: not run in this environment.

Follow-up validation:

- Run SemiMonth and calendar tests after rebuilding extensions.
- Cover every valid anchor, positive/zero/negative n, leap/century
  years, all month lengths, NaT, resolutions, and intraday timestamps.
- Benchmark SemiMonthBegin/SemiMonthEnd array application.

## Batch 24: row-wise apply label cache

Sources inspected:

- Exported `7d7fb8bf02` and the row-apply subset of `1246018d48`.
- Prior pandas 3.0.3 port `f4f5eeb92c`.
- pandas 3.0.1 Series lookup, FrameColumnApply generator, CoW refs, and
  manager set-values call sites.

Adaptation checks:

- Cache is enabled only when columns are unique.
- Cache identity is tied to the exact row Series index.
- Missing/unhashable/callable keys fall back to normal Series behavior.
- Every reused ndarray row updates `_row_apply_values`.
- Ref reset is requested only after a Series result is shallow-copied
  and consumed on the next reused row.
- EA rows continue through `_ixs` and do not use the ndarray cache.
- `SingleBlockManager.set_values` remains documented as apply-only.

Executed:

- `git diff --cached --check`
- Conflict-marker and cache-protocol inspection.
- `python -m py_compile` on all four modified Python/test files.
- `python -m compileall -q pandas`
- Focused apply pytest invocation stopped during conftest import because
  `dateutil` is missing.

Not executed:

- Runtime apply/CoW assertions: pytest did not reach collection.
- ASV: not run in this environment.

Follow-up validation:

- Run frame apply and Copy-on-Write test suites.
- Cover unique/duplicate/unhashable/callable labels, mutation, returned
  scalar/Series/view, EA rows, mixed dtype, and Series subclasses.
- Benchmark `DataFrame.apply(axis=1)` label-heavy workloads.

## Batch 25: dense object-integer value_counts

Sources inspected:

- Object-integer value_counts subset of exported `1246018d48`.
- Prior pandas 3.0.3 port `01e840cd05`.
- pandas 3.0.1 value_counts construction, sorting, and fallback paths.

Adaptation checks:

- Fast-path guards constrain semantics to one existing default result
  mode.
- Exact Python ints use overflow-aware conversion.
- Bool, NA, non-integer, and out-of-int64 values force fallback.
- Dense range allocation is bounded by input length.
- First-seen key offsets plus stable count sort preserve tie order.
- Fast-path result bypasses the second Series sort.
- Tests cover selected path and sparse-range fallback.

Executed:

- `git diff --cached --check`
- Static type/range/sort/fallback inspection.
- `python -m py_compile pandas/core/algorithms.py
  pandas/tests/base/test_value_counts.py`
- `python -m compileall -q pandas`
- Focused value_counts pytest invocation stopped during conftest import
  because `dateutil` is missing.

Not executed:

- Cython compilation: Cython is not installed.
- Runtime assertions: pytest did not reach collection.
- ASV: not run in this environment.

Follow-up validation:

- Run base/Series/Index value_counts tests after rebuilding extensions.
- Cover bool, missing, NumPy ints, subclasses, min/max int64,
  out-of-range ints, ties, sparse/dense ranges, all sort/ascending/
  normalize/dropna combinations.
- Benchmark large dense-range object integer value_counts.

## Batch 26: scalar putmask and bool-object take

Sources inspected:

- Putmask and bool-take subsets of exported `1246018d48`.
- Prior pandas 3.0.3 port `999c26f3c5`.
- pandas 3.0.1 take dispatch/wrapper API and prior migrated take batch.

Adaptation checks:

- Scalar putmask exits before list-like length/repetition logic.
- Bool take requires 2-D bool input, object output, and floating NaN
  fill.
- All-mask, partial-mask, and no-mask paths initialize every output
  position.
- Axis 0 and axis 1 assignments preserve output shape.
- F-contiguous transpose/axis restoration remains outside the helper.
- Generic fallback keeps the pandas 3.0.1 four-argument wrapper
  signature; no 3.0.3 `allow_fill` API is imported.

Executed:

- `git diff --cached --check`
- Conflict-marker and dispatch inspection.
- `python -m py_compile` on implementation and test files.
- `python -m compileall -q pandas`
- Focused putmask/take pytest invocation stopped during conftest import
  because `dateutil` is missing.

Not executed:

- Runtime assertions: pytest did not reach collection.
- ASV: not run in this environment.

Follow-up validation:

- Run array-algorithm putmask and take suites.
- Cover C/F order, both axes, all/partial/no `-1`, `allow_fill`
  true/false, empty dimensions, scalar/list-like new values, and CoW
  callers.
- Benchmark scalar putmask and 2-D bool-to-object take.

## Batch 27: Arrow integer scalar fillna

Sources inspected:

- Arrow fillna subset of exported `1246018d48`.
- Prior pandas 3.0.3 port `9ea1009920`.
- pandas 3.0.1 Arrow mask utility, boxing, and fill fallback.

Adaptation checks:

- Fast path requires one chunk, integer type, and a valid boxed scalar.
- Existing `limit` and array-like value handling occur before the helper.
- Null-free input returns the existing Arrow array.
- Primitive data is copied before replacement.
- Arrow validity uses `False` for missing positions.
- Rebuilt array has no validity buffer and retains the exact Arrow type.
- Every unsupported case returns `None` to the existing safe kernel.
- Unrelated duration helpers from the conflict were not imported.

Executed:

- `git diff --cached --check`
- Conflict-marker and narrow-diff inspection.
- `python -m py_compile` on implementation and test.
- `python -m compileall -q pandas`
- PyArrow availability check: PyArrow is not installed.
- Focused pytest invocation stopped earlier at conftest import because
  `dateutil` is missing.

Not executed:

- PyArrow runtime test: PyArrow is unavailable.
- ASV: not run in this environment.

Follow-up validation:

- Run Arrow extension fillna tests with supported minimum/current
  PyArrow versions.
- Cover all integer widths/signs, slices/offsets, null-free/all-null,
  multi-chunk, invalid/overflowing scalar, limit, non-integer types, and
  buffer lifetime.
- Benchmark single-chunk integer scalar fillna.

## Batch 28: nullable integer value_counts

Sources inspected:

- Masked integer value_counts subset of exported `1246018d48`.
- Prior pandas 3.0.3 port `d108e65555`.
- pandas 3.0.1 BaseMaskedArray dtype construction and stable sorting.

Adaptation checks:

- Fast path requires large BaseMasked integer arrays and default
  descending, non-normalized sorting.
- Negative values, sparse ranges, high NA density, and short inputs
  return to the hashtable.
- Dense allocation is bounded relative to valid count and intp range.
- NA is appended with a masked key only when `dropna=False`.
- Keys are cast to the original dtype's NumPy storage before nullable
  Index construction.
- Counts remain int64-backed nullable integers.
- First occurrence is retained for equal counts before stable sorting.
- Added regression coverage for Int8 dtype and tie order.

Executed:

- `git diff --cached --check`
- Static mask/range/dtype/order inspection.
- `python -m py_compile` on implementation and test.
- `python -m compileall -q pandas`
- Focused nullable integer pytest invocation stopped during conftest
  import because `dateutil` is missing.

Not executed:

- Runtime assertions: pytest did not reach collection.
- ASV: not run in this environment.

Follow-up validation:

- Run nullable integer function/value_counts tests.
- Cover all signed/unsigned widths, min/max values, negative fallback,
  empty/all-NA/high-NA, sparse/dense ranges, ties, dropna, ascending,
  sort=False, and normalize.
- Benchmark dense nullable integer value_counts at threshold boundaries.

## Batch 29: FrameRowApply Series reuse

Sources inspected:

- Axis=0 apply subset of exported `1246018d48`.
- Prior pandas 3.0.3 port `8afaf5bca6`.
- pandas 3.0.1 block layout, Series generator, result wrapping, and the
  Batch 24 CoW reset protocol.

Adaptation checks:

- Reuse requires exactly one non-extension block.
- Column values and name are replaced before every yield.
- Shallow Series results request a ref reset before the next reuse.
- Mixed and EA inputs remain on per-column `_ixs`.
- Fast result construction requires every Series result to share the
  exact original index object.
- Result columns are replaced by the existing result index only when
  lengths match.
- Tests cover Series-valued axis=0 output and reuse overwrite safety.

Executed:

- `git diff --cached --check`
- Static generator/result/CoW inspection.
- `python -m py_compile pandas/core/apply.py
  pandas/tests/apply/test_frame_apply.py`
- `python -m compileall -q pandas`
- Focused apply pytest invocation stopped during conftest import because
  `dateutil` is missing.

Not executed:

- Runtime apply/CoW assertions: pytest did not reach collection.
- ASV: not run in this environment.

Follow-up validation:

- Run frame apply and CoW tests.
- Cover scalar/Series/view returns, mutation, same/equal/different
  indexes, empty frames, mixed blocks, EA blocks, subclasses, and
  duplicate labels.
- Benchmark homogeneous `DataFrame.apply(axis=0)`.

## Batch 30: broad astype/fillna audit

Sources inspected:

- Remaining astype/fillna portions of exported `1246018d48`.
- Prior pandas 3.0.3 audit `429781002a`.
- pandas 3.0.1 generic astype, Block.fillna, manager apply, CoW, and
  warning paths.

Decision:

- Do not port the broad pandas 2.x DataFrame astype helper with stale
  `copy=True/False` branching.
- Do not bypass BlockManager apply/warning handling with the broad dict
  fillna helper.
- Retain existing scalar Block.fillna coverage.
- Evaluate later narrow 3.0.3 astype/fillna ports independently with
  focused CoW tests.

Executed:

- Static source and API comparison.
- `git diff --check`

Not executed:

- Runtime CoW/warning tests: pytest is blocked by missing `dateutil`.
- ASV: not run in this environment.

Follow-up validation:

- Before accepting narrow implementations, run CoW tests for
  `astype(copy=True/False)`, referenced blocks, mutation isolation,
  duplicate/MultiIndex columns, dict/Series fill values, inplace mode,
  and warnings.
- Benchmark homogeneous numeric astype and object-block dict fillna only
  after correctness validation.
