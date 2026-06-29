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
