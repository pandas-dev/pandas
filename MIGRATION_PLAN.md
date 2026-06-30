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
| 9-38 | `2f489b472c` through `4ceff7983d` | remaining audited batches | completed with pandas 3.0.1 adaptations |

## Original private commit inventory

The authoritative 89 candidates now have terminal pandas 3.0.1 dispositions.

| # | Export hash | Subject | Category | Batch | Decision | Source files |
| --- | --- | --- | --- | --- | --- | --- |
| 1 | `26f530f608` | C 实现swisstable原型 | swisstable | B1-existing-audit | migrated in Batch 3 | `asv_bench/benchmarks/swisstable.py`, `_libs/meson.build`, `_libs/swisstable*`, `core/algorithms.py`, `core/config_init.py`, tests |
| 2 | `3cbe555493` | feat(swisstable): 添加 SwissTable C++17 核心实现 | swisstable | B1-existing-audit | migrated in Batch 3 | `_libs/swisstable/swisstable.hpp` |
| 3 | `063703e782` | perf(swisstable): 添加批处理优化方法和 resize | swisstable | B1-existing-audit | migrated in Batch 3 | `_libs/swisstable/swisstable.hpp` |
| 4 | `55b5cb4acb` | feat(swisstable): 添加类型别名和实例化框架 | swisstable | B1-existing-audit | migrated in Batch 3 | `.gitignore`, `_libs/swisstable/swisstable_aliases.hpp`, instances headers |
| 5 | `0bd2e202c9` | feat(swisstable): 添加 pyx 文件结构和 SwissInt64Map 类 | swisstable | B1-existing-audit | migrated in Batch 3 | `_libs/swisstable/swisstable.pyx` |
| 6 | `bb5778c2b4` | feat(swisstable): 添加 int64 辅助函数 | swisstable | B1-existing-audit | migrated in Batch 3 | `_libs/swisstable/swisstable.pyx` |
| 7 | `e9a2ff724a` | feat(swisstable): 添加 SwissUInt64Map 类 | swisstable | B1-existing-audit | migrated in Batch 3 | `_libs/swisstable/swisstable.pyx` |
| 8 | `752d2193df` | feat(swisstable): 添加 uint64 辅助函数和 SwissInt32Map 类 | swisstable | B1-existing-audit | migrated in Batch 3 | `_libs/swisstable/swisstable.pyx` |
| 9 | `66103d5ee8` | feat(swisstable): 添加 SwissUInt32Map 类 | swisstable | B1-existing-audit | migrated in Batch 3 | `_libs/swisstable/swisstable.pyx` |
| 10 | `0f92566a5c` | feat(swisstable): 添加 SwissFloat64Map 类 | swisstable | B1-existing-audit | migrated in Batch 3 | `_libs/swisstable/swisstable.pyx` |
| 11 | `4d37bc846c` | feat(swisstable): 添加 SwissFloat32Map 类 | swisstable | B1-existing-audit | migrated in Batch 3 | `_libs/swisstable/swisstable.pyx` |
| 12 | `e7ec6a4d43` | feat(swisstable): 添加 SwissInt16Map 类 | swisstable | B1-existing-audit | migrated in Batch 3 | `_libs/swisstable/swisstable.pyx` |
| 13 | `39b3922be7` | feat(swisstable): 添加 SwissUInt16Map 类 | swisstable | B1-existing-audit | migrated in Batch 3 | `_libs/swisstable/swisstable.pyx` |
| 14 | `cc3c7cf8f8` | feat(swisstable): 添加 SwissInt8Map 类 | swisstable | B1-existing-audit | migrated in Batch 3 | `_libs/swisstable/swisstable.pyx` |
| 15 | `9028ac986d` | feat(swisstable): 添加 SwissUInt8Map 类 | swisstable | B1-existing-audit | migrated in Batch 3 | `_libs/swisstable/swisstable.pyx` |
| 16 | `67c9d4e13b` | feat(swisstable): 添加 SwissComplex64Map 类 | swisstable | B1-existing-audit | migrated in Batch 3 | `_libs/swisstable/swisstable.pyx` |
| 17 | `ad3b462701` | feat(swisstable): 添加 SwissComplex128Map 类 | swisstable | B1-existing-audit | migrated in Batch 3 | `_libs/swisstable/swisstable.pyx` |
| 18 | `0744c0555f` | optimize BaseMultiIndexCodesEngine init | index | B2-index-join | migrated in `index/search primitives` sub-batch | `_libs/index.pyx` |
| 19 | `1e811cf41e` | optimize take_1d | algorithms | B1-existing-audit | no direct migration: already covered by the pandas 3.0.1 target state | `_libs/algos_take_helper.pxi.in` |
| 20 | `02d2113331` | optimize take_1d | index/algorithms | B1-existing-audit | no direct migration: this reverts row 19, and its index change is superseded by rows 18/21 | `_libs/algos_take_helper.pxi.in`, `_libs/index.pyx` |
| 21 | `df0c75b9fd` | optimize BaseMultiIndexCodesEngine init | index | B2-index-join | migrated in `index/search primitives` sub-batch | `_libs/index.pyx` |
| 22 | `33da51b84c` | optimize take_1d | algorithms | B1-existing-audit | no direct migration: already covered by the pandas 3.0.1 target state | `_libs/algos_take_helper.pxi.in` |
| 23 | `66863407cc` | roll_sum 是 rolling 窗口求和的高频路径，核心操作模式为： | low-level/window | B7-low-level-window | migrated in roll_sum batch | `_libs/window/aggregations.pyx` |
| 24 | `2777612697` | kh put likely | low-level/window | B7-low-level-window | migrated in khash micro-optimization sub-batch | `_libs/include/pandas/vendored/klib/khash.h` |
| 25 | `598e709a9b` | group_cummin_max性能优化 | groupby | B1-existing-audit | migrated in Batch 5 | `_libs/groupby.pyx` |
| 26 | `e38ed47c0e` | group_all_any 性能优化 | groupby | B1-existing-audit | migrated in Batch 5 | `_libs/groupby.pyx` |
| 27 | `d019bc4d8e` | group_cumprod性能优化 | groupby | B1-existing-audit | migrated in Batch 5 | `_libs/groupby.pyx` |
| 28 | `63b8520f8c` | group_shift_indexer性能优化 | groupby | B1-existing-audit | migrated in Batch 5 | `_libs/groupby.pyx` |
| 29 | `171fe464a6` | take_xxx 消减引用计数增减 | algorithms | B1-existing-audit | migrated in Batch 2b | `_libs/algos.pyx`, `_libs/algos_take_helper.pxi.in` |
| 30 | `9e640a5d2c` | add_overflowsafe: 连续场景直接使用裸指针访问替代原本的PyArray_MultiIter迭代器 | tslibs | B5-tslibs | migrated in add_overflowsafe batch | `_libs/tslibs/np_datetime.pyx` |
| 31 | `a60133b265` | array_equivalent_object: 连续场景直接使用裸指针访问替代原本的PyArray_MultiIter迭代器 | lib | B3-lib-object | migrated in lib pointer-helper batch | `_libs/lib.pyx` |
| 32 | `eaaffa04b4` | is_monotonic优化 | index | B2-index-join | no direct migration: reverted by row 65 after a MergeJoin failure | `_libs/algos.pyx` |
| 33 | `fefbc7b5a4` | cython construct_1d_object_array_from_listlike impl | lib | B3-lib-object | migrated in object construction batch | `_libs/lib.pyx`, `_libs/lib.pyi`, `core/dtypes/cast.py` |
| 34 | `53b373262d` | cython _searchsorted_left impl | index | B2-index-join | migrated: typed left/right hooks and DatetimeEngine call site use the Cython binary search | `_libs/index.pyx`, `_libs/index_class_helper.pxi.in` |
| 35 | `82d3f0d7c5` | group_lastx性能优化 | groupby | B1-existing-audit | migrated in Batch 5 | `_libs/groupby.pyx` |
| 36 | `91019df988` | skiplist优化 | low-level/window | B7-low-level-window | migrated in skiplist header batch | `_libs/include/pandas/skiplist.h` |
| 37 | `eaa0ec8472` | group_mean 性能优化 | groupby | B1-existing-audit | migrated in Batch 5 | `_libs/groupby.pyx` |
| 38 | `c60600ce36` | group_min_max性能优化 | groupby | B1-existing-audit | migrated in Batch 5 | `_libs/groupby.pyx` |
| 39 | `829e42a602` | group_prod性能优化 | groupby | B1-existing-audit | migrated in Batch 5 | `_libs/groupby.pyx` |
| 40 | `4c4a5096dc` | 优化 kth_smallest_c 算法性能，使用指针操作替代数组索引，减少内存访问开销 | algorithms | B6-algorithms | migrated in algos Cython batch | `_libs/algos.pyx` |
| 41 | `12c293f745` | 默认启用swisstable | swisstable | B1-existing-audit | migrated in Batch 3 | `core/config_init.py` |
| 42 | `8c1c10d14f` | swisstable 同步蓝区代码 | swisstable | B1-existing-audit | migrated in Batch 3 | hashing/hashtable/swisstable files |
| 43 | `7d2d23951e` | fix debug build | swisstable | B1-existing-audit | migrated in Batch 3 | `_libs/swisstable.pxd` |
| 44 | `4423bd227e` | factorize: fix int na_value bug | swisstable | B1-existing-audit | migrated in Batch 3 | `_libs/swisstable_class_helper.pxi.in` |
| 45 | `72efd58130` | fix object engine performance in _get_loc_duplicates | index | B2-index-join | no direct migration: pandas3 already uses `_bin_search/_bin_search_right` | `_libs/index.pyx` |
| 46 | `a985978612` | fix stride case in swisstable | swisstable | B1-existing-audit | migrated in Batch 3 | `_libs/swisstable.pyx`, `_libs/swisstable_class_helper.pxi.in` |
| 47 | `a3739ab1e4` | 添加携程 benchmark | benchmarks | B9-benchmarks | migrated in xiecheng benchmark batch | `asv_bench/benchmarks/xiecheng.py` |
| 48 | `42b376312f` | add dataframe construction bench for xiecheng | benchmarks | B9-benchmarks | migrated in xiecheng benchmark batch | `asv_bench/benchmarks/xiecheng.py` |
| 49 | `baca9d6de2` | fix xiecheng StringCategorical time_to_categorical | benchmarks | B9-benchmarks | migrated in xiecheng benchmark batch | `asv_bench/benchmarks/xiecheng.py` |
| 50 | `403df2a143` | join: 优化 groupsort_indexer 和 take | join | B2-index-join | migrated where applicable: `groupsort_indexer` unroll retained; join take helper is obsolete after row 57 writes `sort=False` output directly | `_libs/algos.pyx`, `_libs/join.pyx` |
| 51 | `56996122c8` | add swisstable Factorizer implementation | swisstable | B1-existing-audit | migrated in Batch 3 | hashtable/swisstable/merge files |
| 52 | `120b7a79a1` | Fix benchmark script arguments and test script build process | script | B1-existing-audit | no direct migration: canceled by row 53 | `scripts/benchmark_swiss_prefetch.py`, `scripts/test_prefetch_server.sh` |
| 53 | `6c0d804544` | Revert "Fix benchmark script arguments and test script build process" | script | B1-existing-audit | no direct migration: revert marker for row 52 | deletes the scripts from row 52 |
| 54 | `108a06f167` | swisstable get labels batch | swisstable | B1-existing-audit | migrated in Batch 3 | `_libs/swisstable*` |
| 55 | `e951e5f721` | 优化 ints_to_pydatetime | tslibs | B5-tslibs | migrated in vectorized/fields tslibs batch | `_libs/tslibs/vectorized.pyx` |
| 56 | `a22b61a327` | 优化 get_date_name_field | tslibs | B5-tslibs | migrated in Batch 21 because pandas 3.0.1 lacks the 3.0.3 upstream equivalent | `_libs/tslibs/fields.pyx` |
| 57 | `07ccc6e28e` | faster many:many join with sort=False | join | B2-index-join | migrated in Batch 31 with direct original-left-order output | `_libs/join.pyx` |
| 58 | `a041e1e5b3` | group_sum优化 | groupby | B1-existing-audit | migrated with pandas3-safe mask/skipna loop specialization while retaining Kahan, overflow, object, and min_count semantics | `_libs/groupby.pyx` |
| 59 | `fcdd01e0a4` | update xiecheng bench | benchmarks | B9-benchmarks | migrated in xiecheng benchmark batch | `asv_bench/benchmarks/xiecheng.py` |
| 60 | `1246018d48` | 性能优化：优化 apply、astype、fillna、take 和 value_counts 核心执行路径 | python-layer | B8-python-layer | migrated in B8a-B8h: apply, astype, dict/Arrow fillna, value_counts, putmask, and take paths are adapted to pandas3 | `lib.pyx/.pyi`, `core/algorithms.py`, apply/putmask/take/arrow/generic/internals/series/tests |
| 61 | `51a2b98159` | perf: Use stable sort in safe_sort for better ARM performance | algorithms | B6-algorithms | migrated in stable-sort batch | `core/algorithms.py` |
| 62 | `31a63abb20` | perf: use stable sort for all argsort calls in factorize | algorithms | B6-algorithms | migrated in stable-sort batch | `core/algorithms.py` |
| 63 | `1c4c300a36` | perf: Add sort kind in safe_sort | algorithms | B6-algorithms | migrated in stable-sort batch | `core/algorithms.py` |
| 64 | `7d7fb8bf02` | PERF: optimize row-wise apply string access and label overhead | python-layer | B8-python-layer | migrated in row-wise apply cache batch | `core/apply.py`, `core/internals/managers.py`, `core/series.py`, tests |
| 65 | `0b958ce8fd` | is_monotonic 携程 MergeJoin用例失败 退回开源版本 | join/index | B2-index-join | no direct migration: pandas3 already retains the restored upstream behavior | `_libs/algos.pyx` |
| 66 | `c064f8cd15` | eq_NA_compat 优化： 避免引用增减 | lib | B3-lib-object | migrated in lib pointer-helper batch | `_libs/lib.pyx` |
| 67 | `defb42e92e` | fix array_equivalent_object build | lib | B3-lib-object | migrated with array_equivalent_object pointer casts | `_libs/lib.pyx` |
| 68 | `286eba1dc7` | groupby_idmin_max 优化 | groupby | B1-existing-audit | migrated with max/min and mask/skipna loop specialization plus a separate NA-poison state that fixes the exported sentinel ambiguity | `_libs/groupby.pyx` |
| 69 | `52f9c792ec` | join indexer: object optimize | join | B2-index-join | migrated in object join indexer batch | `_libs/join.pyx` |
| 70 | `1b285b9697` | range _concat impl in lib.pyx | index/lib | B2-index-join | migrated with pandas3 fallback and repeated-range semantics preserved | `_libs/lib.pyx`, `_libs/lib.pyi`, `core/indexes/range.py` |
| 71 | `42dc387d2c` | khash 局部缓存 keys/flags 指针 | low-level/window | B7-low-level-window | migrated in khash micro-optimization sub-batch | `_libs/include/pandas/vendored/klib/khash.h` |
| 72 | `454f5e27f1` | optimize fast_zip | lib | B3-lib-object | migrated in lib pointer-helper batch | `_libs/lib.pyx` |
| 73 | `09c956697b` | quantile 性能优化 | groupby | B1-existing-audit | adapted in Batch 31 with the stronger O(n) `kth_smallest_c` quickselect path missing from pandas 3.0.1 | `_libs/groupby.pyx` |
| 74 | `9686250639` | PERF: add early-day fast path and contiguous 1-D direct loop in shift_months(day_opt=None) | tslibs | B5-tslibs | migrated in offsets shift batch | `_libs/tslibs/offsets.pyx`, tests |
| 75 | `15dba60291` | PERF: dense get_dummies helper | reshape | B4-reshape | migrated in dense get_dummies batch | `_libs/reshape.pyi`, `_libs/reshape.pyx`, `core/reshape/encoding.py`, tests |
| 76 | `8232deb8a0` | PERF: split numeric unstack mask/value writes | reshape | B4-reshape | migrated in Batch 14 because pandas 3.0.1 still passes `new_mask` | `_libs/reshape.pyx` |
| 77 | `706991b28f` | PERF: 2D no-limit pad fast path | algorithms | B6-algorithms | migrated in Batch 2a | `_libs/algos.pyx` |
| 78 | `d72753640f` | PERF: flat list/tuple object-array construction fast path | lib/cast | B3-lib-object | migrated in object construction batch | `core/dtypes/cast.py`, tests |
| 79 | `3f77b3b691` | PERF: add_overflowsafe left-contiguous/right-scalar paths | tslibs | B5-tslibs | migrated in add_overflowsafe batch | `_libs/tslibs/np_datetime.pyx` |
| 80 | `0e2677777a` | PERF: checknull uses C `isnan()` | algorithms | B6-algorithms | migrated with pandas3 signature adaptation | `_libs/missing.pyx`, tests |
| 81 | `42c76d772f` | PERF: `_shift_bdays` avoids full date rebuild where possible | tslibs | B5-tslibs | migrated in offsets shift batch | `_libs/tslibs/offsets.pyx`, tests |
| 82 | `82081460b3` | PERF: add_overflowsafe scalar/contiguous cached flags | tslibs | B5-tslibs | migrated in add_overflowsafe batch | `_libs/tslibs/np_datetime.pyx`, tests |
| 83 | `d1f4ed4720` | PERF: SemiMonthBegin/SemiMonthEnd narrow date rebuild | tslibs | B5-tslibs | migrated relevant SemiMonth changes; intentionally rejected the unrelated cast.py rollback that would undo row 78 and delete its tests | `ccalendar.pyx`, `offsets.pyx`, tests |
| 84 | `38f97b58aa` | PERF: add pad_inplace no-limit fast path | algorithms | B6-algorithms | migrated in Batch 2a | `_libs/algos.pyx` |
| 85 | `6fd0359851` | PERF: add pad_inplace loop unroll | algorithms | B6-algorithms | migrated in Batch 2a | `_libs/algos.pyx` |
| 86 | `844af97539` | PERF: maybe_convert_objects homogeneous-block fast path | lib | B1-existing-audit | migrated in Batch 4 | `_libs/lib.pyx`, tests |
| 87 | `aae84a43d5` | PERF: avoid redundant Python-to-C complex conversion | lib | B1-existing-audit | migrated in Batch 4 | `_libs/lib.pyx` |
| 88 | `eba6d76f8c` | PERF: nancorr high-validity mask fast path | algorithms | B1-existing-audit | migrated in Batch 1 | `_libs/algos.pyx` |
| 89 | `2e964bd967` | PERF: take_2d_axis1 broader contiguous fast path | algorithms | B1-existing-audit | migrated in Batch 2b | `_libs/algos_take_helper.pxi.in` |

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
- Initial pandas 3.0.1 result: the prior port applied without textual
  conflicts, but a runtime reachability audit found that this was not
  sufficient to complete the migration.
- Compatibility checks cover the generated integer maps/factorizers,
  explicit float/complex implementations, `.pxd` declarations, C++17
  Meson target, `compute.use_swisstable`, factorize/unique/safe_sort,
  and merge factorizer selection.
- Runtime audit result:
  - `isin`, `duplicated`, and `value_counts_arraylike` still called
    khash even when `compute.use_swisstable=True`; they now dispatch to
    the matching numeric SwissTable helpers.
  - pandas 3.0.1 merge passes `uses_mask` to Factorizers and requires
    `hash_inner_join`; both interfaces were absent from the prior port.
  - the original `uniques_array` abstraction was omitted from the
    normal pandas Factorizers and merge still accessed
    `rizer.uniques.to_array()` directly; the exported abstraction is
    now restored for both khash and SwissTable Factorizers.
  - integer, float, and complex `ismember` helpers now handle empty and
    non-contiguous arrays without passing strided storage to contiguous
    batch routines.
  - the C++ header now includes the standard headers required by MSVC,
    and integer special methods accept the full `uint64` key range.
- Export discrepancy: the final C++ implementation removed the early C
  prototype's `delete`, `clear`, and `load_factor` Python APIs, but
  `test_swisstable.py` retained assertions for them. These stale tests
  are strict xfails; no incomplete tombstone deletion support is
  reintroduced.
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

## Batch 17: object join indexers

- Original private commits: `52f9c792ec` and audited
  `07ccc6e28e`.
- Prior pandas 3.0.3 port: `6e766c4650`.
- pandas 3.0.1 result: object-specific monotonic unique/left/inner/outer
  indexer helpers apply directly.
- Fast paths require the object fused specialization and C-contiguous
  arrays; all numeric and non-contiguous inputs retain generic pandas
  3.0.1 loops.
- Initial audit incorrectly treated the existing post-hoc
  `sort=False` reorder as equivalent to the direct-output optimization.
  Batch 31 corrects this by porting the stronger 3.0.3 implementation.
- Risk: Python rich-comparison exceptions, duplicate-run advancement,
  empty sides, and result reference ownership need native tests.

## Batch 18: skiplist allocation and traversal

- Original private commit: `91019df988`.
- Prior pandas 3.0.3 port: `8ece4c5407`.
- pandas 3.0.1 result: adapted after a header conflict caused by the
  target's older math compatibility helpers.
- Each node, next-pointer array, and width array now share one aligned
  allocation; node destruction therefore frees one block.
- Portable likely/unlikely/prefetch macros and cached traversal state
  reduce branch and member-load overhead in get/rank/insert/remove.
- pandas 3.0.1 `PANDAS_NAN` and `Log2` helpers are retained instead of
  importing the 3.0.3 C17 `NAN/log2` cleanup.
- Risk: recursive ref-counted destruction, allocation failure, duplicate
  ranks, and prefetch safety need native/debug builds.

## Batch 19: Xiecheng ASV workloads

- Original private commits: `a3739ab1e4`, `42b376312f`,
  `baca9d6de2`, and `fcdd01e0a4`.
- Prior pandas 3.0.3 port: `473b0e68ff`.
- pandas 3.0.1 result: benchmark-only suite applies without product-code
  dependencies.
- The final workload state covers construction, GroupBy, merge/join,
  time-series, strings/categoricals, pivot/crosstab, row apply, and
  vectorized operations.
- Later fixes for constructor input reuse, object-string columns, and
  categorical conversion into new columns are retained; duplicate
  benchmark method naming is resolved.
- Risk: parameter sizes are large and require an appropriately
  provisioned ASV machine.

## Batch 20: add_overflowsafe fast paths

- Original private commits: `9e640a5d2c`, `3f77b3b691`, and
  `82081460b3`.
- Prior pandas 3.0.3 port: `e48041d9cc`.
- pandas 3.0.1 result: implementation applies directly; test conflict
  was resolved without importing an unrelated 3.0.3 warning filter.
- Raw pointers cover C-contiguous left plus scalar right and equal-shape
  C-contiguous operands. Broadcast and non-contiguous inputs retain
  `PyArray_MultiIter`.
- Existing `@cython.overflowcheck(True)` and NaT sentinel propagation
  remain active.
- Risk: signed-overflow translation, scalar broadcasting, shape checks,
  NaT, and arbitrary strides need a native build.

## Batch 21: datetime object and name construction

- Original private commits: `e951e5f721` and `a22b61a327`.
- Prior pandas 3.0.3 port: `766dc263de`.
- pandas 3.0.1 result: CPython datetime constructors apply directly;
  date-name caching required an additional target-specific migration.
- `ints_to_pydatetime` initializes the datetime C API and calls
  `datetime_new`, `date_new`, and `time_new`.
- pandas 3.0.3 already contained upstream equivalent date-name caching,
  but pandas 3.0.1 does not. Locale names are now capitalized once when
  the names array is built, not once per result element.
- This is a documented expansion beyond the prior 3.0.3 port, based on
  the authoritative export diff.
- Risk: timezone/fold semantics, resolution conversion, locale casing,
  NaT, and CPython C-API initialization need native tests.

## Batch 22: month and business-day shifts

- Original private commits: `9686250639` and `42c76d772f`.
- Prior pandas 3.0.3 port: `4d00dfadf4`.
- pandas 3.0.1 result: re-adapted around the pre-`_DayOpt` offsets
  architecture rather than taking the 3.0.3 conflict side.
- `shift_months(day_opt=None)` uses a direct contiguous 1-D pointer loop
  and skips month-length lookup for days at most 28; arbitrary shapes,
  strides, and named day options retain MultiIter/string-helper paths.
- BusinessDay derives weekday from whole epoch days and adds whole-day
  periods to the original value, preserving the intraday remainder.
- Risk: negative epochs, weekend roll rules, negative periods,
  sub-day resolutions, leap/month ends, NaT, and non-contiguous arrays
  require native tests.

## Batch 23: SemiMonth offset fast paths

- Original private commit: `d1f4ed4720`.
- Prior pandas 3.0.3 port: `de87482916`.
- pandas 3.0.1 result: applicable ccalendar/offsets portions apply
  directly.
- Branch-light days-in-month helpers avoid a second leap-year table
  lookup. SemiMonthBegin avoids month-end checks; SemiMonthEnd computes
  them only for source month-end detection or a destination day of 31.
- The export's unrelated `cast.py` rollback is intentionally excluded
  because it would undo the migrated object-array construction fast
  path and its regression tests.
- Risk: scalar/array parity across anchors, signs, leap years, month
  ends, NaT, and intraday components requires native tests.

## Batch 24: row-wise apply label cache

- Original private commits: `7d7fb8bf02` and the row-apply subset of
  `1246018d48`.
- Prior pandas 3.0.3 port: `f4f5eeb92c`.
- pandas 3.0.1 result: adapted around the earlier
  `FrameColumnApply.series_generator` while excluding later EA row
  builder changes.
- Unique columns get a label-to-position cache on the reused row Series;
  duplicate columns retain normal Series lookup. Callable keys still go
  through `apply_if_callable`.
- The current row values are cached after manager replacement, and CoW
  refs are reset only after an applied function returned a Series that
  required a shallow copy.
- `SingleBlockManager.set_values` avoids rebuilding unchanged
  `BlockPlacement`.
- Risk: Series subclass attributes, returned views, mutation inside
  apply, duplicate labels, EA rows, and CoW reference state need runtime
  tests.

## Batch 25: dense object-integer value_counts

- Original private commit: the object-int64 value_counts subset of
  `1246018d48`.
- Prior pandas 3.0.3 port: `01e840cd05`.
- pandas 3.0.1 result: directly applied while retaining target Index
  reconstruction and the migrated stable-sort controls.
- Large object arrays use dense int64 counting only for descending,
  sorted, non-normalized, `dropna=False` requests and only when the
  integer range is no larger than the input length.
- Bool, missing, non-integer, int64-overflow, sparse-range, other sort
  modes, and short inputs fall back to the existing hashtable path.
- Stable count sorting preserves first-appearance order for ties.
- Risk: Python/NumPy integer subclasses, extreme bounds, tie ordering,
  and object Index dtype require native runtime tests.

## Batch 26: scalar putmask and bool-object take

- Original private commit: the putmask and bool-to-object take subsets
  of `1246018d48`.
- Prior pandas 3.0.3 port: `999c26f3c5`.
- pandas 3.0.1 result: putmask applies directly; take is adapted to the
  target's pre-3.0.3 four-argument dispatch wrappers.
- Scalar `putmask_without_repeat` delegates immediately to NumPy.
- Two-dimensional bool input with object output and NaN fill uses
  vectorized row/column assignment; every other dtype, dimension, fill
  value, and output mode retains existing dispatch.
- Risk: F-order axis flipping, all/partial/no fill masks, `allow_fill`,
  shape, and object output ownership need runtime tests.

## Batch 27: Arrow integer scalar fillna

- Original private commit: the Arrow fillna subset of `1246018d48`.
- Prior pandas 3.0.3 port: `9ea1009920`.
- pandas 3.0.1 result: narrowly reimplemented after rejecting a conflict
  side that included unrelated 3.0.3 duration helpers.
- Valid scalar fills on single-chunk integer arrays copy primitive data,
  replace invalid positions from the Arrow validity mask, and rebuild a
  null-free Arrow array.
- Invalid scalar, multi-chunk, non-integer, missing-buffer, limit, and
  unsupported cases retain existing `_safe_fill_null`/EA fallback.
- Risk: sliced chunks, signed/unsigned widths, scalar overflow, buffer
  lifetime, Arrow version behavior, and CoW need PyArrow runtime tests.

## Batch 28: nullable integer value_counts

- Original private commit: the masked nullable-integer value_counts
  subset of `1246018d48`.
- Prior pandas 3.0.3 port: `d108e65555`.
- pandas 3.0.1 result: adapted with additional dtype and stable tie-order
  corrections found during review.
- Large one-dimensional nonnegative masked integer arrays use bounded
  `bincount` only when range density and NA density are acceptable.
- Keys are cast back to the original nullable integer storage dtype;
  counts remain nullable Int64. First-appearance order is preserved
  before the existing stable count sort.
- Negative, sparse, high-NA, small, non-integer, alternate sort, and
  normalized requests retain existing EA/hashtable behavior.
- Risk: all nullable widths/signs, NA placement, tie ordering, and range
  thresholds need runtime and performance validation.

## Batch 29: FrameRowApply Series reuse

- Original private commit: the axis=0 FrameRowApply subset of
  `1246018d48`.
- Prior pandas 3.0.3 port: `8afaf5bca6`.
- pandas 3.0.1 result: directly applied on top of the shared conditional
  CoW reset protocol from Batch 24.
- Homogeneous non-extension blocks reuse one column Series and replace
  its values/name. Mixed-block and ExtensionArray frames retain `_ixs`
  per column.
- Same-index Series results use a values dictionary construction path;
  all other result forms retain existing wrapping and error behavior.
- Risk: returned views, mutation, index identity/equality, subclasses,
  zero columns, mixed/EA blocks, and CoW refs need runtime tests.

## Batch 30: broad astype/fillna audit

- Original private commit: remaining DataFrame astype/fillna portions of
  `1246018d48`.
- Prior pandas 3.0.3 audit: `429781002a`.
- No mechanical port of the broad pandas 2.x helpers: their explicit
  copy branching conflicts with pandas 3 Copy-on-Write lazy-copy
  semantics, and dict fillna bypasses current manager/warning flow.
- pandas 3.0.1 already has scalar `Block.fillna` Cython coverage for
  relevant float/object paths.
- Later, narrower 3.0.3 migration commits for homogeneous astype and
  object-dict fillna remain separate batches and will be re-audited
  against 3.0.1 rather than treating this decision as a permanent block.
- Risk: any implementation must prove CoW isolation, duplicate/MultiIndex
  column behavior, warning semantics, and referenced-block handling.

## Batch 31: complete index, join, and quantile hot paths

- Original private commits: `53b373262d`, superseded revert
  `02d2113331`, join/groupsort `403df2a143`, many-to-many
  `07ccc6e28e`, and quantile `09c956697b`.
- Prior pandas 3.0.3 completion commit: `77aa7184b3`.
- DatetimeEngine now routes monotonic lookup through the typed
  Int64Engine binary search migrated in Batch 7.
- pandas 3.0.1 lacked two optimizations that the prior completion commit
  considered upstream in 3.0.3:
  - Group quantile now uses the 3.0.3 O(n) quickselect implementation
    from upstream `8352ab938e`, which is stronger and more semantically
    robust than the private O(n log n) pointer/argsort rewrite.
  - Many-to-many `sort=False` joins now generate indexers directly in
    original left order using upstream `b8692948be`, eliminating the
    post-hoc take/reorder targeted by the private join helper.
- `02d2113331` remains a superseded revert; later take and typed index
  migrations define the final state.
- Risk: quantile interpolation/NA/datetimelike behavior and duplicate
  join ordering/cardinality need native tests.

## Batch 32: remaining GroupBy reductions

- Original private commits: `a041e1e5b3` and `286eba1dc7`.
- Prior pandas 3.0.3 port: `4d716d4cec`.
- pandas 3.0.1 result: directly applied around the now-retained
  quickselect quantile implementation.
- `group_sum` specializes mask and skipna loop modes while retaining
  Kahan compensation, object handling, min_count, and result masks.
- idxmin/idxmax specialize operation/mask/skipna combinations and use a
  separate poison bitmap so `skipna=False` keeps `-1` after any NA
  without conflating unseen and poisoned state.
- Risk: object/numeric/datetimelike fused paths, masks, NA poisoning,
  min_count, negative labels, and tie-first behavior need native tests.

## Batch 33: homogeneous numeric astype

- Original private commit: the astype subset of `1246018d48`.
- Prior pandas 3.0.3 port: `fb67599b4d`.
- pandas 3.0.1 result: narrow CoW-aware helpers apply after resolving a
  shared-document insertion conflict.
- Homogeneous NumPy/BaseMasked numeric conversions build arrays directly
  without Series concat. All-Arrow numeric frames cast each immutable
  ChunkedArray directly.
- Every NumPy/masked conversion requests new storage; Arrow conversion
  returns immutable cast arrays. Existing same-dtype lazy-copy shortcut
  and all mixed/non-numeric/errors-ignore fallbacks remain.
- This supersedes only the narrow astype part of Batch 30's broad audit;
  it does not restore pandas 2.x public `copy` branching.
- Risk: CoW isolation, nullable masks, duplicate columns, subclasses,
  error modes, Arrow casts, and empty/mixed frames need runtime tests.

## Batch 34: object-block dict fillna

- Original private commit: the dict-object fillna subset of
  `1246018d48`.
- Prior pandas 3.0.3 port: `8d730e0d0f`.
- pandas 3.0.1 result: the narrow manager/block implementation applies
  directly and completes the safe subset deferred by Batch 30.
- A single all-column object NumpyBlock can fill selected scalar-mapped
  columns with one mask/broadcast operation and per-column limit.
- Block `_get_refs_and_copy` enforces always-on CoW isolation for
  inplace and non-inplace calls.
- Duplicate/MultiIndex columns, multiple/non-object blocks, non-scalar
  mappings, unsupported locations, and axis=1 retain the existing
  Series/manager path and warning behavior.
- Risk: inplace/view refs, limit orientation, no-op mappings, dtype
  preservation, subclass finalization, and fallback warnings need
  runtime tests.

## Terminal audit summary

- Original private optimization candidates: 89.
- Migrated or pandas 3.0.1-adapted: 81.
- No direct code migration required: 8.
  - `1e811cf41e`, `33da51b84c`: take implementations already covered
    by the target take state.
  - `02d2113331`: source-side revert/superseded index step.
  - `eaaffa04b4`, `0b958ce8fd`: failed monotonic experiment and its
    correctness rollback.
  - `72efd58130`: ObjectEngine tuple-safe equivalent already exists.
  - `120b7a79a1`, `6c0d804544`: temporary benchmark scripts and their
    matching revert.
- Partially resolved: 0.
- Blocked: 0.

Important pandas 3.0.1 differences from the prior 3.0.3 audit:

- `8232deb8a0` numeric unstack applies because 3.0.1 still passes
  `new_mask`.
- `a22b61a327` date-name capitalization caching is absent from 3.0.1
  and was migrated.
- `07ccc6e28e` direct `sort=False` join output is absent from 3.0.1 and
  was backported with the stronger 3.0.3 implementation.
- `09c956697b` was adapted with the stronger O(n) 3.0.3 quickselect
  implementation.
- pad commits `706991b28f`, `38f97b58aa`, and `6fd0359851` required
  explicit 3.0.1 migration.
- Nullable integer value_counts additionally preserves original key
  dtype and first-appearance tie order.
- The prior take commit mixed unrelated pad changes; the 3.0.1 history
  splits them into reviewable commits.

No ASV result is recorded. Native extension builds and runtime tests
remain required before performance or release acceptance.

## Post-migration fix: SwissTable Cython type allocation

- Affected migration: `01e19fd632` (prior pandas 3.0.3 reference
  `f6d861f8a6`).
- pandas 3.0.1 enables `CYTHON_USE_TYPE_SPECS=1` only for C targets,
  while the migrated SwissTable extension is compiled as C++.
- This makes `HashTable` a heap type and Swiss map subclasses static
  types, which CPython 3.14 rejects during `pandas._libs.swisstable`
  initialization.
- Resolution: apply the Cython type-spec macro consistently to C and
  C++ targets and retain the existing class hierarchy and algorithms.
- Risk: low and build-scoped. The change affects Cython extension type
  creation but does not change SwissTable data structures or behavior.
