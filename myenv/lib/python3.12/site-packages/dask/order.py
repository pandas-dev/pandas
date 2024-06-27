from __future__ import annotations

r""" Static order of nodes in dask graph

Dask makes decisions on what tasks to prioritize both

*  Dynamically at runtime
*  Statically before runtime

Dynamically we prefer to run tasks that were just made available.  However when
several tasks become available at the same time we have an opportunity to break
ties in an intelligent way

        d
        |
    b   c
     \ /
      a

For example after we finish ``a`` we can choose to run either ``b`` or ``c``
next.  Making small decisions like this can greatly affect our performance,
especially because the order in which we run tasks affects the order in which
we can release memory, which operationally we find to have a large affect on
many computation.  We want to run tasks in such a way that we keep only a small
amount of data in memory at any given time.


Static Ordering
---------------

And so we create a total ordering over all nodes to serve as a tie breaker.  We
represent this ordering with a dictionary mapping keys to integer values.
Lower scores have higher priority.  These scores correspond to the order in
which a sequential scheduler would visit each node.

    {'a': 0,
     'c': 1,
     'd': 2,
     'b': 3}

There are several ways in which we might order our keys.  This is a nuanced
process that has to take into account many different kinds of workflows, and
operate efficiently in linear time.  We strongly recommend that readers look at
the docstrings of tests in dask/tests/test_order.py.  These tests usually have
graph types laid out very carefully to show the kinds of situations that often
arise, and the order we would like to be determined.

"""
import copy
from collections import defaultdict, deque, namedtuple
from collections.abc import Iterable, Mapping, MutableMapping
from typing import Any, Callable, Literal, NamedTuple, overload

from dask.core import get_dependencies, get_deps, getcycle, istask, reverse_dict
from dask.typing import Key


class Order(NamedTuple):
    priority: int
    critical_path: float | int


@overload
def order(
    dsk: Mapping[Key, Any],
    dependencies: Mapping[Key, set[Key]] | None = None,
    *,
    return_stats: Literal[True],
) -> dict[Key, Order]:
    ...


@overload
def order(
    dsk: Mapping[Key, Any],
    dependencies: Mapping[Key, set[Key]] | None = None,
    *,
    return_stats: Literal[False] = False,
) -> dict[Key, int]:
    ...


def order(
    dsk: Mapping[Key, Any],
    dependencies: Mapping[Key, set[Key]] | None = None,
    *,
    return_stats: bool = False,
) -> dict[Key, Order] | dict[Key, int]:
    """Order nodes in dask graph

    This produces an ordering over our tasks that we use to break ties when
    executing.  We do this ahead of time to reduce a bit of stress on the
    scheduler and also to assist in static analysis.

    This currently traverses the graph as a single-threaded scheduler would
    traverse it.

    Examples
    --------
    >>> inc = lambda x: x + 1
    >>> add = lambda x, y: x + y
    >>> dsk = {'a': 1, 'b': 2, 'c': (inc, 'a'), 'd': (add, 'b', 'c')}
    >>> order(dsk)
    {'a': 0, 'c': 1, 'b': 2, 'd': 3}
    """
    if not dsk:
        return {}  # type: ignore

    dsk = dict(dsk)
    expected_len = len(dsk)

    if dependencies is None:
        dependencies_are_copy = True
        dependencies = {k: get_dependencies(dsk, k) for k in dsk}
    else:
        # Below we're removing items from the sets in this dict
        # We need a deepcopy for that but we only want to do this if necessary
        # since this only happens for special cases.
        dependencies_are_copy = False
        dependencies = dict(dependencies)

    dependents = reverse_dict(dependencies)

    leaf_nodes = {k for k, v in dependents.items() if not v}
    root_nodes = {k for k, v in dependencies.items() if not v}

    result: dict[Key, Order | int] = {}

    # Normalize the graph by removing leaf nodes that are not actual tasks, see
    # for instance da.store where the task is merely an alias
    # to multiple keys, i.e. [key1, key2, ...,]

    # Similarly, we are removing root nodes that are pure data tasks. Those task
    # are embedded in the run_spec of a task and are not runnable. We have to
    # assign a priority but their priority has no impact on performance.
    # The removal of those tasks typically transforms the graph topology in a
    # way that is simpler to handle
    all_tasks = False
    n_removed_leaves = 0
    requires_data_task = defaultdict(set)
    while not all_tasks:
        all_tasks = True
        for leaf in list(leaf_nodes):
            if leaf in root_nodes:
                continue
            if (
                not istask(dsk[leaf])
                # Having a linear chain is fine
                and len(dependencies[leaf]) > 1
            ):
                all_tasks = False
                # Put non-tasks at the very end since they are merely aliases
                # and have no impact on performance at all
                prio = len(dsk) - 1 - n_removed_leaves
                if return_stats:
                    result[leaf] = Order(prio, -1)
                else:
                    result[leaf] = prio
                n_removed_leaves += 1
                leaf_nodes.remove(leaf)
                del dsk[leaf]
                del dependents[leaf]
                for dep in dependencies[leaf]:
                    dependents[dep].remove(leaf)
                    if not dependents[dep]:
                        leaf_nodes.add(dep)

        for root in list(root_nodes):
            if root in leaf_nodes:
                continue
            if not istask(dsk[root]) and len(dependents[root]) > 1:
                if not dependencies_are_copy:
                    dependencies_are_copy = True
                    dependencies = copy.deepcopy(dependencies)
                root_nodes.remove(root)
                for dep in dependents[root]:
                    requires_data_task[dep].add(root)
                    dependencies[dep].remove(root)
                    if not dependencies[dep]:
                        root_nodes.add(dep)
                del dsk[root]
                del dependencies[root]
                del dependents[root]

    num_needed, total_dependencies = ndependencies(dependencies, dependents)
    if len(total_dependencies) != len(dsk):
        cycle = getcycle(dsk, None)
        raise RuntimeError(
            "Cycle detected between the following keys:\n  -> %s"
            % "\n  -> ".join(str(x) for x in cycle)
        )
    assert dependencies is not None
    roots_connected, max_dependents = _connecting_to_roots(dependencies, dependents)
    leafs_connected, _ = _connecting_to_roots(dependents, dependencies)
    i = 0

    runnable_hull = set()
    reachable_hull = set()

    runnable = []

    known_runnable_paths: dict[Key, list[list[Key]]] = {}
    known_runnable_paths_pop = known_runnable_paths.pop

    crit_path_counter = 0
    scrit_path: set[Key] = set()
    _crit_path_counter_offset: int | float = 0

    _sort_keys_cache: dict[Key, tuple[int, int, int, int, str]] = {}

    def sort_key(x: Key) -> tuple[int, int, int, int, str]:
        try:
            return _sort_keys_cache[x]
        except KeyError:
            assert dependencies is not None
            _sort_keys_cache[x] = rv = (
                total_dependencies[x],
                len(dependencies[x]),
                len(roots_connected[x]),
                -max_dependents[x],
                # Converting to str is actually faster than creating some
                # wrapper class and comparisons that come this far are
                # relatively rare so we prefer fast init over fast comparison
                str(x),
            )
            return rv

    def add_to_result(item: Key) -> None:
        nonlocal crit_path_counter
        # Earlier versions recursed into this method but this could cause
        # recursion depth errors. This is the only reason for the while loop
        next_items = [item]
        nonlocal i
        while next_items:
            item = next_items.pop()
            runnable_hull.discard(item)
            reachable_hull.discard(item)
            leaf_nodes.discard(item)
            if item in result:
                continue

            while requires_data_task[item]:
                add_to_result(requires_data_task[item].pop())

            if return_stats:
                result[item] = Order(i, crit_path_counter - _crit_path_counter_offset)
            else:
                result[item] = i
            i += 1
            # Note: This is a `set` and therefore this introduces a certain
            # randomness. However, this randomness should not have any impact on
            # the final result since the `process_runnable` should produce
            # equivalent results regardless of the order in which runnable is
            # populated (not identical but equivalent)
            for dep in dependents.get(item, ()):
                num_needed[dep] -= 1
                reachable_hull.add(dep)
                if not num_needed[dep]:
                    if len(dependents[item]) == 1:
                        next_items.append(dep)
                    else:
                        runnable.append(dep)

    def _with_offset(func: Callable[..., None]) -> Callable[..., None]:
        # This decorator is only used to reduce indentation levels. The offset
        # is purely cosmetical and used for some visualizations and I haven't
        # settled on how to implement this best so I didn't want to have large
        # indentations that make things harder to read
        nonlocal _crit_path_counter_offset

        def wrapper(*args: Any, **kwargs: Any) -> None:
            nonlocal _crit_path_counter_offset
            _crit_path_counter_offset = 0.5
            try:
                func(*args, **kwargs)
            finally:
                _crit_path_counter_offset = 0

        return wrapper

    @_with_offset
    def process_runnables() -> None:
        """Compute all currently runnable paths and either cache or execute them

        This is designed to ensure we are running tasks that are free to execute
        (e.g. the result of a splitter task) not too eagerly. If we executed
        such free tasks too early we'd be walking the graph in a too wide /
        breadth first fashion that is not optimal. If instead we were to only
        execute them once they are needed for a final result, this can cause
        very high memory pressure since valuable reducers are executed too
        late.

        The strategy here is to take all runnable tasks and walk forwards until
        we hit a reducer node (i.e. a node with more than one dependency). We
        will remember/cache the path to this reducer node.
        If this path leads to a leaf or if we find enough runnable paths for a
        reducer to be runnable, we will execute the path.

        If instead of a reducer a splitter is encountered that is runnable, we
        will follow its splitter paths individually and apply the same logic to
        each branch.
        """
        while runnable:
            candidates = runnable.copy()
            runnable.clear()
            while candidates:
                key = candidates.pop()
                if key in runnable_hull or key in result:
                    continue
                if key in leaf_nodes:
                    add_to_result(key)
                    continue
                path = [key]
                branches = deque([path])
                while branches:
                    path = branches.popleft()
                    while True:
                        current = path[-1]
                        runnable_hull.add(current)
                        deps_downstream = dependents[current]
                        deps_upstream = dependencies[current]  # type: ignore
                        if not deps_downstream:
                            # FIXME: The fact that it is possible for
                            # num_needed[current] == 0 means we're doing some
                            # work twice
                            if num_needed[current] <= 1:
                                for k in path:
                                    add_to_result(k)
                        elif len(path) == 1 or len(deps_upstream) == 1:
                            if len(deps_downstream) > 1:
                                for d in sorted(deps_downstream, key=sort_key):
                                    # This ensures we're only considering splitters
                                    # that are genuinely splitting and not
                                    # interleaving
                                    if len(dependencies[d]) == 1:  # type: ignore
                                        branch = path.copy()
                                        branch.append(d)
                                        branches.append(branch)
                                break
                            runnable_hull.update(deps_downstream)
                            path.extend(sorted(deps_downstream, key=sort_key))
                            continue
                        elif current in known_runnable_paths:
                            known_runnable_paths[current].append(path)
                            if (
                                len(known_runnable_paths[current])
                                >= num_needed[current]
                            ):
                                pruned_branches: deque[list[Key]] = deque()
                                for path in known_runnable_paths_pop(current):
                                    if path[-2] not in result:
                                        pruned_branches.append(path)
                                if len(pruned_branches) < num_needed[current]:
                                    known_runnable_paths[current] = list(
                                        pruned_branches
                                    )
                                else:
                                    while pruned_branches:
                                        path = pruned_branches.popleft()
                                        for k in path:
                                            if num_needed[k]:
                                                pruned_branches.append(path)
                                                break
                                            add_to_result(k)
                        else:
                            if (
                                len(dependencies[current]) > 1  # type: ignore
                                and num_needed[current] <= 1
                            ):
                                for k in path:
                                    add_to_result(k)
                            else:
                                known_runnable_paths[current] = [path]
                        break

    # Pick strategy
    # Note: We're trying to be smart here by picking a strategy on how to
    # determine the critical path. This is not always clear and we may want to
    # consider just calculating both orderings and picking the one with less
    # pressure. The only concern to this would be performance but at time of
    # writing, the most expensive part of ordering is the prep work (mostly
    # connected roots + sort_key) which can be reused for multiple orderings.

    def use_longest_path() -> bool:
        size = 0
        # Heavy reducer / splitter topologies often benefit from a very
        # traditional critical path that expresses the longest chain of
        # tasks.
        if abs(len(root_nodes) - len(leaf_nodes)) / len(root_nodes) < 0.8:
            # If the graph stays about the same, we are checking for symmetry
            # and choose a "quickest path first" approach if the graph appears
            # to be asymmetrical
            for r in root_nodes:
                if not size:
                    size = len(leafs_connected[r])
                elif size != len(leafs_connected[r]):
                    return False

        return True

    longest_path = use_longest_path()

    def get_target() -> Key:
        raise NotImplementedError()

    if not longest_path:

        def _build_get_target() -> Callable[[], Key]:
            occurrences: defaultdict[Key, int] = defaultdict(int)
            for t in leaf_nodes:
                for r in roots_connected[t]:
                    occurrences[r] += 1
            occurences_grouped = defaultdict(set)
            for root, occ in occurrences.items():
                occurences_grouped[occ].add(root)
            occurences_grouped_sorted = {}
            for k, v in occurences_grouped.items():
                occurences_grouped_sorted[k] = sorted(v, key=sort_key, reverse=True)
            del occurences_grouped, occurrences

            def pick_seed() -> Key | None:
                while occurences_grouped_sorted:
                    key = max(occurences_grouped_sorted)
                    picked_root = occurences_grouped_sorted[key][-1]
                    if picked_root in result:
                        occurences_grouped_sorted[key].pop()
                        if not occurences_grouped_sorted[key]:
                            del occurences_grouped_sorted[key]
                        continue
                    return picked_root
                return None

            def get_target() -> Key:
                candidates = leaf_nodes
                skey: Callable = sort_key

                if runnable_hull:
                    skey = lambda k: (num_needed[k], sort_key(k))
                    candidates = runnable_hull & candidates
                elif reachable_hull:
                    skey = lambda k: (num_needed[k], sort_key(k))
                    candidates = reachable_hull & candidates

                if not candidates:
                    if seed := pick_seed():
                        candidates = leafs_connected[seed]
                    else:
                        candidates = runnable_hull or reachable_hull
                # FIXME: This can be very expensive
                return min(candidates, key=skey)

            return get_target

        get_target = _build_get_target()
    else:
        leaf_nodes_sorted = sorted(leaf_nodes, key=sort_key, reverse=False)
        get_target = leaf_nodes_sorted.pop
        del leaf_nodes_sorted

    # *************************************************************************
    # CORE ALGORITHM STARTS HERE
    #
    # 0. Nomenclature
    #
    #   - roots: Nodes that have no dependencies (i.e. typically data producers)
    #   - leafs: Nodes that have no dependents (i.e. user requested keys)
    #   - critical_path: The strategic path through the graph.
    #   - walking forwards: Starting from a root node we walk the graph as if we
    #     were to compute the individual nodes, i.e. along dependents
    #   - walking backwards: Starting from a leaf node we walk the graph in
    #     reverse direction, i.e. along dependencies
    #   - runnable: Nodes that are ready to be computed because all their
    #     dependencies are in result
    #   - runnable_hull: Nodes that could be reached and executed without
    #     "walking back". This typically means that these are tasks than can be
    #     executed without loading additional data/executing additional root
    #     nodes
    #  -  reachable_hull: Nodes that are touching the result, i.e. all nodes in
    #     reachable_hull have at least one dependency in result
    #
    # A. Build the critical path
    #
    #   To build the critical path we will use a provided `get_target` function
    #   that returns a node that is anywhere in the graph, typically a leaf
    #   node. This node is not required to be runnable. We will walk the graph
    #   backwards, i.e. from leafs to roots and append nodes to the graph as we
    #   go. The critical path is a
    #   linear path in the graph. While this is a viable strategy, it is not
    #   required for the critical path to be a classical "longest path" but it
    #   can define any route through the graph that should be considered as top
    #   priority.
    #
    #   1. Determine the target node by calling `get_target`` and append the
    #      target to the critical path stack
    #   2. Take the _most valuable_ (max given a `sort_key`) of its dependents
    #      and append it to the critical path stack. This key is the new target.
    #   3. Repeat step 2 until we reach a node that has no dependencies and is
    #      therefore runnable
    #
    # B. Walk the critical path
    #
    #   Only the first element of the critical path is an actually runnable node
    #   and this is where we're starting the sort. Strategically, this is the
    #   most important goal to achieve but since not all of the nodes are
    #   immediately runnable we have to walk back and compute other nodes first
    #   before we can unlock the critical path. This typically requires us also
    #   to load more data / run more root tasks.
    #   While walking the critical path we will also unlock non-critical tasks
    #   that could be run but are not contributing to our strategic goal. Under
    #   certain circumstances, those runnable tasks are allowed to be run right
    #   away to reduce memory pressure. This is described in more detail in
    #   `process_runnable`.
    #   Given this, the algorithm is as follows:
    #
    #   1. Pop the first element of the critical path
    #   2a. If the node is already in the result, continue
    #   2b. If the node is not runnable, we will put it back on the stack and
    #       put all its dependencies on the stack and continue with step 1. This
    #       is what we refer to as "walking back"
    #   2c. Else, we add the node to the result
    #   3.  If we previously had to walk back we will consider running
    #       non-critical tasks (by calling process_runnables)
    #   4a. If critical path is not empty, repeat step 1
    #   4b. Go back to A.) and build a new critical path given a new target that
    #       accounts for the already computed nodes.
    #
    # *************************************************************************

    critical_path: list[Key] = []
    cpath_append = critical_path.append
    scpath_add = scrit_path.add

    def path_append(item: Key) -> None:
        cpath_append(item)
        scpath_add(item)

    scpath_update = scrit_path.update
    cpath_extend = critical_path.extend

    def path_extend(items: Iterable[Key]) -> None:
        cpath_extend(items)
        scpath_update(items)

    cpath_pop = critical_path.pop
    scpath_discard = scrit_path.discard

    def path_pop() -> Key:
        item = cpath_pop()
        scpath_discard(item)
        return item

    while len(result) < expected_len:
        crit_path_counter += 1
        assert not critical_path
        assert not scrit_path

        # A. Build the critical path
        target = get_target()
        next_deps = dependencies[target]
        path_append(target)

        while next_deps:
            item = max(next_deps, key=sort_key)
            path_append(item)
            next_deps = dependencies[item]
            path_extend(next_deps)

        # B. Walk the critical path

        walked_back = False
        while critical_path:
            item = path_pop()
            if item in result:
                continue
            if num_needed[item]:
                if item in known_runnable_paths:
                    for path in known_runnable_paths_pop(item):
                        path_extend(reversed(path))
                    continue
                path_append(item)
                deps = dependencies[item].difference(result)
                unknown: list[Key] = []
                known: list[Key] = []
                k_append = known.append
                uk_append = unknown.append
                for d in sorted(deps, key=sort_key):
                    if d in known_runnable_paths:
                        k_append(d)
                    else:
                        uk_append(d)
                if len(unknown) > 1:
                    walked_back = True

                for d in unknown:
                    path_append(d)
                for d in known:
                    for path in known_runnable_paths_pop(d):
                        path_extend(reversed(path))

                del deps
                continue
            else:
                if walked_back and len(runnable) < len(critical_path):
                    process_runnables()
                add_to_result(item)
        process_runnables()

    assert len(result) == expected_len
    return result  # type: ignore


def _connecting_to_roots(
    dependencies: Mapping[Key, set[Key]], dependents: Mapping[Key, set[Key]]
) -> tuple[dict[Key, set[Key]], dict[Key, int]]:
    """Determine for every node which root nodes are connected to it (i.e.
    ancestors). If arguments of dependencies and dependents are switched, this
    can also be used to determine which leaf nodes are connected to which node
    (i.e. descendants).

    Also computes a weight that is defined as (cheaper to compute here)

            `max(len(dependents[k]) for k in connected_roots[key])`

    """
    result = {}
    current = []
    num_needed = {k: len(v) for k, v in dependencies.items() if v}
    max_dependents = {}
    roots = set()
    for k, v in dependencies.items():
        if not v:
            # Note: Hashing the full keys is relatively expensive. Hashing
            # integers would be much faster so this could be sped up by just
            # introducing a counter here. However, the order algorithm is also
            # sometimes interested in the actual keys and the only way to
            # benefit from the speedup of using integers would be to convert
            # this back on demand which makes the code very hard to read.
            roots.add(k)
            result[k] = {k}
            deps = dependents[k]
            max_dependents[k] = len(deps)
            for child in deps:
                num_needed[child] -= 1
                if not num_needed[child]:
                    current.append(child)
    while current:
        key = current.pop()
        for parent in dependents[key]:
            num_needed[parent] -= 1
            if not num_needed[parent]:
                current.append(parent)
        # At some point, all the roots are the same, particularly for dense
        # graphs. We don't want to create new sets over and over again
        new_set = None
        identical_sets = True
        result_first = None

        for child in dependencies[key]:
            r_child = result[child]
            if not result_first:
                result_first = r_child
                max_dependents[key] = max_dependents[child]
            # This clause is written such that it can circuit break early
            elif not (  # type: ignore[unreachable]
                identical_sets
                and (result_first is r_child or r_child.issubset(result_first))
            ):
                identical_sets = False
                if not new_set:
                    new_set = result_first.copy()
                max_dependents[key] = max(max_dependents[child], max_dependents[key])
                new_set.update(r_child)

        assert new_set is not None or result_first is not None
        result[key] = new_set or result_first

    # The order algo doesn't care about this but this makes it easier to
    # understand and shouldn't take that much time
    empty_set: set[Key] = set()
    for r in roots:
        result[r] = empty_set
    return result, max_dependents


def ndependencies(
    dependencies: Mapping[Key, set[Key]], dependents: Mapping[Key, set[Key]]
) -> tuple[dict[Key, int], dict[Key, int]]:
    """Number of total data elements on which this key depends

    For each key we return the number of tasks that must be run for us to run
    this task.

    Examples
    --------
    >>> inc = lambda x: x + 1
    >>> dsk = {'a': 1, 'b': (inc, 'a'), 'c': (inc, 'b')}
    >>> dependencies, dependents = get_deps(dsk)
    >>> num_dependencies, total_dependencies = ndependencies(dependencies, dependents)
    >>> sorted(total_dependencies.items())
    [('a', 1), ('b', 2), ('c', 3)]

    Returns
    -------
    num_dependencies: Dict[key, int]
    total_dependencies: Dict[key, int]
    """
    num_needed = {}
    result = {}
    for k, v in dependencies.items():
        num_needed[k] = len(v)
        if not v:
            result[k] = 1

    num_dependencies = num_needed.copy()
    current: list[Key] = []
    current_pop = current.pop
    current_append = current.append

    for key in result:
        for parent in dependents[key]:
            num_needed[parent] -= 1
            if not num_needed[parent]:
                current_append(parent)
    while current:
        key = current_pop()
        result[key] = 1 + sum(result[child] for child in dependencies[key])
        for parent in dependents[key]:
            num_needed[parent] -= 1
            if not num_needed[parent]:
                current_append(parent)
    return num_dependencies, result


OrderInfo = namedtuple(
    "OrderInfo",
    (
        "order",
        "age",
        "num_data_when_run",
        "num_data_when_released",
        "num_dependencies_freed",
    ),
)


def diagnostics(
    dsk: MutableMapping[Key, Any],
    o: Mapping[Key, int] | None = None,
    dependencies: MutableMapping[Key, set[Key]] | None = None,
) -> tuple[dict[Key, OrderInfo], list[int]]:
    """Simulate runtime metrics as though running tasks one at a time in order.

    These diagnostics can help reveal behaviors of and issues with ``order``.

    Returns a dict of `namedtuple("OrderInfo")` and a list of the number of outputs held over time.

    OrderInfo fields:
    - order : the order in which the node is run.
    - age : how long the output of a node is held.
    - num_data_when_run : the number of outputs held in memory when a node is run.
    - num_data_when_released : the number of outputs held in memory when the output is released.
    - num_dependencies_freed : the number of dependencies freed by running the node.
    """
    if dependencies is None:
        dependencies, dependents = get_deps(dsk)
    else:
        dependents = reverse_dict(dependencies)
    assert dependencies is not None
    if o is None:
        o = order(dsk, dependencies=dependencies, return_stats=False)

    pressure = []
    num_in_memory = 0
    age = {}
    runpressure = {}
    releasepressure = {}
    freed = {}
    num_needed = {key: len(val) for key, val in dependents.items()}
    for i, key in enumerate(sorted(dsk, key=o.__getitem__)):
        pressure.append(num_in_memory)
        runpressure[key] = num_in_memory
        released = 0
        for dep in dependencies[key]:
            num_needed[dep] -= 1
            if num_needed[dep] == 0:
                age[dep] = i - o[dep]
                releasepressure[dep] = num_in_memory
                released += 1
        freed[key] = released
        if dependents[key]:
            num_in_memory -= released - 1
        else:
            age[key] = 0
            releasepressure[key] = num_in_memory
            num_in_memory -= released

    rv = {
        key: OrderInfo(
            val, age[key], runpressure[key], releasepressure[key], freed[key]
        )
        for key, val in o.items()
    }
    return rv, pressure


def _f() -> None:
    ...


def _convert_task(task: Any) -> Any:
    if istask(task):
        assert callable(task[0])
        new_spec: list[Any] = []
        for el in task[1:]:
            if isinstance(el, (str, int)):
                new_spec.append(el)
            elif isinstance(el, tuple):
                if istask(el):
                    new_spec.append(_convert_task(el))
                else:
                    new_spec.append(el)
            elif isinstance(el, list):
                new_spec.append([_convert_task(e) for e in el])
        return (_f, *new_spec)
    elif isinstance(task, tuple):
        return (_f, task)
    else:
        return (_f, *task)


def sanitize_dsk(dsk: MutableMapping[Key, Any]) -> dict:
    """Take a dask graph and replace callables with a dummy function and remove
    payload data like numpy arrays, dataframes, etc.
    """
    new = {}
    for key, values in dsk.items():
        new_key = key
        new[new_key] = _convert_task(values)
    if get_deps(new) != get_deps(dsk):
        # The switch statement in _convert likely dropped some keys
        raise RuntimeError("Sanitization failed to preserve topology.")
    return new
