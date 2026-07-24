"""Helpers for manipulations with graphs."""

from __future__ import annotations

from collections.abc import Iterator, Set as AbstractSet
from typing import TypeVar

T = TypeVar("T")


def strongly_connected_components(
    vertices: AbstractSet[T], edges: dict[T, list[T]]
) -> Iterator[set[T]]:
    """Compute Strongly Connected Components of a directed graph.

    Args:
      vertices: the labels for the vertices
      edges: for each vertex, gives the target vertices of its outgoing edges

    Returns:
      An iterator yielding strongly connected components, each
      represented as a set of vertices.  Each input vertex will occur
      exactly once; vertices not part of a SCC are returned as
      singleton sets.

    From https://code.activestate.com/recipes/578507/.
    """
    identified: set[T] = set()
    stack: list[T] = []
    index: dict[T, int] = {}
    boundaries: list[int] = []

    def dfs(v: T) -> Iterator[set[T]]:
        index[v] = len(stack)
        stack.append(v)
        boundaries.append(index[v])

        for w in edges[v]:
            if w not in index:
                yield from dfs(w)
            elif w not in identified:
                while index[w] < boundaries[-1]:
                    boundaries.pop()

        if boundaries[-1] == index[v]:
            boundaries.pop()
            scc = set(stack[index[v] :])
            del stack[index[v] :]
            identified.update(scc)
            yield scc

    for v in vertices:
        if v not in index:
            yield from dfs(v)


def prepare_sccs(
    sccs: list[set[T]], edges: dict[T, list[T]]
) -> dict[AbstractSet[T], set[AbstractSet[T]]]:
    """Use original edges to organize SCCs in a graph by dependencies between them."""
    sccsmap = {}
    for scc in sccs:
        scc_frozen = frozenset(scc)
        for v in scc:
            sccsmap[v] = scc_frozen
    data: dict[AbstractSet[T], set[AbstractSet[T]]] = {}
    for scc in sccs:
        deps: set[AbstractSet[T]] = set()
        for v in scc:
            deps.update(sccsmap[x] for x in edges[v])
        data[frozenset(scc)] = deps
    return data


class topsort(Iterator[set[T]]):  # noqa: N801
    """Topological sort using Kahn's algorithm.

    Uses in-degree counters and a reverse adjacency list, so the total work
    is O(V + E).

    Implemented as a class rather than a generator for better mypyc
    compilation.

    Args:
      data: A map from vertices to all vertices that it has an edge
            connecting it to. NOTE: dependency sets in this data
            structure are modified in place to remove self-dependencies.
            Orphans are handled internally and are not added to `data`.

    Returns:
      An iterator yielding sets of vertices that have an equivalent
      ordering.

    Example:
      Suppose the input has the following structure:

        {A: {B, C}, B: {D}, C: {D}}

      The algorithm treats orphan dependencies as if normalized to:

        {A: {B, C}, B: {D}, C: {D}, D: {}}

      It will yield the following values:

        {D}
        {B, C}
        {A}
    """

    def __init__(self, data: dict[T, set[T]]) -> None:
        # Single pass: remove self-deps, build reverse adjacency list,
        # compute in-degree counts, detect orphans, and find initial ready set.
        in_degree: dict[T, int] = {}
        rev: dict[T, list[T]] = {}
        ready: set[T] = set()
        for item, deps in data.items():
            deps.discard(item)  # Ignore self dependencies.
            deg = len(deps)
            in_degree[item] = deg
            if deg == 0:
                ready.add(item)
            if item not in rev:
                rev[item] = []
            for dep in deps:
                if dep in rev:
                    rev[dep].append(item)
                else:
                    rev[dep] = [item]
                    if dep not in data:
                        # Orphan: appears as dependency but has no entry in data.
                        in_degree[dep] = 0
                        ready.add(dep)

        self.in_degree = in_degree
        self.rev = rev
        self.ready = ready
        self.remaining = len(in_degree) - len(ready)

    def __iter__(self) -> Iterator[set[T]]:
        return self

    def __next__(self) -> set[T]:
        ready = self.ready
        if not ready:
            assert self.remaining == 0, (
                f"A cyclic dependency exists amongst "
                f"{[k for k, deg in self.in_degree.items() if deg > 0]!r}"
            )
            raise StopIteration
        in_degree = self.in_degree
        rev = self.rev
        new_ready: set[T] = set()
        for item in ready:
            for dependent in rev[item]:
                new_deg = in_degree[dependent] - 1
                in_degree[dependent] = new_deg
                if new_deg == 0:
                    new_ready.add(dependent)
        self.remaining -= len(new_ready)
        self.ready = new_ready
        return ready
